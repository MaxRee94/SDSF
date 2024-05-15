#pragma once
#include "dispersal.h"


class Dynamics {
public:
	Dynamics() = default;
	Dynamics(
		int _timestep, float _cell_width, float _self_ignition_factor, float _rainfall, float _seed_bearing_threshold,
		float _growth_rate_multiplier, float _unsuppressed_flammability, float _min_suppressed_flammability, float _max_suppressed_flammability,
		float _radius_suppr_flamm_min, float radius_range_suppr_flamm, float _max_dbh, float _saturation_threshold, float _fire_resistance_argmin,
		float _fire_resistance_argmax, float _fire_resistance_stretch, float _background_mortality, map<string, map<string, float>> _strategy_distribution_params,
		float _resource_grid_relative_size, float _mutation_rate, int _verbosity
	) :
		timestep(_timestep), cell_width(_cell_width), unsuppressed_flammability(_unsuppressed_flammability),
		self_ignition_factor(_self_ignition_factor), rainfall(_rainfall), seed_bearing_threshold(_max_dbh* _seed_bearing_threshold),
		growth_rate_multiplier(_growth_rate_multiplier), radius_suppr_flamm_min(_max_dbh* _radius_suppr_flamm_min),
		flamm_delta_radius((_cell_width* _min_suppressed_flammability - _cell_width * _max_suppressed_flammability) / (_max_dbh * radius_range_suppr_flamm)),
		max_suppressed_flammability(_cell_width* _max_suppressed_flammability),
		min_suppressed_flammability(_cell_width* _min_suppressed_flammability),
		max_dbh(_max_dbh), verbosity(_verbosity), saturation_threshold(1.0f / _saturation_threshold), fire_resistance_argmin(_fire_resistance_argmin),
		fire_resistance_argmax(_fire_resistance_argmax), fire_resistance_stretch(_fire_resistance_stretch), background_mortality(_background_mortality),
		strategy_distribution_params(_strategy_distribution_params), resource_grid_relative_size(_resource_grid_relative_size), mutation_rate(_mutation_rate)
	{
		time = 0;
		help::init_RNG();
		pop = &state.population;
		grid = &state.grid;
	};
	void init_state(int gridsize, float radius_q1, float radius_q2, float _seed_mass) {
		state = State(
			gridsize, cell_width, max_dbh, radius_q1, radius_q2, seed_bearing_threshold, _seed_mass, saturation_threshold, strategy_distribution_params,
			mutation_rate
		);
		linear_disperser = Disperser();
		wind_disperser = WindDispersal();
		animal_dispersal = AnimalDispersal();
		neighbor_offsets = state.neighbor_offsets;
	}
	void update() {
		// Prepare next iteration
		time++;
		printf("Time: %i\n", time);
		if (verbosity > 0) printf("Resetting state distr... \n");
		grid->reset_state_distr();

		// Do simulation
		Timer timer; timer.start();
		if (verbosity > 0) printf("Beginning dispersal... \n");
		if (time > 0) disperse();
		timer.stop();
		if (verbosity > 0) printf("Dispersal took %f seconds. Beginning burn... \n", timer.elapsedSeconds());
		timer.start();
		burn();
		timer.stop();
		if (verbosity > 0) printf("Percolation took %f seconds. Beginning growth... \n", timer.elapsedSeconds());
		grow();
		timer.start();
		induce_background_mortality();
		timer.stop();
		if (verbosity > 0) printf("Inducing background mortality took %f seconds. \n", timer.elapsedSeconds());
		timer.start();
		int pre_thinning_popsize = pop->size();
		//thin_crowds();
		timer.stop();
		if (verbosity > 0) printf("Thinning crowds by removing %i trees took %f seconds.\n", pre_thinning_popsize - pop->size(), timer.elapsedSeconds());

		// Do post-simulation cleanup and data reporting
		printf("Repopulating grid... \n");
		state.repopulate_grid(verbosity);
		if (verbosity > 0) printf("Redoing grid count... \n");
		grid->redo_count();
		report_state();
	}
	void report_state() {
		int grid_memory_size = 0;
		for (int i = 0; i < grid->no_cells; i++) {
			int element_size = grid->distribution[i].trees.capacity();
			grid_memory_size += element_size;
		}

		printf("// STATE at t=%i: Tree cover: %f, #trees: %s \n", time, grid->get_tree_cover(), help::readable_number(pop->size()).c_str());
		/*printf("----- MEMORY REPORT: ------\n- Trees memory size: %i \n", (int)((pop->members.size() * sizeof(Tree)) >> 10));
		printf("- Kernels memory size: %i \n", (int)(pop->kernels_individual.size() * sizeof(pair<int, Kernel>)) >> 10);
		printf("- Crops memory size: %i \n", (pop->crops.size() * sizeof(pair<int, Crop>)) >> 10);
		printf("- Grid memory size: %i \n", (grid_memory_size >> 10));
		printf("------ END REPORT ------\n");*/
		if (verbosity == 2) for (auto& [id, tree] : pop->members) if (id % 500 == 0) printf("Radius of tree %i : %f \n", id, tree.radius);
	}
	void free() {
		delete[] neighbor_offsets;
		grid->free();
		pop->free();
		resource_grid.free();
	}
	vector<float> get_ordered_fire_ignition_times() {
		int i = 0;
		float fire_count = round((self_ignition_factor * rainfall * (float)grid->no_savanna_cells * grid->cell_area) / (float)1e6 );
		if (fire_count < 1) fire_count = help::get_rand_float(0, 1) < fire_count; // If fire_count is less than 1, we stochastically decide whether to ignite a fire or not.
		vector<float> fire_ignition_times = {};
		while (i < fire_count) {
			float t_start = help::get_rand_float(0.0, 1.0) + time;
			fire_ignition_times.push_back(t_start);
			i++;
		}
		std::sort(fire_ignition_times.begin(), fire_ignition_times.end());
		return fire_ignition_times;
	}
	void grow() {
		for (auto& [id, tree] : pop->members) {
			bool became_reproductive = tree.grow(seed_bearing_threshold);
			if (became_reproductive) {
				pop->add_reproduction_system(tree);
			}
		}
	}
	/*void thin_crowds() {
		PairSet population_sorted_by_height;
		pop->sort_by_trait("height", population_sorted_by_height);
		for (auto& [id, _] : population_sorted_by_height) {
			Tree* tree = pop->get(id);
			float shade = grid->compute_shade_on_individual_tree(tree);
			if (shade > 1.0f) {
				kill_tree(tree);
			}
		}
	}*/
	void set_global_linear_kernel(float lin_diffuse_q1, float lin_diffuse_q2, float min, float max) {
		global_kernels["linear"] = Kernel(1, lin_diffuse_q1, lin_diffuse_q2, min, max);
		pop->add_kernel("linear", global_kernels["linear"]);
		cout << "Global kernel created (Linear diffusion). " << endl;
	}
	void set_global_wind_kernel(float wspeed_gmean, float wspeed_stdev, float wind_direction, float wind_direction_stdev) {
		global_kernels["wind"] = Kernel(1, grid->width_r * 2.0f, wspeed_gmean, wspeed_stdev, wind_direction, wind_direction_stdev);
		pop->add_kernel("wind", global_kernels["wind"]);
		cout << "Global kernel created (Wind dispersal). " << endl;
	}
	void set_global_animal_kernel(map<string, map<string, float>>& animal_kernel_params) {
		global_kernels["animal"] = Kernel(1, "animal");
		init_resource_grid(animal_kernel_params);
		animal_dispersal.animals = Animals(&state, animal_kernel_params);
		animal_dispersal.animals.initialize_population();
		pop->add_kernel("animal", global_kernels["animal"]);
		cout << "Global kernel created (Dispersal by animals). \n";
	}
	void set_global_kernels(map<string, map<string, float>> nonanimal_kernel_params, map<string, map<string, float>> animal_kernel_params) {
		map<string, float> params;
		params = nonanimal_kernel_params["linear"];
		set_global_linear_kernel(params["lin_diffuse_q1"], params["lin_diffuse_q2"], params["min"], params["max"]);
		params = nonanimal_kernel_params["wind"];
		set_global_wind_kernel(params["wspeed_gmean"], params["wspeed_stdev"], params["wind_direction"], params["wind_direction_stdev"]);
		set_global_animal_kernel(animal_kernel_params);
	}
	void init_resource_grid(map<string, map<string, float>>& animal_kernel_params) {
		int resource_grid_no_cells_x = round((float)grid->width * resource_grid_relative_size);
		float resource_grid_cell_width = grid->width_r / (float)resource_grid_no_cells_x;
		vector<string> species = {};
		for (auto& [_species, _] : animal_kernel_params) species.push_back(_species);
		resource_grid = ResourceGrid(&state, resource_grid_no_cells_x, resource_grid_cell_width, species);
	}
	bool global_kernel_exists(string type) {
		return global_kernels.find(type) != global_kernels.end();
	}
	bool ensure_kernel_exists(int id) {
		if (pop->get_kernel(id)->id == -1) {
			string tree_dispersal_vector = pop->get_crop(id)->strategy.vector;
			if (global_kernel_exists(tree_dispersal_vector))
				pop->add_kernel(tree_dispersal_vector, global_kernels[tree_dispersal_vector]);
			else {
				printf("\n\n-------------- ERROR: No global kernel found for tree dispersal vector %s. \n\n", tree_dispersal_vector.c_str());
				//exit(1); Let's not exit the program for now
				return false;
			}
		}
		return true;
	}
	void disperse_wind_seeds_and_init_fruits(int& no_seed_bearing_trees, int& total_no_seeds, int& wind_trees) {
		int pre_dispersal_popsize = pop->size();
		int no_wind_seeds = 0;
		Timer timer; timer.start();
		for (auto& [id, tree] : pop->members) {
			// Get crop and kernel
			if (tree.life_phase < 2) continue;
			no_seed_bearing_trees++;
			if (id == -1 || tree.id == -1) {
				pop->remove(id);
				continue;
			}
			Crop* crop = pop->get_crop(id);
			if (crop->id == -1) {
				pop->remove(id);
				continue;
			}
			bool kernel_exists = ensure_kernel_exists(id);
			if (!kernel_exists) {
				pop->remove(id);
				continue;
			}
			crop->update(tree);
			total_no_seeds += crop->no_seeds;

			// Add fruit crop or disperse seeds, depending on dispersal vector type
			if (pop->get_kernel(id)->type == "animal") {
				resource_grid.add_crop(tree.position, crop);
			}
			else if (pop->get_kernel(id)->type == "wind") {
				pop->get_kernel(id)->update(tree.height);
				wind_disperser.disperse_crop(crop, &state);
				wind_trees++;
				no_wind_seeds += crop->no_seeds;
			}
			else {
				linear_disperser.disperse_crop(crop, &state);
			}
		}
		timer.stop(); printf(
			"-- Dispersing %s wind-dispersed seeds and initializing %s fruits took %f seconds. \n",
			help::readable_number(no_wind_seeds).c_str(), help::readable_number(resource_grid.total_no_fruits).c_str(), timer.elapsedSeconds()
		);
	}
	void disperse_animal_seeds() {
		int no_animal_seeds = 0;
		Timer timer; timer.start();
		if (resource_grid.has_fruits) {
			no_animal_seeds = animal_dispersal.disperse(&state, &resource_grid, 1);
		}
		timer.stop(); printf("-- Dispersing %s animal seeds took %f seconds. \n", help::readable_number(no_animal_seeds).c_str(), timer.elapsedSeconds());
	}
	void recruit() {
		Timer timer; timer.start();
		int pre_recruitment_popsize = pop->size();
		for (int i = 0; i < grid->no_cells; i++) {
			Cell* cell = &grid->distribution[i];
			if (cell->seedling_present) {
				Tree* tree = pop->add(
					grid->get_random_location_within_cell(cell->idx), &pop->get_crop(cell->largest_stem.second)->strategy
				);
				cell->add_tree(tree);
			}
		}

		int no_recruits = pop->size() - pre_recruitment_popsize;
		timer.stop(); printf("-- Recruitment of %s trees took %f seconds. \n", help::readable_number(no_recruits).c_str(), timer.elapsedSeconds());
	}
	void disperse() {
		resource_grid.reset();
		int pre_dispersal_popsize = pop->size();
		seeds_dispersed = 0;
		int no_seed_bearing_trees = 0;
		int no_wind_trees = 0;
		disperse_wind_seeds_and_init_fruits(no_seed_bearing_trees, seeds_dispersed, no_wind_trees);
		//disperse_animal_seeds();
		recruit();

		if (verbosity > 0) printf("-- Fraction of trees that are seed-bearing: %f, #seeds (all): %i\n", (float)no_seed_bearing_trees / (float)pop->size(), seeds_dispersed);
		if (verbosity > 0) printf("-- Proportion wind dispersed trees: %f \n", no_wind_trees / (float)no_seed_bearing_trees);
	}
	void induce_background_mortality() {
		for (auto& [id, tree] : pop->members) {
			if (help::get_rand_float(0, 1) < background_mortality) {
				state.population.remove(&tree);
			}
		}
	}
	int* get_resource_grid_colors(string species, string type) {
		return resource_grid.get_color_distribution(species, type);
	}
	void burn() {
		if (verbosity == 2) printf("Updated tree flammabilities.\n");
		vector<float> fire_ignition_times = get_ordered_fire_ignition_times();
		int no_burned_cells = 0;
		int re_ignitions = 0;
		fire_spatial_extent = 0;
		for (int i = 0; i < fire_ignition_times.size(); i++) {
			Cell* sav_cell = grid->get_random_savanna_cell();
			int _no_burned_cells = percolate(sav_cell, fire_ignition_times[i]);
			no_burned_cells += _no_burned_cells;
			if (_no_burned_cells <= 1 && re_ignitions < 5) { // If fire did not spread beyond ignition point, re-do percolation.
				re_ignitions++;
				i--; continue;
			}
		}
		fire_spatial_extent = ((float)no_burned_cells * grid->cell_area) / (float)fire_ignition_times.size();
		if (verbosity > 0) {
			printf("-- Cells burned: %i \n", no_burned_cells);
			printf("-- Number of fires: %i \n", (int)fire_ignition_times.size());
		}
	}
	float get_forest_flammability(Cell* cell, float fire_free_interval) {
		float fuel_load = cell->get_fuel_load();
		return fire_free_interval * unsuppressed_flammability * fuel_load; // We assume flammability is directly proportional to fuel load and fire-free interval.
	}
	float get_savanna_flammability(float fire_free_interval) {
		return fire_free_interval * unsuppressed_flammability; // We assume grass flammability is directly proportional to fire-free interval.
	}
	float get_cell_flammability(Cell* cell, float fire_free_interval) {
		if (cell->state == 1) {
			return get_forest_flammability(cell, fire_free_interval);
		}
		else return get_savanna_flammability(fire_free_interval);
	}
	bool tree_dies(Tree* tree, float fire_free_interval) {
		// if (verbosity == 2) printf("stem diameter: %f cm, bark thickness: %f mm, survival probability: %f \n", dbh, bark_thickness, survival_probability);
		// COMMENT: We currently assume topkill always implies death, but resprouting should also be possible. (TODO: make death dependent on fire-free interval)
		return !tree->survives_fire(fire_resistance_argmin, fire_resistance_argmax, fire_resistance_stretch);
	}
	void kill_tree(Tree* tree, float time_last_fire, queue<Cell*>& queue) {
		if (verbosity == 2) printf("Burning tree %i ... \n", tree->id);
		grid->burn_tree_domain(tree, queue, time_last_fire);
		pop->remove(tree);
	}
	void kill_tree(Tree* tree) {
		if (verbosity == 2) printf("Killing tree %i ... \n", tree->id);
		grid->kill_tree_domain(tree, false);
		pop->remove(tree);
	}
	void induce_tree_mortality(Cell* cell, float fire_free_interval, queue<Cell*>& queue) {
		for (auto tree_id : cell->trees) {
			Tree* tree = pop->get(tree_id);
			if (tree->last_mortality_check == time) continue; // Skip mortality evaluation if this was already done in the current timestep.
			if (tree_dies(tree, fire_free_interval)) {
				kill_tree(tree, cell->time_last_fire, queue);
			}
			else tree->last_mortality_check = time;
		}
	}
	void kill_saplings(Cell* cell) {
		for (auto tree_id : cell->trees) {
			pop->remove(tree_id);
			cell->remove_tree(tree_id);
		}
	}
	inline bool cell_will_ignite(Cell* cell, float t_start) {
		if (t_start - cell->time_last_fire < 10e-4) {
			return false; // Do not ignite cells which have already been burned by the current fire.
		}
		return help::get_rand_float(0.0, 1.0) < get_cell_flammability(cell, min(t_start - cell->time_last_fire, 1));
	}
	inline void burn_cell(Cell* cell, float t_start, queue<Cell*>& queue) {
		cell->time_last_fire = t_start;
		grid->state_distribution[grid->pos_2_idx(cell->pos)] = -5;
		if (cell->state == 1) {
			induce_tree_mortality(cell, t_start - cell->time_last_fire, queue);
		}
		else {
			kill_saplings(cell);
		}
	}
	int percolate(Cell* cell, float t_start) {
		std::queue<Cell*> queue;
		burn_cell(cell, t_start, queue);
		queue.push(cell);
		int no_burned_cells = 1;
		if (verbosity == 2) printf("Percolating fire...\n");
		pop_size = state.population.size();
		while (!queue.empty()) {
			Cell* cell = queue.front();
			queue.pop();

			// Percolate to neighbors
			for (int i = 0; i < 8; i++) {
				Cell* neighbor = grid->get_cell_at_position(cell->pos + neighbor_offsets[i]);
				if (cell_will_ignite(neighbor, t_start)) {
					burn_cell(neighbor, t_start, queue);
					queue.push(neighbor);
					no_burned_cells++;
				}
			}
		}
		return no_burned_cells;
	}
	float* get_firefree_intervals() {
		float* histo = new float[grid->no_cells];
		for (int i = 0; i < grid->no_cells; i++) {
			float interval = time - grid->distribution[i].time_last_fire;
			histo[i] = interval;
		}
		return histo;
	}
	float unsuppressed_flammability = 0;
	float min_suppressed_flammability = 0;
	float max_suppressed_flammability = 0;
	float self_ignition_factor = 0;
	float rainfall = 0;
	float seed_bearing_threshold = 0;
	float growth_rate_multiplier = 0;
	float radius_suppr_flamm_min = 0;
	float flamm_delta_radius = 0;
	float max_dbh = 0;
	float cell_width = 0;
	float fire_spatial_extent = 0;
	float saturation_threshold = 0;
	float fire_resistance_argmin = 0;
	float fire_resistance_argmax = 0;
	float fire_resistance_stretch = 0;
	float background_mortality = 0;
	float resource_grid_relative_size = 0;
	float mutation_rate = 0;
	int timestep = 0;
	int time = 0;
	int pop_size = 0;
	int verbosity = 0;
	int seeds_dispersed = 0;
	State state;
	Population* pop = 0;
	Grid* grid = 0;
	Disperser linear_disperser;
	WindDispersal wind_disperser;
	AnimalDispersal animal_dispersal;
	ResourceGrid resource_grid;
	pair<int, int>* neighbor_offsets = 0;
	map<string, Kernel> global_kernels;
	map<string, map<string, float>> strategy_distribution_params;
	Animals animals;
};

