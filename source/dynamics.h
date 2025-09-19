#pragma once
#include "dispersal.h"


class Dynamics {
public:
	Dynamics() = default;
	Dynamics(Dynamics&&) = default;
	Dynamics& operator=(Dynamics&&) = default;
	Dynamics(
		int _timestep, float _cell_width, float _self_ignition_factor, float _rainfall, float _seed_bearing_threshold,
		float _growth_rate_multiplier, float _unsuppressed_flammability, float _max_dbh, float _saturation_threshold, map<string, float> fire_resistance_params,
		float _background_mortality, map<string, map<string, float>> _strategy_distribution_params,
		int _resource_grid_width, float _mutation_rate, float _STR, int _verbosity, int random_seed, int firefreq_random_seed, float __enforce_no_recruits,
		int _animal_group_size
	) :
		timestep(_timestep), cell_width(_cell_width), unsuppressed_flammability(_unsuppressed_flammability),
		self_ignition_factor(_self_ignition_factor), rainfall(_rainfall), seed_bearing_threshold(_max_dbh * _seed_bearing_threshold),
		growth_rate_multiplier(_growth_rate_multiplier), 
		max_dbh(_max_dbh), seedling_discard_dbh(0.05 * _max_dbh), verbosity(_verbosity), saturation_threshold(1.0f / _saturation_threshold),
		fire_resistance_argmin(fire_resistance_params["argmin"]), fire_resistance_argmax(fire_resistance_params["argmax"]), 
		fire_resistance_stretch(fire_resistance_params["stretch"]),
		background_mortality(_background_mortality), strategy_distribution_params(_strategy_distribution_params),
		resource_grid_width(_resource_grid_width), mutation_rate(_mutation_rate), STR(_STR), _enforce_no_recruits(__enforce_no_recruits), 
		animal_group_size(_animal_group_size)
	{
		time = 0;
		help::init_RNG(random_seed);
		random_generator = default_random_engine(firefreq_random_seed);
	};
	void init_state(int gridsize, float dbh_q1, float dbh_q2, float growth_multiplier_stdev, float growth_multiplier_min, float growth_multiplier_max,
		float minimum_patch_size) {
		state = State(
			gridsize, cell_width, max_dbh, dbh_q1, dbh_q2, seed_bearing_threshold, saturation_threshold, strategy_distribution_params,
			mutation_rate, growth_multiplier_stdev, growth_multiplier_min, growth_multiplier_max, minimum_patch_size
		);
		linear_disperser = Disperser();
		wind_disperser = WindDispersal();
		animal_dispersal = AnimalDispersal();
		neighbor_offsets = state.grid.neighbor_offsets;
		pop = &state.population;
		grid = &state.grid;
		fire_free_interval_averages = make_shared<float[]>(grid->no_cells);
		for (int i = 0; i < grid->no_cells; i++) fire_free_interval_averages[i] = 0;
	}
	bool invalid_tree_ids() {
		for (auto& [id, tree] : pop->members) {
			if (id == -1 || tree.id == -1) {
				printf("Key: %i, Tree id: %i \n", id, tree.id);
				return true;
			}
		}
		return false;
	}
	void update() {
		// Prepare next iteration
		time++;
		printf("\nTime: %i\n", time);
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

		if (verbosity > 0) printf("Burns took %f seconds. Beginning growth... \n", timer.elapsedSeconds());
		timer.start();
		grow();
		timer.stop();
		if (verbosity > 0) printf("Growth took %f seconds.\n", timer.elapsedSeconds());
		induce_background_mortality();
		if (verbosity > 0) printf("Induced background mortality. Repopulating grid...\n");

		update_forest_patch_detection();

		// Do post-simulation cleanup and data reporting
		state.repopulate_grid(verbosity);
		if (verbosity > 1) printf("Redoing grid count... \n");
		grid->redo_count();
		report_state();

		update_firefree_interval_averages();
	}
	void update_forest_patch_detection() {
		// Obtain forest patches, report their sizes, their perimeter lengths, and the sizes of the neighboring savanna areas.
		if (verbosity > 0) printf("Updating forest patch detection... \n");
		state.grid.get_patches(verbosity);
	}
	void report_state() {
		int grid_memory_size = 0;
		for (int i = 0; i < grid->no_cells; i++) {
			int element_size = grid->distribution[i].trees.capacity();
			grid_memory_size += element_size;
		}

		printf("Tree cover: %f, Number of trees: %s \n", grid->get_tree_cover(), help::readable_number(pop->size()).c_str());
		if (verbosity == 2) for (auto& [id, tree] : pop->members) if (id % 500 == 0) printf("Radius of tree %i : %f \n", id, tree.radius);
	}
	void free() {
		pop->free();
	}
	void update_firefree_interval_averages() {
		shared_ptr<float[]> fire_free_intervals = get_firefree_intervals("current_iteration");
		for (int i = 0; i < grid->no_cells; i++) {
			if (grid->distribution[i].state == 1) continue; // Skip forest cells
			float cur_interval = fire_free_intervals[i];
			float avg_interval = fire_free_interval_averages[i];
			fire_free_interval_averages[i] = (avg_interval * ((float)time - 1) + cur_interval) / (float)time;
		}
	}
	void grow() {
		vector<int> tree_deletion_schedule = {};
		map<float, int> increment_counts;
		for (auto& [id, tree] : pop->members) {
			float shade = state.compute_shade_on_individual_tree(&tree);
			tree.shade = shade;
			auto [became_reproductive, dies_due_to_light_limitation] = tree.grow(seed_bearing_threshold, shade);
			if (dies_due_to_light_limitation) {
				tree_deletion_schedule.push_back(id);
			}
		}
		for (int id : tree_deletion_schedule) {
			pop->remove(id);
		}
		printf("-- No trees dead due to light limitation: %i \n", tree_deletion_schedule.size());
	}
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
		animal_dispersal.animals = Animals(& state, animal_kernel_params, animal_group_size);
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
		float resource_grid_cell_width = grid->width_r / (float)resource_grid_width;
		vector<string> species = {};
		for (auto& [_species, _] : animal_kernel_params) species.push_back(_species);
		if (resource_grid.selection_probabilities.probabilities != nullptr) printf("not nullptr");
		resource_grid = ResourceGrid(&state, resource_grid_width, resource_grid_cell_width, species, animal_kernel_params);
		printf("resource grid's probmodel obj id: %i \n", resource_grid.selection_probabilities.id);
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
				if (id != -1) printf("No global kernel found for tree dispersal vector '%s'.\n", tree_dispersal_vector.c_str());
				//exit(1); Let's not exit the program for now
				return false;
			}
		}
		return true;
	}
	void disperse_wind_seeds_and_init_fruits(int& no_seed_bearing_trees, int& no_wind_seedlings, int& wind_seeds_dispersed, int& animal_seeds_dispersed, int& wind_trees) {
		int pre_dispersal_popsize = pop->size();
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
			bool kernel_exists = ensure_kernel_exists(tree.id);
			if (!kernel_exists) {
				pop->remove(id);
				continue;
			}
			crop->update(tree, STR);

			// Add fruit crop or disperse seeds, depending on dispersal vector type
			if (pop->get_kernel(id)->type == "animal") {
				animal_seeds_dispersed += crop->no_seeds;
				resource_grid.has_fruits = true;
			}
			else if (pop->get_kernel(id)->type == "wind") {
				wind_seeds_dispersed += crop->no_seeds;
				int enforce_no_recruits = -1;
				if (_enforce_no_recruits >= 0) enforce_no_recruits = (float)crop->no_seeds * _enforce_no_recruits; // Enforce a certain fraction of the number of produced seeds to be recruited.
				pop->get_kernel(id)->update(tree.height);
				wind_disperser.disperse_crop(
					crop, &state, no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees,
					no_germination_attempts, no_cases_seedling_competition_and_shading, no_cases_oldstem_competition_and_shading,
					enforce_no_recruits, no_wind_seedlings
				);
				wind_trees++;
			}
			else {
				int enforce_no_recruits = -1;
				if (_enforce_no_recruits >= 0) enforce_no_recruits = (float)crop->no_seeds * _enforce_no_recruits; // Enforce a certain fraction of the number of produced seeds to be recruited.
				linear_disperser.disperse_crop(
					crop, &state, no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees,
					no_germination_attempts, no_cases_seedling_competition_and_shading, no_cases_oldstem_competition_and_shading,
					enforce_no_recruits, no_wind_seedlings
				);
			}
		}
		timer.stop(); printf(
			"-- Dispersing %s wind-dispersed seeds and initializing %s fruits took %f seconds. \n",
			help::readable_number(wind_seeds_dispersed).c_str(), help::readable_number(resource_grid.total_no_fruits).c_str(), timer.elapsedSeconds()
		);
	}
	void disperse_animal_seeds(int no_seeds_to_disperse, int& no_recruits) {
		Timer timer; timer.start();
		int enforce_no_recruits = -1;
		if (enforce_no_recruits >= 0) enforce_no_recruits = (float)no_seeds_to_disperse * enforce_no_recruits; // Enforce a certain fraction of the number of produced seeds to be recruited.)
		if (resource_grid.has_fruits) {
			no_recruits = animal_dispersal.disperse(
				&state, &resource_grid, no_seeds_to_disperse, fraction_time_spent_moving, no_seedlings_dead_due_to_shade, no_seedling_competitions, 
				no_competitions_with_older_trees, no_germination_attempts, no_cases_seedling_competition_and_shading, no_cases_oldstem_competition_and_shading, 
				enforce_no_recruits, 1
			);
		}
		timer.stop(); printf("-- Dispersing %s animal seeds took %f seconds. \n", help::readable_number(no_seeds_to_disperse).c_str(), timer.elapsedSeconds());
	}
	void recruit() {
		Timer timer; timer.start();
		int pre_recruitment_popsize = pop->size();
		for (int i = 0; i < grid->no_cells; i++) {
			Cell* cell = &grid->distribution[i];
			if (cell->seedling_present) {
				Tree* tree = pop->add(
					grid->get_real_cell_position(cell),
					&pop->get_crop(cell->stem.second)->strategy
				);
				cell->insert_sapling(tree, grid->cell_area, grid->cell_halfdiagonal_sqrt);
				grid->state_distribution[i] = -7;
			}
		}

		no_recruits = pop->size() - pre_recruitment_popsize;
		if (time == 1) initial_no_effective_dispersals = no_recruits; // The number of recruits is really the number of effective dispersals, since some seedlings may be burned right after germinating.
		timer.stop();
		printf("-- Recruitment of %s trees took %f seconds. \n", help::readable_number(no_recruits).c_str(), timer.elapsedSeconds());
	}
	int get_no_recruits(string type) {
		if (type == "wind") return no_wind_seedlings;
		else if (type == "animal") return no_animal_seedlings;
		else return (no_recruits);
	}
	void disperse() {
		resource_grid.reset();
		int pre_dispersal_popsize = pop->size();
		int animal_seeds_dispersed = 0;
		int wind_seeds_dispersed = 0;
		int no_seed_bearing_trees = 0;
		int no_wind_trees = 0;
		no_wind_seedlings = 0;
		no_cases_oldstem_competition_and_shading = 0;
		no_cases_seedling_competition_and_shading = 0;
		no_seedlings_dead_due_to_shade = 0;
		no_competitions_with_older_trees = 0;
		no_germination_attempts = 0;
		no_seedling_competitions = 0;
		no_animal_seedlings = 0;
		disperse_wind_seeds_and_init_fruits(no_seed_bearing_trees, no_wind_seedlings, wind_seeds_dispersed, animal_seeds_dispersed, no_wind_trees);
		disperse_animal_seeds(animal_seeds_dispersed, no_animal_seedlings);
		recruit();

		if (verbosity > 0) {
			printf(
				"-- Fraction of trees that are seed-bearing: %f, #seeds (all): %s\n",
				(float)no_seed_bearing_trees / (float)pop->size(), help::readable_number(seeds_produced).c_str()
			);
		}
		if (verbosity > 0) printf("-- Proportion wind dispersed trees: %f \n", no_wind_trees / (float)no_seed_bearing_trees);
	}
	void induce_background_mortality() {
		vector<int> tree_deletion_schedule = {};
		for (auto& [id, tree] : pop->members) {
			if (help::get_rand_float(0, 1) < background_mortality) {
				tree_deletion_schedule.push_back(id);
			}
		}
		for (int id : tree_deletion_schedule) {
			pop->remove(id);
		}
		printf("-- Number of trees dead due to background mortality: %i \n", tree_deletion_schedule.size());
	}
	shared_ptr<int[]> get_resource_grid_colors(string species, string type) {
		return resource_grid.get_color_distribution(species, type);
	}
	int get_no_fires() {
		// Get number of fires to ignite
		if (self_ignition_factor == -1) {
			return 1;
		}
		std::binomial_distribution<int> no_fires_distribution(grid->no_savanna_cells, self_ignition_factor / 1e6);
		return no_fires_distribution(random_generator);
	}
	void burn() {
		if (verbosity == 2) printf("Updated tree flammabilities.\n");
		int no_fires = get_no_fires();
		int no_ash_cells = 0;
		int popsize_before_burns = pop->size();
		int re_ignitions = 0;
		no_fire_induced_topkills = 0;
		no_fire_induced_nonseedling_topkills = 0;
		fires.clear();
		for (int i = 0; i < no_fires; i++) {
			Cell* cell = grid->get_random_cell();
			if (cell->time_last_fire == time) {
				no_fires--;
				continue;
			}
			no_fires++;
			auto [_no_ash_cells, _no_grassy_ash_cells] = percolate(cell, time, no_fire_induced_topkills, no_fire_induced_nonseedling_topkills);
			no_ash_cells += _no_ash_cells;
			fires.push_back((float)_no_ash_cells * grid->cell_area);
		}
		no_fire_induced_deaths = popsize_before_burns - pop->size();
		printf("no fire induced deaths (time = %i): %i \n", time, no_fire_induced_deaths);
		printf("no fire induced topkills (time = %i): %i \n", time, no_fire_induced_topkills);
		cout.precision(2);
		printf(
			"-- Fires: %i, Topkills: %s, Kills: %s \n",
			no_fires, help::readable_number(no_fire_induced_topkills).c_str(), help::readable_number(no_fire_induced_deaths).c_str()
		);
		cout <<
			"-- Fraction of domain burned: " << (float)no_ash_cells / (float)grid->no_cells << ", Area burned: " <<
			scientific << (float)no_ash_cells * grid->cell_area << " / " << grid->area << " m^2 \n";
		cout << fixed;
	}
	float get_forest_flammability(Cell* cell, bool grass_has_recovered) {
		float fuel_load = cell->get_fuel_load();
		return grass_has_recovered * unsuppressed_flammability * fuel_load; // We assume forest flammability is directly proportional to fuel load
	}
	float get_savanna_flammability(bool grass_has_recovered) {
		return grass_has_recovered * unsuppressed_flammability; // We assume grass flammability is directly proportional to fire-free interval.
	}
	float get_cell_flammability(Cell* cell, bool grass_has_recovered) {
		if (cell->state == 1) {
			return get_forest_flammability(cell, grass_has_recovered);
		}
		else return get_savanna_flammability(grass_has_recovered);
	}
	bool tree_is_topkilled(Tree* tree) {
		// if (verbosity == 2) printf("stem diameter: %f cm, bark thickness: %f mm, survival probability: %f \n", dbh, bark_thickness, survival_probability);
		// COMMENT: We currently assume topkill always implies death, but resprouting should also be possible. (TODO: make death dependent on fire-free interval)

		if (tree->dbh < seedling_discard_dbh) return true; // We assume that seedlings with a dbh below this 'discard'-value are always killed by fire.
		return !tree->survives_fire(fire_resistance_argmin, fire_resistance_argmax, fire_resistance_stretch);
	}
	void kill_tree(Tree* tree, float time_last_fire, queue<Cell*>& queue, Cell* cell) {
		if (verbosity > 1) printf("Burning tree %i ... \n", tree->id);
		Cell* stem_cell = grid->burn_tree_domain(tree, queue, time_last_fire, true, true, cell->idx);
		if (tree->life_phase == 2 || tree->life_phase == 0) {
			// A tree which has been burned once is allowed to resprout, in line with findings of Hoffmann et al. (2012).
			tree->resprout(seed_bearing_threshold);
			stem_cell->resprout_present = true;
		}
		else {
			// Resprouts are not allowed to resprout again, since previous work appears to indicate that forest species are not able to consistently recover from repeated burning (Fensham et al 2003)
			bool removed = pop->remove(tree->id);
			if (!removed) {
				printf("Tree %i could not be removed from the population. \n", tree->id);
			}
		}
	}
	void kill_tree(Tree* tree) {
		if (verbosity == 2) printf("Removing tree %i ... \n", tree->id);
		grid->kill_tree_domain(tree, false);
		pop->remove(tree);
	}
	void induce_tree_mortality(Cell* cell, queue<Cell*>& queue, int& no_trees_topkilled, int& no_fire_induced_nonseedling_topkills) {
		int tree_id = cell->stem.second;
		//printf("idx: %i , no trees: %i\n", idx, cell->trees.size());
		if (tree_id == 0) return; // If no tree stem is present in this cell, skip mortality evaluation.

		Tree* tree = pop->get(tree_id);
		vector<int> trees = cell->trees;
		if (tree->id == -1) {
			//printf("\n\n ------- ERROR: Tree %i has been removed from the population but is still present in cell %i, %i. \n", tree_id, cell->pos.first, cell->pos.second);
			//printf("Trees in cell before starting this mortality loop: ");
			//help::print_vector(&trees);
			//bool present = state.check_grid_for_tree_presence(tree_id, 0);
			return;
		}
		if (tree->last_mortality_check == time) return; // Skip mortality evaluation if this was already done in the current timestep.
		if (tree_is_topkilled(tree)) {
			no_trees_topkilled++;
			if (tree->age > -1) no_fire_induced_nonseedling_topkills++;
			kill_tree(tree, cell->time_last_fire, queue, cell);
		}
		else tree->last_mortality_check = time;
	}
	inline bool cell_will_ignite(Cell* cell, float t_start) {
		if (t_start - cell->time_last_fire < 10e-4) {
			return false; // Do not ignite cells which have already been burned by the current fire.
		}
		return help::get_rand_float(0.0, 1.0) < get_cell_flammability(cell, min(t_start - cell->time_last_fire, 1));
	}
	inline void burn_cell(Cell* cell, float t_start, queue<Cell*>& queue, int& no_trees_topkilled, int& no_fire_induced_nonseedling_topkills) {
		cell->time_last_fire = t_start;
		grid->state_distribution[grid->pos_2_idx(cell->pos)] = -5;
		induce_tree_mortality(cell, queue, no_trees_topkilled, no_fire_induced_nonseedling_topkills);
	}
	pair<int, int> percolate(Cell* cell, float t_start, int& no_trees_topkilled, int& no_fire_induced_nonseedling_topkills) {
		std::queue<Cell*> queue;
		burn_cell(cell, t_start, queue, no_trees_topkilled, no_fire_induced_nonseedling_topkills);
		queue.push(cell);
		int no_ash_cells = 1;
		int no_grassy_ash_cells = 1;
		if (verbosity == 2) printf("Percolating fire...\n");
		pop_size = state.population.size();
		while (!queue.empty()) {
			Cell* cell = queue.front();
			queue.pop();

			// Percolate to neighbors
			for (int i = 0; i < 8; i++) {
				Cell* neighbor = grid->get_cell_at_position(cell->pos + neighbor_offsets[i]);
				if (cell_will_ignite(neighbor, t_start)) {
					burn_cell(neighbor, t_start, queue, no_trees_topkilled, no_fire_induced_nonseedling_topkills);
					queue.push(neighbor);
					no_ash_cells++;
					if (neighbor->state == 0) no_grassy_ash_cells++;
				}
			}
		}
		return pair<int, int>(no_ash_cells, no_grassy_ash_cells);
	}
	shared_ptr<float[]> get_firefree_intervals(string type = "current_iteration") {
		if (type == "current_iteration") {
			shared_ptr<float[]> intervals = make_shared<float[]>(grid->no_cells);
			for (int i = 0; i < grid->no_cells; i++) {
				if (grid->distribution[i].state == 1) continue; // Skip forest cells
				float interval = time + 1 - grid->distribution[i].time_last_fire;
				intervals[i] = interval;
			}
			return intervals;
		}
		else if (type == "average") {
			return fire_free_interval_averages;
		}
	}
	float unsuppressed_flammability = 0;
	shared_ptr<float[]> fire_free_interval_averages;
	float min_suppressed_flammability = 0;
	float max_suppressed_flammability = 0;
	float self_ignition_factor = 0;
	float rainfall = 0;
	float seed_bearing_threshold = 0;
	float growth_rate_multiplier = 0;
	float radius_suppr_flamm_min = 0;
	float flamm_delta_radius = 0;
	float max_dbh = 0;
	float seedling_discard_dbh = 0;
	float cell_width = 0;
	vector<float> fires;
	float saturation_threshold = 0;
	float fire_resistance_argmin = 0;
	float fire_resistance_argmax = 0;
	float fire_resistance_stretch = 0;
	float background_mortality = 0;
	float fraction_time_spent_moving = 0;
	float mutation_rate = 0;
	float STR = 0;
	float _enforce_no_recruits = 0;
	int no_cases_seedling_competition_and_shading = 0;
	int no_cases_oldstem_competition_and_shading = 0;
	int no_germination_attempts = 0;
	int no_competitions_with_older_trees = 0;
	int no_seedlings_dead_due_to_shade = 0;
	int no_seedling_competitions = 0;
	int no_wind_seedlings = 0;
	int no_animal_seedlings = 0;
	int no_recruits = 0;
	int timestep = 0;
	int time = 0;
	int pop_size = 0;
	int verbosity = 0;
	int seeds_produced = 0;
	int resource_grid_width = 0;
	int animal_group_size = 0;
	int no_fire_induced_deaths = 0;
	int no_fire_induced_topkills = 0;
	int no_fire_induced_nonseedling_topkills = 0;
	int initial_no_effective_dispersals = 0;
	State state;
	Population* pop = 0;
	Grid* grid = 0;
	Disperser linear_disperser;
	WindDispersal wind_disperser;
	AnimalDispersal animal_dispersal;
	default_random_engine random_generator;
	ResourceGrid resource_grid;
	shared_ptr<pair<int, int>[]> neighbor_offsets = 0;
	map<string, Kernel> global_kernels;
	map<string, map<string, float>> strategy_distribution_params;
	Animals animals;
};

