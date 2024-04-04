#pragma once
#include "dispersal.h"


class Dynamics {
public:
	Dynamics() = default;
	Dynamics(
		int _timestep, float _cellsize, float _self_ignition_factor, float _rainfall, float _seed_bearing_threshold, float _mass_budget_factor,
		float _growth_rate_multiplier, float _unsuppressed_flammability, float _min_suppressed_flammability, float _max_suppressed_flammability,
		float _radius_suppr_flamm_min, float radius_range_suppr_flamm, float _max_radius, int _verbosity
	) :
		timestep(_timestep), cellsize(_cellsize), unsuppressed_flammability(_unsuppressed_flammability),
		self_ignition_factor(_self_ignition_factor), rainfall(_rainfall), seed_bearing_threshold(_seed_bearing_threshold),
		mass_budget_factor(_mass_budget_factor), growth_rate_multiplier(_growth_rate_multiplier),
		radius_suppr_flamm_min(_max_radius * _radius_suppr_flamm_min),
		flamm_d_radius((_cellsize * _min_suppressed_flammability - _cellsize * _max_suppressed_flammability) / (_max_radius * radius_range_suppr_flamm)),
		max_suppressed_flammability(_cellsize* _max_suppressed_flammability),
		min_suppressed_flammability(_cellsize* _min_suppressed_flammability),
		max_radius(_max_radius), verbosity(_verbosity)
	{
		time = 0;
		help::init_RNG();
		pop = &state.population;
		grid = &state.grid;
	};
	void init_state(int gridsize, float radius_q1, float radius_q2, float _seed_mass) {
		state = State(gridsize, cellsize, max_radius, radius_q1, radius_q2, seed_bearing_threshold, mass_budget_factor, _seed_mass);
		disperser = Disperser(&state);
		neighbor_offsets = state.neighbor_offsets;
	}
	void update() {
		// Prepare next iteration
		time++;
		printf("Time: %i\n", time);
		grid->reset_state_distr();

		// Do simulation
		disperse();
		grow();
		burn();

		// Do post-simulation cleanup and data reporting
		state.repopulate_grid(0);
		if (verbosity > 0) report_state();
	}
	void report_state() {
		printf("- Tree cover: %f, #trees: %i \n", grid->get_tree_cover(), pop->size());
		if (time == 1 && global_kernel.id != -1 && verbosity > 0) {
			Crop* crop = pop->get_crop(1);
			printf("total biomass/tree: %f \n", crop->mass);
			printf("seed mass: %f \n", crop->seed_mass);
			printf("no seeds: %i \n", crop->no_seeds);
		}
		if (verbosity == 2) for (auto& [id, tree] : pop->members) if (id % 500 == 0) printf("Radius of tree %i : %f \n", id, tree.radius);
	}
	vector<float> get_ordered_fire_ignition_times() {
		int i = 0;
		int fire_count = round(self_ignition_factor * rainfall * (float)grid->area / (float)1e6);
		vector<float> fire_ignition_times = {};
		while (i < fire_count) {
			float t_start = help::get_rand_float(0.0, 1.0) + time;
			fire_ignition_times.push_back(t_start);
			i++;
		}
		std::sort(fire_ignition_times.begin(), fire_ignition_times.end());
		return fire_ignition_times;
	}
	void grow_tree(Tree &tree) {
		// TEMP: constant growth rate. TODO: Make dependent on radius and life phase (resprout or seedling)
		tree.radius_tmin1 = tree.radius;
		tree.radius = min(tree.radius + sqrtf(tree.radius) * growth_rate_multiplier, pop->max_radius);
	}
	void grow() {
		for (auto& [id, tree] : pop->members) {
			grow_tree(tree);
		}
	}
	void set_global_kernel(float lin_diffuse_q1, float lin_diffuse_q2, float min, float max) {
		global_kernel = Kernel(1, lin_diffuse_q1, lin_diffuse_q2, min, max);
		cout << "Global kernel created. " << endl;
	}
	void disperse() {
		int x = 0;
		int j = 0;
		for (auto& [id, tree] : pop->members) {
			if (tree.radius < (seed_bearing_threshold * pop->max_radius)) continue;
			x++;
			Crop* crop = pop->get_crop(id);
			if (crop->kernel == nullptr) {
				if (global_kernel.id != -1) crop->kernel = pop->add_kernel(id, global_kernel);
			}
			crop->update(tree);
			for (int i = 0; i < crop->no_seeds; i++) {
				disperser.disperse(crop);
				j++;
			}
		}
		if (verbosity > 0) printf("Number of seed bearing trees: %i, #seeds dispersed: %i \n", x, j);
	}
	void burn() {
		update_tree_flammabilities();
		if (verbosity == 2) printf("Updated tree flammabilities.\n");
		vector<float> fire_ignition_times = get_ordered_fire_ignition_times();
		int no_burned_cells = 0;
		int re_ignitions = 0;
		for (int i = 0; i < fire_ignition_times.size(); i++) {
			Cell* sav_cell = grid->get_random_savanna_cell();
			int _no_burned_cells = percolate(sav_cell, fire_ignition_times[i]);
			no_burned_cells += _no_burned_cells;
			if (_no_burned_cells <= 1 && re_ignitions < 5) { // If fire did not spread beyond ignition point, re-do percolation.
				re_ignitions++;
				i--;
			}
		}
		if (verbosity > 0) {
			printf("Cells burned: %i \n", no_burned_cells);
			printf("Number of fires: %i \n", fire_ignition_times.size());
		}
	}
	float get_suppressed_flammability(Tree* tree, float fire_free_interval) {
		float cumulative_radius = tree->radius;
		float d_radius = max(cumulative_radius - radius_suppr_flamm_min, 0);
		return max(max_suppressed_flammability + d_radius * flamm_d_radius, min_suppressed_flammability);
	}
	float get_unsuppressed_flammability(float fire_free_interval) {
		return fire_free_interval * unsuppressed_flammability;
	}
	float get_flammability(Cell* cell, float fire_free_interval) {
		if (cell->state == 1) {
			return get_suppressed_flammability(pop->get(cell->tree), fire_free_interval);
		}
		else return get_unsuppressed_flammability(fire_free_interval);
	}
	void update_tree_flammabilities() {
		vector<Tree*> trees = {};
		int j = 0;
		for (int i = 0; i < grid->no_cells; i++) {
			if (grid->distribution[i].state == 0) continue;

			// TODO: Once functionality to store multiple trees per cell is implemented, add a check here to only compute flammability
			// in case there is only one tree in the cell.
			pop->get(grid->distribution[i].tree)->flammability = get_suppressed_flammability(pop->get(grid->distribution[i].tree), 1);
			if (verbosity == 2 && i % 1000 == 0) printf("Flammability value: %f \n", pop->get(grid->distribution[i].tree)->flammability);
		}
	}
	bool tree_dies(Cell* cell, float fire_free_interval) {
		return 1; // TEMP: trees always die if ignited (TODO: make dependent on fire-free interval and fire resistance)
	}
	void store_tree_colors(bool reset=true) {
		if (reset) grid->reset_state_distr();
		vector<Tree*> trees = {};
		int no_internally_stored_cells = 0;
		int no_grid_counted_cells = 0;
		for (int i = 0; i < grid->no_cells; i++) {
			if (grid->distribution[i].state == 0) continue;
			Cell cell = grid->distribution[i];
			if (cell.tree != 0) no_grid_counted_cells++;
			bool new_tree = find(trees.begin(), trees.end(), pop->get(cell.tree)) == trees.end();
			if (new_tree) {
				trees.push_back(pop->get(cell.tree));
			}
			int j = 0;
			while (j < trees.size()) {
				if (pop->get(cell.tree) == trees[j]) {
					grid->state_distribution[i] = (float)trees[j]->id;
				}
				j++;
			}
		}
	}
	void kill_tree(Cell* cell, queue<Cell*>& queue) {
		Tree* tree = pop->get(cell->tree);
		if (verbosity == 2) printf("Killing tree %i ... \n", tree->id);

		// Obtain cells within radius of killed tree that are not part of neighboring trees.
		vector<Cell*> burned_cells = {};
		grid->get_cells_within_radius(tree, &burned_cells); // Get all cells within killed tree radius
		vector<Tree*> neighbors = state.get_tree_neighbors(tree);
		for (auto& neighbor : neighbors) {
			grid->get_cells_within_radius(neighbor, &burned_cells, true); // Exclude cells that fall within radius of neighbors
		}

		// Burn cells and add to queue
		for (Cell* _cell : burned_cells) {
			queue.push(_cell);
			grid->set_to_savanna(_cell->pos, cell->time_last_fire);
			grid->state_distribution[grid->pos_2_idx(_cell->pos)] = -5;
		}

		// Change tree occupancy for overlapped cells to that of neighbors.
		for (auto& neighbor : neighbors) {
			grid->populate_tree_domain(neighbor, true);
		}

		state.population.remove(tree);
	}
	void induce_mortality(Cell* cell, float fire_free_interval, queue<Cell*>& queue) {
		if (tree_dies(cell, fire_free_interval)) {
			kill_tree(cell, queue);
		}
	}
	inline bool cell_will_ignite(Cell* cell, float t_start) {
		if (t_start - cell->time_last_fire < 10e-4) {
			return false; // Do not ignite cells which have already been burned by the current fire.
		}
		return help::get_rand_float(0.0, 1.0) < get_flammability(cell, min(t_start - cell->time_last_fire, 1));
	}
	inline void burn_cell(Cell* cell, float t_start, queue<Cell*> &queue) {
		cell->time_last_fire = t_start;
		grid->state_distribution[grid->pos_2_idx(cell->pos)] = -5;
		if (cell->state == 1) {
			induce_mortality(cell, t_start - cell->time_last_fire, queue);
		}
	}
	int percolate(Cell* cell, float t_start) {
		std::queue<Cell*> queue;
		queue.push(cell);
		burn_cell(cell, t_start, queue);
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
					queue.push(neighbor);
					burn_cell(neighbor, t_start, queue);
					no_burned_cells++;
				}
			}
		}
		return no_burned_cells;
	}
	float unsuppressed_flammability = 0;
	float min_suppressed_flammability = 0;
	float max_suppressed_flammability = 0;
	float mass_budget_factor = 0;
	float self_ignition_factor = 0;
	float rainfall = 0;
	float seed_bearing_threshold = 0;
	float growth_rate_multiplier = 0;
	float radius_suppr_flamm_min = 0;
	float flamm_d_radius = 0;
	float max_radius = 0;
	float cellsize = 0;
	int timestep = 0;
	int time = 0;
	int pop_size = 0;
	int verbosity = 0;
	State state;
	Population* pop = 0;
	Grid* grid = 0;
	Disperser disperser;
	pair<int, int>* neighbor_offsets = 0;
	Kernel global_kernel;
};

