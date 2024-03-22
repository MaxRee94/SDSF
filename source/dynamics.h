#pragma once
#include "state.h"


class Dynamics {
public:
	Dynamics() = default;
	Dynamics(int _timestep, float __unsuppressed_flammability, float __suppressed_flammability, float _self_ignition_factor, float _rainfall) :
		timestep(_timestep), _unsuppressed_flammability(__unsuppressed_flammability), _suppressed_flammability(__suppressed_flammability),
		self_ignition_factor(_self_ignition_factor), rainfall(_rainfall)
	{
		time = 0;
		help::init_RNG();
		init_neighbor_indices();
	};
	void init_neighbor_indices() {
		neighbor_indices = new pair<int, int>[8];
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if (i == 0 && j == 0) continue;
				neighbor_indices[i * 3 + j] = pair<int, int>(i, j);
			}
		}
	}
	void init_state(int gridsize, float cellsize, float max_radius) {
		state = State(gridsize, cellsize, max_radius);
	}
	void update(bool verbose = false) {
		time++;
		printf("Updating... Time: %i\n", time);
		//printf("Pop size: %i \n", state.population.size());
		state.repopulate_grid(verbose);
		simulate_fires(verbose);
	}
	vector<float> get_ordered_fire_ignition_times() {
		int i = 0;
		int fire_count = round(self_ignition_factor * rainfall);
		vector<float> fire_ignition_times = {};
		while (i < fire_count) {
			float t_start = help::get_rand_float(0.0, 1.0) + time;
			fire_ignition_times.push_back(t_start);
			i++;
		}
		std::sort(fire_ignition_times.begin(), fire_ignition_times.end());
		return fire_ignition_times;
	}
	void simulate_fires(bool verbose = false) {
		vector<float> fire_ignition_times = get_ordered_fire_ignition_times();
		for (int i = 0; i < fire_ignition_times.size(); i++) {
			Cell* sav_cell = state.grid.get_random_savanna_cell();
			percolate(sav_cell, fire_ignition_times[i], verbose);
		}
	}
	float get_suppressed_flammability(Cell* cell, float fire_free_interval) {
		if (_suppressed_flammability > -999) return _suppressed_flammability;
	}
	float get_unsuppressed_flammability(Cell* cell, float fire_free_interval) {
		if (_unsuppressed_flammability > -999) return _unsuppressed_flammability;
	}
	float get_flammability(Cell* cell, float fire_free_interval) {
		if (cell->state == 0) return get_unsuppressed_flammability(cell, fire_free_interval);
		else return get_suppressed_flammability(cell, fire_free_interval);
	}
	bool tree_dies(Cell* cell, float fire_free_interval) {
		return 1; // TEMP: trees always die if ignited (TODO: make dependent on fire-free interval and fire resistance)
	}
	void kill_tree(Tree* tree) {
		state.grid.populate_tree_domain(tree, false);
		vector<Tree*> neighbors = state.population.get_neighbors(tree);
		for (auto neighbor : neighbors) state.grid.populate_tree_domain(neighbor, true);
		state.population.remove(tree);
	}
	void induce_mortality(Cell* cell, float fire_free_interval) {
		if (tree_dies(cell, fire_free_interval)) {
			kill_tree(cell->tree);
		}
	}
	inline bool cell_will_ignite(Cell* cell, float t_start) {
		if (t_start - cell->time_last_fire < 10e-4) {
			return false; // Do not ignite cells which have already been burned by the current fire.
		}
		return help::get_rand_float(0.0, 1.0) < get_flammability(cell, t_start - cell->time_last_fire);
	}
	inline void burn_cell(Cell* cell, float t_start) {
		cell->time_last_fire = t_start;
		if (cell->state == 1) induce_mortality(cell, t_start - cell->time_last_fire);
	}
	void percolate(Cell* cell, float t_start, bool verbose=false) {
		std::queue<Cell*> queue;
		queue.push(cell);
		burn_cell(cell, t_start);
		int no_burned_cells = 1;
		if (verbose) printf("Percolating fire...\n");
		while (!queue.empty()) {
			Cell* cell = queue.front();
			queue.pop();

			// Percolate to neighbors
			for (int i = 0; i < 8; i++) {
				Cell* neighbor = state.grid.get_cell_at_position(cell->pos + neighbor_indices[i]);
				if (cell_will_ignite(neighbor, t_start)) {
					queue.push(neighbor);
					burn_cell(neighbor, t_start);
					no_burned_cells++;
				}
			}
		}
		if (verbose) printf("Cells burned: %i \n", no_burned_cells);
	}
	float _unsuppressed_flammability = 0;
	float _suppressed_flammability = 0;
	float self_ignition_factor = 0;
	float rainfall = 0;
	int timestep = 0;
	int time = 0;
	pair<int, int>* neighbor_indices = 0;
	State state;
};

