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
	void init_state(int gridsize, float cellsize, float mean_radius) {
		state = State(gridsize, cellsize, mean_radius);
	}
	void update() {
		time++;
		printf("Updating... Time: %i\n", time);
		state.repopulate_grid();
		simulate_fires();
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
	void simulate_fires() {
		vector<float> fire_ignition_times = get_ordered_fire_ignition_times();
		for (int i = 0; i < fire_ignition_times.size(); i++) {
			Cell* sav_cell = state.grid.get_random_savanna_cell();
			percolate(sav_cell, fire_ignition_times[i]);
		}
	}
	float get_tree_flammability(Cell* cell, float fire_free_interval) {
		if (_unsuppressed_flammability > -999) return _unsuppressed_flammability;
	}
	float get_unsuppressed_flammability(Cell* cell, float fire_free_interval) {
		if (_unsuppressed_flammability > -999) {
			if (cell->state == 0) return _unsuppressed_flammability;
			else return get_tree_flammability(cell, fire_free_interval);
		}
	}
	bool tree_dies(Cell* cell, float fire_free_interval) {
		return 1; // TEMP: all trees die when ignited (TODO: make dependent on fire-free interval and fire resistance)
	}
	void induce_mortality(Cell* cell, float fire_free_interval) {
		if (tree_dies(cell, fire_free_interval)) {
			Tree* tree = cell->tree;
			state.population.remove(tree);
			state.grid.populate_tree_domain(tree, false);
		}
	}
	inline bool cell_will_ignite(Cell* cell, float t_start) {
		if (t_start - cell->time_last_fire < 1e-4) {
			return false; // Do not ignite cells which have already been burned by the current fire.
		}
		return help::get_rand_float(0.0, 1.0) < get_unsuppressed_flammability(cell, t_start - cell->time_last_fire);
	}
	inline void burn_cell(Cell* cell, float t_start) {
		cell->time_last_fire = t_start;
		if (cell->state == 1) induce_mortality(cell, t_start - cell->time_last_fire);
	}
	void percolate(Cell* cell, float t_start) {
		std::queue<Cell*> queue;
		queue.push(cell);
		burn_cell(cell, t_start);
		int i = 0;
		printf("Percolating fire...\n");
		while (!queue.empty()) {
			i++;
			Cell* cell = queue.front();
			queue.pop();
			if (i % 10 == 0) {
				printf("Iteration %i, queue length: %i \n", i, queue.size());
			}

			// Percolate to neighbors
			for (int i = 0; i < 8; i++) {
				Cell* neighbor = state.grid.get_cell_at_position(cell->pos + neighbor_indices[i]);
				if (cell_will_ignite(neighbor, t_start)) {
					queue.push(neighbor);
					burn_cell(neighbor, t_start);
				}
			}
		}
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

