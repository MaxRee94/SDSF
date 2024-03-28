#pragma once
#include "state.h"


class Dynamics {
public:
	Dynamics() = default;
	Dynamics(int _timestep, float __unsuppressed_flammability, float __suppressed_flammability, float _self_ignition_factor, float _rainfall, int _verbosity) :
		timestep(_timestep), _unsuppressed_flammability(__unsuppressed_flammability), _suppressed_flammability(__suppressed_flammability),
		self_ignition_factor(_self_ignition_factor), rainfall(_rainfall)
	{
		time = 0;
		help::init_RNG();
		init_neighbor_indices();
		verbosity = _verbosity;
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
	void init_state(int gridsize, float cellsize, float max_radius, float radius_q1, float radius_q2) {
		state = State(gridsize, cellsize, max_radius, radius_q1, radius_q2);
	}
	void update() {
		// TEMP
		printf("Checking all cells...\n");
		for (int i = 0; i < state.grid.no_cells; i++) {
			Cell cell = state.grid.distribution[i];
			printf(
				" -- Checking cell %i, %i...\n", cell.pos.first, cell.pos.second
			);
			if (cell.state == 1) {
				if (cell.tree->radius > 10e6) {
					printf("Dummy printout\n");
				}
			}
			if (i > 150) break;
		}

		time++;
		printf("Time: %i\n", time);
		//printf("Pop size: %i \n", state.population.size());
		simulate_fires();
		//state.population.grow();
		state.repopulate_grid(verbosity);
		if (verbosity > 0) report_statistics();
	}
	void report_statistics() {
		printf("- Tree cover: %f \n", state.grid.get_tree_cover());
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
		int no_burned_cells = 0;
		for (int i = 0; i < fire_ignition_times.size(); i++) {
			Cell* sav_cell = state.grid.get_random_savanna_cell();
			no_burned_cells += percolate(sav_cell, fire_ignition_times[i]);
		}
		if (verbosity > 0) {
			printf("Cells burned: %i \n", no_burned_cells);
			printf("Number of fires: %i \n", fire_ignition_times.size());
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
	void kill_tree(Cell* cell) {
		Tree* tree = cell->tree;
		if (!state.population.is_population_member(tree)) {
			printf("Tree is not a population member.\n");
			state.grid.populate_tree_domain(tree, false);
			return;
		}
		vector<Tree*> neighbors = state.population.get_neighbors(tree);
		state.grid.burn_tree(tree, neighbors, cell->time_last_fire);

		for (int i = 0; i < state.grid.no_cells; i++) {
			Cell cell = state.grid.distribution[i];
			if (cell.tree == tree) {
				pair<int, int> pos = state.grid.idx_2_pos(i);
				printf(
					"Cell %i, %i, still references tree to be deleted.\n",
					pos.first, pos.second
				);
			}
		}
		//state.population.remove(tree);
	}
	void induce_mortality(Cell* cell, float fire_free_interval) {
		if (tree_dies(cell, fire_free_interval)) {
			// TEMP
			printf(
				"Checking cell (in induce mortality, before killing tree) %i, %i...\n", cell->pos.first, cell->pos.second
			);
			if (cell->state == 1 && cell->tree->radius > 10e6) {
				cout << "dummy printout\n";
			}
			kill_tree(cell);
		}
	}
	inline bool cell_will_ignite(Cell* cell, float t_start) {
		if (t_start - cell->time_last_fire < 10e-4) {
			return false; // Do not ignite cells which have already been burned by the current fire.
		}
		return help::get_rand_float(0.0, 1.0) < get_flammability(cell, t_start - cell->time_last_fire);
	}
	inline void burn_cell(Cell* cell, float t_start) {
		// TEMP
		printf(
			"Checking cell (in burn cell, before inducing mortality) %i, %i...\n", cell->pos.first, cell->pos.second
		);
		if (cell->state == 1 && cell->tree->radius > 10e6) {
			cout << "dummy printout\n";
		}
		cell->time_last_fire = t_start;
		if (cell->state == 1) {
			induce_mortality(cell, t_start - cell->time_last_fire);
		}
	}
	int percolate(Cell* cell, float t_start) {
		std::queue<Cell*> queue;
		queue.push(cell);
		burn_cell(cell, t_start);
		int no_burned_cells = 1;
		if (verbosity == 2) printf("Percolating fire...\n");

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
			cout << "no burned cells: " << no_burned_cells << endl;
		}
		return no_burned_cells;
	}
	float _unsuppressed_flammability = 0;
	float _suppressed_flammability = 0;
	float self_ignition_factor = 0;
	float rainfall = 0;
	int timestep = 0;
	int time = 0;
	int verbosity = 0;
	pair<int, int>* neighbor_indices = 0;
	State state;
};

