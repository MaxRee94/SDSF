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
		verbosity = _verbosity;
		pop = &state.population;
		grid = &state.grid;
	};
	void init_state(int gridsize, float cellsize, float max_radius, float radius_q1, float radius_q2) {
		state = State(gridsize, cellsize, max_radius, radius_q1, radius_q2);
		for (int i = 0; i < 8; i++) neighbor_offsets.push_back(state.neighbor_offsets[i]);
	}
	void update() {
		grid->reset_state_distr();
		if (verbosity == 2) {
			cout << endl;
			for (auto& [id, tree] : state.population.members) {
				printf("*** (update func) - tree (ptr %i) radius: %f\n", &tree, tree.radius);
			}
			vector<Tree*> trees = {};
			for (int i = 0; i < grid->no_cells; i++) {
				Cell cell = grid->distribution[i];
				if (cell.state == 0) continue;
				if (find(trees.begin(), trees.end(), pop->get(cell.tree)) == trees.end()) {
					printf("*** (update func) - Unique tree in grid with id %i and radius %f \n",
						cell.tree, pop->get(cell.tree)->radius);
					trees.push_back(pop->get(cell.tree));
				}
			}
			printf("*** (update func) - no trees: %i \n\n", trees.size());
		}

		time++;
		printf("Time: %i\n", time);
		simulate_fires();
		state.repopulate_grid(0);
		store_tree_colors(false);
		if (verbosity > 0) report_statistics();
	}
	void report_statistics() {
		printf("- Tree cover: %f \n", grid->get_tree_cover());
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
					grid->state_distribution[i] = (float)trees[j]->id * INV_RAND_MAX * (float)pop->size() + 1;
				}
				j++;
			}
			
		}
		if (verbosity > 0) printf("grid counted cells: %i, total intern cells: %i \n", no_grid_counted_cells, no_internally_stored_cells);
	}
	void kill_tree(Cell* cell, queue<Cell*>& queue) {
		Tree* tree = pop->get(cell->tree);
		if (verbosity == 2) printf("Killing tree %i ... \n", tree->id);

		// Obtain cells within radius of killed tree that are not part of neighboring trees.
		vector<Cell*> burned_cells = {};
		grid->get_cells_within_radius(tree, &burned_cells); // Get all cells within killed tree radius
		vector<Cell*> orig_burned_cells = burned_cells;
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
		if (verbosity == 2) printf("new pop size (after removal of ptr %i): %i \n", tree, state.population.size());
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
		return help::get_rand_float(0.0, 1.0) < get_flammability(cell, t_start - cell->time_last_fire);
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
			vector<pair<int, int>> _neighbor_offsets = neighbor_offsets;
			for (int i = 0; i < 8; i++) {
				Cell* neighbor = state.grid.get_cell_at_position(cell->pos + state.neighbor_offsets[i]);
				if (cell_will_ignite(neighbor, t_start)) {
					queue.push(neighbor);
					burn_cell(neighbor, t_start, queue);
					no_burned_cells++;
				}
			}
		}
		return no_burned_cells;
	}
	float _unsuppressed_flammability = 0;
	float _suppressed_flammability = 0;
	float self_ignition_factor = 0;
	float rainfall = 0;
	int timestep = 0;
	int time = 0;
	int pop_size = 0;
	int verbosity = 0;
	State state;
	Population* pop = 0;
	Grid* grid = 0;
	vector<pair<int, int>> neighbor_offsets = {};
};

