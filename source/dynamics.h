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
	bool check_faulty_refs(Tree* __tree = 0) {
		// TEMP
		int faulty_references = 0;
		/*for (int i = 0; i < state.grid.no_cells; i++) {
			Cell _cell = state.grid.distribution[i];
			if (_cell.state == 1 && _cell.tree->id == tree->id) {
				state.grid.set_to_savanna(_cell.pos);
			}
		}*/

		//TEMP
		int* distribution = new int[state.grid.no_cells];
		Tree* _tree = 0;
		for (auto& tree : state.population.members) {
			if (__tree != nullptr) tree = *__tree;
			for (auto& _cell : tree.cells) {
				state.grid.set_to_savanna(_cell->pos, -1);
			}
			for (int i = 0; i < state.grid.no_cells; i++) {
				Cell _cell = state.grid.distribution[i];
				if (_cell.state == 1 && _cell.tree->id == tree.id && _cell.tree->radius == tree.radius) {
					_tree = &tree;
					pair<int, int> pos = state.grid.idx_2_pos(i);
					if (faulty_references % 20 == 0) {
						/*printf(
							"Cell %i, %i still references tree to be deleted.\n",
							pos.first, pos.second
						);*/
					}
					faulty_references++;
					distribution[i] = 2;
					/*for (auto& c : _cell.tree->cells) {
						distribution[state.grid.pos_2_idx(c->pos)] += 4;
					}
					if (faulty_references == 1) printf("Referenced tree pos: (%f, %f), base pos: (%f, %f)", _cell.tree->position.first, _cell.tree->position.second, tree->position.first, tree->position.second);*/
				}
			}
			if (_tree != nullptr) break;
			if (__tree != nullptr) break;
		}
		if (faulty_references) printf("--- References to burned tree detected. No tree cells: %i, radius: %f, no faulty references: %i \n", _tree->cells.size(), _tree->radius, faulty_references);
		if (faulty_references) {
			for (auto& _cell : _tree->cells) {
				distribution[state.grid.pos_2_idx(_cell->pos)] += 1;
			}
			state.grid.set_state_distribution(distribution);
		}
		delete[] distribution;
		state.repopulate_grid(0);
		return faulty_references > 0;
	}
	void update() {
		if (check_faulty_refs()) return;
		time++;
		printf("Time: %i\n", time);
		//state.repopulate_grid(verbosity);
		simulate_fires();
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
		//vector<Cell*> burned_cells = {};
		//state.grid.burn_tree(tree, neighbors, cell->time_last_fire, burned_cells);
		//state.grid.populate_tree_domain(tree, false, cell->time_last_fire);
		for (auto& _cell : tree->cells) {
			state.grid.set_to_savanna(_cell->pos, cell->time_last_fire);
		}
		vector<Tree*> neighbors = state.population.get_neighbors(tree);
		for (auto& neighbor : neighbors) {
			state.grid.populate_tree_domain(neighbor, true);
		}

		if (check_faulty_refs(tree)) return;

		/*if (faulty_references) {
			printf("No burned cells: %i, no faulty references: %i \n", burned_cells.size(), faulty_references);
			printf("Burned cells:\n");
			for (auto& _cell : burned_cells) {
				pair<int, int> pos = _cell->pos;
				state.grid.cap(pos);
				cout << "(" + to_string(pos.first) + ", " + to_string(pos.second) + "), ";
			}
			cout << endl;
		}*/
		state.population.remove(tree);
	}
	void induce_mortality(Cell* cell, float fire_free_interval) {
		if (tree_dies(cell, fire_free_interval)) {
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
		//check_faulty_refs();
		bool fault = false;
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
					//fault = check_faulty_refs();
					//cout << "checking..\n";
					if (fault) {
						printf("Fault at burned cell %i \n", no_burned_cells);
					}
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
	int verbosity = 0;
	pair<int, int>* neighbor_indices = 0;
	State state;
};

