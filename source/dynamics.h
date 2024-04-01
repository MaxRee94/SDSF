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
	};
	void init_state(int gridsize, float cellsize, float max_radius, float radius_q1, float radius_q2) {
		state = State(gridsize, cellsize, max_radius, radius_q1, radius_q2);
	}
	bool check_faulty_refs(Tree* tree = 0) {
		int faulty_references = 0;
		int* distribution = new int[state.grid.no_cells];
		pair<float, float> tree_pos = tree->position;

		for (int i = 0; i < state.grid.no_cells; i++) {
			Cell _cell = state.grid.distribution[i];
			if (_cell.state == 0 && _cell.tree == nullptr) continue;
			float dist = state.get_dist(
				state.grid.get_real_position(_cell.pos), tree_pos
			);
			int cell_tree_id = _cell.tree->id;
			int tree_id = tree->id;
			float radius = tree->radius;
			bool cell_outside_tree_radius_references_tree = cell_tree_id == tree_id && (int)dist > radius;
			if (cell_outside_tree_radius_references_tree) {
				if (faulty_references == 0) printf(
					"- Cell's tree ptr: %i, tree's ptr: %i, cell tree id %i, tree id %i \n",
					_cell.tree, tree, _cell.tree->id, tree->id
				);
				faulty_references++;
				state.get_dist(
					state.grid.get_real_position(_cell.pos), tree_pos, true
				);
				distribution[i] = 2;
			}
		}

		if (faulty_references) printf("--- References to burned tree detected. No tree cells: %i, radius: %f, no faulty references: %i, tree ptr: %i \n", tree->cells.size(), tree->radius, faulty_references, tree);
		if (faulty_references) {
			for (auto& _cell : tree->cells) {
				int idx = state.grid.pos_2_idx(_cell->pos);
				if (state.grid.distribution[idx].state == 0) continue;
				int grid_id = state.grid.distribution[idx].tree->id;
				int tree_id = tree->id;
				if (grid_id != tree_id) distribution[idx] = 4;
				else distribution[idx] += 1;
			}
			/*if (state.population.size() == 0) {
				state.grid.set_state_distribution(distribution);
				stop_early = true;
			}*/
			state.grid.set_state_distribution(distribution);
			stop_early = true;
		}
		delete[] distribution;
		return faulty_references > 0;
	}
	void update() {
		state.population.store_positions();
		stop_early = false;
		state.grid.reset_state_distr();
		//state.get_dist(pair<float, float>(4, 5), pair<float, float>(140, 140));

		//if (check_faulty_refs()) return;
		if (verbosity == 2) {
			cout << endl;
			for (auto& tree : state.population.members) {
				printf("*** (update func) - tree (ptr %i) radius: %f, no cells: %i \n", &tree, tree.radius, tree.cells.size());
			}
			vector<Tree*> trees = {};
			for (int i = 0; i < state.grid.no_cells; i++) {
				Cell cell = state.grid.distribution[i];
				if (cell.state == 0) continue;
				if (find(trees.begin(), trees.end(), cell.tree) == trees.end()) {
					printf("*** (update func) - Unique tree in grid with ptr %i, id %i and radius %f \n",
						cell.tree, cell.tree->id, cell.tree->radius);
					trees.push_back(cell.tree);
					//check_faulty_refs(cell.tree);
				}
			}
			printf("*** (update func) - no trees: %i \n\n", trees.size());
		}

		time++;
		printf("\nTime: %i\n", time);
		//state.repopulate_grid(verbosity);
		simulate_fires();
		//state.repopulate_grid(verbosity);
		if (verbosity > 0) report_statistics();
		state.population.check_positions();
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
	void save_2_trees(int method) {
		state.grid.reset_state_distr();
		vector<Tree*> trees = {};
		int no_internally_stored_cells = 0;
		int no_grid_counted_cells = 0;
		for (int i = 0; i < state.grid.no_cells; i++) {
			if (state.grid.distribution[i].state == 0) continue;
			if (method == 1) {
				if (find(state.population.members[0].cells.begin(), state.population.members[0].cells.end(), &state.grid.distribution[i]) !=
					state.population.members[0].cells.end()
					)
					state.grid.state_distribution[i] = 1;
				if (find(state.population.members[1].cells.begin(), state.population.members[1].cells.end(), &state.grid.distribution[i]) !=
					state.population.members[1].cells.end()
					)
					state.grid.state_distribution[i] += 2;
			}
			else if (method == 2) {
				Cell cell = state.grid.distribution[i];
				if (cell.tree != nullptr) no_grid_counted_cells++;
				bool new_tree = find(trees.begin(), trees.end(), cell.tree) == trees.end();
				if (new_tree) {
					trees.push_back(cell.tree);
					no_internally_stored_cells += cell.tree->cells.size();
				}
				int j = 0;
				while (j < trees.size()) {
					if (cell.tree == trees[j]) state.grid.state_distribution[i] = j + 1;
					if (find(state.population.members.begin(), state.population.members.end(), *cell.tree) == state.population.members.end()) state.grid.state_distribution[i] = 999;
					j++;
				}
			}
		}
		printf("grid counted cells: %i, total intern cells: %i \n", no_grid_counted_cells, no_internally_stored_cells);
	}
	void kill_tree(Cell* cell, queue<Cell*>& queue) {
		Tree* tree = cell->tree;

		if (verbosity == 2) printf("before emptying tree with ptr %i\n", tree);
		//if (check_faulty_refs(tree)) return;
		pop_size = state.population.size();

		// Obtain cells within radius of killed tree that are not part of neighboring trees.
		vector<Cell*> burned_cells = {};
		state.grid.get_cells_within_radius(tree, &burned_cells); // Get all cells within killed tree radius
		vector<Cell*> orig_burned_cells = burned_cells;
		vector<Tree*> neighbors = state.get_tree_neighbors(tree);
		for (auto& neighbor : neighbors) {
			state.grid.get_cells_within_radius(neighbor, &burned_cells, true); // Exclude cells that fall within radius of neighbors
		}

		// Burn cells and add to queue
		for (Cell* _cell : burned_cells) {
			queue.push(_cell);
			state.grid.state_distribution[state.grid.pos_2_idx(_cell->pos)] = 5;
			state.grid.set_to_savanna(_cell->pos, cell->time_last_fire);
		}

		// Fill cells which may not have belonged to neighbors previously (because of dominance by killed tree), but do fall within their radii.
		for (auto& neighbor : neighbors) {
			if (verbosity == 2) printf("Re-filling tree with ptr %i and no cells %i\n", neighbor, neighbor->cells.size());
			state.grid.populate_tree_domain(neighbor, true);
		}

		printf("No burned cells: %i \n", burned_cells.size());
		if (burned_cells.size() == 0) printf("No cells in ptr %i : %i\n", tree, tree->cells.size());

		// Check whether all grid cell references to killed tree have been removed
		for (int i = 0; i < state.grid.no_cells; i++) {
			Cell _cell = state.grid.distribution[i];
			if (_cell.state == 0) continue;
			if (_cell.tree == tree) {
				if (!stop_early) {
					state.grid.reset_state_distr();
					for (auto& __cell : orig_burned_cells) {
						int idx = state.grid.pos_2_idx(__cell->pos);
						state.grid.state_distribution[idx] = 1;
					}
					for (auto& __cell : burned_cells) {
						int idx = state.grid.pos_2_idx(__cell->pos);
						state.grid.state_distribution[idx] = 6;
					}
				}
				stop_early = true;
				state.grid.state_distribution[i] = 999;
				bool is_orig_burned_cell = find(orig_burned_cells.begin(), orig_burned_cells.end(), &_cell) != orig_burned_cells.end();
				bool is_burned_cell = find(burned_cells.begin(), burned_cells.end(), &_cell) != burned_cells.end();
				printf("Cell (%i, %i) is %s a burned cell, and %s an original burned cell\n", _cell.pos.first, _cell.pos.second, (is_burned_cell ? "" : "not"), (is_orig_burned_cell ? "" : "not"));
				printf("Tree ptr %i has position %f, %f \n", tree, tree->position.first, tree->position.second);
			}
		}

		state.population.remove(tree);
		//state.repopulate_grid(0);
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
		//cout << "--- Checking before percolating.\n";
		//check_faulty_refs();
		//cout << "--- Checked before percolating.\n";
		bool fault = false;
		pop_size = state.population.size();
		while (!queue.empty()) {
			Cell* cell = queue.front();
			queue.pop();

			// Percolate to neighbors
			for (int i = 0; i < 8; i++) {
				Cell* neighbor = state.grid.get_cell_at_position(cell->pos + state.neighbor_offsets[i]);
				if (cell_will_ignite(neighbor, t_start)) {

					// Check trees on grid
					vector<Tree*> trees = {};
					for (int i = 0; i < state.grid.no_cells; i++) {
						Cell _cell = state.grid.distribution[i];
						if (_cell.state == 0) continue;
						if (find(trees.begin(), trees.end(), _cell.tree) == trees.end() &&
							_cell.tree->cells.size() == 0) {
							//printf
							//	("*** (percolate func) - Unique tree in grid with ptr %i, id %i, radius %f and no cells %i \n",
							//	_cell.tree, _cell.tree->id, _cell.tree->radius, _cell.tree->cells.size()
							//);
							//printf("double-check: no cells of ptr %i : %i, id: %i \n", &state.population.members[0], state.population.members[0].cells.size(), state.population.members[0].id);
							trees.push_back(_cell.tree);
						}
					}
					//if (trees.size() > 0) cout << endl;

					queue.push(neighbor);
					burn_cell(neighbor, t_start, queue);
					if (stop_early) return no_burned_cells;
					/*if (state.population.size() < pop_size) {
						pop_size = state.population.size();
						fault = check_faulty_refs(&state.population.removed_tree);
					}*/
					no_burned_cells++;
					state.grid.state_distribution[state.grid.pos_2_idx(neighbor->pos)] = 5;
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
	bool stop_early = false;
	int time = 0;
	int pop_size = 0;
	int verbosity = 0;
	State state;
};

