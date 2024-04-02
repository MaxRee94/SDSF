#pragma once
#include "agents.h"
#include "grid_agent.forward.h"

class Cell {
public:
	Cell() = default;
	int state = 0;
	int tree = 0;
	pair<int, int> pos;
	float time_last_fire = -999.0;
	bool operator==(const Cell& cell) const
	{
		return pos == cell.pos;
	}
};


class Grid {
public:
	Grid() = default;
	Grid(int _size, float _cellsize) {
		size = _size;
		cellsize = _cellsize;
		size_r = (float)size * cellsize;
		no_cells = size * size;
		no_savanna_cells = no_cells;
		init_grid_cells();
		reset_state_distr();
	}
	void init_grid_cells() {
		distribution = new Cell[no_cells];
		for (int i = 0; i < no_cells; i++) {
			pair<int, int> pos = idx_2_pos(i);
			distribution[i].pos = pos;
		}
		state_distribution = new int[no_cells];
	}
	int pos_2_idx(pair<int, int> pos) {
		return size * pos.second + pos.first;
	}
	pair<float, float> get_random_real_position() {
		float x = help::get_rand_float(0, size_r);
		float y = help::get_rand_float(0, size_r);
		pair<float, float> position = pair(x, y);
		return position;
	}
	pair<int, int> get_random_grid_position() {
		int x = help::get_rand_uint(0, size);
		int y = help::get_rand_uint(0, size);
		pair<int, int> position(x, y);
		return position;
	}
	Cell* get_random_savanna_cell() {
		int i = 0;
		int fetch_attempt_limit = 1e6;
		while (i < 1e6) {
			pair<int, int> pos = get_random_grid_position();
			Cell* cell = get_cell_at_position(pos);
			if (cell->state == 0) return cell;
		}
		throw("Runtime error: Could not find savanna cell after %i attempts.\n", fetch_attempt_limit);
	}
	void reset() {
		for (int i = 0; i < no_cells; i++) {
			distribution[i].state = 0;
			distribution[i].tree = 0;
		}
		no_forest_cells = 0;
		no_savanna_cells = no_cells;
	}
	void reset_state_distr() {
		for (int i = 0; i < no_cells; i++) {
			state_distribution[i] = 0;
		}
	}
	void redo_count() {
		no_savanna_cells = 0;
		no_forest_cells = 0;
		for (int i = 0; i < no_cells; i++) {
			no_savanna_cells += !distribution[i].state;
			no_forest_cells += distribution[i].state;
		}
	}
	pair<int, int> idx_2_pos(int idx) {
		int x = idx % size;
		int y = idx / size;
		return pair<int, int>(x, y);
	}
	float get_tree_cover() {
		return (float)no_forest_cells / (float)(no_cells);
	}
	Cell* get_cell_at_position(pair<int, int> pos) {
		cap(pos);
		return &distribution[pos.second * size + pos.first];
	}
	Cell* get_cell_at_position(pair<float, float> _pos) {
		pair<int, int> pos = get_gridbased_position(_pos);
		return get_cell_at_position(pos);
	}
	void get_cells_within_radius(Tree* tree, vector<Cell*>* cells, bool remove_cells = false) {
		pair<float, float> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = tree->radius / cellsize;
		for (int x = tree_center_gb.first - radius_gb; x <= tree_center_gb.first + radius_gb; x++) {
			for (int y = tree_center_gb.second - radius_gb; y <= tree_center_gb.second + radius_gb; y++) {
				if (help::get_dist(pair<float, float>(x, y), tree_center_gb) < radius_gb) {
					Cell* cell = get_cell_at_position(pair<int, int>(x, y));
					if (remove_cells) {
						auto it = find(cells->begin(), cells->end(), cell);
						if (it != cells->end()) cells->erase(it);
					}
					else cells->push_back(cell);
				}
			}
		}
	}
	void populate_tree_domain(Tree* tree, bool add_tree = true, vector<Cell*>* burned_cells = 0, float time_last_fire = -1) {
		pair<float, float> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = tree->radius / cellsize;
		if (!add_tree) radius_gb = round(radius_gb * 1.5);
		int i = 0;
		for (int x = tree_center_gb.first - radius_gb; x <= tree_center_gb.first + radius_gb; x++) {
			for (int y = tree_center_gb.second - radius_gb; y <= tree_center_gb.second + radius_gb; y++) {
				if (help::get_dist(pair<float, float>(x, y), tree_center_gb) < radius_gb) {
					Cell* cell = get_cell_at_position(pair<int, int>(x, y));
					if (add_tree) {
						i++;
						set_to_forest(pair<int, int>(x, y), tree);
						if (burned_cells != nullptr) {
							auto it = find(burned_cells->begin(), burned_cells->end(), cell);
							if (it != burned_cells->end()) burned_cells->erase(it);
						}
					}
					else {
						set_to_savanna(pair<int, int>(x, y), time_last_fire);
						burned_cells->push_back(cell);
					}
				}
			}
		}
		//printf("%s %i cells for tree ptr %i. Stored no cells %i\n", (add_tree ? "added" : "removed"), i, tree, tree->cells.size());
		//printf("Populated %i cells; tree contains %i cells, radius is %f \n", i, tree->cells.size(), tree->radius);
		/*if (tree->cells.size() == 0 && tree->radius > 2.0) {
			printf("WEIRD ----- Tree at %i, %i has no cells but should.\n", tree->position.first, tree->position.second);
		}*/
	}
	void burn_tree(Tree* tree, vector<Tree*> neighbors, float time_last_fire, vector<Cell*>* burned_cells=0) {
		pair<float, float> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = (tree->radius * 1.1) / cellsize;
		for (int x = tree_center_gb.first - radius_gb; x < tree_center_gb.first + radius_gb; x++) {
			for (int y = tree_center_gb.second - radius_gb; y < tree_center_gb.second + radius_gb; y++) {
				pair<int, int> cell_pos_gb = pair<int, int>(x, y);
				if (help::get_dist(cell_pos_gb, tree_center_gb) < radius_gb) {
					pair<float, float> cell_pos_real = get_real_position(cell_pos_gb);
					Tree* overlapping_neighbor = 0;
					for (auto& neighbor : neighbors) {
						if (neighbor->is_within_radius(cell_pos_real)) overlapping_neighbor = neighbor;
					}
					Cell* cell = get_cell_at_position(cell_pos_gb);
					if (overlapping_neighbor != nullptr) {
						cell->tree = overlapping_neighbor->id;
					}
					else {
						set_to_savanna(cell_pos_gb, time_last_fire);
						if (burned_cells != nullptr) burned_cells->push_back(cell);
					}
				}
			}
		}
	}
	int* get_state_distribution(bool collect = true) {
		if (collect) {
			for (int i = 0; i < no_cells; i++) {
				if (distribution[i].state == 1) state_distribution[i] = distribution[i].tree;
			}
		}
		return state_distribution;
	}
	void set_state_distribution(int* distr) {
		for (int i = 0; i < no_cells; i++) {
			if (distr[i] >= 0 && distr[i] < 10)
				state_distribution[i] = distr[i];
			else {
				state_distribution[i] = 0;
			}
		}
	}
	void set_to_forest(int idx, Tree* tree) {
		if (distribution[idx].state == 0) {
			no_savanna_cells--;
			no_forest_cells++;
		}
		distribution[idx].state = 1;
		distribution[idx].tree = tree->id;
		distribution[idx].time_last_fire = -1.0;
	}
	void set_to_savanna(int idx, float _time_last_fire = -1) {
		if (distribution[idx].state == 1) {
			no_savanna_cells++;
			no_forest_cells--;
		}
		//else cout << "was already savanna\n";
		distribution[idx].state = 0;
		distribution[idx].tree = 0;
		if (_time_last_fire > -1) distribution[idx].time_last_fire = _time_last_fire;
	}
	void set_to_forest(pair<int, int> position_grid, Tree* tree) {
		cap(position_grid);
		set_to_forest(position_grid.second * size + position_grid.first, tree);
	}
	void set_to_savanna(pair<int, int> position_grid, float time_last_fire = -1) {
		cap(position_grid);
		set_to_savanna(pos_2_idx(position_grid), time_last_fire);
	}
	void cap(pair<int, int> &position_grid) {
		if (position_grid.first < 0) position_grid.first = size + position_grid.first;
		if (position_grid.second < 0) position_grid.second = size + position_grid.second;
		position_grid.first %= size;
		position_grid.second %= size;
	}
	pair<int, int> get_gridbased_position(pair<float, float> position) {
		return pair<int, int>(position.first / cellsize, position.second / cellsize);
	}
	pair<float, float> get_real_position(pair<int, int> position) {
		return pair<int, int>(position.first * cellsize, position.second * cellsize);
	}
	int size = 0;
	int no_cells = 0;
	float size_r = 0;
	float cellsize = 1.5;
	Cell* distribution = 0;
	int* state_distribution = 0;
	int no_savanna_cells = 0;
	int no_forest_cells = 0;
};
