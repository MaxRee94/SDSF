#pragma once
#include "agents.h"
#include "grid_agent.forward.h"


class Cell {
public:
	Cell() = default;
	int state = 0;
	int idx = 0;
	float time_last_fire = 0;
	map<int, int> trees;
	float LAI = 0; // Tree leaf area index
	pair<int, int> pos;
	float get_grass_LAI() {
		return 0.241f * (LAI * LAI) - 1.709f * LAI + 2.899f; // Relationship between grass and tree LAI from Hoffman et al. (2012), figure 2b.
	}
	float get_fuel_load() {
		float grass_LAI = get_grass_LAI();
		return 0.344946533f * grass_LAI; // Normalize to range [0, 1] by multiplying with inverse of max value (2.899).
	}
	bool operator==(const Cell& cell) const
	{
		return pos == cell.pos;
	}
};


class Grid {
public:
	Grid() = default;
	Grid(int _width, float _cellsize) {
		width = _width;
		cellsize = _cellsize;
		width_r = (float)width * cellsize;
		no_cells = width * width;
		no_savanna_cells = no_cells;
		init_grid_cells();
		reset_state_distr();
		area = no_cells * cellsize * cellsize;
	}
	void init_grid_cells() {
		distribution = new Cell[no_cells];
		for (int i = 0; i < no_cells; i++) {
			pair<int, int> pos = idx_2_pos(i);
			distribution[i].pos = pos;
			distribution[i].idx = i;
		}
		state_distribution = new int[no_cells];
	}
	int pos_2_idx(pair<float, float> pos) {
		return width * (pos.second / cellsize) + (pos.first / cellsize);
	}
	int pos_2_idx(pair<int, int> pos) {
		return width * pos.second + pos.first;
	}
	pair<float, float> get_random_real_position() {
		float x = help::get_rand_float(0, width_r);
		float y = help::get_rand_float(0, width_r);
		pair<float, float> position = pair(x, y);
		return position;
	}
	pair<int, int> get_random_grid_position() {
		int x = help::get_rand_uint(0, width);
		int y = help::get_rand_uint(0, width);
		pair<int, int> position(x, y);
		return position;
	}
	Cell* get_random_forest_cell() {
		int i = 0;
		int fetch_attempt_limit = 1e6;
		while (i < 1e6) {
			pair<int, int> pos = get_random_grid_position();
			Cell* cell = get_cell_at_position(pos);
			if (cell->state == 1) return cell;
		}
		throw("Runtime error: Could not find forest cell after %i attempts.\n", fetch_attempt_limit);
	}
	pair<float, float> get_random_location_within_cell(pair<int, int> &gridbased_location) {
		pair<float, float> real_position((float)gridbased_location.first * cellsize, (float)gridbased_location.second * cellsize);
		float x = help::get_rand_float(real_position.first, real_position.first + cellsize);
		float y = help::get_rand_float(real_position.second, real_position.second + cellsize);
		return pair<float, float>(x, y);
	}
	pair<float, float> get_random_location_within_cell(int idx) {
		pair<int, int> pos = idx_2_pos(idx);
		return get_random_location_within_cell(pos);
	}
	Cell* get_random_cell() {
		pair<int, int> pos = get_random_grid_position();
		Cell* cell = get_cell_at_position(pos);
		return cell;
	}
	pair<float, float> get_real_cell_position(Cell* cell) {
		return pair<float, float>(cell->pos.first * cellsize, cell->pos.second * cellsize);
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
	virtual void reset() {
		for (int i = 0; i < no_cells; i++) {
			distribution[i].state = 0;
			distribution[i].trees.clear();
			distribution[i].LAI = 0;
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
		int x = idx % width;
		int y = idx / width;
		return pair<int, int>(x, y);
	}
	float get_tree_cover() {
		tree_cover = (float)no_forest_cells / (float)(no_cells);
		return tree_cover;
	}
	Cell* get_cell_at_position(pair<int, int> pos) {
		cap(pos);
		return &distribution[pos.second * width + pos.first];
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
	void populate_tree_domain(Tree* tree) {
		pair<float, float> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = tree->radius / cellsize;
		for (float x = tree_center_gb.first - radius_gb; x <= tree_center_gb.first + radius_gb; x+=1) {
			for (float y = tree_center_gb.second - radius_gb; y <= tree_center_gb.second + radius_gb; y+=1) {
				pair<float, float> position(x * cellsize, y * cellsize);
				if (tree->is_within_radius(position)) {
					set_to_forest(pair<float, float>(x, y), tree);
				}
			}
		}
	}
	void burn_tree_domain(Tree* tree, queue<Cell*> &queue, float time_last_fire = -1, bool store_tree_death_in_color_distribution = true) {
		pair<float, float> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = round((tree->radius * 1.5) / cellsize);
		for (float x = tree_center_gb.first - radius_gb; x <= tree_center_gb.first + radius_gb; x += 1) {
			for (float y = tree_center_gb.second - radius_gb; y <= tree_center_gb.second + radius_gb; y += 1) {
				pair<float, float> position(x * cellsize, y * cellsize);
				if (tree->is_within_radius(position)) {
					Cell* cell = get_cell_at_position(position);
					if (distribution[cell->idx].state == 0) continue;
					
					// Remove tree id from cell->trees.
					int map_size = cell->trees.size();
					auto it = cell->trees.find(tree->id);
					if (it != cell->trees.end()) cell->trees.erase(it);

					// Update LAI
					cell->LAI -= tree->LAI;

					// Set cell to savanna if it is no longer occupied by trees larger than the cell itself.
					if (cell->LAI < (cellsize * 0.2)) {
						queue.push(cell);
						set_to_savanna(cell->idx, time_last_fire);
						if (store_tree_death_in_color_distribution) state_distribution[cell->idx] = -6;
						continue;
					}
					state_distribution[cell->idx] = -5;
				}
			}
		}
	}
	void kill_tree_domain(Tree* tree, bool store_tree_death_in_color_distribution = true) {
		queue<Cell*> dummy;
		burn_tree_domain(tree, dummy, -1, store_tree_death_in_color_distribution);
	}
	int* get_state_distribution(bool collect = true) {
		if (collect) {
			for (int i = 0; i < no_cells; i++) {
				//if (distribution[i].state == 1) state_distribution[i] = distribution[i].trees.begin()->first % 100 + 1;
				if (distribution[i].state == 1) state_distribution[i] = min(distribution[i].LAI, 9) * 10 + 1;
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
		distribution[idx].trees[tree->id] = tree->id;
		if (tree->id == 0) cout << "pushing back tree id 0\n";
		distribution[idx].LAI += tree->LAI;
		distribution[idx].time_last_fire = 0;
	}
	void set_to_savanna(int idx, float _time_last_fire = -1) {
		no_savanna_cells += (distribution[idx].state == 1);
		no_forest_cells -= (distribution[idx].state == 1);

		distribution[idx].state = 0;
		//distribution[idx].LAI = 0;
		if (_time_last_fire != -1) distribution[idx].time_last_fire = _time_last_fire;
	}
	float get_tree_cover_within_bb(pair<int, int> bb_min, pair<int, int> bb_max) {
		int no_forest_cells = 0;
		int no_cells = 0;
		for (int x = bb_min.first; x < bb_max.first; x++) {
			for (int y = bb_min.second; y < bb_max.second; y++) {
				Cell* cell = get_cell_at_position(pair<int, int>(x, y));
				no_forest_cells += cell->state;
				no_cells++;
			}
		}
		return (float)no_forest_cells / (float)no_cells;
	}
	void set_to_forest(pair<int, int> position_grid, Tree* tree) {
		cap(position_grid);
		set_to_forest(position_grid.second * width + position_grid.first, tree);
	}
	void set_to_savanna(pair<int, int> position_grid, float time_last_fire = -1) {
		cap(position_grid);
		set_to_savanna(pos_2_idx(position_grid), time_last_fire);
	}
	void cap(pair<int, int> &position_grid) {
		if (position_grid.first < 0) position_grid.first = width + (position_grid.first % width);
		if (position_grid.second < 0) position_grid.second = width + (position_grid.second % width);
		position_grid.first %= width;
		position_grid.second %= width;
	}
	pair<int, int> get_gridbased_position(pair<float, float> position) {
		return pair<int, int>(position.first / cellsize, position.second / cellsize);
	}
	pair<float, float> get_real_position(pair<int, int> position) {
		return pair<float, float>((float)position.first * cellsize, (float)position.second * cellsize);
	}
	pair<float, float> get_real_position(int idx) {
		pair<int, int> position = idx_2_pos(idx);
		return pair<float, float>((float)position.first * cellsize, (float)position.second * cellsize);
	}
	int width = 0;
	int no_cells = 0;
	float width_r = 0;
	float tree_cover = 0;
	float cellsize = 0;
	Cell* distribution = 0;
	int* state_distribution = 0;
	int no_savanna_cells = 0;
	int no_forest_cells = 0;
	float area = 0;
};
