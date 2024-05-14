#pragma once
#include "agents.h"
#include "grid_agent.forward.h"


class Cell {
public:
	Cell() = default;
	int state = 0;
	int idx = 0;
	float time_last_fire = 0;
	vector<int> trees;
	float LAI = 0; // Cumulative Leaf Area Index (LAI) of all trees in the cell.
	float grass_LAI = 0; // Cumulative Leaf Area Index (LAI) of all trees in the cell.
	pair<int, int> pos;
	bool seedling_present = false;
	pair<float, int> largest_stem;	// < float: dbh of largest tree or seedling that has its stem in this cell,
									//	 int:   Largest tree id, or id of parent tree if the largest stem belongs to a seedling >
	bool cell_is_occupied_by_larger_stem(pair<float, int> &tree_proxy) {
		if (tree_proxy.first < largest_stem.first) {
			return true;
		}
		return false;
	}
	bool seedling_is_shaded_out() {
		float shade = get_shading();
		if (help::get_rand_float(0, 1) < shade) {
			return true;
		}
		return false;
	}
	void set_largest_stem(float dbh, int id) {
		largest_stem = pair<float, int>(dbh, id);
	}
	void reset_largest_stem() {
		largest_stem = pair<float, int>(0, 0);
	}
	void update_grass_LAI() {
		grass_LAI = get_grass_LAI();
	}
	bool is_hospitable(pair<float, int> tree_proxy) {
		// If tree_proxy is larger than the current largest stem in the cell, it is assumed to be able to outcompete the other tree.
		// If not, tree_proxy is assumed to be shaded out or otherwise outcompeted by the existing larger tree.
		bool hospitable = !cell_is_occupied_by_larger_stem(tree_proxy);
		hospitable = hospitable && !seedling_is_shaded_out();
		return hospitable;
	}
	void add_tree_if_not_exists(Tree* tree) {
		if (find(trees.begin(), trees.end(), tree->id) == trees.end()) {
			add_tree(tree);
		}
	}
	void add_tree(Tree* tree) {
		trees.push_back(tree->id);
	}
	void remove_tree(int id) {
		//printf("\nRemoving tree %i from cell %i.\n", id, idx);
		//help::print_vector(&trees);
		help::remove_from_vec(&trees, id);
		//help::print_vector(&trees);
	}
	float get_grass_LAI(float _LAI) {
		return 0.241f * (_LAI * _LAI) - 1.709f * _LAI + 2.899f; // Relationship between grass- and tree LAI from Hoffman et al. (2012), figure 2b.
	}
	float get_grass_LAI() {
		return get_grass_LAI(LAI);
	}
	float get_fuel_load(float _grass_LAI) {
		return 0.344946533f * _grass_LAI; // Normalize to range [0, 1] by multiplying with inverse of max value (2.899).
	}
	float get_fuel_load() {
		return get_fuel_load(grass_LAI);
	}
	float get_LAI_above_tree(Tree* tree, Population* population = 0) {
		if (population == nullptr) return LAI;
		float LAI_above_tree = 0;
		float given_tree_height = tree->height;
		for (int tree_id : trees) {
			if (tree_id == tree->id) continue;
			Tree* _tree = population->get(tree_id);
			float height = _tree->height;
			if (height < given_tree_height) LAI_above_tree += _tree->LAI;
		}
		return LAI_above_tree;
	}
	float get_tree_shading(Tree* tree, Population* population = 0) {
		float LAI_above_tree = get_LAI_above_tree(tree, population);
		//float PAR = -0.056299 * LAI_above_tree + 0.58102; // Percentage Photosynthetically Active Radiation (PAR) in W/m^2 that reaches the given tree.
															 // From Mukherjee (2016), figure 4.
		//PAR = max(0.0f, min(PAR, 1.0f)); // Cap PAR to range [0, 1].
		return LAI_above_tree / 5.0f; // Max LAI for forests is 5.0, so we normalize the shading to range [0, 1] by dividing by 5.
	}
	float get_shading() {
		Tree tree; tree.LAI = LAI;
		return get_tree_shading(&tree);
	}
	bool operator==(const Cell& cell) const
	{
		return pos == cell.pos;
	}
};


class Grid {
public:
	Grid() = default;
	Grid(int _width, float _cell_width) {
		printf("\nInitializing grid with width %i and cell width %f.\n", _width, _cell_width);
		width = _width;
		cell_width = _cell_width;
		width_r = (float)width * cell_width;
		no_cells = width * width;
		no_savanna_cells = no_cells;
		init_grid_cells();
		reset_state_distr();
		area = no_cells * cell_width * cell_width;
		cell_area = cell_width * cell_width;
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
	void free() {
		delete[] distribution;
		delete[] state_distribution;
	}
	int pos_2_idx(pair<float, float> pos) {
		return width * (pos.second / cell_width) + (pos.first / cell_width);
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
		pair<float, float> real_position((float)gridbased_location.first * cell_width, (float)gridbased_location.second * cell_width);
		float x = help::get_rand_float(real_position.first, real_position.first + cell_width);
		float y = help::get_rand_float(real_position.second, real_position.second + cell_width);
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
		return pair<float, float>(cell->pos.first * cell_width, cell->pos.second * cell_width);
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
	void update_grass_LAIs() {
		for (int i = 0; i < no_cells; i++) {
			distribution[i].update_grass_LAI();
		}
	}
	virtual void reset() {
		for (int i = 0; i < no_cells; i++) {
			distribution[i].state = 0;
			distribution[i].trees.clear();
			distribution[i].LAI = 0;
			distribution[i].grass_LAI = 0;
			distribution[i].largest_stem = pair<float, int>(0, -1);
			distribution[i].seedling_present = false;
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
		int radius_gb = tree->radius / cell_width;
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
		pair<int, int> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = tree->radius / cell_width;
		for (float x = tree_center_gb.first - radius_gb; x <= tree_center_gb.first + radius_gb; x+=1) {
			for (float y = tree_center_gb.second - radius_gb; y <= tree_center_gb.second + radius_gb; y+=1) {
				pair<float, float> position(x * cell_width, y * cell_width);
				if (tree->is_within_radius(position)) {
					set_to_forest(pair<float, float>(x, y), tree);
				}
			}
		}
		cap(tree_center_gb);
		distribution[pos_2_idx(tree_center_gb)].set_largest_stem(tree->dbh, tree->id);
		distribution[pos_2_idx(tree_center_gb)].add_tree_if_not_exists(tree);
	}
	float compute_shade_on_individual_tree(Tree* tree) {
		float shade = 0;
		float no_cells = 0;
		pair<int, int> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = tree->radius / cell_width;
		for (float x = tree_center_gb.first - radius_gb; x <= tree_center_gb.first + radius_gb; x += 1) {
			for (float y = tree_center_gb.second - radius_gb; y <= tree_center_gb.second + radius_gb; y += 1) {
				pair<float, float> position(x * cell_width, y * cell_width);
				if (tree->is_within_radius(position)) {
					Cell* cell = get_cell_at_position(pair<int, int>(x, y));
					shade += cell->get_tree_shading(tree);
					no_cells += 1;
				}
			}
		}
		return shade / no_cells;
	}
	void burn_tree_domain(Tree* tree, queue<Cell*> &queue, float time_last_fire = -1, bool store_tree_death_in_color_distribution = true) {
		pair<int, int> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = round((tree->radius * 1.5) / cell_width);
		for (float x = tree_center_gb.first - radius_gb; x <= tree_center_gb.first + radius_gb; x += 1) {
			for (float y = tree_center_gb.second - radius_gb; y <= tree_center_gb.second + radius_gb; y += 1) {
				pair<float, float> position(x * cell_width, y * cell_width);
				if (tree->is_within_radius(position)) {
					Cell* cell = get_cell_at_position(position);
					if (distribution[cell->idx].state == 0) continue;
					
					// Remove tree id from cell->trees.
					cell->remove_tree(tree->id);

					// Update LAI
					cell->LAI -= tree->LAI * cell_area;

					// Set cell to savanna if the cumulative leaf area is less than half of the area of the cell (LAI * cell_area < 0.5 * cell_area, i.e., LAI < 0.5)).
					if (tree->LAI < 0.5) {
						queue.push(cell);
						set_to_savanna(cell->idx, time_last_fire);
						if (store_tree_death_in_color_distribution) state_distribution[cell->idx] = -6;
						continue;
					}
					state_distribution[cell->idx] = -5;
				}
			}
		}
		cap(tree_center_gb);
		distribution[pos_2_idx(tree_center_gb)].reset_largest_stem();
	}
	void kill_tree_domain(Tree* tree, bool store_tree_death_in_color_distribution = true) {
		queue<Cell*> dummy;
		burn_tree_domain(tree, dummy, -1, store_tree_death_in_color_distribution);
	}
	void update_grass_LAIs_for_individual_tree(Tree* tree) {
		pair<int, int> tree_center_gb = get_gridbased_position(tree->position);
		int radius_gb = tree->radius / cell_width;
		for (float x = tree_center_gb.first - radius_gb; x <= tree_center_gb.first + radius_gb; x += 1) {
			for (float y = tree_center_gb.second - radius_gb; y <= tree_center_gb.second + radius_gb; y += 1) {
				distribution[pos_2_idx(pair<int, int>(x, y))].update_grass_LAI();
			}
		}
	}
	int* get_state_distribution(bool collect = true) {
		if (collect) {
			for (int i = 0; i < no_cells; i++) {
				if (distribution[i].state == 1) state_distribution[i] = max(99.0f - (distribution[i].LAI * 19.0f), 1);
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
		distribution[idx].add_tree(tree);
		distribution[idx].LAI += tree->LAI * cell_area;
	}
	void set_to_savanna(int idx, float _time_last_fire = -1) {
		no_savanna_cells += (distribution[idx].state == 1);
		no_forest_cells -= (distribution[idx].state == 1);

		distribution[idx].state = 0;
		distribution[idx].LAI = 0;
		if (_time_last_fire != -1) distribution[idx].time_last_fire = _time_last_fire;
	}
	float get_LAI_within_bb(pair<int, int> bb_min, pair<int, int> bb_max, float bb_area) {
		float cumulative_LAI = 0;
		for (int x = bb_min.first; x < bb_max.first; x++) {
			for (int y = bb_min.second; y < bb_max.second; y++) {
				Cell* cell = get_cell_at_position(pair<int, int>(x, y));
				cumulative_LAI += cell->LAI;
			}
		}
		float crown_area = get_tree_cover_within_bb(bb_min, bb_max) * bb_area;
		return cumulative_LAI / crown_area;
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
		return pair<int, int>(position.first / cell_width, position.second / cell_width);
	}
	pair<float, float> get_real_position(pair<int, int> position) {
		return pair<float, float>((float)position.first * cell_width, (float)position.second * cell_width);
	}
	pair<float, float> get_real_position(int idx) {
		pair<int, int> position = idx_2_pos(idx);
		return pair<float, float>((float)position.first * cell_width, (float)position.second * cell_width);
	}
	int width = 0;
	int no_cells = 0;
	float width_r = 0;
	float tree_cover = 0;
	float cell_width = 0;
	Cell* distribution = 0;
	int* state_distribution = 0;
	int no_savanna_cells = 0;
	int no_forest_cells = 0;
	float area = 0;
	float cell_area = 0;
};
