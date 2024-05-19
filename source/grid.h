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
	float grass_LAI = 0;
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
		float shade = LAI / 5.0f; // Normalize to range [0, 1] by dividing by 5.0 (max LAI is 5.0).
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
	bool tree_is_present(Tree* tree) {
		return find(trees.begin(), trees.end(), tree->id) != trees.end();
	}
	bool tree_is_present(int id) {
		Tree dummy; dummy.id = id;
		return tree_is_present(&dummy);
	}
	bool is_smaller_than_cell(Tree* tree, float cell_halfdiagonal_sqrt) {
		return tree->radius < cell_halfdiagonal_sqrt;
	}
	void add_tree_if_not_present(Tree* tree, float cell_area, bool smaller_than_cell = true) {
		if (!tree_is_present(tree)) {
			add_tree(tree, smaller_than_cell, cell_area);
		}
	}
	bool remove_tree_if_smaller_than_cell(Tree* tree, float cell_area, float cell_halfdiagonal_sqrt, bool smaller_than_cell = true) {
		if (is_smaller_than_cell(tree, cell_halfdiagonal_sqrt)) {
			remove_tree(tree);
		}
		return true;
	}
	void add_tree(Tree* tree, bool sapling = false, float cell_area = 0, float cell_halfdiagonal_sqrt = 0) {
		trees.push_back(tree->id);
		if (sapling || is_smaller_than_cell(tree, cell_halfdiagonal_sqrt))
			add_LAI_of_tree_smaller_than_cell(tree, cell_area);
		else LAI += tree->LAI;
	}
	void remove_tree(Tree* tree, bool sapling = false, float cell_area = 0) {
		help::remove_from_vec(&trees, tree->id);
		if (sapling) remove_LAI_of_tree_smaller_than_cell(tree, cell_area);
		else LAI -= tree->LAI;
	}
	float get_grass_LAI(float _LAI) {
		_LAI -= int(_LAI / 3.0f) * 3.0f;	// Effectively _LAI modulo 3; Done to avoid re-intersecting the y=0 line at about LAI=4.3 
											// (grass LAI would then (incorrectly) start rising again).
		float _grass_LAI = max(0, 0.241f * (_LAI * _LAI) - 1.709f * _LAI + 2.899f); // Relationship between grass- and tree LAI from 
																					// Hoffman et al. (2012), figure 2b.									
		return _grass_LAI;
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
	float get_LAI_of_crown_intersection_and_above(Tree* tree, Population* population = 0) {
		if (population == nullptr) {
			return LAI;
		}
		float shade = 0;
		float crown_reach = tree->height - tree->lowest_branch;
		for (int tree_id : trees) {
			if (tree_id == tree->id) {
				shade += LAI; // Add the tree's own LAI to the shade (since it will cast shade on itself).
				continue;
			}
			Tree* neighbor = population->get(tree_id);
			if (neighbor->height >= tree->height) shade += neighbor->LAI;
			else if (neighbor->height > tree->lowest_branch) { // Implies crown intersection; a portion of the neighbor's crown will cast shade on our tree.
				float neighbor_crown_reach = neighbor->height - neighbor->lowest_branch;
				float crown_intersection = neighbor->height - tree->lowest_branch;
				float LAI_intersection = neighbor->LAI * (crown_intersection / neighbor_crown_reach); // Fraction of the smaller tree's LAI that is above
																									  // the height of the lowest branch.															
				float crown_reach_above_neighbor = (tree->height - neighbor->height);
				float shade_contribution = LAI_intersection * (crown_reach_above_neighbor / crown_reach);
				shade += shade_contribution;
				if (isnan(shade_contribution)) {
					printf(" -- (in cell) shade_contribution is nan. LAI_intersection: %f, crown_reach_above_neighbor: %f, crown_reach: %f\n", LAI_intersection, crown_reach_above_neighbor, crown_reach);
				}
			}
		}
		return shade;
	}
	float get_shading_on_tree(Tree* tree, Population* population = 0) {
		return get_LAI_of_crown_intersection_and_above(tree, population);
	}
	float get_LAI() {
		return LAI;
	}
	void set_LAI(float _LAI) {
		LAI = _LAI;
	}
	void insert_stem(Tree* tree, float cell_area) {
		set_largest_stem(tree->dbh, tree->id);
		add_tree_if_not_present(tree, cell_area);
	}
	void remove_stem(Tree* tree, float cell_area, float cell_halfdiagonal_sqrt) {
		reset_largest_stem();
		remove_tree_if_smaller_than_cell(tree, cell_area, cell_halfdiagonal_sqrt);
	}
	void reset() {
		state = 0;
		trees.clear();
		LAI = 0;
		grass_LAI = 0;
		largest_stem = pair<float, int>(0, -1);
		seedling_present = false;
	}
	bool operator==(const Cell& cell) const
	{
		return pos == cell.pos;
	}
private:
	float get_leaf_area_over_cell_area(Tree* tree, float cell_area) {
		return (tree->LAI * tree->crown_area) / cell_area;
	}
	void add_LAI_of_tree_smaller_than_cell(Tree* tree, float cell_area) {
		LAI += get_leaf_area_over_cell_area(tree, cell_area);	// Obtain the leaf area of the tree and divide by the area of the cell to get the tree's
																// contribution to the cell's LAI.
	}
	void remove_LAI_of_tree_smaller_than_cell(Tree* tree, float cell_area) {
		LAI -= get_leaf_area_over_cell_area(tree, cell_area);	// Same as above but now we subtract the tree's contribution to the cell's LAI.
	}
	float LAI = 0; // Cumulative Leaf Area Index (LAI) of all trees in the cell.
};


class TreeDomainIterator {
public:
	TreeDomainIterator() = default;
	TreeDomainIterator(float _cell_width, Tree* _tree) {
		tree = _tree;
		cell_width = _cell_width;
		tree_center_gb = pair<float, float>(round(tree->position.first / cell_width), round(tree->position.second / cell_width));
		radius_gb = round(tree->radius / cell_width);
		x = tree_center_gb.first - radius_gb;
		y = tree_center_gb.second - radius_gb;
	}
	bool can_increment_x() {
		x++;
		if (x > tree_center_gb.first + radius_gb) {
			return false;
		}
		return true;
	}
	bool can_increment_y() {
		y++;
		if (y > tree_center_gb.second + radius_gb) {
			y = tree_center_gb.second - radius_gb; // Reset to lower bound.
			return false;
		}
		return true;
	}
	bool next() {
		if (begin) {
			begin = false; // If this is the first time 'next' is called, we leave the initialized x, y coordinates unchanged.
		}
		else if (!can_increment_y()) {
			if (!can_increment_x()) {
				return false;
			}
		}
		real_cell_position.first = (float)x * cell_width;
		real_cell_position.second = (float)y * cell_width;
		gb_cell_position.first = x;
		gb_cell_position.second = y;
		return true;
	}
	int x = 0;
	int y = 0;
	bool begin = true;
	Tree* tree;
	pair<float, float> tree_center_gb;
	pair<float, float> real_cell_position;
	pair<int, int> gb_cell_position;
	float radius_gb;
	float cell_width;
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
		cell_halfdiagonal_sqrt = help::get_dist(pair<float, float>(0, 0), pair<float, float>(0.5f * cell_width, 0.5f * cell_width));
		cell_area_half = cell_area * 0.5f;
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
			distribution[i].reset();
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
	int get_capped_center_idx(pair<float, float> &tree_center_gridbased) {
		pair<int, int> center = tree_center_gridbased;
		cap(center);
		return pos_2_idx(center);
	}
	bool populate_tree_domain(Tree* tree) {
		TreeDomainIterator it(cell_width, tree);
		while (it.next()) {
			if (tree->crown_area < cell_area_half) break; // Do not populate cells with trees that are smaller than half the cell area.
			if (tree->radius_spans(it.real_cell_position)) {
				set_to_forest(it.gb_cell_position, tree);
			}
		}
		int center_idx = get_capped_center_idx(it.tree_center_gb);
		distribution[center_idx].insert_stem(tree, cell_area);
		return true;
	}
	void burn_tree_domain(Tree* tree, queue<Cell*> &queue, float time_last_fire = -1, bool store_tree_death_in_color_distribution = true) {
		TreeDomainIterator it(cell_width, tree);
		while (it.next()) {
			if (tree->radius_spans(it.real_cell_position)) {
				Cell* cell = get_cell_at_position(it.gb_cell_position);
					
				// Remove tree id from cell->trees.
				cell->remove_tree(tree);

				// Set cell to savanna if the cumulative leaf area is less than half of the area of the cell
				// (leaf area < 0.5 * cell_area   <==>   (LAI * cell_area) < 0.5 * cell_area   <==>   LAI < 0.5)).
				if (cell->get_LAI() < 0.5f) { 
					queue.push(cell);
					set_to_savanna(cell->idx, time_last_fire);
					if (store_tree_death_in_color_distribution) state_distribution[cell->idx] = -6;
					continue;
				}
				state_distribution[cell->idx] = -5;
			}
		}
		int center_idx = get_capped_center_idx(it.tree_center_gb);
		distribution[center_idx].remove_stem(tree, cell_area, cell_halfdiagonal_sqrt);
	}
	void kill_tree_domain(Tree* tree, bool store_tree_death_in_color_distribution = true) {
		queue<Cell*> dummy;
		burn_tree_domain(tree, dummy, -1, store_tree_death_in_color_distribution);
	}
	void update_grass_LAIs_for_individual_tree(Tree* tree) {
		TreeDomainIterator it(cell_width, tree);
		while (it.next()) {
			distribution[pos_2_idx(it.gb_cell_position)].update_grass_LAI();
		}
	}
	int* get_state_distribution(bool collect = true) {
		if (collect) {
			for (int i = 0; i < no_cells; i++) {
				if (distribution[i].state == 1) state_distribution[i] = max(99.0f - (distribution[i].get_LAI() * 19.0f), 1);
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
	}
	void set_to_savanna(int idx, float _time_last_fire = -1) {
		no_savanna_cells += (distribution[idx].state == 1);
		no_forest_cells -= (distribution[idx].state == 1);

		distribution[idx].state = 0;
		if (_time_last_fire != -1) distribution[idx].time_last_fire = _time_last_fire;
	}
	float get_LAI_within_bb(pair<int, int> bb_min, pair<int, int> bb_max, float bb_area) {
		float cumulative_LAI = 0;
		for (int x = bb_min.first; x < bb_max.first; x++) {
			for (int y = bb_min.second; y < bb_max.second; y++) {
				Cell* cell = get_cell_at_position(pair<int, int>(x, y));
				cumulative_LAI += cell->get_LAI();
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
		return pair<int, int>(round(position.first / cell_width), round(position.second / cell_width));
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
	float cell_area_half = 0;
	float cell_halfdiagonal_sqrt = 0;
};
