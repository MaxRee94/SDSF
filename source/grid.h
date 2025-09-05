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
	pair<int, int> pos;
	bool seedling_present = false;
	bool resprout_present = false;
	pair<float, int> stem = pair<float, int>(0, 0);	// < float: dbh of largest tree or seedling that has its stem in this cell,
													//	 int:   Largest tree id, or id of parent tree if the largest stem belongs to a seedling >
	bool cell_is_occupied_by_larger_stem(pair<float, int> &tree_proxy) {
		if (tree_proxy.first < stem.first) {
			return true;
		}
		return false;
	}
	bool seedling_is_shaded_out() {
		float shade = LAI * 0.18f; // Normalize to range [0, 1] by dividing by 5.0 (max LAI is 5.0).
		if (help::get_rand_float(0, 1) < shade) {
			return true;
		}
		return false;
	}
	void set_stem(float dbh, int id) {
		stem = pair<float, int>(dbh, id);
	}
	void reset_stem() {
		stem = pair<float, int>(0, 0);
	}
	void update_grass_LAI(float tree_LAI_local_neighborhood) {
		grass_LAI = compute_grass_LAI(tree_LAI_local_neighborhood);
	}
	bool is_hospitable(
		pair<float, int> tree_proxy, int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions,
		int& no_competitions_with_older_trees, int& no_cases_seedling_competition_and_shading, int& no_cases_oldstem_competition_and_shading
	) {
		// If tree_proxy is larger than the current largest stem in the cell, it is assumed to be able to outcompete the other tree.
		// If not, tree_proxy is assumed to be shaded out or otherwise outcompeted by the existing larger tree.
		//no_seedling_competitions += seedling_present;
		bool outcompeted = false;
		bool shaded_out = false;
		bool outcompeted_by_sapling = false;
		bool outcompeted_by_older_tree = false;
		if (cell_is_occupied_by_larger_stem(tree_proxy)) {
			outcompeted_by_sapling = seedling_present; // If there is an existing sapling in the cell, and its stem is larger than the new sapling's stem, the new sapling is outcompeted.
			outcompeted_by_older_tree = !seedling_present; // If there is an existing older stem in the cell, and it is larger than the new sapling's stem, the new sapling is outcompeted.
			no_competitions_with_older_trees += outcompeted_by_older_tree;
			no_seedling_competitions += outcompeted_by_sapling;
			outcompeted = true;
		}
		if (seedling_is_shaded_out()) {
			no_seedlings_dead_due_to_shade += !outcompeted; // This counter indicates how many seedlings died due to shading exclusively (i.e. without suffering from competition).
			//no_seedling_competitions -= seedling_present; // If the new sapling dies due to shade, the number of seedling competitions should remain the same.
			shaded_out = true;
		}

		// Determine viability of the cell for the new sapling.
		bool viable = !(outcompeted || shaded_out);

		// Track seed fate
		no_cases_seedling_competition_and_shading += outcompeted_by_sapling && shaded_out;
		no_cases_oldstem_competition_and_shading += outcompeted_by_older_tree && shaded_out;
		no_seedling_competitions += viable && seedling_present; // If the new sapling is viable, it will replace (outcompete) the existing sapling.

		return viable;
	}
	float query_grass_LAI() {
		return grass_LAI;
	}
	bool tree_is_present(Tree* tree) {
		return find(trees.begin(), trees.end(), tree->id) != trees.end();
	}
	bool tree_is_present(int id) {
		Tree dummy; dummy.id = id;
		return tree_is_present(&dummy);
	}
	bool is_sapling(Tree* tree, float cell_halfdiagonal_sqrt) {
		return tree->radius < cell_halfdiagonal_sqrt;
	}
	void add_tree_if_not_present(Tree* tree, float cell_area, bool _is_sapling = true, float cell_halfdiagonal_sqrt = 0) {
		if (!tree_is_present(tree)) {
			add_tree(tree, _is_sapling, cell_area);
		}
	}
	bool remove_tree_if_sapling(Tree* tree, float cell_area, float cell_halfdiagonal_sqrt, bool sapling = true) {
		if (is_sapling(tree, cell_halfdiagonal_sqrt)) {
			remove_tree(tree);
		}
		return true;
	}
	void add_tree(Tree* tree, float cell_area = 0, float cell_halfdiagonal_sqrt = 0) {
		trees.push_back(tree->id);
		if (is_sapling(tree, cell_halfdiagonal_sqrt))
			add_LAI_of_tree_sapling(tree, cell_area);
		else LAI += tree->LAI;
	}
	void remove_tree(Tree* tree, float cell_area = 0, float cell_halfdiagonal_sqrt = 0) {
		help::remove_from_vec(&trees, tree->id);
		if (is_sapling(tree, cell_halfdiagonal_sqrt)) remove_LAI_of_tree_sapling(tree, cell_area);
		else LAI -= tree->LAI;
	}
	float compute_grass_LAI(float tree_LAI) {
		tree_LAI = min(3.0f, tree_LAI);			// Done to avoid re-intersecting the y=0 line at about LAI=4.3 
												// (grass LAI would then (incorrectly) start rising again).
		float _grass_LAI = max(0, 0.241f * (tree_LAI * tree_LAI) - 1.709f * tree_LAI + 2.899f);		// Relationship between grass- and tree LAI from 
																									// Hoffman et al. (2012), figure 2b.									
		return _grass_LAI;
	}
	float get_fuel_load() {
		return 0.344946533f * grass_LAI; // Normalize to range [0, 1] by multiplying with inverse of max value (2.899).
	}
	float get_LAI_of_crown_intersection_and_above(Tree* tree, Population* population = 0) {
		if (population == nullptr) {
			return get_LAI();
		}
		float shade = 0;
		float crown_reach = tree->height - tree->lowest_branch;
		for (int tree_id : trees) {
			if (tree_id == tree->id) {
				shade += tree->LAI; // Add the tree's own LAI to the shade (since it will cast shade on itself).
				continue;
			}
			Tree* neighbor = population->get(tree_id);
			if (neighbor->height >= tree->height) shade += neighbor->LAI;
			else if (neighbor->height > tree->lowest_branch) { // Implies crown intersection; a portion of the neighbor's crown will cast shade on our tree.
				float neighbor_crown_reach = neighbor->height - neighbor->lowest_branch;
				float crown_intersection = neighbor->height - tree->lowest_branch;
				float LAI_intersection = neighbor->LAI * (crown_intersection / neighbor_crown_reach); // Fraction of the smaller tree's LAI that is above
																									  // the height of the lowest branch.															
				float crown_reach_overlap = crown_reach - (tree->height - neighbor->height);
				float shade_contribution = LAI_intersection * (crown_reach_overlap / crown_reach); // Fraction of our tree's crown that is affected by the shade cast by the smaller neighbor.
				shade += shade_contribution;
				if (shade_contribution < 0) printf("---------------------- ERROR: shade contribution negative. Shade contribution: %f\n", shade_contribution);
			}
		}
		return shade;
	}
	float get_LAI_of_taller_trees(Tree* tree, Population* population) {
		if (population == nullptr) {
			return get_LAI();
		}
		float LAI_taller_trees = 0;
		float crown_reach = tree->height - tree->lowest_branch;
		for (int tree_id : trees) {
			Tree* neighbor = population->get(tree_id);
			if (neighbor->height > tree->height) LAI_taller_trees += neighbor->LAI;
		}
		return LAI_taller_trees;
	}
	float get_shading_on_tree(Tree* tree, Population* population = 0) {
		return get_LAI_of_taller_trees(tree, population) + tree->LAI;
	}
	float get_LAI() {
		return LAI;
	}
	void insert_sapling(Tree* tree, float cell_area, float cell_halfdiagonal_sqrt) {
		insert_stem(tree, cell_area, cell_halfdiagonal_sqrt);
	}
	void set_LAI(float _LAI) {
		printf(" --------- WARNING: Setting cell LAI manually. This should only be done for testing purposes.\n");
		LAI = _LAI;
	}
	void insert_stem(Tree* tree, float cell_area, float cell_halfdiagonal_sqrt) {
		set_stem(tree->dbh, tree->id);
		add_tree_if_not_present(tree, cell_area, cell_halfdiagonal_sqrt);
	}
	void remove_stem(Tree* tree, float cell_area, float cell_halfdiagonal_sqrt) {
		reset_stem();
		remove_tree(tree, cell_area);
	}
	void reset() {
		state = 0;
		trees.clear();
		LAI = 0;
		grass_LAI = 0;
		stem = pair<float, int>(0, 0);
		seedling_present = false;
		resprout_present = false;
	}
	bool operator==(const Cell& cell) const
	{
		return pos == cell.pos;
	}
private:
	float get_leaf_area_over_cell_area(Tree* tree, float cell_area) {
		return (tree->LAI * tree->crown_area) / cell_area;
	}
	void add_LAI_of_tree_sapling(Tree* tree, float cell_area) {
		float _LAI = LAI;
		LAI += get_leaf_area_over_cell_area(tree, cell_area);	// Obtain the leaf area of the tree and divide by the area of the cell to get the tree's
																// contribution to the cell's LAI.
	}
	void remove_LAI_of_tree_sapling(Tree* tree, float cell_area) {
		LAI -= get_leaf_area_over_cell_area(tree, cell_area);	// Same as above but now we subtract the tree's contribution to the cell's LAI.
	}
	float LAI = 0; // Cumulative Leaf Area Index (LAI) of all trees in the cell.
	float grass_LAI = 0;
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
		y_lowbound = y;
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
			y = y_lowbound; // Reset to lower bound.
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
	void update_real_cell_position() {
		real_cell_position.first = tree_center_gb.first * cell_width;
		real_cell_position.second = tree_center_gb.second * cell_width;
	}
	int get_no_tree_cells() {
		return tree->crown_area / (cell_width * cell_width);
	}
	int x = 0;
	int y = 0;
	int y_lowbound = 0;
	int grid_width = 0;
	bool begin = true;
	Tree* tree;
	pair<float, float> tree_center_gb;
	pair<float, float> real_cell_position;
	pair<int, int> gb_cell_position;
	pair<int, int> grid_bb_min = pair<int, int>(-(2 << 28), -(2 << 28));
	pair<int, int> grid_bb_max = pair<int, int>((2 << 28), (2 << 28));
	float radius_gb;
	float cell_width;
};


class ForestCluster {
public:
	ForestCluster() = default;
	ForestCluster(
		int centr_x, int centr_y, int centr_idx, vector<int> _cells, vector<pair<int, int>> _perimeter,
		float _cell_width, int _id
	) {
		centroid_idx = centr_idx;
		centroid_x = centr_x;
		centroid_y = centr_y;
		cells = _cells;
		perimeter = _perimeter;
		cell_width = _cell_width;
		id = _id;
		perimeter_length = (float)perimeter.size() * cell_width;
		area = cells.size() * cell_width * cell_width;
	};
	vector<int> cells;
	vector<pair<int, int>> perimeter; // Indices of forest-savanna pairs on the edge of the cluster.
	vector<int> centroid;
	int centroid_idx = -1;
	int centroid_x = -1;
	int centroid_y = -1;
	int id = -1;
	float cell_width = 0;
	float perimeter_length = 0;
	float area = 0;
};


class Grid {
public:
	Grid() = default;
	Grid(int _width, float _cell_width) {
		printf("\nInitializing grid with width %i and cell width %f.\n", _width, _cell_width);
		width = _width;
		cell_width = _cell_width;
		cell_width_inv = 1.0f / cell_width;
		cell_half_width = cell_width * 0.5f;
		width_r = (float)width * cell_width;
		no_cells = width * width;
		no_savanna_cells = no_cells;
		init_grid_cells();
		init_neighbor_offsets();
		reset_state_distr();
		area = no_cells * cell_width * cell_width;
		cell_area = cell_width * cell_width;
		cell_area_inv = 1.0f / cell_area;
		cell_halfdiagonal_sqrt = help::get_dist(pair<float, float>(0, 0), pair<float, float>(0.5f * cell_width, 0.5f * cell_width));
		cell_area_half = cell_area * 0.5f;
	}
	void init_grid_cells() {
		distribution = make_shared<Cell[]>(no_cells);
		for (int i = 0; i < no_cells; i++) {
			pair<int, int> pos = idx_2_pos(i);
			distribution[i].pos = pos;
			distribution[i].idx = i;
		}
		state_distribution = make_shared<int[]>(no_cells);
	}
	void init_neighbor_offsets() {
		neighbor_offsets = make_shared<pair<int, int>[]>(8);
		/*neighbor_offsets = new pair<int, int>[4];
		neighbor_offsets[0] = pair<int, int>(-1, 0);
		neighbor_offsets[1] = pair<int, int>(0, -1);
		neighbor_offsets[2] = pair<int, int>(1, 0);
		neighbor_offsets[3] = pair<int, int>(0, 1);
		return;*/
		int q = 0;
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if (i == 0 && j == 0) continue;
				neighbor_offsets[q] = pair<int, int>(i, j);
				q++;
			}
		}
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
	double get_tree_cover() {
		tree_cover = (double)no_forest_cells / (double)(no_cells);
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
				add_tree_to_cell(it.gb_cell_position, tree);
			}
		}
		int center_idx = get_capped_center_idx(it.tree_center_gb);
		distribution[center_idx].insert_stem(tree, cell_area, cell_halfdiagonal_sqrt);
		return true;
	}
	pair<float, float> get_random_position_within_crown(
		Tree* tree, bool success,
		pair<int, int> grid_bb_min = pair<int, int>(-2 << 28, -2 << 28),
		pair<int, int> grid_bb_max = pair<int, int>(2 << 28, 2 << 28)
	) {
		TreeDomainIterator it(cell_width, tree);
		vector<pair<int, int>> cells;
		while (it.next()) {
			cap(it.gb_cell_position);
			if (it.gb_cell_position.first > grid_bb_max.first || it.gb_cell_position.second > grid_bb_max.second) continue; // Skip cells outside the bounding box.
			if (it.gb_cell_position.first < grid_bb_min.first || it.gb_cell_position.second < grid_bb_min.second) continue; // Skip cells outside the bounding box.
			it.update_real_cell_position();
			if (tree->radius_spans(it.real_cell_position)) {
				cells.push_back(it.gb_cell_position);
			}
		}
		success = cells.size() > 0;
		if (!success) return pair<float, float>(-1, -1);
		pair<int, int> random_cell = cells[help::get_rand_int(0, cells.size() - 1)];
		return get_random_location_within_cell(random_cell);
	}
	pair<float, float> get_random_position_within_crown(Tree* tree) {
		bool success = true;
		return get_random_position_within_crown(tree, success);
	}
	Cell* burn_tree_domain(Tree* tree, queue<Cell*> &queue, float time_last_fire = -1, bool store_tree_death_in_color_distribution = false,
		bool store_burn_events = true, int ignition_cell_idx = -1) {
		TreeDomainIterator it(cell_width, tree);
		while (it.next()) {
			if (tree->radius_spans(it.real_cell_position)) {
				Cell* cell = get_cell_at_position(it.gb_cell_position);
					
				// Remove tree id from cell->trees.
				cell->remove_tree(tree);

				// Set cell to savanna if the cumulative leaf area is less than half of the area of the cell
				// (leaf area < 0.5 * cell_area   <==>   (LAI * cell_area) < 0.5 * cell_area   <==>   LAI < 0.5).
				if (cell->get_LAI() < 1.0f) { 
					if (cell->idx != ignition_cell_idx) queue.push(cell); // The ignition cell (responsible for setting the tree on fire) is already in the queue.
					set_to_savanna(cell->idx, time_last_fire);
					if (store_tree_death_in_color_distribution) state_distribution[cell->idx] = -6;
					continue;
				}
				if (store_burn_events) state_distribution[cell->idx] = -5;
			}
		}
		int center_idx = get_capped_center_idx(it.tree_center_gb);
		distribution[center_idx].remove_stem(tree, cell_area, cell_halfdiagonal_sqrt);

		return &distribution[center_idx];
	}
	void kill_tree_domain(Tree* tree, bool store_tree_death_in_color_distribution = true) {
		queue<Cell*> dummy;
		burn_tree_domain(tree, dummy, -1, store_tree_death_in_color_distribution, false);
	}
	float get_cumulative_onering_LAI_for_cell(Cell* cell) {
		float LAI_sum = 0;
		for (int i = 0; i < 8; i++) {
			pair<int, int> pos = cell->pos + neighbor_offsets[i];
			Cell* neighbor = get_cell_at_position(pos);
			LAI_sum += neighbor->get_LAI();
		}
		return LAI_sum;
	}
	float get_tree_LAI_of_local_neighborhood(Cell* cell, bool debug=false) {
		float LAI = get_cumulative_onering_LAI_for_cell(cell);
		LAI += cell->get_LAI();
		LAI /= 9.0f; // Get the average LAI of the cell and its neighbors.
		
		if (debug) {
			float neighbor_LAI_sum = 0;
			for (int i = 0; i < 8; i++) {
				neighbor_LAI_sum += get_tree_LAI_of_local_neighborhood(get_cell_at_position(cell->pos + neighbor_offsets[i]), false);
			}
			float mean_neighbor_LAI = neighbor_LAI_sum / 8.0f;
			/*if (mean_neighbor_LAI > 1.0f && LAI < 0.5f) {
				printf("Average smoothed neighbor LAI: %f, smoothed cell LAI: %f, unsmoothed cell LAI: %f, average unsmoothed neighbor LAI: %f\n",
					mean_neighbor_LAI, LAI, cell->get_LAI(), get_cumulative_onering_LAI_for_cell(cell) / 8.0f);
			}*/
		}
		
		return LAI;
		//return cell->get_LAI();
	}
	bool is_forest(float LAI) {
		return LAI > 1.0f;
	}
	bool is_forest(int x, int y) {
		Cell* cell = get_cell_at_position(pair<int, int>(x, y));
		return is_forest(cell->get_LAI());
	}
	void update_grass_LAI(Cell* cell) {
		float tree_LAI_local_neighborhood = get_tree_LAI_of_local_neighborhood(cell);
		cell->update_grass_LAI(tree_LAI_local_neighborhood);
	}
	void update_grass_LAIs() {
		for (int i = 0; i < no_cells; i++) {
			update_grass_LAI(&distribution[i]);
		}
	}
	void update_grass_LAIs_for_individual_tree(Tree* tree) {
		TreeDomainIterator it(cell_width, tree);
		while (it.next()) {
			Cell* cell = get_cell_at_position(it.gb_cell_position);
			update_grass_LAI(cell);
		}
	}
	shared_ptr<int[]> get_state_distribution(int collect = 0) {
		if (collect > 0) {
			for (int i = 0; i < no_cells; i++) {
				if (collect == 1) {
					if (distribution[i].state == 1) state_distribution[i] = max(99.0f - (distribution[i].get_LAI() * 19.0f), 1);
				}
				else if (collect == 2) {
					state_distribution[i] = distribution[i].query_grass_LAI() * 33;
				}
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
	void add_tree_to_cell(int idx, Tree* tree) {
		distribution[idx].add_tree(tree);

		// If the cell was not previously classified as forest, check if it should be reclassified now.
		if (is_forest(distribution[idx].get_LAI())) {
			if (distribution[idx].state == 0) {
				// The cell was previously classified as savanna, but is now forest. Update the counters accordingly.
				no_savanna_cells--;
				no_forest_cells++;
			}
			// Set the cell state to forest.
			distribution[idx].state = 1;
		}
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
	void add_tree_to_cell(pair<int, int> position_grid, Tree* tree) {
		cap(position_grid);
		add_tree_to_cell(position_grid.second * width + position_grid.first, tree);
	}
	void set_to_savanna(pair<int, int> position_grid, float time_last_fire = -1) {
		cap(position_grid);
		set_to_savanna(pos_2_idx(position_grid), time_last_fire);
	}
	float get_cumulative_fuel_load() {
		float cumulative_load = 0;
		for (int i = 0; i < no_cells; i++) {
			cumulative_load += distribution[i].get_fuel_load();
		}
		return cumulative_load;
	}
	void get_forest_clusters() {
		std::vector<std::vector<bool>> visited(width, std::vector<bool>(width, false));

		// Directions for 4-neighbor connectivity
		const std::vector<std::pair<int, int>> directions = {
			{1,0}, {-1,0}, {0,1}, {0,-1}
		};

		for (int x = 0; x < width; ++x) {
			for (int y = 0; y < width; ++y) {
				if (!visited[x][y] && is_forest(x, y)) {
					// Start BFS for a new cluster
					std::queue<std::pair<int, int>> q;
					std::vector<int> cell_indices;
					vector<pair<int, int>> perimeter;

					q.push({ x, y });
					visited[x][y] = true;

					long long sum_x = 0, sum_y = 0;
					int count = 0;

					while (!q.empty()) {
						auto [cx, cy] = q.front();
						q.pop();

						// Add to cluster
						pair<int, int> pos(cx, cy);
						int index = pos_2_idx(pos);
						cell_indices.push_back(index);

						sum_x += cx;
						sum_y += cy;
						count++;

						// Explore neighbors
						for (auto [dx, dy] : directions) {
							int nx = cx + dx, ny = cy + dy;

							// Correct for periodic boundaries
							pair<int, int> neighbor(nx, ny);
							cap(neighbor);
							nx = neighbor.first;
							ny = neighbor.second;

							if (!visited[nx][ny] && is_forest(nx, ny)) {
								visited[nx][ny] = true;
								q.push({ nx, ny });
							}
							else if (!is_forest(nx, ny)) {
								// If the coordinate is valid but it is a savanna cell, we've encountered a forest-savanna edge.
								int savanna_cell_idx = pos_2_idx(pair(nx, ny));
								perimeter.push_back(pair<int, int>(index, savanna_cell_idx));
							}
						}
					}

					// Compute centroid
					double centroid_x = static_cast<double>(sum_x) / count;
					double centroid_y = static_cast<double>(sum_y) / count;
					int centroid_idx = pos_2_idx(pair<int, int>(round(centroid_x), round(centroid_y)));

					// Create ForestCluster
					clusters.push_back(ForestCluster(centroid_x, centroid_y, centroid_idx, cell_indices, perimeter, cell_width, clusters.size()));
				}
			}
		}
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
	float cell_width_inv = 0;
	float cell_half_width = 0;
	shared_ptr<Cell[]> distribution = 0;
	shared_ptr<int[]> state_distribution = 0;
	int no_savanna_cells = 0;
	int no_forest_cells = 0;
	float area = 0;
	float cell_area = 0;
	float cell_area_inv = 0;
	float cell_area_half = 0;
	float cell_halfdiagonal_sqrt = 0;
	vector<ForestCluster> clusters;
	shared_ptr<pair<int, int>[]> neighbor_offsets = 0;
};

