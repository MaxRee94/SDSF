#pragma once
#include "diaspora.h"



class ResourceCell {
public:
	ResourceCell() = default;
	ResourceCell(pair<int, int> _position, int _idx): pos(_position), idx(_idx) {}
	void reset() {
		fruits.clear();
		trees.clear();
	}
	bool extract_random_fruit(Fruit &fruit, int tree_id) {
		if (fruits.no_fruits() == 0) {
			return false;
		}
		bool success = fruits.get(fruit, tree_id);
		return success;
	}
	float get_fruit_abundance_index() {
		if (fruits.no_fruits() == 0) return 0.0f;
		return log10(fruits.no_fruits()); // Fruit Abundance Index, as used by Morales et al 2013.
	}
	void add_tree(int tree_id) {
		trees.push_back(tree_id);
	}
	pair<int, int> pos;
	pair<int, int> grid_bb_min;
	pair<int, int> grid_bb_max;
	pair<int, int> normal_grid_bb_min;
	pair<int, int> normal_grid_bb_max;
	vector<int> trees;
	int idx = 0;
	Fruits fruits;
};


class ResourceGrid : public Grid {
public:
	ResourceGrid() = default;
	~ResourceGrid() {
		free();
	}
	ResourceGrid(ResourceGrid&&) = default;
	ResourceGrid& operator=(ResourceGrid&&) = default;
	ResourceGrid(State* _state, int _width, float _cell_width, vector<string> _species, map<string, map<string, float>>& _animal_kernel_params): Grid(_width, _cell_width) {
		width = _width;
		state = _state;
		grid = &state->grid;
		width_r = (float)width * cell_width;
		size = width * width;
		lookup_table_size = size * size;
		cells = make_shared<ResourceCell[]>(size);
		selection_probabilities = DiscreteProbabilityModel(size);
		species = _species;
		animal_kernel_params = _animal_kernel_params;
		init_property_distributions(species);
		init_cells();
		init_neighbor_offsets();
	}
	void free() {
		printf("calling free on resource grid whose probmodel has id %i \n", selection_probabilities.id);
		cells = nullptr;
		delete_c();
		delete_f();
		delete[] d;
		d = nullptr;
		delete[] cover;
		cover = nullptr;
		delete[] fruit_abundance;
		fruit_abundance = nullptr;
		delete[] dist_aggregate;
		dist_aggregate = nullptr;
		delete[] color_distribution;
		color_distribution = nullptr;
		delete[] visits;
		visits = nullptr;
		delete[] neighbor_offsets;
		neighbor_offsets = nullptr;
		delete_lookup_table();
		Grid::free();
	}
	void delete_c() {
		for (auto it = c.begin(); it != c.end(); it++) {
			delete[] it->second;
			it->second = nullptr;
		}
	}
	void delete_f() {
		for (auto it = f.begin(); it != f.end(); it++) {
			delete[] it->second;
			it->second = nullptr;
		}
	}
	void delete_lookup_table() {
		for (auto it = dist_lookup_table.begin(); it != dist_lookup_table.end(); it++) {
			delete[] it->second;
			it->second = nullptr;
		}
	}
	void init_property_distributions(vector<string> &species) {
		d = new float[size];
		cover = new float[size];
		fruit_abundance = new float[size];
		dist_aggregate = new float[size];
		color_distribution = new int[size];
		visits = new int[size];
		for (int i = 0; i < size; i++) visits[i] = 0;
		for (int i = 0; i < size; i++) dist_aggregate[i] = 0;
		for (int i = 0; i < species.size(); i++) {
			c[species[i]] = new float[size];
			f[species[i]] = new float[size];
			dist_lookup_table[species[i]] = new float[lookup_table_size];
		}
	}
	void reset() {
		for (int i = 0; i < size; i++) {
			cells[i].reset();
		}
		has_fruits = false;
		total_no_fruits = 0;
	}
	void add_tree(Tree* tree) {
		ResourceCell* cell = get_resource_cell_at_position(tree->position);
		cell->add_tree(tree->id);
	}
	float get_tree_cover_within_resourcegrid_bb(ResourceCell* rcell, pair<int, int> bb_min, pair<int, int> bb_max, vector<int>& trees) {
		int no_forest_cells = 0;
		int no_cells = 0;
		for (int x = bb_min.first; x < bb_max.first; x++) {
			for (int y = bb_min.second; y < bb_max.second; y++) {
				Cell* cell = grid->get_cell_at_position(pair<int, int>(x, y));
				for (int tree_id : cell->trees) {
					if (!help::is_in(&trees, tree_id)) {
						trees.push_back(tree_id);
					}
					if (state->population.get(tree_id)->life_phase == 2) {
						rcell->fruits.add_fruits(
							state->population.get_crop(tree_id),
							state->population.get(tree_id)->crown_area * cell_area_inv
						);
					}
					//printf("Tree %i, fruit abundance: %i \n", tree_id, rcell->fruits[tree_id]);
				}
				no_forest_cells += cell->state;
				no_cells++;
			}
		}
		return (float)no_forest_cells / (float)no_cells;
	}
	void compute_cover() {
		for (int i = 0; i < size; i++) {
			ResourceCell* cell = &cells[i];
			if (cell->trees.size() == 0) cover[i] = 0.0f;
			cover[i] = get_tree_cover_within_resourcegrid_bb(cell, cell->grid_bb_min, cell->grid_bb_max, cell->trees);
			cover[i] = asin(sqrt(cover[i]));
		}
	}
	void compute_fruit_abundance() {
		for (int i = 0; i < size; i++) {
			ResourceCell* cell = &cells[i];
			fruit_abundance[i] = cell->get_fruit_abundance_index();
		}
	}
	void init_neighbor_offsets() {
		neighbor_offsets = new pair<float, float>[8];
		int q = 0;
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if (i == 0 && j == 0) continue;
				neighbor_offsets[q] = pair<float, float>(i, j);
				q++;
			}
		}
	}
	void get_random_stategrid_location(pair<float, float>& location) {
		get_random_location_within_cell(location);
		location = state->grid.get_gridbased_position(location);
		location = state->grid.cell_width * location; // Convert to real position, at the origin of a stategrid-cell.
	}
	bool extract_fruit(pair<float, float> pos, Fruit &fruit, int tree_id) {
		ResourceCell* cell = get_resource_cell_at_position(pos);
		bool success = cell->extract_random_fruit(fruit, tree_id);
		total_no_fruits -= success;
		return success;
	}
	pair<int, int> _idx_2_pos(int idx) {
		int x = idx % width;
		int y = idx / width;
		return pair<int, int>(x, y);
	}
	int get_random_forested_location(ResourceCell* cell, pair<float, float>& location) {
		// Attempt to find a fruit-producing tree. If none is found after 10 attempts, choose a non-fruit-producing tree.
		int tree_id = cell->trees[help::get_rand_int(0, cell->trees.size() - 1)];
		int i = 0;
		while (state->population.get(tree_id)->life_phase != 2 && i < 10) {
			tree_id = cell->trees[help::get_rand_int(0, cell->trees.size() - 1)];
			i++;
		}
		
		// Get random cell within the crown of the tree
		bool success = true;
		location = state->grid.get_random_position_within_crown(
			state->population.get(tree_id), success, cell->grid_bb_min, cell->grid_bb_max
		);
		if (!success) return get_random_forested_location(cell, location);

		return tree_id;
	}
	void get_random_location_within_cell(ResourceCell* cell, pair<float, float>& location) {
		location.first = help::get_rand_float((float)cell->pos.first * cell_width, (float)(cell->pos.first + 1) * cell_width);
		location.second = help::get_rand_float((float)cell->pos.second * cell_width, (float)(cell->pos.second + 1) * cell_width);
	}
	void get_random_location_within_cell(pair<float, float>& deposition_location) {
		ResourceCell* cell = get_resource_cell_at_position(deposition_location);
		get_random_location_within_cell(cell, deposition_location);
	}
	ResourceCell* get_resource_cell_at_position(pair<int, int> pos) {
		cap(pos);
		return &cells[pos.second * width + pos.first];
	}
	ResourceCell* get_resource_cell_at_position(pair<float, float> _pos) {
		pair<int, int> pos = get_rc_gridbased_position(_pos);
		//pos = pos - pair<int, int>(1, 1);
		return get_resource_cell_at_position(pos);
	}
	pair<float, float> get_real_cell_position(ResourceCell* cell) {
		return pair<float, float>(cell->pos.first * cell_width, cell->pos.second * cell_width);
	}
	pair<int, int> get_rc_gridbased_position(pair<float, float> position) {
		return pair<int, int>(position.first * cell_width_inv, position.second* cell_width_inv);
	}
	pair<float, float> get_rc_real_position(pair<int, int> position) {
		return pair<float, float>((float)position.first * cell_width, (float)position.second * cell_width);
	}
	pair<float, float> get_rc_real_position(int idx) {
		pair<int, int> position = idx_2_pos(idx);
		return pair<float, float>((float)position.first * cell_width, (float)position.second * cell_width);
	}
	ResourceCell* get_random_resource_cell() {
		int idx = help::get_rand_int(0, size - 1);
		return &cells[idx];
	}
	float get_resourcegrid_dist(pair<float, float> a, pair<float, float>* b, bool verbose = false) {
		// This function yields the minimum distance between a and b, taking into account periodic boundary conditions.
		// WARNING: Position b may be modified in the process to a position outside of the grid.
		vector<float> dists = { help::get_manhattan_dist(a, *b)};
		float min_dist = dists[0];
		int min_idx = 0;
		for (int i = 0; i < 8; i++) {
			float dist = help::get_manhattan_dist(
				neighbor_offsets[i] * width_r + *b, a
			);

			if (dist < min_dist && dist > 0) {
				min_dist = dist;
				min_idx = i + 1;
			}
			dists.push_back(dist);
		}
		float dist;
		if (min_idx != 0) {
			dist = help::get_dist(neighbor_offsets[min_idx - 1] * width_r + *b, a);
			*b = neighbor_offsets[min_idx - 1] * width_r + *b;
		}
		else dist = help::get_dist(a, *b);
		if (verbose) printf("normal dist: %f, periodic dist: %f \n", dists[0], dist);
		return dist;
	}
	float get_resourcegrid_dist(pair<float, float> a, pair<float, float> b, bool verbose = false) {
		// This function yields the minimum distance between a and b, taking into account periodic boundary conditions.
		// Position b is not modified in the process.
		pair<float, float> _b = b;
		return get_resourcegrid_dist(a, &_b, verbose);
	}
	pair<float, float> get_shortest_trajectory(pair<float, float> a, pair<float, float> _b, float &distance) {
		pair<float, float> b = _b;
		distance = get_resourcegrid_dist(a, &b); // Modify position b if needed to ensure it is at the shortest distance from a.
		//printf("new b: %f, %f. Trajectory: %f, %f \n", b.first, b.second, (b-a).first, (b-a).second);
		return b - a;
	}
	int* get_color_distribution(string species, string collect = "distance", int verbosity = 0) {
		if (collect == "distance") {
			for (int i = 0; i < no_cells; i++) {
				float color = dist_aggregate[i] * 10000;
				color_distribution[i] = color;
			}
		}
		if (collect == "distance_single") {
			for (int i = 0; i < no_cells; i++) {
				float color = d[i] * 10000;
				color_distribution[i] = color;
			}
		}
		else if (collect == "fruits") {
			for (int i = 0; i < no_cells; i++) {
				color_distribution[i] = f[species][i] * 100;
			}
		}
		else if (collect == "cover") {
			for (int i = 0; i < no_cells; i++) {
				color_distribution[i] = c[species][i] * 100000;
			}
		}
		else if (collect == "visits") {
			float sum_recipr = 1.0f / visits_sum;
			int max_visits = 0;
			for (int i = 0; i < no_cells; i++) {
				if (visits[i] > max_visits) max_visits = visits[i];
				color_distribution[i] = (float)visits[i] * sum_recipr * 100000;
			}
		}
		else if (collect == "k") {
			for (int i = 0; i < no_cells; i++) {
				float color = selection_probabilities.probabilities[i] * 1000000;
				color_distribution[i] = color;
			}
		}
		if (verbosity > 0) printf("collected %s for species %s \n", collect.c_str(), species.c_str());
		return color_distribution;
	}
	float* get_lookup_table(string species) {
		return dist_lookup_table[species];
	}
	ResourceCell* select_random_cell() {
		int idx = help::get_rand_int(0, size - 1);
		return &cells[idx];
	}
	void precompute_dist_lookup_table(string species) {
		float a_d = animal_kernel_params[species]["a_d"];
		float b_d = animal_kernel_params[species]["b_d"];
		float a_d_recipr = 1.0f / a_d;
		for (int i = 0; i < size; i++) {
			pair<float, float> curpos = get_rc_real_position(cells[i].pos);
			for (int j = 0; j < size; j++) {
				pair<float, float> target_pos = get_rc_real_position(cells[j].pos);
				float dist = get_resourcegrid_dist(curpos, target_pos);
				float val = tanh(pow((-dist * a_d_recipr), b_d));
				dist_lookup_table[species][size * (cells[i].pos.first + cells[i].pos.second * width) + cells[j].pos.first + cells[j].pos.second * width] = val;
			}
			if (i % width == 0) printf("Finished row %i / %i\n", i / width, width);
		}
		printf("Computed dist lookup table for species %s \n", species.c_str());
	}
	void set_dist_lookup_table(float* lookup_table, string species) {
		for (int i = 0; i < width * width * width * width; i++) {
			dist_lookup_table[species][i] = lookup_table[i];
		}
	}
	void compute_d(pair<int, int>& curpos, string species) {
		for (int i = 0; i < size; i++) {
			//pair<float, float> cell_pos = get_real_position(cells[i].pos);
			//float dist = get_resourcegrid_dist(get_real_position(curpos), cell_pos);
			//d[i] = tanh(pow((-dist * 0.5), 0.4));
			//continue;
			
			pair<int, int> target_pos = cells[i].pos;
		    int lookup_idx = size * (curpos.first + curpos.second * width) + target_pos.first + target_pos.second * width;
			/*if (lookup_idx < 0 || lookup_idx > size * size) {
				printf("cur pos: %i, %i, target pos: %i, %i \n", curpos.first, curpos.second, target_pos.first, target_pos.second);
				printf("Idx: %i (limits: %i - %i) \n", lookup_idx, 0, size * size);
			}*/
			d[i] = dist_lookup_table[species][lookup_idx];
		}
	}
	void compute_c(string species, float a_c, float b_c) {
		float* _c = c[species];
		float a_c_recipr = 1.0f / a_c;
		for (int i = 0; i < size; i++) {
			_c[i] = tanh(pow((cover[i] * a_c_recipr), b_c));
		}
	}
	void compute_f(string species, float a_f, float b_f) {
		float* _f = f[species];
		float a_f_recipr = 1.0f / a_f;
		for (int i = 0; i < size; i++) {
			_f[i] = tanh(pow((fruit_abundance[i] * a_f_recipr), b_f));
		}
	}
	void update_cover_probabilities(string species, map<string, float>& species_params) {
		compute_cover();
		compute_c(species, species_params["a_c"], species_params["b_c"]);
	}
	void update_fruit_probabilities(string species, map<string, float>& species_params) {
		compute_fruit_abundance();
		compute_f(species, species_params["a_f"], species_params["b_f"]);
	}
	void reset_color_arrays() {
		for (int i = 0; i < size; i++) visits[i] = 0;
		for (int i = 0; i < size; i++) dist_aggregate[i] = 0;
		visits_sum = 0;
	}
	ResourceCell* select_cell(string species, pair<float, float> cur_position, bool fruit_agnostic_selection = false) {
		pair<int, int> gridbased_curpos = get_rc_gridbased_position(cur_position);
		cap(gridbased_curpos);
		compute_d(gridbased_curpos, species);
		compute_k(species, fruit_agnostic_selection);
		int idx = selection_probabilities.sample();
		visits[idx] += 1;
		visits_sum += 1;
		return &cells[idx];
	}
	State* state = 0;
	Grid* grid = 0;
	shared_ptr<ResourceCell[]> cells = 0;
	map<string, float*> c;
	map<string, float*> f;
	map<string, float*> dist_lookup_table;
	vector<string> species;
	map<string, map<string, float>> animal_kernel_params;	
	float* dist_aggregate = 0;
	float* cover = 0;
	float* fruit_abundance = 0;
	float* d = 0;
	float visits_sum = 0;
	int* visits = 0;
	int* color_distribution = 0;
	int iteration = -1;
	int total_no_fruits = 0;
	int lookup_table_size = 0;
	int size = 0;
	DiscreteProbabilityModel selection_probabilities;
	pair<float, float>* neighbor_offsets = 0;
	bool has_fruits = false;

private:
	void init_cells() {
		int no_gridcells_along_x_per_resource_cell = round((float)grid->width / (float)width);
		for (int i = 0; i < size; i++) {
			cells[i] = ResourceCell(idx_2_pos(i), i);
			cells[i].grid_bb_min = no_gridcells_along_x_per_resource_cell * cells[i].pos;
			cells[i].grid_bb_max = cells[i].grid_bb_min + pair<int, int>(no_gridcells_along_x_per_resource_cell - 1, no_gridcells_along_x_per_resource_cell - 1);
		}
	}
	void compute_k(string species, bool try_fruit_agnostic_selection = false) {
		float* _c = c[species];
		float* _f = f[species];
		float sum = 0.0f;
		for (int i = 0; i < size; i++) {
			if (try_fruit_agnostic_selection) {
				selection_probabilities.probabilities[i] = d[i] * _c[i];
			}
			else {
				selection_probabilities.probabilities[i] = d[i] * _c[i] * _f[i];
				sum += selection_probabilities.probabilities[i];
			}
		}
		selection_probabilities.normalize(sum);
		selection_probabilities.build_cdf();
	}
};

