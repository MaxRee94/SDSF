#pragma once
#include "diaspora.h"



class ResourceCell {
public:
	ResourceCell() = default;
	ResourceCell(pair<int, int> _position, int _idx): pos(_position), idx(_idx) {}
	void reset() {
		fruits.clear();
	}
	bool extract_random_fruit(Fruit &fruit) {
		if (fruits.no_fruits() <= 0) return false;
		int rand_idx = help::get_rand_int(0, fruits.no_fruits() - 1);
		bool success = fruits.get(rand_idx, fruit);
		fruits.remove(rand_idx);
		return success;
	}
	float get_fruit_abundance_index() {
		if (fruits.no_fruits() == 0) return 0.0f;
		return log10(fruits.no_fruits()); // Fruit Abundance Index, as used by Morales et al 2013.
	}
	pair<int, int> pos;
	pair<int, int> grid_bb_min;
	pair<int, int> grid_bb_max;
	int idx = 0;
	Fruits fruits;
};


class CoarseCell {
public:
	CoarseCell() = default;
	CoarseCell(pair<int, int> _position, int _idx): pos(_position), idx(_idx) {}
	pair<int, int> pos;
	pair<int, int> grid_bb_min;
	pair<int, int> grid_bb_max;
	int idx = 0;
	Fruits fruits;
};


class ResourceGrid : public Grid {
public:
	ResourceGrid() = default;
	ResourceGrid(State* _state, int _width, float _cell_width, vector<string>& species): Grid(_width, _cell_width) {
		width = _width;
		state = _state;
		grid = &state->grid;
		width_r = (float)width * cell_width;
		no_coarse_cells_along_x = round(sqrtf(width));
		no_coarse_cells = no_coarse_cells_along_x * no_coarse_cells_along_x;
		size = width * width;
		coarse_cells = new CoarseCell[no_coarse_cells];
		cells = new ResourceCell[size];
		selection_probabilities = DiscreteProbabilityModel(size);
		coarse_selection_probabilities = DiscreteProbabilityModel(no_coarse_cells);
		init_property_distributions(species);
		init_cells();
		init_coarse_cells();
		init_neighbor_offsets();
		printf("Initialized resource grid with width %i, where each resource cell has a discrete width of %i and a real width %f, %i coarse cells along x of grid width %i and real width %f\n",
			width, state->grid.width / width, cell_width, no_coarse_cells_along_x, (int)((float)width / (float)no_coarse_cells_along_x), coarse_cell_width_r
		);
	}
	void free() {
		delete[] cells;
		delete[] coarse_cells;
		delete_c();
		delete_f();
		delete_c_coarse();
		delete_f_coarse();
		delete[] d;
		delete[] cover;
		delete[] fruit_abundance;
		delete[] dist_aggregate;
		delete[] coarse_dist_aggregate;
		delete[] color_distribution;
		delete[] visits;
		selection_probabilities.free();
		coarse_selection_probabilities.free();
		Grid::free();
	}
	void delete_c() {
		for (auto it = c.begin(); it != c.end(); it++) {
			delete[] it->second;
		}
	}
	void delete_f() {
		for (auto it = f.begin(); it != f.end(); it++) {
			delete[] it->second;
		}
	}
	void delete_c_coarse() {
		for (auto it = c_coarse.begin(); it != c_coarse.end(); it++) {
			delete[] it->second;
		}
	}
	void delete_f_coarse() {
		for (auto it = f_coarse.begin(); it != f_coarse.end(); it++) {
			delete[] it->second;
		}
	}
	void init_property_distributions(vector<string> &species) {
		d = new float[size];
		d_coarse = new float[no_coarse_cells];
		cover = new float[size];
		fruit_abundance = new float[size];
		dist_aggregate = new float[size];
		coarse_dist_aggregate = new float[no_coarse_cells];
		color_distribution = new int[size];
		visits = new int[size];
		for (int i = 0; i < size; i++) visits[i] = 0;
		for (int i = 0; i < size; i++) dist_aggregate[i] = 0;
		for (int i = 0; i < no_coarse_cells; i++) coarse_dist_aggregate[i] = 0;
		for (int i = 0; i < species.size(); i++) {
			c[species[i]] = new float[size];
			f[species[i]] = new float[size];
			c_coarse[species[i]] = new float[no_coarse_cells];
			f_coarse[species[i]] = new float[no_coarse_cells];
		}
	}
	void reset() {
		for (int i = 0; i < size; i++) {
			cells[i].reset();
		}
		has_fruits = false;
		total_no_fruits = 0;
	}
	void add_crop(pair<float, float> position, Crop* crop) {
		ResourceCell* cell = get_resource_cell_at_position(position);
		cell->fruits.add_fruits(crop);
		total_no_fruits += crop->no_diaspora;
		has_fruits = total_no_fruits > 0;
	}
	float get_tree_cover_within_resourcegrid_bb(pair<int, int> bb_min, pair<int, int> bb_max) {
		int no_forest_cells = 0;
		int no_cells = 0;
		for (int x = bb_min.first; x < bb_max.first; x++) {
			for (int y = bb_min.second; y < bb_max.second; y++) {
				Cell* cell = grid->get_cell_at_position(pair<int, int>(x, y));
				no_forest_cells += cell->state;
				no_cells++;
			}
		}
		return (float)no_forest_cells / (float)no_cells;
	}
	void compute_cover_and_fruit_abundance() {
		for (int i = 0; i < size; i++) {
			ResourceCell* cell = &cells[i];
			cover[i] = get_tree_cover_within_resourcegrid_bb(cell->grid_bb_min, cell->grid_bb_max);
			cover[i] = asin(sqrt(cover[i]));
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
	void update_fruit_abundance(ResourceCell* cell, string species, map<string, float> &species_params) {
		fruit_abundance[cell->idx] = cell->get_fruit_abundance_index();
		f[species][cell->idx] = tanh(pow((fruit_abundance[cell->idx] / species_params["a_f"]), species_params["b_f"]));
	}
	void update_fruit_abundance(pair<float, float> position, string species, map<string, float> &species_params) {
		ResourceCell* cell = get_resource_cell_at_position(position);
		update_fruit_abundance(cell, species, species_params);
	}
	bool extract_fruit(pair<int, int> pos, Fruit &fruit) {
		ResourceCell* cell = get_resource_cell_at_position(pos);
		bool success = cell->extract_random_fruit(fruit);
		total_no_fruits -= success;
		return success;
	}
	pair<int, int> _idx_2_pos(int idx) {
		int x = idx % width;
		int y = idx / width;
		return pair<int, int>(x, y);
	}
	void get_random_location_within_cell(ResourceCell* cell, pair<float, float>& deposition_location) {
		deposition_location.first = help::get_rand_float((float)cell->pos.first * cell_width, (float)(cell->pos.first + 1) * cell_width);
		deposition_location.second = help::get_rand_float((float)cell->pos.second * cell_width, (float)(cell->pos.second + 1) * cell_width);
	}
	void get_random_location_within_cell(pair<float, float>& deposition_location) {
		ResourceCell* cell = get_resource_cell_at_position(deposition_location);
		get_random_location_within_cell(cell, deposition_location);
	}
	CoarseCell* get_coarse_cell_at_position(pair<int, int> pos) {
		coarse_cap(pos);
		return &coarse_cells[pos.second * no_coarse_cells_along_x + pos.first];
	}
	CoarseCell* get_coarse_cell_at_position(pair<float, float> _pos) {
		pair<int, int> pos = get_gridbased_coarsecell_position(_pos);
		return get_coarse_cell_at_position(pos);
	}
	ResourceCell* get_resource_cell_at_position(pair<int, int> pos) {
		cap(pos);
		return &cells[pos.second * width + pos.first];
	}
	ResourceCell* get_resource_cell_at_position(pair<float, float> _pos) {
		pair<int, int> pos = get_gridbased_position(_pos);
		return get_resource_cell_at_position(pos);
	}
	pair<float, float> get_real_cell_position(ResourceCell* cell) {
		return pair<float, float>(cell->pos.first * cell_width, cell->pos.second * cell_width);
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
		if (collect == "distance_single_coarse") {
			for (int i = 0; i < no_cells; i++) {
				pair<int, int> coarse_pos = normal_to_coarse_position(cells[i].pos);
				CoarseCell* coarse_cell = get_coarse_cell_at_position(coarse_pos);
				float color = d_coarse[coarse_cell->idx] * 10000;
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
				color_distribution[i] = c[species][i] * 100;
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
		else if (collect == "k_coarse") {
			for (int i = 0; i < no_cells; i++) {
				pair<int, int> coarse_pos = normal_to_coarse_position(cells[i].pos);
				CoarseCell* coarse_cell = get_coarse_cell_at_position(coarse_pos);
				float color = coarse_selection_probabilities.probabilities[coarse_cell->idx] * 1000000;
				color_distribution[i] = color;
			}
		}
		if (verbosity > 0) printf("collected %s for species %s \n", collect.c_str(), species.c_str());
		return color_distribution;
	}
	ResourceCell* select_random_cell() {
		int idx = help::get_rand_int(0, size - 1);
		return &cells[idx];
	}
	float compute_mean_across_coarse_cell(CoarseCell& coarse_cell, string species, string type) {
		float* source_distribution = 0;
		if (type == "c") source_distribution = c[species];
		else if (type == "f") source_distribution = f[species];
		float sum = 0.0f;
		for (int _x = 0; _x < coarse_cell_width; _x++) {
			for (int _y = 0; _y < coarse_cell_width; _y++) {
				int resourcegrid_x = coarse_cell.grid_bb_min.first + _x;
				int resourcegrid_y = coarse_cell.grid_bb_min.second + _y;
				int resourcegrid_idx = pos_2_idx(pair<int, int>(resourcegrid_x, resourcegrid_y));
				sum += source_distribution[resourcegrid_idx];
			}
		}
		return sum / (float)(coarse_cell_width * coarse_cell_width);
	}
	void compute_d(pair<float, float>& cur_position, CoarseCell coarse_cell, float a_d, float b_d) {
		float a_d_recipr = 1.0f / a_d;
		for (int i = 0; i < size; i++) {
			if (!resource_cell_is_within_coarse_cell(&cells[i], coarse_cell)) {
				d[i] = 0.0f;
				continue;
			}
			pair<float, float> cell_pos = get_real_position(cells[i].pos);
			float dist = get_resourcegrid_dist(cur_position, cell_pos);
			d[i] = tanh(pow((-dist * a_d_recipr), b_d));
			dist_aggregate[i] += d[i];
		}
	}
	void compute_coarse_d(pair<float, float>& cur_position, float a_d, float b_d) {
		float a_d_recipr = 1.0f / a_d;
		//cur_position = pair<float, float>(help::get_rand_float(0, 1000), help::get_rand_float(0, 1000));
		for (int i = 0; i < no_coarse_cells; i++) {
			pair<float, float> cell_pos = get_real_coarsecell_position(coarse_cells[i].pos);
			//cell_pos = cell_pos + 0.5f * pair<float, float>(coarse_cell_width, coarse_cell_width);
			float dist = get_resourcegrid_dist(cur_position, cell_pos);
			d_coarse[i] = tanh(pow((-dist * a_d_recipr), b_d));
			coarse_dist_aggregate[i] += d_coarse[i];
		}
	}
	void compute_c(string species, float a_c, float b_c) {
		float* _c = c[species];
		float a_c_recipr = 1.0f / a_c;
		for (int i = 0; i < size; i++) {
			_c[i] = tanh(pow((cover[i] * a_c_recipr), b_c));
		}
	}
	void update_coarse_distribution(string species, float a_c, float b_c, string type) {
		float* coarse_distr = 0;
		if (type == "c") coarse_distr = c_coarse[species];
		else if (type == "f") coarse_distr = f_coarse[species];
		for (int i = 0; i < no_coarse_cells; i++) {
			CoarseCell &coarse_cell = coarse_cells[i];
			coarse_distr[i] = compute_mean_across_coarse_cell(coarse_cell, species, type);
		}
	}
	void compute_f(string species, float a_f, float b_f) {
		float* _f = f[species];
		float a_f_recipr = 1.0f / a_f;
		for (int i = 0; i < size; i++) {
			_f[i] = tanh(pow((fruit_abundance[i] * a_f_recipr), b_f));
		}
	}
	void update_cover_and_fruit_probabilities(string species, map<string, float>& species_params) {
		compute_c(species, species_params["a_c"], species_params["b_c"]);
		compute_f(species, species_params["a_f"], species_params["b_f"]);
		update_coarse_distribution(species, species_params["a_f"], species_params["b_f"], "c");
		update_coarse_distribution(species, species_params["a_f"], species_params["b_f"], "f");
	}
	void update_coarse_probability_distribution(string species, map<string, float> &species_params, pair<float, float> &cur_position) {
		pair<int, int> gb_pos = get_gridbased_position(cur_position);
		cap(gb_pos);
		visits[pos_2_idx(gb_pos)] += 1;
		visits_sum += 1;
		compute_coarse_d(cur_position, species_params["a_d"], species_params["b_d"]);
		compute_coarse_k(species);
	}
	void reset_color_arrays() {
		for (int i = 0; i < size; i++) visits[i] = 0;
		for (int i = 0; i < size; i++) dist_aggregate[i] = 0;
		visits_sum = 0;
	}
	pair<int, int> idx_2_coarse_pos(int idx) {
		int x = idx % no_coarse_cells_along_x;
		int y = idx / no_coarse_cells_along_x;
		return pair<int, int>(x, y);
	}
	int coarse_pos_2_coarse_idx(pair<int, int> pos) {
		return no_coarse_cells_along_x * pos.second + pos.first;
	}
	pair<float, float> get_real_coarsecell_position(pair<int, int> position) {
		return pair<float, float>((float)position.first * coarse_cell_width_r, (float)position.second * coarse_cell_width_r);
	}
	pair<int, int> get_gridbased_coarsecell_position(pair<float, float> position) {
		return pair<int, int>(position.first / coarse_cell_width_r, position.second / coarse_cell_width_r);
	}
	CoarseCell* select_coarse_cell() {
		int coarse_idx = coarse_selection_probabilities.sample();
		return &coarse_cells[coarse_idx];
	}
	bool resource_cell_is_within_coarse_cell(ResourceCell* cell, CoarseCell &coarse_cell) {
		if (cell->pos.first >= coarse_cell.grid_bb_min.first && cell->pos.second >= coarse_cell.grid_bb_min.second &&
			cell->pos.first <= coarse_cell.grid_bb_max.first && cell->pos.second <= coarse_cell.grid_bb_max.second) {
			//printf("cell with pos %i, %i is within coarse cell with pos %i, %i \n", cell->pos.first, cell->pos.second, coarse_cell->grid_bb_min.first, coarse_cell->grid_bb_min.second);
			return true;
		}
		return false;
	}
	void update_selection_probabilities_within_coarse_cell(CoarseCell &coarse_cell, string species) {
		float sum = 0.0f;
		float* _c = c[species];
		float* _f = f[species];
		for (int i = 0; i < size; i++) {
			ResourceCell* cell = &cells[i];
			if (resource_cell_is_within_coarse_cell(cell, coarse_cell)) {
				selection_probabilities.probabilities[i] = d[i] * _c[i] * _f[i];
				sum += selection_probabilities.probabilities[i];
			}
			else selection_probabilities.probabilities[i] = 0.0f;
		}
		selection_probabilities.normalize(sum);
		selection_probabilities.build_cdf();
	}
	void add_random_offset(CoarseCell &coarse_cell) {
		// Offset the coarse cell randomly to create a falloff effect.
		pair<int, int> coarse_cell_offset = pair<int, int>(
			help::get_rand_int((-0.3f * coarse_cell_width), (0.3f * coarse_cell_width)),
			help::get_rand_int((-0.3f * coarse_cell_width), (0.3f * coarse_cell_width))
		);
		coarse_cell.grid_bb_min = coarse_cell.grid_bb_min + coarse_cell_offset;
		coarse_cell.grid_bb_max = coarse_cell.grid_bb_max + coarse_cell_offset;
	}
	ResourceCell* select_cell(string species, map<string, float>& species_params, pair<float, float> cur_position) {
		CoarseCell coarse_cell = *select_coarse_cell();
		add_random_offset(coarse_cell);
		compute_d(cur_position, coarse_cell, species_params["a_d"], species_params["b_d"]);
		update_selection_probabilities_within_coarse_cell(coarse_cell, species);
		int idx = selection_probabilities.sample();
		return &cells[idx];
	}
	void coarse_cap(pair<int, int>& position_grid) {
		if (position_grid.first < 0) position_grid.first = no_coarse_cells_along_x + (position_grid.first % no_coarse_cells_along_x);
		if (position_grid.second < 0) position_grid.second = no_coarse_cells_along_x + (position_grid.second % no_coarse_cells_along_x);
		position_grid.first %= no_coarse_cells_along_x;
		position_grid.second %= no_coarse_cells_along_x;
	}
	pair<int, int> normal_to_coarse_position(pair<int, int> position) {
		return pair<int, int>(position.first / coarse_cell_width, position.second / coarse_cell_width);
	}
	State* state = 0;
	Grid* grid = 0;
	ResourceCell* cells = 0;
	CoarseCell* coarse_cells = 0;
	map<string, float*> c;
	map<string, float*> c_coarse;
	map<string, float*> f;
	map<string, float*> f_coarse;
	float* dist_aggregate = 0;
	float* coarse_dist_aggregate = 0;
	float* cover = 0;
	float* fruit_abundance = 0;
	float* d = 0;
	float* d_coarse = 0;
	float visits_sum = 0;
	int* visits = 0;
	int* color_distribution = 0;
	int iteration = -1;
	int total_no_fruits = 0;
	int size = 0;
	int no_coarse_cells = 0;
	int no_coarse_cells_along_x = 0;
	int coarse_cell_width = 0;
	int coarse_cell_width_r = 0;
	DiscreteProbabilityModel selection_probabilities;
	DiscreteProbabilityModel coarse_selection_probabilities;
	pair<float, float>* neighbor_offsets = 0;
	bool has_fruits = false;

private:
	void init_cells() {
		int no_gridcells_along_x_per_resource_cell = (float)grid->width / (float)width;
		for (int i = 0; i < size; i++) {
			cells[i] = ResourceCell(idx_2_pos(i), i);
			cells[i].grid_bb_min = no_gridcells_along_x_per_resource_cell * cells[i].pos;
			cells[i].grid_bb_max = cells[i].grid_bb_min + pair<int, int>(no_gridcells_along_x_per_resource_cell, no_gridcells_along_x_per_resource_cell);
		}
	}
	void init_coarse_cells() {
		coarse_cell_width = (float)width / (float)no_coarse_cells_along_x;		// Number of gridcells along x per coarse cell.
		coarse_cell_width_r = (float)coarse_cell_width * cell_width;			// Real width of a coarse cell.
		for (int i = 0; i < no_coarse_cells; i++) {
			coarse_cells[i] = CoarseCell(idx_2_coarse_pos(i), i);
			coarse_cells[i].grid_bb_min = coarse_cell_width * coarse_cells[i].pos;
			coarse_cells[i].grid_bb_max = coarse_cells[i].grid_bb_min + pair<int, int>(coarse_cell_width, coarse_cell_width);
		}
	}
	void compute_k(string species) {
		float* _c = c[species];
		float* _f = f[species];
		float sum = 0.0f;
		for (int i = 0; i < size; i++) {
			selection_probabilities.probabilities[i] = d[i] * _c[i] * _f[i];
			sum += selection_probabilities.probabilities[i];
		}
		selection_probabilities.normalize(sum);
		selection_probabilities.build_cdf();
	}
	void compute_coarse_k(string species) {
		float* _c = c_coarse[species];
		float* _f = f_coarse[species];
		float sum = 0.0f;
		for (int i = 0; i < no_coarse_cells; i++) {
			coarse_selection_probabilities.probabilities[i] = d_coarse[i] * _c[i] * _f[i];
			sum += coarse_selection_probabilities.probabilities[i];
		}
		coarse_selection_probabilities.normalize(sum);
		coarse_selection_probabilities.build_cdf();
	}
};

