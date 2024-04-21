#pragma once
#include "diaspora.h"



class ResourceCell {
public:
	ResourceCell() = default;
	ResourceCell(pair<int, int> _position, int _idx): pos(_position), idx(_idx) {}
	void reset() {
		trees.clear();
		tree_cover = 0.0f;
		fruits.clear();
	}
	bool extract_random_fruit(Fruit &fruit) {
		if (fruits.no_fruits() == 0) return false;
		int rand_idx = help::get_rand_int(0, fruits.no_fruits() - 1);
		fruits.get(rand_idx, fruit);
		fruits.remove(rand_idx);
		return true;
	}
	float get_fruit_abundance_index() {
		if (fruits.no_fruits() == 0) return 0.0f;
		return log10(fruits.no_fruits()); // Fruit Abundance Index, as used by Morales et al 2013.
	}
	pair<int, int> pos;
	pair<int, int> grid_bb_min;
	pair<int, int> grid_bb_max;
	int idx = 0;
	vector<int> trees;
	float tree_cover = 0.0f;
	Fruits fruits;
};


class ResourceGrid : public Grid {
public:
	ResourceGrid() = default;
	ResourceGrid(State* _state, int _width, float _cellsize, vector<string>& species): Grid(_width, _cellsize) {
		width = _width;
		state = _state;
		grid = &state->grid;
		pop = &state->population;
		width_r = (float)width * cellsize;
		size = width * width;
		cells = new ResourceCell[size];
		init_property_distributions(species);
		init_cells();
		init_neighbor_offsets();
	}
	void init_property_distributions(vector<string> &species) {
		d = new float[size];
		cover = new float[size];
		fruit_abundance = new float[size];
		k = new float[size];
		color_distribution = new float[size];
		for (int i = 0; i < species.size(); i++) {
			c[species[i]] = new float[size];
			f[species[i]] = new float[size];
		}
	}
	void reset() {
		for (int i = 0; i < size; i++) {
			cells[i].reset();
		}
		has_fruits = false;
	}
	void add_crop(Tree &tree, Crop* crop) {
		ResourceCell* cell = get_resource_cell_at_position(tree.position);
		cell->fruits.add_fruits(crop);
		total_no_fruits += crop->no_diaspora;
		cell->trees.push_back(tree.id);
		has_fruits = true;
	}
	float get_tree_cover_within_bb(pair<int, int> bb_min, pair<int, int> bb_max) {
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
			cover[i] = get_tree_cover_within_bb(cell->grid_bb_min, cell->grid_bb_max);
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
	bool extract_fruit(pair<int, int> pos, Fruit &fruit) {
		ResourceCell* cell = get_resource_cell_at_position(pos);
		bool success = cell->extract_random_fruit(fruit);
		fruit_abundance[cell->idx] = cell->get_fruit_abundance_index();
		return success;
	}
	pair<int, int> _idx_2_pos(int idx) {
		int x = idx % width;
		int y = idx / width;
		return pair<int, int>(x, y);
	}
	void get_random_location_within_cell(pair<float, float> &deposition_location) {
		ResourceCell* cell = get_resource_cell_at_position(deposition_location);
		deposition_location.first = help::get_rand_float((float)cell->pos.first * cellsize, (float)(cell->pos.first + 1) * cellsize);
		deposition_location.second = help::get_rand_float((float)cell->pos.second * cellsize, (float)(cell->pos.second + 1) * cellsize);
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
		return pair<float, float>(cell->pos.first * cellsize, cell->pos.second * cellsize);
	}
	float get_dist(pair<float, float> a, pair<float, float> b, bool manhattan = true, bool verbose = false) {
		vector<float> dists = { help::get_manhattan_dist(a, b)};
		float min_dist = 1e9f;
		int min_idx = 0;
		for (int i = 0; i < 8; i++) {
			float dist = help::get_manhattan_dist(
				neighbor_offsets[i] * width_r + b, a
			);
			if (dist < min_dist) {
				min_dist = dist;
				min_idx = i;
			}
			dists.push_back(dist);
		}
		if (verbose) printf("normal dist: %f, periodic dist: %f \n", dists[0], dists[min_idx]);
		return help::get_dist(a, neighbor_offsets[min_idx] * width_r + b);
	}
	float* get_color_distribution(string collect = "fruits") {
		if (collect == "fruits") {
			for (int i = 0; i < no_cells; i++) {
				color_distribution[i] = cells[i].fruits.no_fruits();
			}
		}
		return color_distribution;
	}
	void compute_d(pair<float, float> &cur_position, float a_d, float b_d) {
		float a_d_recipr = 1.0f / a_d;
		for (int i = 0; i < size; i++) {
			pair<float, float> cell_pos = get_real_position(_idx_2_pos(i));
			float dist = get_dist(cur_position, cell_pos) * 0.05f;
			d[i] = 1.0f - tanh(pow((-dist * a_d_recipr), b_d));
		}
	}
	void compute_c(string species, float a_c, float b_c) {
		float* _c = c[species];
		float a_c_recipr = 1.0f / a_c;
		for (int i = 0; i < size; i++) {
			_c[i] = tanh(pow((cover[i] * a_c_recipr), b_c));
			//printf("cover: %f, pow: %f, c: %f \n", cover[i], pow((cover[i] * a_c_recipr), b_c), _c[i]);
		}
	}
	void compute_f(string species, float a_f, float b_f) {
		float* _f = f[species];
		float a_f_recipr = 1.0f / a_f;
		for (int i = 0; i < size; i++) {
			_f[i] = tanh(pow((fruit_abundance[i] * a_f_recipr), b_f));
			//printf("fai: %f, pow: %f, f: %f \n", fruit_abundance[i], pow((fruit_abundance[i] * a_f_recipr), b_f), _f[i]);
		}
	}
	void update_probability_distribution(string species, map<string, float> &species_params, pair<float, float> &cur_position) {
		compute_d(cur_position, species_params["a_d"], species_params["b_d"]);
		compute_c(species, species_params["a_c"], species_params["b_c"]);
		compute_f(species, species_params["a_f"], species_params["b_f"]);
		compute_k(species);
	}
	State* state = 0;
	Grid* grid = 0;
	Population* pop = 0;
	ResourceCell* cells = 0;
	map<string, float*> c;
	map<string, float*> f;
	float* cover = 0;
	float* fruit_abundance = 0;
	float* d = 0;
	float* k = 0;
	pair<float, float>* neighbor_offsets = 0;
	int size = 0;
	float* color_distribution = 0;
	bool has_fruits = false;
	int total_no_fruits = 0;
private:
	void init_cells() {
		int no_gridcells_along_x_per_resource_cell = (float)grid->width / (float)width;
		for (int i = 0; i < size; i++) {
			cells[i] = ResourceCell(idx_2_pos(i), i);
			cells[i].grid_bb_min = no_gridcells_along_x_per_resource_cell * cells[i].pos;
			cells[i].grid_bb_max = cells[i].grid_bb_min + pair<int, int>(no_gridcells_along_x_per_resource_cell, no_gridcells_along_x_per_resource_cell);
		}
	}
	void compute_k(string species) {
		float* _c = c[species];
		float* _f = f[species];
		float sum = 0.0f;
		for (int i = 0; i < size; i++) {
			k[i] = d[i] * _c[i] * _f[i];
			sum += k[i];
		}
		float sum_recipr = 1.0f / sum;
		for (int i = 0; i < size; i++) {
			k[i] *= sum_recipr;
			printf("prob at %d: %f \n", i, k[i]);
		}
	}
};

