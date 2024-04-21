#pragma once
#include "diaspora.h"



class ResourceCell {
public:
	ResourceCell() = default;
	ResourceCell(pair<int, int> _position, int _idx): pos(_position), idx(_idx) {}
	void reset() {
		trees.clear();
		tree_cover = 0.0f;
		no_fruits = 0;
	}
	pair<int, int> pos;
	int idx = 0;
	vector<int> trees;
	float tree_cover = 0.0f;
	int no_fruits = 0.0f;
};


class ResourceGrid : public Grid {
public:
	ResourceGrid() = default;
	ResourceGrid(State* _state, int _width, float _cellsize): Grid(_width, _cellsize) {
		width = _width;
		state = _state;
		grid = &state->grid;
		width_r = (float)width * cellsize;
		size = width * width;
		cells = new ResourceCell[size];
		init_cells();
		init_neighbor_offsets();
	}
	void reset() {
		for (int i = 0; i < size; i++) {
			cells[i].reset();
		}
	}
	void init_cells() {
		for (int i = 0; i < size; i++) {
			cells[i] = ResourceCell(idx_2_pos(i), i);
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
	float get_dist(pair<float, float> a, pair<float, float> b, bool verbose = false) {
		vector<float> dists = { help::get_dist(a, b) };
		for (int i = 0; i < 8; i++) {
			float dist = help::get_dist(
				neighbor_offsets[i] * width_r, a
			);
			dists.push_back(dist);
		}
		if (verbose) printf("normal dist: %f, periodic dist: %f \n", dists[0], help::get_min(&dists));
		return help::get_min(&dists);
	}
	State* state = 0;
	Grid* grid = 0;
	int width = 0;
	int width_r = 0;
	int size = 0;
	float cellsize = 0.0f;
	ResourceCell* cells = 0;
	pair<float, float>* neighbor_offsets = 0;
};

