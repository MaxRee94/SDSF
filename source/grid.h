#pragma once
#include "agents.h"


class Cell {
public:
	Cell() = default;
	int state = 0;
	Tree* tree = 0;
};


class Grid {
public:
	Grid() = default;
	Grid(int _size, float _cellsize) {
		size = _size;
		cellsize = _cellsize;
		distribution = new Cell[size * size];
		no_savanna_cells = size * size;
		size_r = (float)size * cellsize;
		no_cells = size * size;
	}
	void reset() {
		for (int i = 0; i < size * size; i++) {
			distribution[i].state = 0;
		}
		no_forest_cells = 0;
		no_savanna_cells = size * size;
	}
	void redo_count() {
		no_savanna_cells = 0;
		no_forest_cells = 0;
		for (int i = 0; i < size * size; i++) {
			no_savanna_cells += !distribution[i].state;
			no_forest_cells += distribution[i].state;
		}
	}
	float get_tree_cover() {
		return (float)no_forest_cells / (float)no_savanna_cells;
	}
	Cell* get_cell_at_position(pair<float, float> pos) {
		return &distribution[(int)(pos.second * size + pos.first)];
	}
	void populate_tree_domain(Tree* tree) {
		pair<float, float> tree_center_gb = get_gridbased_position(tree);
		int radius_gb = tree->radius / cellsize;
		int radius_gb_half = radius_gb / 2.0;
		for (int x = tree_center_gb.first - radius_gb; x < tree_center_gb.first + radius_gb; x++) {
			for (int y = tree_center_gb.second - radius_gb; y < tree_center_gb.second + radius_gb; y++) {
				if (help::get_dist(pair<float, float>(x, y), tree_center_gb) < radius_gb) {
					set_to_forest(pair<int, int>(x, y), tree);
				}
			}
		}
	}
	int* get_state_distribution() {
		if (state_distribution == 0) state_distribution = new int[size * size];
		for (int i = 0; i < size * size; i++) {
			state_distribution[i] = distribution[i].state;
		}
		return state_distribution;
	}
	void set_to_forest(int idx, Tree* tree) {
		if (distribution[idx].state == 0) {
			no_savanna_cells--;
			no_forest_cells++;
		}
		distribution[idx].state = 1;
		distribution[idx].tree = tree;
	}
	void set_to_savanna(int idx, Tree* tree) {
		if (distribution[idx].state == 1) {
			no_savanna_cells++;
			no_forest_cells--;
		}
		distribution[idx].state = 0;
		distribution[idx].tree = 0;
	}
	void set_to_forest(pair<int, int> position_grid, Tree* tree) {
		cap(position_grid);
		set_to_forest(position_grid.second * size + position_grid.first, tree);
	}
	void set_to_savanna(pair<int, int> position_grid, Tree* tree) {
		cap(position_grid);
		set_to_savanna(position_grid.second * size + position_grid.first, tree);
	}
	void cap(pair<int, int> &position_grid) {
		if (position_grid.first < 0) position_grid.first = size + position_grid.first;
		if (position_grid.second < 0) position_grid.second = size + position_grid.second;
		position_grid.first %= size;
		position_grid.second %= size;
	}
	pair<int, int> get_gridbased_position(Tree* tree) {
		return pair<int, int>(tree->position.first / cellsize, tree->position.second / cellsize);
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
