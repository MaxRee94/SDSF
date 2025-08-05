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
	Grid(int _size) {
		size = _size;
		distribution = new Cell[size * size];
		no_savanna_cells = size;
		size_r = (float)size * cell_width;
	}
	void reset() {
		for (int i = 0; i < size * size; i++) {
			distribution[i].state = 0;
		}
		no_forest_cells = 0;
		no_savanna_cells = size * size;
	}
	void set_to_savanna(int idx) {
		if (distribution[idx].state == 1) {
			no_forest_cells -= 1;
			no_savanna_cells += 1;
		}
		distribution[idx].state = 0;
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
		pair<float, float> tree_center = get_gridbased_position(tree);
		cout << "tree center: " << tree_center.first << ", " << tree_center.second << endl;
		for (int x = tree_center.first - (tree->radius / 2); x < tree_center.first + (tree->radius / 2); x++) {
			for (int y = tree_center.second - (tree->radius / 2); y < tree_center.second + (tree->radius / 2); y++) {
				if (help::get_dist(pair<float, float>(x, y), tree_center) < tree->radius) {
					set_to_forest(pair<int, int>(x, y), tree);
				}
			}
		}
	}
	int* get_state_distribution() {
		if (state_distribution = 0) state_distribution = new int[size * size];
		for (int i = 0; i < size * size; i++) {
			state_distribution[i] = distribution[i].state;
		}
		return state_distribution;
	}
	void set_to_forest(int idx, Tree* tree) {
		if (distribution[idx].state == 0) {
			no_savanna_cells -= 1;
			no_forest_cells += 1;
		}
		distribution[idx].state = 1;
		distribution[idx].tree = tree;
	}
	void set_to_forest(pair<int, int> position_grid, Tree* tree) {
		set_to_forest(position_grid.second * size + position_grid.first, tree);
	}
	pair<int, int> get_gridbased_position(Tree* tree) {
		return pair<int, int>(tree->position.first / cell_width, tree->position.second / cell_width);
	}
	int size = 1000;
	float size_r = 0;
	float cell_width = 1.5;
	Cell* distribution = 0;
	int* state_distribution = 0;
	int no_savanna_cells = 0;
	int no_forest_cells = 0;
};
