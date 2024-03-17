#pragma once
#include <iostream>
#include "agents.h"
#include "grid.h"


class State {
public:
	State() {
		grid = Grid();
		population = Population();
	}
	State(int gridsize, float cellsize, float mean_radius) {
		grid = Grid(gridsize, cellsize);
		population = Population(mean_radius);
	}
	void populate_grid() {
		cout << "Populating  grid..." << endl;
	}
	void set_tree_cover(float _tree_cover) {
		help::init_RNG();
		grid.reset();
		while (grid.get_tree_cover() < _tree_cover) {
			float x = help::get_rand_float(0, grid.size_r);
			float y = help::get_rand_float(0, grid.size_r);
			pair<float, float> position = pair(x, y);
			Tree* tree = population.add(position);
			grid.populate_tree_domain(tree);
			if (population.size() % 1000 == 0) {
				printf("Current tree cover: %f, current population size: %i\n", grid.get_tree_cover(), population.size());
			}
			continue;
		}
		printf("Final tree cover: %f\n", grid.get_tree_cover());
	}
	float TEMP_RADIUS = 200.0;
	Grid grid;
	Population population;
};

