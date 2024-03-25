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
	State(int gridsize, float cellsize, float max_radius) {
		population = Population(max_radius, cellsize);
		grid = Grid(gridsize, cellsize);
	}
	void repopulate_grid(bool verbose = false) {
		grid.reset();
		for (auto& tree : population.members) {
			grid.populate_tree_domain(&tree);
		}
		if (verbose) cout << "Repopulated grid." << endl;
	}
	void set_tree_cover(float _tree_cover) {
		help::init_RNG();
		grid.reset();
		while (grid.get_tree_cover() < _tree_cover) {
			pair<int, int> position = grid.get_random_real_position();
			Tree* tree = population.add(position);
			grid.populate_tree_domain(tree);
			if (population.size() % 1000 == 0) {
				printf("Current tree cover: %f, current population size: %i\n", grid.get_tree_cover(), population.size());
			}
			continue;
		}
		printf("Final tree cover: %f\n", grid.get_tree_cover());
	}
	Grid grid;
	Population population;
};

