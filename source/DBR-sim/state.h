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
	State(int _gridsize) {
		grid = Grid(_gridsize);
		population = Population();
	}
	void populate_grid() {
		cout << "Populating  grid..." << endl;
	}
	void set_tree_cover(float _tree_cover) {
		help::init_RNG();
		grid.reset();
		while (grid.get_tree_cover() < _tree_cover) {
			cout << "start\n";
			float x = help::get_rand_float(0, grid.size);
			float y = help::get_rand_float(0, grid.size);
			cout << "made randnums\n";
			pair<float, float> position = pair(x, y);
			cout << "made pos\n";
			Tree tree(position);
			cout << "made tree\n";
			tree.radius = TEMP_RADIUS;
			cout << "set tree radius.\n";
			grid.set_to_forest(position, &tree);
			cout << "Done All.\n\n";
			/*cout << "no forest cells: " << grid.no_forest_cells << endl;
			cout << "no savanna cells: " << grid.no_savanna_cells << endl;
			cout << "breaking? " << ((grid.get_tree_cover() > _tree_cover) ? "yes" : "no") << endl;
			printf("Current tree cover: %f\n", grid.get_tree_cover());
			cout << "current population size: " << population.members.size() << endl;*/
			continue;

			population.add(tree);
			grid.populate_tree_domain(&tree);
		}
		
	}
	float TEMP_RADIUS = 6.0;
	Grid grid;
	Population population;
};

