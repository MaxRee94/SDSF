#pragma once
#include <iostream>
#include "agents.h"
#include "grid.h"

/**
 * @class State
 * @brief Represents the overall state of the DBR simulation.
 * 
 * Contains the grid and population, and provides methods for state manipulation
 * such as populating the grid and setting tree cover.
 */
class State {
public:
	/**
	 * @brief Default constructor.
	 * 
	 * Initializes with default grid size and empty population.
	 */
	State() {
		grid = Grid();
		population = Population();
	}
	
	/**
	 * @brief Constructor with custom grid size.
	 * @param _gridsize The size of the grid to initialize.
	 */
	State(int _gridsize) {
		grid = Grid(_gridsize);
		population = Population();
	}
	
	/**
	 * @brief Populate the grid with trees.
	 * 
	 * Prints status message indicating grid population is starting.
	 */
	void populate_grid() {
		cout << "Populating  grid..." << endl;
	}
	
	/**
	 * @brief Set the target tree cover for the simulation.
	 * 
	 * Continuously adds trees to the grid until the target tree cover is reached.
	 * Uses random positions and a temporary radius for new trees.
	 * 
	 * @param _tree_cover The target tree cover as a fraction (0.0 to 1.0).
	 */
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
	
	float TEMP_RADIUS = 6.0;   ///< Temporary radius used for tree creation
	Grid grid;                 ///< The spatial grid for the simulation
	Population population;    ///< The population of trees
};