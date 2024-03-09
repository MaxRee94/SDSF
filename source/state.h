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
		cout << "Populating grid..." << endl;
	}

	Grid grid;
	Population population;
};

