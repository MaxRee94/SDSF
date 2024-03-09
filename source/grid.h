#pragma once
#include "helpers.h"



class Grid {
public:
	Grid() = default;
	Grid(int _gridsize) {
		gridsize = _gridsize;
		distribution = new int[gridsize * gridsize];
	}
	int gridsize = 1000;
	float cellsize = 1.5;
	int* distribution = 0;
};
