#pragma once
#include "state.h"


class Dynamics {
public:
	Dynamics() = default;
	Dynamics(int _timestep) : timestep(_timestep) {
	};
	void init_state(int gridsize, float cellsize, float mean_radius) {
		state = State(gridsize, cellsize, mean_radius);
	}
	void update() {
		printf("updating\n");
	}
	int timestep = 0;
	State state;
};

