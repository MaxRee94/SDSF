#pragma once
#include "state.h"
#include "grid_agent.forward.h"


using namespace help;


class Disperser {
public:
	Disperser() = default;
	Disperser(State* _state) {
		state = _state;
		Population* pop = &state->population;
	}
	void disperse(Crop* crop) {
		pair<float, float> direction = get_random_direction();
		float distance = crop->kernel->sample();
		pair<float, float> deposition_point = crop->origin + distance * direction;
		Strategy strategy = *(state->population.get_strat(crop->id));
		state->population.add(deposition_point, strategy, 0.1); // TEMP: Arbitrary starting radius of 0.1. TODO: replace with 0 once growth curve is implemented.
	}
	State* state = 0;
};


class Anemochory : Disperser {
public:
	Anemochory() = default;
	Anemochory(State* _state) : Disperser(_state) {

	}
};


class Zoochory : Disperser {
public:
	Zoochory() = default;
	Zoochory(State* _state) : Disperser(_state) {

	}
};

