#pragma once
#include "animals.h"



class Disperser {
public:
	Disperser() = default;
	Disperser(State* _state) {
		state = _state;
		_pop = &state->population;
		_grid = &state->grid;
	}
	virtual float get_dist(Crop* crop) {
		return crop->kernel->get_ld_dist();
	}
	void disperse(Crop* crop) {
		// Compute seed deposition location
		pair<float, float> direction = get_random_direction();
		float distance = get_dist(crop);
		pair<float, float> seed_deposition_location = crop->origin + distance * direction;

		// Create seed and germinate if viable
		Seed seed(crop->strategy, seed_deposition_location, state);
		seed.germinate_if_location_is_viable();
	}
	Population* _pop = 0;
	Grid* _grid = 0;
	State* state = 0;
};


class WindDispersal : public Disperser {
public:
	WindDispersal() = default;
	WindDispersal(State* _state) : Disperser(_state) {};
	float get_dist(Crop* crop) {
		return crop->kernel->get_wind_dispersed_dist();
	}
};


class AnimalDispersal : public Disperser {
public:
	AnimalDispersal() = default;
	AnimalDispersal(State* _state) : Disperser(_state) {};
	void disperse(Animals &animals, Fruits &fruits) {
		animals.initialize_population();
		animals.place();
		vector<Seed> seeds;
		animals.disperse(fruits, seeds);
		for (Seed& seed : seeds) {
			seed.germinate_if_location_is_viable();
		}
	}
	Animals animals;
};

