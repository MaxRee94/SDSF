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
	void disperse_seed(Crop* crop) {
		// Compute seed deposition location
		pair<float, float> direction = get_random_direction();
		float distance = get_dist(crop);
		pair<float, float> seed_deposition_location = crop->origin + distance * direction;

		// Create seed and germinate if viable
		Seed seed(*_pop->get_strat(crop->id), seed_deposition_location);
		seed.germinate_if_location_is_viable(state);
	}
	void disperse_crop(Crop* crop) {
		for (int i = 0; i < crop->no_diaspora; i++) {
			disperse_seed(crop);
		}
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
	void disperse(Animals &animals, ResourceGrid* resource_grid, int verbosity = 0) {
		resource_grid->compute_cover_and_fruit_abundance();
		animals.initialize_population(resource_grid);
		animals.place();
		vector<Seed> seeds;
		animals.disperse(seeds);
		if (verbosity > 0) printf("Animal popsize: %i (%f), No seeds dispersed: %i\n", animals.popsize(), animals.total_no_animals, seeds.size());
		for (int i = 0; i < seeds.size(); i++) {
			seeds[i].germinate_if_location_is_viable(state);
		}
	}
	Animals animals;
};

