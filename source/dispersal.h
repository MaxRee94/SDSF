#pragma once
#include "animals.h"



class Disperser {
public:
	Disperser() = default;
	virtual float get_dist(Crop* crop, State* state) {
		return state->population.get_kernel(crop->id)->get_ld_dist();
	}
	void disperse_seed(Crop* crop, State* state) {
		// Compute seed deposition location
		pair<float, float> direction = get_random_direction();
		float distance = get_dist(crop, state);
		pair<float, float> seed_deposition_location = crop->origin + distance * direction;

		// Create seed and germinate if viable
		Seed seed(*(state->population.get_strat(crop->id)), seed_deposition_location);
		seed.germinate_if_location_is_viable(state);
	}
	void disperse_crop(Crop* crop, State* state) {
		for (int i = 0; i < crop->no_diaspora; i++) {
			disperse_seed(crop, state);
		}
	}
};


class WindDispersal : public Disperser {
public:
	WindDispersal() : Disperser() {};
	float get_dist(Crop* crop, State* state) {
		return state->population.get_kernel(crop->id)->get_wind_dispersed_dist();
	}
};


class AnimalDispersal : public Disperser {
public:
	AnimalDispersal() : Disperser() {};
	void disperse(State* state, ResourceGrid* resource_grid, int verbosity = 0) {
		resource_grid->compute_cover_and_fruit_abundance();
		animals.place(state);
		vector<Seed> seeds;
		animals.disperse(seeds, state, resource_grid);
		if (verbosity > 0) printf("Animal popsize: %i (%f), No seeds dispersed: %i\n", animals.popsize(), animals.total_no_animals, seeds.size());
		for (int i = 0; i < seeds.size(); i++) {
			bool germinated = seeds[i].germinate_if_location_is_viable(state);
			if (germinated) {
				int id = state->population.no_created_trees - 1;
				//printf("Germinated tree id: %i (key: %i)\n", state->population.get(id)->id, id);
			}
		}
	}
	Animals animals;
};

