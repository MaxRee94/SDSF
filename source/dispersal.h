#pragma once
#include "animals.h"



class Disperser {
public:
	Disperser() = default;
	virtual float get_dist(Kernel* kernel) {
		return kernel->get_ld_dist();
	}
	virtual pair<float, float> get_direction(Kernel* kernel) {
		return get_random_direction();
	}
	void disperse_seed(Crop* crop, State* state) {
		// Compute seed deposition location
		Kernel* kernel = state->population.get_kernel(crop->id);
		pair<float, float> direction = get_direction(kernel);
		float distance = get_dist(kernel);
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
	pair<float, float> get_direction(Kernel* kernel) {
		if (kernel->wind_direction_stdev < 360) return get_normal_distributed_direction(kernel->wind_direction, kernel->wind_direction_stdev);
		else return get_random_direction();
	}
	float get_dist(Kernel* kernel) {
		return kernel->get_wind_dispersed_dist();
	}
};


class AnimalDispersal : public Disperser {
public:
	AnimalDispersal() : Disperser() {};
	void disperse(State* state, ResourceGrid* resource_grid, int verbosity = 0) {
		Timer timer; timer.start();
		resource_grid->compute_cover_and_fruit_abundance();
		timer.stop(); printf("Time taken to compute cover and fruit abundance: %f\n", timer.elapsedSeconds());
		animals.place(state);
		int no_seeds_dispersed = 0;
		animals.disperse(no_seeds_dispersed, state, resource_grid);
		if (verbosity > 0) printf("Animal popsize: %i, No seeds dispersed: %i\n", animals.popsize(), no_seeds_dispersed);
	}
	Animals animals;
};

