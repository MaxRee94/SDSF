#pragma once
#include "animals.h"



class Disperser {
public:
	Disperser() = default;
	virtual float get_dist(Kernel* kernel) {
		return kernel->get_ld_dist();
	}
	virtual void get_direction(Kernel* kernel, pair<float, float> &direction) {
		return get_random_unit_vector(direction);
	}
	void compute_deposition_location(Crop* crop, State* state, pair<float, float> &deposition_location) {
		// Compute seed deposition location
		Kernel* kernel = state->population.get_kernel(crop->id);
		pair<float, float> direction;
		get_direction(kernel, direction);
		float distance = get_dist(kernel);
		deposition_location = crop->origin + distance * direction;
	}
	bool seed_can_germinate(Crop* crop) {
		return true;
	}
	void germinate_seed(Crop* crop, State* state, pair<float, float>& deposition_location) {
		// Create seed and germinate if location has suitable conditions (existing LAI not too high).
		Seed seed(crop->strategy, deposition_location);
		seed.germinate_if_location_is_viable(state);
	}
	void attempt_seed_germination(Crop* crop, State* state, pair<float, float>& deposition_location) {
		for (int i = 0; i < crop->strategy.no_seeds_per_diaspore; i++) {
			if (seed_can_germinate(crop)) {
				germinate_seed(crop, state, deposition_location);
			}
		}
	}
	void disperse_diaspore(Crop* crop, State* state) {
		pair<float, float> deposition_location;
		compute_deposition_location(crop, state, deposition_location);
		//deposition_location = state->grid.get_random_real_position();
		attempt_seed_germination(crop, state, deposition_location);
	}
	void disperse_crop(Crop* crop, State* state) {
		for (int i = 0; i < crop->no_diaspora; i++) {
			disperse_diaspore(crop, state);
		}
	}
};


class WindDispersal : public Disperser {
public:
	WindDispersal() : Disperser() {};
	void get_direction(Kernel* kernel, pair<float, float>& direction) {
		if (kernel->wind_direction_stdev < 360) get_normal_distributed_direction(direction, kernel->wind_direction, kernel->wind_direction_stdev);
		else {
			get_random_unit_vector(direction);
		}
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

