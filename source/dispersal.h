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
		Kernel* kernel = state->population.get_kernel(crop->id);
		pair<float, float> direction;
		get_direction(kernel, direction);
		float distance = get_dist(kernel);
		pair<float, float> release_location = state->grid.get_random_position_within_crown(state->population.get(crop->id));
		deposition_location = crop->origin + distance * direction;
	}
	bool germinate_seed(
		Crop* crop, State* state, pair<float, float>& deposition_location, int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions,
		int& no_competitions_with_older_trees, int& no_germination_attempts
	) {
		// Create seed and germinate if location has suitable conditions (existing LAI not too high).
		Seed seed(crop->strategy, deposition_location);
		return seed.germinate_if_location_is_viable(
			state, no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts
		);
	}
	void germinate_seeds_in_diaspore(
		Crop* crop, State* state, pair<float, float>& deposition_location, int& no_recruits,
		int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions, int& no_competitions_with_older_trees, int& no_germination_attempts
	) {
		for (int i = 0; i < crop->strategy.no_seeds_per_diaspore; i++) {
			no_recruits += germinate_seed(
				crop, state, deposition_location, no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts
			);
		}
	}
	void disperse_diaspore(
		Crop* crop, State* state, int& no_recruits, int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions, int& no_competitions_with_older_trees,
		int& no_germination_attempts
	) {
		pair<float, float> deposition_location;
		compute_deposition_location(crop, state, deposition_location);
		germinate_seeds_in_diaspore(
			crop, state, deposition_location, no_recruits, no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts
		);
	}
	int disperse_crop(
		Crop* crop, State* state, int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions, int& no_competitions_with_older_trees, int& no_germination_attempts
	) {
		int no_recruits = 0;
		for (int i = 0; i < crop->no_diaspora; i++) {
			disperse_diaspore(crop, state, no_recruits, no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts);
		}
		return no_recruits;
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
	int disperse(State* state, ResourceGrid* resource_grid, int no_seeds_to_disperse, float& fraction_time_spent_moving, int& no_seedlings_dead_due_to_shade,
		int& no_seedling_competitions, int& no_competitions_with_older_trees, int& no_germination_attempts, int verbosity = 0
	) {
		animals.place(state);
		int no_seeds_dispersed = 0;
		animals.disperse(
			no_seeds_dispersed, no_seeds_to_disperse, state, resource_grid, fraction_time_spent_moving, no_seedlings_dead_due_to_shade,
			no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts
		);
		return no_seeds_dispersed;
	}
	Animals animals;
};

