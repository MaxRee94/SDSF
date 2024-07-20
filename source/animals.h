#pragma once
#include "animal_resources.h"
#include "grid_agent.forward.h"


class Animal {
public:
	Animal() = default;
	Animal(map<string, float> _traits, string _species) {
		traits = _traits;
		traits["a_c_recipr"] = 1.0f / traits["a_c"];
		traits["a_f_recipr"] = 1.0f / traits["a_f"];
		species = _species;
		recipr_speed = 1.0f / traits["speed"];
		init_prob_distributions();
	}
	void init_prob_distributions() {
		gut_passage_time_distribution = GammaProbModel(traits["gut_passage_time_shape"], traits["gut_passage_time_scale"]);
		rest_time_distribution = GammaProbModel(traits["rest_time_shape"], traits["rest_time_scale"]);
	}
	void reset() {
		curtime = 0;
		stomach_content.clear();
		total_no_seeds_consumed = 0;
	}
	void update(
		int& no_seeds_dispersed, int _iteration, State* state, ResourceGrid* resource_grid, int& no_seeds_eaten,
		float& time_spent_resting, float& distance_travelled, float& average_gut_passage_time, float& average_gut_time,
		int& no_seeds_defecated, float& time_spent_moving, int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions,
		int& no_competitions_with_older_trees, int& no_germination_attempts
	) {
		iteration = _iteration;

		pair<int, int> move_time_interval(curtime, -1);
		move(resource_grid, distance_travelled);
		move_time_interval.second = curtime;
		time_spent_moving += travel_time;

		digest(
			resource_grid, no_seeds_dispersed, no_seeds_defecated, state, move_time_interval,
			no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts
		);

		pair<int, int> rest_time_interval(curtime, -1);
		rest(time_spent_resting);
		eat(resource_grid, rest_time_interval.first, no_seeds_eaten);
		rest_time_interval.second = curtime;

		digest(
			resource_grid, no_seeds_dispersed, no_seeds_defecated, state, rest_time_interval,
			no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts
		);
	}
	void rest(float& time_spent_resting) {
		moving = false;
		float rest_time = rest_time_distribution.get_gamma_sample() * 60.0f; // Multiply by 60 to convert minutes to seconds.
		curtime += rest_time;
		time_spent_resting += rest_time;
	}
	float get_biomass_appetite() {
		float appetite = 0;
		for (int i = 0; i < animal_group_size; i++) {
			// Check if the current animal will engage in fruit consumption.
			if (help::get_rand_float(0, 1) > traits["probability_fruit_consumption"]) continue;

			// Add the mass of the average feeding bout to the total appetite.
			appetite += traits["average_fruit_mass_consumed"];
		}
		return appetite;
	}
	void eat(ResourceGrid* resource_grid, float begin_time, int& no_seeds_eaten) {
		float biomass_appetite = get_biomass_appetite();
		vector<pair<float, int>> consumed_seed_gpts_plus_cropids = {};
		float fruit_mass = resource_grid->state->population.get_crop(last_tree_visited)->strategy.diaspore_mass;
		while (biomass_appetite >= fruit_mass) {
			Fruit fruit;

			bool fruit_available = resource_grid->extract_fruit(position, fruit, last_tree_visited);

			if (!fruit_available) break;

			eat_fruit(fruit, consumed_seed_gpts_plus_cropids);

			no_seeds_eaten += fruit.strategy.no_seeds_per_diaspore;
			biomass_appetite -= fruit_mass;
		}
		set_defecation_times(begin_time, consumed_seed_gpts_plus_cropids);
	}
	void set_defecation_times(float begin_time, vector<pair<float, int>>& consumed_seed_gpts_plus_cropids) {
		float time_stepsize = (curtime - begin_time) / (float)consumed_seed_gpts_plus_cropids.size();
		for (int i = consumed_seed_gpts_plus_cropids.size() - 1; i >= 0; i--) {
			float ingestion_time = curtime - (float)(i + 1) * time_stepsize;
			float gut_passage_time = consumed_seed_gpts_plus_cropids[i].first;
			int crop_id = consumed_seed_gpts_plus_cropids[i].second;
			int defecation_time = round(gut_passage_time + ingestion_time);
			if (stomach_content[defecation_time].size() == 0) stomach_content[defecation_time] = { crop_id };
			else stomach_content[defecation_time].push_back(crop_id);
		}
	}
	void eat_fruit(Fruit &fruit, vector<pair<float, int>>& consumed_seed_gpts_plus_cropids) {
		vector<Seed> seeds;
		for (int i = 0; i < fruit.strategy.no_seeds_per_diaspore; i++) {
			float gut_passage_time = gut_passage_time_distribution.get_gamma_sample();
			
			// Multiply GPT by 60 to convert minutes to seconds.
			gut_passage_time *= 60.0f;

			consumed_seed_gpts_plus_cropids.push_back(pair<float, int>(gut_passage_time, fruit.id));
			total_no_seeds_consumed++;
		}
	}
	void digest(
		ResourceGrid* resource_grid, int& no_seeds_dispersed, int& no_seeds_defecated,
		State* state, pair<int, int>& time_interval, int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions,
		int& no_competitions_with_older_trees, int& no_germination_attempts
	) {
		vector<int> deletion_schedule;
		for (int defecation_time = time_interval.first; defecation_time <= time_interval.second; defecation_time++) {
			if (stomach_content[defecation_time].size() == 0) continue; // Skip if there are no seeds with this defecation time in the stomach.

			for (int crop_id : stomach_content[defecation_time]) {
				if (iteration < 5) continue; // Do not disperse in the first 5 iterations (after Morales et al 2013)

				// Create seed and calculate time since defecation.
				Seed seed(state->population.get_crop(crop_id)->strategy);
				float time_since_defecation = curtime - defecation_time;

				// Defecate seed.
				defecate(
					seed, time_since_defecation, no_seeds_dispersed, resource_grid, no_seedlings_dead_due_to_shade,
					no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts
				);

				no_seeds_defecated++;
			}
			stomach_content[defecation_time].clear(); // Remove all seeds with this defecation time from the stomach.
		}
	}
	pair<float, float> select_destination(ResourceGrid* resource_grid) {
		ResourceCell* cell = resource_grid->select_cell(species, position);
		pair<float, float> destination;
		if (cell->trees.size() == 0) return select_destination(resource_grid);
		last_tree_visited = resource_grid->get_random_forested_location(cell, destination);
		return destination;
	}
	void move(ResourceGrid* resource_grid, float& distance_travelled) {
		moving = true;
		pair<float, float> destination = select_destination(resource_grid);
		float distance;
		trajectory = resource_grid->get_shortest_trajectory(position, destination, distance);
		position = destination;
		distance_travelled += distance;
		travel_time = distance * recipr_speed;
		curtime += travel_time;
	}
	pair<float, float> get_backwards_traced_location(float time_since_defecation) {
		return position - (time_since_defecation / travel_time) * trajectory;
	}
	void defecate(
		Seed &seed, float time_since_defecation, int &no_seeds_dispersed, ResourceGrid* resource_grid,
		int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions, int& no_competitions_with_older_trees,
		int& no_germination_attempts
	) {
		pair<float, float> seed_deposition_location;
		if (moving) {
			seed_deposition_location = get_backwards_traced_location(time_since_defecation);
		}
		else {
			seed_deposition_location = position;
		}

		int prev_seed_dispersed_number = no_seeds_dispersed;
		seed.deposition_location = seed_deposition_location;
		bool germination = seed.germinate_if_location_is_viable(
			resource_grid->state, no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees, no_germination_attempts
		);
		no_seeds_dispersed += germination;
	}
	map<string, float> traits;
	unordered_map<int, vector<int>> stomach_content; // { defecation_time: crop_ids }
	pair<float, float> position;
	pair<float, float> trajectory;
	string species;
	GammaProbModel gut_passage_time_distribution;
	GammaProbModel rest_time_distribution;
	float curtime = 0;
	float recipr_speed = 0;
	float travel_time = 0;
	int last_tree_visited = -1;
	int iteration = -1;
	int animal_group_size = 20;
	int verbosity = 0;
	int total_no_seeds_consumed = 0;
	bool moving = false;
};


class Animals {
public:
	Animals() = default;
	Animals(State* state, map<string, map<string, float>> _animal_kernel_params) {
		animal_kernel_params = _animal_kernel_params;
		total_no_animals = round(animal_kernel_params["population"]["density"] * (state->grid.area / 1e6));
		no_iterations = animal_kernel_params["population"]["no_iterations"];
		animal_kernel_params.erase("population");
	}
	void print_animal_kernel_params(map<string, map<string, float>>& _animal_kernel_params) {
		for (auto& [species, params] : animal_kernel_params) {
			printf("species: %s\n", species.c_str());
			for (auto& [param, value] : params) {
				printf("%s: %f\n", param.c_str(), value);
			}
		}
	}
	void initialize_population() {
		float cumulative_population_fraction = 0;
		for (auto& [species, params] : animal_kernel_params) {
			int popsize = round(total_no_animals * params["population_fraction"]);
			cumulative_population_fraction += params["population_fraction"];
			vector<Animal> species_population;
			create_species_population(popsize, species_population, params, species);
			total_animal_population[species] = species_population;
		}
		if (
			(cumulative_population_fraction < 0.95f) || (cumulative_population_fraction > 1.05f)
			) {
			printf("\nERROR: Cumulative population fraction %f outside of range (0.95 - 1.05).\n", cumulative_population_fraction);
			throw;
		}
	}
	void create_species_population(
		int popsize, vector<Animal>& species_population, map<string, float> traits, string species
	) {
		for (int i = 0; i < popsize; i++) {
			Animal animal(traits, species);
			species_population.push_back(animal);
		}
	}
	void place(State* state) {
		for (auto& [species, species_population] : total_animal_population) {
			for (auto& animal : species_population) {
				animal.reset();
				animal.position = state->grid.get_random_real_position();
			}
		}
	}
	void disperse(int& no_seeds_dispersed, int no_seeds_to_disperse, State* state, ResourceGrid* resource_grid,
		float& fraction_time_spent_moving, int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions,
		int& no_competitions_with_older_trees, int& no_germination_attempts
	) {
		place(state);
		resource_grid->reset_color_arrays();
		int iteration = 0;

		int no_seeds_eaten = 0;
		int no_seeds_defecated = 0;
		float time_spent_resting = 0;
		float time_spent_moving = 0;
		while (no_seeds_defecated < no_seeds_to_disperse) {
			int prev_no_seeds_dispersed = no_seeds_dispersed;
			int prev_no_seeds_defecated = no_seeds_defecated;
			float average_gut_passage_time = 0;
			int cur_no_seeds_defecated = 0;
			float distance_travelled = 0;
			float average_gut_time = 0;
			int cur_no_seeds_eaten = 0;
			float average_curtime = 0;
			for (auto& [species, species_population] : total_animal_population) {
				if (iteration == 0) resource_grid->update_cover_probabilities(species, species_population[0].traits);
				resource_grid->update_fruit_probabilities(species, species_population[0].traits);
				for (auto& animal : species_population) {
					float _average_gut_time = 0;
					float _average_gut_passage_time = 0;
					animal.update(
						no_seeds_dispersed, iteration, state, resource_grid, cur_no_seeds_eaten, time_spent_resting,
						distance_travelled, _average_gut_passage_time, _average_gut_time, no_seeds_defecated, time_spent_moving,
						no_seedlings_dead_due_to_shade, no_seedling_competitions, no_competitions_with_older_trees,
						no_germination_attempts
					);
					average_curtime += animal.curtime;
					average_gut_time += _average_gut_time;
					average_gut_passage_time += _average_gut_passage_time;
				}
			}
			if (iteration % 60 == 0 && iteration > 0) {
				printf("-- Animal dispersal progress: %f %%\n", (float)no_seeds_defecated * 100.0f / (float)no_seeds_to_disperse);
			}
			no_seeds_eaten += cur_no_seeds_eaten;
			iteration++;
		}
		printf("-- Number of iterations spent dispersing fruits: %d\n", iteration);
		fraction_time_spent_moving = time_spent_moving / (time_spent_moving + time_spent_resting);
	}
	int popsize() {
		int total_popsize = 0;
		for (auto& [species, species_population] : total_animal_population) {
			total_popsize += species_population.size();
		}
		return total_popsize;
	}
	map<string, vector<Animal>> total_animal_population;
	map<string, map<string, float>> animal_kernel_params;
	float total_no_animals = 0;
	int verbose = 0;
	int no_iterations = 0;
};
