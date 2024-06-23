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
	void update(int &no_seeds_dispersed, int _iteration, State* state, ResourceGrid* resource_grid) {
		float begin_time = curtime;
		iteration = _iteration;

		move(resource_grid);
		digest(resource_grid, no_seeds_dispersed);
		rest();
		eat(resource_grid, begin_time);
		digest(resource_grid, no_seeds_dispersed);
	}
	void rest() {
		moving = false;
		float rest_time = rest_time_distribution.get_gamma_sample() * 60.0f; // Multiply by 60 to convert minutes to seconds.
		curtime += rest_time;
	}
	float get_biomass_appetite() {
		return 1000; // We assume that the group of animals will eat 1 kg of fruit in a single session.
	}
	void eat(ResourceGrid* resource_grid, float begin_time) {
		float biomass_appetite = get_biomass_appetite();
		vector<int> consumed_seed_ids = {};
		while (biomass_appetite > 0) {
			Fruit fruit;
			bool fruit_available = resource_grid->extract_fruit(position, fruit);
			if (!fruit_available) break;
			eat_fruit(fruit, consumed_seed_ids);
			biomass_appetite -= fruit.strategy.diaspore_mass;
		}
		resource_grid->update_fruit_abundance(position, species, traits);
		set_ingestion_times(begin_time, consumed_seed_ids);
	}
	void set_ingestion_times(float begin_time, vector<int> &consumed_seed_ids) {
		float time_stepsize = (curtime - begin_time) / (float)consumed_seed_ids.size();
		for (int i = consumed_seed_ids.size() - 1; i >= 0; i--) {
			stomach_content[consumed_seed_ids[i]].second.second = curtime - (float)(i + 1) * time_stepsize;
		}
	}
	void eat_fruit(Fruit &fruit, vector<int> &consumed_seed_ids) {
		vector<Seed> seeds;
		fruit.get_seeds(seeds);
		for (auto& seed : seeds) {
			float gut_passage_time = gut_passage_time_distribution.get_gamma_sample() * 60.0f; // Multiply by 60 to convert minutes to seconds.
			pair<float, float> times(gut_passage_time, curtime);
			stomach_content[total_no_seeds_consumed] = pair<Seed, pair<float, float>>(seed, times);
			consumed_seed_ids.push_back(total_no_seeds_consumed);
			total_no_seeds_consumed++;
		}
	}
	void digest(ResourceGrid* resource_grid, int& no_seeds_dispersed) {
		vector<int> defecation_schedule;
		for (auto& [seed_id, seed_plus_times] : stomach_content) {
			auto [gut_passage_time, ingestion_time] = seed_plus_times.second;
			if (gut_passage_time + ingestion_time <= curtime) {
				if (iteration < 5) continue; // Ignore the first 5 iterations (after Morales et al 2013)
				float defecation_time = gut_passage_time + ingestion_time;
				float time_since_defecation = curtime - defecation_time;
				Seed &to_defecate = seed_plus_times.first;
				defecate(to_defecate, time_since_defecation, no_seeds_dispersed, resource_grid);
			}
		}
		for (auto& seed_id : defecation_schedule) {
			stomach_content.erase(seed_id);
		}
	}
	pair<float, float> select_destination(ResourceGrid* resource_grid) {
		resource_grid->update_coarse_probability_distribution(species, traits, position);
		ResourceCell* cell = resource_grid->select_cell(species, traits, position);
		pair<float, float> destination;
		//ResourceCell* cell = resource_grid->get_random_resource_cell();
		resource_grid->get_random_location_within_cell(cell, destination);
		return destination;
	}
	void move(ResourceGrid* resource_grid) {
		moving = true;
		pair<float, float> destination = select_destination(resource_grid);
		float distance;
		trajectory = resource_grid->get_shortest_trajectory(position, destination, distance);
		position = destination;
		travel_time = distance * recipr_speed;
		curtime += travel_time;
	}
	pair<float, float> get_backwards_traced_location(float time_since_defecation) {
		return position - (time_since_defecation / travel_time) * trajectory;
	}
	void defecate(Seed &seed, float time_since_defecation, int &no_seeds_dispersed, ResourceGrid* resource_grid) {
		pair<float, float> seed_deposition_location;
		if (moving) {
			seed_deposition_location = get_backwards_traced_location(time_since_defecation);
		}
		else {
			seed_deposition_location = position;
		}
		
		for (int i = 0; i < 2; i++) { // Resample once to speed up dispersal.
			resource_grid->get_random_stategrid_location(seed_deposition_location);
			seed.deposition_location = seed_deposition_location;
			bool germination = seed.germinate_if_location_is_viable(resource_grid->state);
			no_seeds_dispersed += germination;
		}
	}
	map<string, float> traits;
	map<int, pair<Seed, pair<float, float>>> stomach_content;
	pair<float, float> position;
	pair<float, float> trajectory;
	string species;
	GammaProbModel gut_passage_time_distribution;
	GammaProbModel rest_time_distribution;
	float curtime = 0;
	float recipr_speed = 0;
	float travel_time = 0;
	int iteration = -1;
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
		for (auto& [species, params] : animal_kernel_params) {
			int popsize = round(total_no_animals * params["population_fraction"]);
			vector<Animal> species_population;
			create_species_population(popsize, species_population, params, species);
			total_animal_population[species] = species_population;
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
	void disperse(int& no_seeds_dispersed, int no_seeds_to_disperse, State* state, ResourceGrid* resource_grid) {
		place(state);
		resource_grid->reset_color_arrays();
		int iteration = 0;
		while (no_seeds_dispersed < no_seeds_to_disperse) {
			for (auto& [species, species_population] : total_animal_population) {
				for (auto& animal : species_population) {
					resource_grid->update_cover_and_fruit_probabilities(species, animal.traits);
					animal.update(no_seeds_dispersed, iteration, state, resource_grid);
				}
			}
			iteration++;
		}
		printf("-- Number of iterations spent dispersing fruits: %d\n", iteration);
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
