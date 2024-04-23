#pragma once
#include "animal_resources.h"
#include "grid_agent.forward.h"


class Animal {
public:
	Animal() = default;
	Animal(map<string, float> _traits, string _species, State* _state, ResourceGrid* _resource_grid) {
		resource_grid = _resource_grid;
		traits = _traits;
		traits["a_c_recipr"] = 1.0f / traits["a_c"];
		traits["a_f_recipr"] = 1.0f / traits["a_f"];
		species = _species;
		state = _state;
		recipr_speed = 1.0f / traits["speed"];
		init_prob_distributions();
	}
	void init_prob_distributions() {
		gut_passage_time_distribution = GammaProbModel(traits["gut_passage_time_shape"], traits["gut_passage_time_scale"]);
		rest_time_distribution = GammaProbModel(traits["rest_time_shape"], traits["rest_time_scale"]);
	}
	void update(vector<Seed> &seeds, int _iteration) {
		float begin_time = curtime;
		iteration = _iteration;
		rest();
		eat(begin_time);
		digest(seeds);
		move();
		digest(seeds);
	}
	void rest() {
		moving = false;
		float rest_time = rest_time_distribution.get_gamma_sample() * 60.0f; // Multiply by 60 to convert minutes to seconds.
		curtime += rest_time;
	}
	float get_biomass_appetite() {
		return 0.2f; // TEMP: replace with value drawn from species-dependent poisson distribution.
	}
	void eat(float begin_time) {
		float biomass_appetite = get_biomass_appetite();
		int no_fruits_consumed = 0;
		while (biomass_appetite > 0) {
			Fruit fruit;
			bool fruit_available = resource_grid->extract_fruit(position, fruit);
			if (!fruit_available) break;
			eat_fruit(fruit);
			biomass_appetite -= fruit.strategy.diaspore_mass;
			no_fruits_consumed++;
		}
		resource_grid->update_fruit_abundance(position, species, traits);
		set_ingestion_times(begin_time, no_fruits_consumed);
	}
	void set_ingestion_times(float begin_time, float no_fruits_consumed) {
		float time_stepsize = (curtime - begin_time) / no_fruits_consumed;
		for (int i = 0; i < stomach_content.size(); i++) {
			stomach_content[stomach_content.size() - 1 - i].second.second = begin_time + (float)(i + 1) * time_stepsize;
		}
	}
	void eat_fruit(Fruit &fruit) {
		vector<Seed> seeds;
		fruit.get_seeds(seeds);
		for (auto& seed : seeds) {
			float gut_passage_time = gut_passage_time_distribution.get_gamma_sample();
			pair<float, float> times(gut_passage_time, curtime);
			stomach_content.push_back(pair<Seed, pair<float, float>>(seed, times));
		}
	}
	void digest(vector<Seed> &seeds_to_disperse) {
		vector<pair<Seed, pair<float, float>>> _stomach_content = {};
		for (int i = 0; i < stomach_content.size(); i++) {
			auto [gut_passage_time, ingestion_time] = stomach_content[i].second;
			if (gut_passage_time + ingestion_time <= curtime) {
				float defecation_time = gut_passage_time + ingestion_time;
				float time_since_defecation = curtime - defecation_time;
				if (iteration > 4) // Ignore the first 5 iterations (after Morales et al 2013)
					defecate(stomach_content[i].first, time_since_defecation, seeds_to_disperse);
			}
			else {
				_stomach_content.push_back(stomach_content[i]);
			}
		}
		stomach_content = _stomach_content;
	}
	pair<float, float> select_destination() {
		resource_grid->update_probability_distribution(species, traits, position);
		ResourceCell* cell = resource_grid->select_cell();
		pair<float, float> destination;
		resource_grid->get_random_location_within_cell(cell, destination);
		return destination;
	}
	void move() {
		moving = true;
		pair<float, float> destination = select_destination();
		float distance = help::get_dist(position, destination);
		trajectory = destination - position;
		position = destination;
		travel_time = distance * recipr_speed;
		curtime += travel_time;
	}
	pair<float, float> get_backwards_traced_location(float time_since_defecation) {
		return position - (time_since_defecation / travel_time) * trajectory;
	}
	void defecate(Seed &seed, float time_since_defecation, vector<Seed>& seeds_to_disperse) {
		pair<float, float> seed_deposition_location;
		if (moving) {
			seed_deposition_location = get_backwards_traced_location(time_since_defecation);
		}
		else {
			seed_deposition_location = position;
		}
		resource_grid->get_random_location_within_cell(seed_deposition_location);
		seed.deposition_location = seed_deposition_location;
		seeds_to_disperse.push_back(seed);
	}
	ResourceGrid* resource_grid = 0;
	map<string, float> traits;
	vector<pair<Seed, pair<float, float>>> stomach_content;
	pair<float, float> position;
	pair<float, float> direction;
	pair<float, float> trajectory;
	string species;
	GammaProbModel gut_passage_time_distribution;
	GammaProbModel rest_time_distribution;
	State* state = 0;
	float curtime = 0;
	float recipr_speed = 0;
	float travel_time = 0;
	int iteration = -1;
	bool moving = false;
};


class Animals {
public:
	Animals() = default;
	Animals(State* _state, map<string, map<string, float>> _animal_kernel_params) {
		state = _state;
		tree_population = &state->population;
		grid = &state->grid;
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
	void initialize_population(ResourceGrid* resource_grid) {
		for (auto& [species, params] : animal_kernel_params) {
			int popsize = round(total_no_animals * params["population_fraction"]);
			vector<Animal> species_population;
			create_species_population(popsize, species_population, params, species, resource_grid);
			total_animal_population[species] = species_population;
		}
	}
	void create_species_population(
		int popsize, vector<Animal>& species_population, map<string, float> traits, string species,
		ResourceGrid* resource_grid
	) {
		for (int i = 0; i < popsize; i++) {
			Animal animal(traits, species, state, resource_grid);
			species_population.push_back(animal);
		}
	}
	void place() {
		for (auto& [species, species_population] : total_animal_population) {
			for (auto& animal : species_population) {
				animal.position = grid->get_random_real_position();
			}
		}
	}
	void disperse(vector<Seed>& seeds) {
		place();
		for (int i = 0; i < no_iterations; i++) {
			for (auto& [species, species_population] : total_animal_population) {
				for (auto& animal : species_population) {
					if (i == 0) animal.resource_grid->update_cover_and_fruit_probabilities(species, animal.traits);
					animal.update(seeds, i);
				}
			}
		}
	}
	int popsize() {
		int total_popsize = 0;
		for (auto& [species, species_population] : total_animal_population) {
			total_popsize += species_population.size();
		}
		return total_popsize;
	}
	Population* tree_population = 0;
	map<string, vector<Animal>> total_animal_population;
	map<string, map<string, float>> animal_kernel_params;
	Grid* grid = 0;
	State* state = 0;
	float total_no_animals = 0;
	int verbose = 0;
	int no_iterations = 0;
};
