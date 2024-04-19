#pragma once
#include "state.h"
#include "diaspora.h"
#include "grid_agent.forward.h"


class Animal {
public:
	Animal() = default;
	Animal(map<string, float> _traits, string _species) {
		traits = _traits;
		species = _species;
	}
	void eat(Fruit& fruit) {
		aborb_seeds(fruit);
	}
	void aborb_seeds(Fruit &fruit) {
		for (int i = 0; i < fruit.strategy.no_seeds_per_diaspore; i++) {
			vector<Seed> seeds; fruit.get_seeds(seeds);
			for (auto& seed : seeds) stomach_content.push_back(seed);
		}
	}
	map<string, float> traits;
	vector<Seed> stomach_content;
	pair<float, float> position;
	string species;
};


class Animals {
public:
	Animals() = default;
	Animals(State* _state, map<string, map<string, float>> &_animal_kernel_params) {
		state = _state;
		tree_population = &state->population;
		grid = &state->grid;
		total_no_animals = animal_kernel_params["population"]["density"] * state->grid.area;
		animal_kernel_params = _animal_kernel_params;
		animal_kernel_params.erase("population");
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
		int popsize, vector<Animal> &species_population, map<string, float> traits, string species
	) {
		for (int i = 0; i < popsize; i++) {
			Animal animal(traits, species);
			species_population.push_back(animal);
		}
	}
	void place() {
		for (auto& [species, species_population] : total_animal_population) {
			for (auto& animal : species_population) {
				animal.position = grid->get_random_real_position(); // TEMP: random position. TODO: Make dependent on fruit availability.
			}
		}
	}
	void disperse(Fruits& fruits, vector<Seed> seeds) {
	}
	Population* tree_population = 0;
	map<string, vector<Animal>> total_animal_population;
	map<string, map<string, float>> animal_kernel_params;
	Grid* grid = 0;
	State* state = 0;
	float total_no_animals = 0;
	int verbose = 0;
};
