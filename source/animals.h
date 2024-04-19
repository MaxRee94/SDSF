#pragma once
#include "state.h"
#include "diaspora.h"
#include "grid_agent.forward.h"


class Animal {
public:
	Animal() = default;
	Animal(map<string, float> _traits, string _species, State* _state) {
		traits = _traits;
		species = _species;
		state = _state;
		inv_speed = 1.0f / traits["speed"];
		init_prob_distributions();
	}
	void init_prob_distributions() {
		gut_passage_time_distribution = GammaProbModel(traits["gut_passage_time_shape"], traits["gut_passage_time_scale"]);
		rest_time_distribution = GammaProbModel(traits["rest_time_shape"], traits["rest_time_scale"]);
	}
	void update(vector<Seed> &seeds) {
		rest();
		digest(seeds);
		move();
		digest(seeds);
	}
	void rest() {
		moving = false;
		float rest_time = rest_time_distribution.get_gamma_sample() * 60.0f;
		curtime += rest_time;
	}
	void eat(Fruit &fruit) {
		vector<Seed> seeds;
		fruit.get_seeds(seeds);
		for (auto& seed : seeds) {
			float gut_passage_time = gut_passage_time_distribution.get_gamma_sample();
			pair<float, float> times(gut_passage_time, curtime);
			stomach_content.push_back(pair<Seed, pair<float, float>>(seed, times));
		}
	}
	void digest(vector<Seed> &seeds_to_disperse) {
		vector<pair<Seed, pair<float, float>>> _stomach_content;
		for (int i = 0; i < stomach_content.size(); i++) {
			auto [gut_passage_time, ingestion_time] = stomach_content[i].second;
			if (gut_passage_time + ingestion_time <= curtime) {
				float defecation_time = gut_passage_time + ingestion_time;
				float time_since_defecation = curtime - defecation_time;
				Seed seed = stomach_content[i].first;
				defecate(seed, time_since_defecation, seeds_to_disperse);
			}
			else {
				_stomach_content.push_back(stomach_content[i]);
			}
		}
		stomach_content = _stomach_content;
	}
	Cell* select_destination() {
		Cell* destination = state->grid.get_random_forest_cell(); // TEMP: replace with landscape-dependent cell selection.
		return destination;
	}
	void move() {
		moving = true;
		Cell* destination = select_destination(); // TEMP. TODO: use coarse grid cells as destinations. Ensure these are populated with fruit abundances beforehand (to get a map of fruit abundance).
		pair<float, float> destination_position = state->grid.get_real_cell_position(destination);
		float distance = help::get_dist(position, destination_position);
		trajectory = destination_position - position;
		travel_time = distance * inv_speed;
		curtime += travel_time;
	}
	pair<float, float> get_backwards_traced_location(float time_since_defecation) {
		return position - (time_since_defecation / travel_time) * trajectory;
	}
	void defecate(Seed &seed, float time_since_defecation, vector<Seed>& seeds_to_disperse) {
		pair<float, float> deposition_location;
		if (moving) {
			deposition_location = get_backwards_traced_location(time_since_defecation);
		}
		seed.deposition_location = deposition_location;
		seeds_to_disperse.push_back(seed);
	}
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
	float inv_speed = 0;
	float travel_time = 0;
	bool moving = false;
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
			Animal animal(traits, species, state);
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
