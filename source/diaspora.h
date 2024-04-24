#pragma once
#include "state.h"



class Seed {
public:
	Seed() = default;
	Seed(Strategy &_strategy
	) :
		strategy(_strategy)
	{};
	Seed(Strategy &_strategy, pair<float, float> _deposition_location
	) :
		strategy(_strategy), deposition_location(_deposition_location)
	{};
	bool germinate_if_location_is_viable(State* state) {
		Cell* cell = state->grid.get_cell_at_position(deposition_location);
		if (!is_outcompeted(cell, state)) {
			germinate(cell, state);
			return true;
		}
		return false;
	}
	void set_deposition_location(pair<float, float> _deposition_location) {
		deposition_location = _deposition_location;
	}
	Strategy strategy;
	pair<float, float> deposition_location;
private:
	bool is_outcompeted(Cell* cell, State* state) {
		// Suppress germination of seeds in areas with increased competition
		float outcompete_likelihood = ((float)cell->trees.size() * state->saturation_threshold);
		if (help::get_rand_float(0, 1) < outcompete_likelihood) {
			return true;
		}
		return false;
	}
	void germinate(Cell* cell, State* state) {
		// TEMP: Arbitrary starting radius of 0.1. TODO: replace with 0 once proper growth curve is implemented.
		Tree* tree = state->population.add(deposition_location, &strategy, 0.1);
		cell->trees[tree->id] = tree->id;
	}
};


class Fruit {
public:
	Fruit() = default;
	Fruit(Strategy &_strategy) : strategy(_strategy) {};
	void get_seeds(vector<Seed>& seeds) {
		for (int i = 0; i < strategy.no_seeds_per_diaspore; i++) {
			Seed seed(strategy);
			seeds.push_back(seed);
		}
	}
	Strategy strategy;
};


class Fruits {
public:
	Fruits() = default;
	void push(Fruit &fruit) {
		fruits.push_back(fruit);
	}
	void add_fruits(Crop* crop) {
		for (int i = 0; i < crop->no_diaspora; i++) {
			Fruit fruit(crop->strategy);
			fruits.push_back(fruit);
		}
	}
	bool are_available() {
		return !fruits.empty();
	}
	int no_fruits() {
		return fruits.size();
	}
	void clear() {
		fruits.clear();
	}
	void get(int idx, Fruit &fruit) {
		fruit = fruits[idx];
	}
	void remove(int idx) {
		fruits.erase(fruits.begin() + idx);
	}
	vector<Fruit> fruits;
};

