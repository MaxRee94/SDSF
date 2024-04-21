#pragma once
#include "state.h"



class Seed {
public:
	Seed() = default;
	Seed(Strategy& _strategy, State* _state
	) :
		strategy(_strategy), state(_state)
	{};
	Seed(Strategy &_strategy, pair<float, float> _deposition_location, State* _state
	) :
		strategy(_strategy), deposition_location(_deposition_location), 
		state(_state)
	{};
	void germinate_if_location_is_viable() {
		Cell* cell = state->grid.get_cell_at_position(deposition_location);
		if (!is_outcompeted(cell)) {
			germinate(cell);
		}
	}
	void set_deposition_location(pair<float, float> _deposition_location) {
		deposition_location = _deposition_location;
	}
	Strategy strategy;
	pair<float, float> deposition_location;
	State* state;
private:
	bool is_outcompeted(Cell* cell) {
		// Suppress germination of seeds in areas with increased competition
		if (help::get_rand_float(0, 1) < ((float)cell->trees.size() * state->saturation_threshold)) {
			return true;
		}
	}
	void germinate(Cell* cell) {
		// TEMP: Arbitrary starting radius of 0.1. TODO: replace with 0 once proper growth curve is implemented.
		Tree* tree = state->population.add(deposition_location, &strategy, 0.1);
		cell->trees[tree->id] = tree->id;
	}
};


class Fruit {
public:
	Fruit() = default;
	Fruit(Strategy &_strategy, State* _state) : strategy(_strategy), state(_state) {};
	void get_seeds(vector<Seed>& seeds) {
		for (int i = 0; i < strategy.no_seeds_per_diaspore; i++) {
			Seed seed(strategy, state);
			seeds.push_back(seed);
		}
	}
	Strategy strategy;
	State* state = 0;
};


class Fruits {
public:
	Fruits() = default;
	Fruits(State* _state) {
		state = _state;
		_pop = &state->population;
		_grid = &state->grid;
	}
	void push(Fruit &fruit) {
		fruits.push_back(fruit);
	}
	void add_fruit(Crop* crop) {
		Fruit fruit(crop->strategy, state);
		fruits.push_back(fruit);
	}
	bool are_available() {
		return !fruits.empty();
	}
	Population* _pop = 0;
	Grid* _grid = 0;
	State* state = 0;
	vector<Fruit> fruits;
};

