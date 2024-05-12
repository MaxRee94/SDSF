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
		if (cell->is_hospitable(pair<float, int>(strategy.growth_rate, strategy.id))) {
			//printf("Germinating seed at location (%f, %f)\n", deposition_location.first, deposition_location.second);
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
	void germinate(Cell* cell, State* state) {
		cell->update_largest_stem(strategy.growth_rate, strategy.id);
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
	bool get(int idx, Fruit &fruit) {
		if (idx >= fruits.size()) return false;
		fruit = fruits[idx];
		return true;
	}
	void remove(int idx) {
		if (idx >= fruits.size()) return;
		fruits.erase(fruits.begin() + idx);
	}
	int size() {
		return fruits.size();
	}
	vector<Fruit> fruits;
};

