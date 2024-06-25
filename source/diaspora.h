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
		if (cell->is_hospitable(pair<float, int>(strategy.seedling_dbh, strategy.id))) {
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
		cell->set_stem(strategy.seedling_dbh, strategy.id);
		cell->seedling_present = true;
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
	void add_fruits(Crop* crop) {
		//printf("Fruit abundance: %d, no diaspora to disperse: %d\n", crop->fruit_abundance, crop->no_diaspora);
		fruits[fruit_types.size()] = pair<Fruit, int>(Fruit(crop->strategy), crop->fruit_abundance);
		fruit_types.push_back(fruits.size());
		_no_fruits += crop->fruit_abundance;
	}
	bool are_available() {
		return _no_fruits > 0;
	}
	void clear() {
		fruits.clear();
		fruit_types.clear();
		_no_fruits = 0;
	}
	bool get(Fruit &fruit, vector<int>& compute_times) {
		auto start = high_resolution_clock::now();
		if (_no_fruits == 0) return false;
		compute_times[12] += help::microseconds_elapsed_since(start);

		// Select fruit type randomly
		start = high_resolution_clock::now();
		int fruit_type = fruit_types[help::get_rand_int(0, fruit_types.size() - 1)];
		compute_times[9] += help::microseconds_elapsed_since(start);
		
		// Extract fruit
		start = high_resolution_clock::now();
		fruit = fruits[fruit_type].first;
		decrement_one(fruit_type);
		compute_times[10] += help::microseconds_elapsed_since(start);

		return true;
	}
	void decrement_one(int fruit_type) {
		fruits[fruit_type].second--;
		_no_fruits--;
		if (fruits[fruit_type].second == 0) {
			fruits.erase(fruit_type);
			fruit_types.erase(find(fruit_types.begin(), fruit_types.end(), fruit_type));
		}
	}
	int no_fruits() {
		return _no_fruits;
	}
	/*int size() {
		int _no_fruits = 0;
		for (auto [type, fruit] : fruits) {
			_no_fruits += fruit.second;
		}
		return _no_fruits;
	}*/
	map<int, pair<Fruit, int>> fruits;
	vector<int> fruit_types;
	Population* population;
private:
	int _no_fruits;
};

