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
	bool germinate_if_location_is_viable(
		State* state, int& no_seedlings_dead_due_to_shade, int& no_seedling_competitions, int& no_competitions_with_older_trees,
		int& no_germination_attempts, int& no_cases_seedling_competition_and_shading, int& no_cases_oldstem_competition_and_shading
	) {
		Cell* cell = state->grid.get_cell_at_position(deposition_location);
		no_germination_attempts++;
		if (cell->is_hospitable(pair<float, int>(strategy.seedling_dbh, strategy.id), no_seedlings_dead_due_to_shade, no_seedling_competitions,
			no_competitions_with_older_trees, no_cases_seedling_competition_and_shading, no_cases_oldstem_competition_and_shading)
		) {
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
	Fruit(Strategy &_strategy, int _id) : strategy(_strategy), id(_id) {};
	void get_seeds(vector<Seed>& seeds) {
		for (int i = 0; i < strategy.no_seeds_per_diaspore; i++) {
			Seed seed(strategy);
			seeds.push_back(seed);
		}
	}
	Strategy strategy;
	int id = -1;
};


class Fruits {
public:
	Fruits() = default;
	void add_fruits(Crop* crop, float abundance) {
		//printf("Fruit abundance: %d, no diaspora to disperse: %d\n", crop->fruit_abundance, crop->no_diaspora);
		fruits[crop->id] = pair<Fruit, int>(Fruit(crop->strategy, crop->id), crop->fruit_abundance);
		//printf("fruits with crop id %d added, fruit id: %i\n", crop->id, fruits[crop->id].first.id);
		if (!help::is_in(&fruit_types, crop->id)) fruit_types.push_back(crop->id);
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
	bool get(Fruit &fruit, int tree_id = -1) {
		if (_no_fruits == 0) return false;

		// Select fruit type randomly
		int fruit_type;
		if (tree_id != -1) {
			fruit_type = tree_id;
		}
		else fruit_type = fruit_types[help::get_rand_int(0, fruit_types.size() - 1)];
		
		// Extract fruit
		fruit = fruits[fruit_type].first;
		if (fruit.id == -1) {
			//printf("---------------------------------------------- Error: fruit id is -1\n\n\n\n\n\n");
			return false;
			for (int i = 0; i < fruit_types.size(); i++) {
				printf("Fruit type: %d\n", fruit_types[i]);
			}
		}
		decrement_one(fruit_type);

		return true;
	}
	void decrement_one(int fruit_type) {
		fruits[fruit_type].second--;
		_no_fruits--;
		if (fruits[fruit_type].second == 0) {
			fruits.erase(fruit_type);
			fruit_types.erase(find(fruit_types.begin(), fruit_types.end(), fruit_type));
			//printf("\n\n\n\n\n\n Erasing fruit type %i.\n", fruit_type);
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

