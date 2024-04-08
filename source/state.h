#pragma once
#include <iostream>
#include "agents.h"
#include "grid.h"


class State {
public:
	State() {
		grid = Grid();
		population = Population();
		init_neighbor_offset();
	}
	State(
		int gridsize, float cellsize, float max_radius, float radius_q1, float radius_q2,
		float _seed_bearing_threshold, float _mass_budget_factor, float _seed_mass
	) {
		seed_bearing_threshold = _seed_bearing_threshold;
		mass_budget_factor = _mass_budget_factor / help::cubed(max_radius); // Normalize by maximum radius
		seed_mass = _seed_mass;
		population = Population(max_radius, cellsize, radius_q1, radius_q2, mass_budget_factor);
		grid = Grid(gridsize, cellsize);
		init_neighbor_offset();
	}
	void init_neighbor_offset() {
		neighbor_offsets = new pair<int, int>[8];
		int q = 0;
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if (i == 0 && j == 0) continue;
				neighbor_offsets[q] = pair<int, int>(i, j);
				q++;
			}
		}
	}
	void repopulate_grid(int verbosity) {
		if (verbosity == 2) cout << "Repopulating grid..." << endl;
		grid.reset();
		for (auto& [id, tree] : population.members) {
			grid.populate_tree_domain(&tree);
		}
		if (verbosity == 2) cout << "Repopulated grid." << endl;
	}
	float get_dist(pair<float, float> a, pair<float, float> b, bool verbose = false) {
		vector<float> dists = {help::get_dist(a, b)};
		for (int i = 0; i < 8; i++) {
			float dist = help::get_dist(
				neighbor_offsets[i] * grid.width * grid.cellsize, a
			);
			dists.push_back(dist);
		}
		if (verbose) printf("normal dist: %f, periodic dist: %f \n", dists[0], help::get_min(&dists));
		return help::get_min(&dists);
	}
	vector<Tree*> get_tree_neighbors(pair<float, float> baseposition, float search_radius, int base_id = -1) {
		vector<Tree*> neighbors;
		for (auto& [id, candidate] : population.members) {
			if (base_id != -1 && id == base_id) {
				continue;
			}
			if (get_dist(candidate.position, baseposition) < search_radius) {
				neighbors.push_back(&candidate);
			}
		}
		return neighbors;
	}
	vector<Tree*> get_tree_neighbors(Tree* base) {
		vector<Tree*> neighbors;
		float search_radius = (round(base->radius * 1.1) + population.max_radius);
		return get_tree_neighbors(base->position, search_radius, base->id);
	}
	void set_tree_cover(float _tree_cover) {
		help::init_RNG();
		grid.reset();

		while (grid.get_tree_cover() < _tree_cover) {
			pair<float, float> position = grid.get_random_real_position();
			
			// TEMP: all strategies are initially the same. TODO: Add option to make these variable according to prob. distribution
			Strategy strategy(seed_mass, 1.0); // TEMP: Seed density of 1 for now
			Tree* tree = population.add(position, strategy, -2);
			grid.populate_tree_domain(tree);

			if (population.size() % 1000 == 0) {
				printf("Current tree cover: %f, current population width: %i\n", grid.get_tree_cover(), population.size());
			}
			continue;
		}
		printf("Final tree cover: %f\n", grid.get_tree_cover());

		// Count no small trees
		int no_small = 0;
		for (auto& [id, tree] : population.members) {
			if (tree.radius < population.max_radius / 2.0) {
				no_small++;
			}
		}
		printf("- Fraction small trees: %f \n", (float)no_small / (float)population.size());

		repopulate_grid(0);
	}
	Grid grid;
	Population population;
	pair<int, int>* neighbor_offsets = 0;
	float seed_bearing_threshold = 0;
	float mass_budget_factor = 0;
	float seed_mass = 0;
	float saturation_threshold = 0.333f;
};

