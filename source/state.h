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
	State(int gridsize, float cellsize, float max_radius, float radius_q1, float radius_q2) {
		population = Population(max_radius, cellsize, radius_q1, radius_q2);
		grid = Grid(gridsize, cellsize);
		init_neighbor_offset();
	}
	void init_neighbor_offset() {
		neighbor_offsets = new pair<int, int>[8];
		for (int i = -1; i < 2; i++) {
			for (int j = -1; j < 2; j++) {
				if (i == 0 && j == 0) continue;
				neighbor_offsets[i * 3 + j] = pair<int, int>(i, j);
			}
		}
	}
	void repopulate_grid(int verbosity) {
		if (verbosity == 2) cout << "Repopulating grid..." << endl;
		grid.reset();
		for (auto& tree : population.members) {
			grid.populate_tree_domain(&tree);
		}
		if (verbosity == 2) cout << "Repopulated grid." << endl;
	}
	float get_dist(pair<float, float> a, pair<float, float> b) {
		vector<float> dists = {help::get_dist(a, b)};
		for (int i = 0; i < 8; i++) {
			float dist = help::get_dist(
				neighbor_offsets[i] * grid.size * grid.cellsize, a
			);
			dists.push_back(dist);
		}
		return help::get_min(&dists);
	}
	vector<Tree*> get_tree_neighbors(pair<float, float> baseposition, float search_radius, Tree* base = 0) {
		vector<Tree*> neighbors;
		for (auto& candidate : population.members) {
			if (base != nullptr && candidate.id == base->id) {
				continue;
			}
			if (get_dist(candidate.position, baseposition) < search_radius) {
				neighbors.push_back(&candidate);
			}
		}
		printf("Neighbors for ptr %i: ", base);
		for (auto& neighbor : neighbors) printf("%i, ", neighbor);
		cout << endl;
		return neighbors;
	}
	vector<Tree*> get_tree_neighbors(Tree* base) {
		vector<Tree*> neighbors;
		float search_radius = (round(base->radius * 1.1) + population.max_radius);
		return get_tree_neighbors(base->position, search_radius, base);
	}
	void set_tree_cover(float _tree_cover) {
		help::init_RNG();
		grid.reset();

		while (grid.get_tree_cover() < _tree_cover) {
			pair<float, float> position = grid.get_random_real_position();
			Tree* tree = population.add(position, "maximum");
			grid.populate_tree_domain(tree);
			printf("radius: %f, ptr: %i \n", population.members.back().radius, tree);
			if (population.members.back().radius == -1) {
				cout << "problem reflected here too.\n";
			}

			if (population.size() % 1000 == 0) {
				printf("Current tree cover: %f, current population size: %i\n", grid.get_tree_cover(), population.size());
			}
			continue;
		}
		for (auto& tree : population.members) {
			printf("*** (tree cover func) - tree (ptr %i) radius: %f, no cells: %i \n", &tree, tree.radius, tree.cells.size());
		}
		printf("Final tree cover: %f\n", grid.get_tree_cover());

		// Count no small trees
		int no_small = 0;
		for (auto& tree : population.members) {
			if (tree.radius < population.max_radius / 2.0) {
				no_small++;
			}
		}
		printf("- Fraction small trees: %f \n", (float)no_small / (float)population.size());
		
		// HOTFIX: Repopulating the grid appears to prevent read access errors (where some trees cannot be read from memory). Cause unknown.
		repopulate_grid(2);
	}
	Grid grid;
	Population population;
	pair<int, int>* neighbor_offsets = 0;
};

