#pragma once
#include <iostream>
#include "agents.h"
#include "grid.h"


class State {
public:
	State() {
		grid = Grid();
		population = Population();
	}
	State(
		int gridsize, float cell_width, float max_dbh, float dbh_q1, float dbh_q2,
		float seed_bearing_threshold, float _saturation_threshold,
		map<string, map<string, float>>& strategy_distribution_params, float mutation_rate
	) {
		saturation_threshold = _saturation_threshold;
		population = Population(
			max_dbh, cell_width, dbh_q1, dbh_q2, strategy_distribution_params, mutation_rate,
			seed_bearing_threshold
		);
		grid = Grid(gridsize, cell_width);
	}
	bool check_grid_for_tree_presence(int tree_id, int verbose = 0) {
		bool presence = false;
		for (int i = 0; i < grid.no_cells; i++) {
			Cell cell = grid.distribution[i];
			if (verbose > 1 && i % 100000 == 0) printf("Checking cell %i for tree presence... \n", i);
			if (cell.tree_is_present(tree_id)) {
				if (verbose > 0) printf("pos: (%i, %i)\n", cell.pos.first, cell.pos.second);
				presence = true;
			}
		}
		return presence;
	}
	void repopulate_grid(int verbosity) {
		if (verbosity == 2) cout << "Repopulating grid..." << endl;
		grid.reset();
		int i = 0;
		bool success;
		for (auto& [id, tree] : population.members) {
			//if (verbosity == 1) printf("Tree id (beginning): %i \n", tree.id);
			//if (i % (population.size()/1000) == 0) printf("i: %i / %i \n", i, population.size());
			if (id == -1 || tree.id == -1) {
				printf("Removing tree %i\n", tree.id);
				population.remove(id); // HOTFIX: Sometimes trees are not initialized properly and need to be removed.
				continue;
			}
			success = grid.populate_tree_domain(&tree);
			if (!success) {
				population.remove(tree.id);
				printf("\n------------- Restarting grid repopulation because tree %i failed to have its domain populated. --------\n", tree.id);
				repopulate_grid(verbosity);
			}
			/*if (verbosity > 0 && !check_grid_for_tree_presence(tree.id)) {
				printf("Repopulation failure: Tree %i with radius %f and position (%f, %f) is not present in the grid. \n", tree.id, tree.radius, tree.position.first, tree.position.second);
				grid.populate_tree_domain(&tree, 1);
			}*/
			i++;
		}
		grid.update_grass_LAIs();
		if (verbosity == 2) cout << "Repopulated grid." << endl;
	}
	float compute_shade_on_individual_tree(Tree* tree) {
		float LAI_shade = 0;
		float no_cells = 0;
		TreeDomainIterator it(grid.cell_width, tree);
		while (it.next()) {
			if (tree->radius_spans(it.real_cell_position)) {
				Cell* cell = grid.get_cell_at_position(it.gb_cell_position);
				float _LAI_shade = cell->get_shading_on_tree(tree, &population);
				LAI_shade += _LAI_shade;
				if (_LAI_shade < tree->LAI) printf(" -- Shade is less than tree LAI. Shade: %f, tree LAI: %f\n", _LAI_shade, tree->LAI);
				no_cells += 1;
			}
		}
		if (no_cells == 0) {
			int center_idx = grid.get_capped_center_idx(it.tree_center_gb);
			return grid.distribution[center_idx].get_shading_on_tree(tree, &population);
		}
		return LAI_shade / no_cells;	// We obtain mean LAI at- or above the height of the lowest branch by dividing by the tree's crown area.
										// We use this as a measure of shading on the tree.
	}
	float get_dist(pair<float, float> a, pair<float, float> b, bool verbose = false) {
		vector<float> dists = { help::get_manhattan_dist(a, b) };
		float min_dist = dists[0];
		int min_idx = 0;
		for (int i = 0; i < 8; i++) {
			float dist = help::get_dist(
				grid.neighbor_offsets[i] * grid.width_r + b, a
			);
			if (dist < min_dist && dist > 0) {
				min_dist = dist;
				min_idx = i + 1;
			}
			dists.push_back(dist);
		}
		float dist;
		if (min_idx != 0) {
			dist = help::get_dist(grid.neighbor_offsets[min_idx - 1] * grid.width_r + b, a);
		}
		else dist = help::get_dist(a, b);
		return dist;
	}
	vector<Tree*> get_tree_neighbors(pair<float, float> baseposition, float search_radius, int& closest_idx, int base_id = -1, bool stop_on_1 = false) {
		vector<Tree*> neighbors;
		float min_dist = INFINITY;
		int idx = 0;
		for (auto& [id, candidate] : population.members) {
			if (base_id != -1 && id == base_id) {
				continue;
			}
			float dist = get_dist(candidate.position, baseposition);
			if (dist < search_radius) {
				neighbors.push_back(&candidate);
				if (dist < min_dist) {
					min_dist = dist;
					closest_idx = idx;
				}
				if (stop_on_1) return neighbors;
				idx++;
			}
		}
		return neighbors;
	}
	vector<Tree*> get_tree_neighbors(Tree* base) {
		vector<Tree*> neighbors;
		float search_radius = (round(base->radius * 1.1) + population.max_dbh);
		int dummy;
		return get_tree_neighbors(base->position, search_radius, dummy, base->id);
	}
	vector<Tree*> get_tree_neighbors(pair<float, float> baseposition, float search_radius) {
		int dummy;
		return get_tree_neighbors(baseposition, search_radius, dummy);
	}
	void thin_crowds(bool recompute_shade = false) {
		PairSet population_sorted_by_height;
		population.sort_by_trait("height", population_sorted_by_height);
		for (auto& [id, height] : population_sorted_by_height) {
			Tree* tree = population.get(id);
			if (tree->life_phase == 2) continue; // We assume that mature trees are no longer vulnerable to death by shading.
			float shade;
			if (recompute_shade) shade = compute_shade_on_individual_tree(tree);
			else shade = tree->shade;
			float death_probability = max(0, (shade - 5.0f)); // We assume that trees always survive if shade <= 5.0, and always die if shade > 6.0.
			death_probability *= death_probability;
			if (help::get_rand_float(0, 1) < death_probability) {
				remove_tree(tree);
			}
		}
	}
	void remove_tree(Tree* tree) {
		grid.kill_tree_domain(tree, false);
		population.remove(tree);
	}
	void set_cover_from_image(float* image, int img_width, int img_height) {
		float pixel_size = grid.width_r / (float)img_width;
		int no_gridcells_along_x_per_pixel = round(grid.width_r / img_width);
		pair<int, int> bb = pair<int, int>(no_gridcells_along_x_per_pixel, no_gridcells_along_x_per_pixel);
		int iteration_budget = 5;
		int j = 0;

		// Create probability model
		DiscreteProbabilityModel probmodel = DiscreteProbabilityModel(grid.no_cells);
		float integral_image_cover;
		probmodel.set_probabilities(image, integral_image_cover);
		probmodel.normalize(integral_image_cover);
		probmodel.build_cdf();
		integral_image_cover /= (float)(img_width * img_height);
		printf("Image cover: %f\n", integral_image_cover);

		// Set tree cover
		int no_crowd_thinning_runs = 0;
		int max_no_crowd_thinning_runs = 10;
		while (grid.get_tree_cover() < integral_image_cover) {
			int idx = probmodel.sample();
			Cell &cell = grid.distribution[idx];
			pair<float, float> position = grid.get_real_cell_position(&cell);
			Tree* tree = population.add(position);
			if (!cell.is_hospitable(pair<float, int>(tree->radius, tree->id))) {
				population.remove(tree->id);
				continue;
			}
			grid.populate_tree_domain(tree);
			grid.update_grass_LAIs_for_individual_tree(tree);

			// Thin crowds when close to target tree cover
			if (grid.tree_cover > 0.95f * integral_image_cover && no_crowd_thinning_runs < max_no_crowd_thinning_runs) {
				no_crowd_thinning_runs++;
				thin_crowds(true);
			}

			if (population.size() % 10000 == 0) printf("Current tree cover: %f, current population size: %i\n", grid.get_tree_cover(), population.size());
		}
		probmodel.free();
		repopulate_grid(0);
		printf("Finished setting tree cover from image.\n");
	}
	void set_tree_cover(float _tree_cover) {
		help::init_RNG();
		grid.reset();

		int wind_trees = 0;
		int animal_trees = 0;
		int no_crowd_thinning_runs = 0;
		int max_no_crowd_thinning_runs = 10;
		while (grid.get_tree_cover() < _tree_cover) {
			pair<float, float> position = grid.get_random_real_position();
			Cell* cell = grid.get_cell_at_position(position);
			Tree* tree = population.add(position);
			if (!cell->is_hospitable(pair<float, int>(tree->radius, tree->id))) {
				population.remove(tree->id);
				continue;
			}
			if (tree->id == -1) population.remove(population.no_created_trees - 1); // HOTFIX: Sometimes trees are not initialized properly and need to be removed.
			grid.populate_tree_domain(tree);
			if (population.get_crop(tree->id)->strategy.vector == "wind") {
				wind_trees++;
			}
			else if (population.get_crop(tree->id)->strategy.vector == "animal") {
				animal_trees++;
			}

			// Thin crowds when close to target tree cover
			if (grid.tree_cover > 0.95f * _tree_cover && no_crowd_thinning_runs < max_no_crowd_thinning_runs) {
				no_crowd_thinning_runs++;
				thin_crowds(true);
			}

			if (population.size() % 1000 == 0) {
				printf("Current tree cover: %f, current population size: %i\n", grid.get_tree_cover(), population.size());
			}
			continue;
		}
		printf("Final tree cover: %f\n", grid.tree_cover);
		printf("Wind trees: %i, Animal trees: %i\n", wind_trees, animal_trees);
		printf("First tree's strategy: \n");
		population.get_crop(10)->strategy.print();

		// Count no small trees
		int no_small = 0;
		for (auto& [id, tree] : population.members) {
			if (tree.dbh < population.max_dbh / 2.0) {
				no_small++;
			}
		}
		printf("- Fraction small trees: %f \n", (float)no_small / (float)population.size());

		initial_tree_cover = grid.tree_cover;
		repopulate_grid(0);
	}
	void get_tree_positions(float* x, float* y) {
		int i = 0;
		for (auto& [id, tree] : population.members) {
			if (tree.life_phase < 2) continue; // Only consider trees in the mature phase (life_phase = 2)
			x[i] = tree.position.first;
			y[i] = tree.position.second;
			i++;
		}
	}
	void get_tree_sizes(float* tree_sizes) {
		int i = 0;
		for (auto& [id, tree] : population.members) {
			tree_sizes[i] = tree.dbh;
			i++;
		}
	}
	Grid grid;
	Population population;
	float initial_tree_cover = 0;
	float saturation_threshold = 0;
};

