#pragma once
#include "dynamics.h"

class Tests {
public:
	Tests() {};
	Tests(int _timestep, float _cell_width, float _self_ignition_factor, float _rainfall, float _seed_bearing_threshold,
		float _growth_rate_multiplier, float _unsuppressed_flammability, float _min_suppressed_flammability, float _max_suppressed_flammability,
		float _radius_suppr_flamm_min, float radius_range_suppr_flamm, float _max_dbh, float _saturation_threshold, float _fire_resistance_argmin,
		float _fire_resistance_argmax, float _fire_resistance_stretch, float _background_mortality, map<string, map<string, float>> _strategy_distribution_params,
		float _resource_grid_relative_size, float _mutation_rate, int _verbosity, int grid_size, float dbh_q1, float dbh_q2
	) {
		dynamics = Dynamics(_timestep, _cell_width, _self_ignition_factor, _rainfall, _seed_bearing_threshold,
			_growth_rate_multiplier, _unsuppressed_flammability, _min_suppressed_flammability, _max_suppressed_flammability,
			_radius_suppr_flamm_min, radius_range_suppr_flamm, _max_dbh, _saturation_threshold, _fire_resistance_argmin,
			_fire_resistance_argmax, _fire_resistance_stretch, _background_mortality, _strategy_distribution_params,
			_resource_grid_relative_size, _mutation_rate, _verbosity
		);
		dynamics.init_state(grid_size, dbh_q1, dbh_q2);
		printf("cell width: %f\n", dynamics.grid->cell_width);
		verbosity = _verbosity;
	}
	bool test_flammability(vector<string>& failed_tests) {
		bool success = true;
		
		if (!success) failed_tests.push_back("Flammability function");
		return success;
	}
	bool test_quadratic_solver(vector<string>& failed_tests) {
		bool success = false;

		// Cases
		vector<float> c1 = { 2.0, 0.5, 1.5, 0.5 };
		map<vector<float>, float> cases;
		cases.insert(pair<vector<float>, float>(c1, 2.0));

		for (auto& _case : cases) {
			float result = help::get_lowest_solution_for_quadratic(_case.first[0], _case.first[1], _case.first[2], _case.first[3]);
			if (result != _case.second) {
				printf("Case %s failed (result obtained: %s, correct answer: %s) \n",
					help::join_as_string(_case.first, ", ").c_str(), to_string(result).c_str(), to_string(_case.second).c_str()
				);
			}
		}

		if (!success) failed_tests.push_back("Linear PDF sampler");
		return success;
	}
	void create_tree_smaller_than_cell(Tree& tree) {
		pair<float, float> position(help::get_rand_float(0, dynamics.grid->width_r), help::get_rand_float(0, dynamics.grid->width_r));
		float cell_width = dynamics.grid->cell_width;
		float max_radius = 0.95f * (cell_width / 2.0);
		float dbh = 5;
		while (true) {
			tree = Tree(1, position, dbh, dynamics.seed_bearing_threshold);
			printf("max radius: %f, tree radius: %f, dbh: %f \n", max_radius, tree.radius, dbh);
			if (tree.radius < max_radius || tree.dbh < 0) break;
			else dbh *= 0.9f;
		}
		if (tree.dbh < 0) { printf("Error: could not create tree smaller than cell. \n"); throw; }
		if (verbosity > 0) printf("Created tree with dbh %f and radius %f (cell width = %f) \n", dbh, tree.radius, cell_width);
	}
	bool test_add_tree_if_not_exists(vector<string>& failed_tests) {
		bool success = true;

		// Setup
		Tree tree;
		create_tree_smaller_than_cell(tree);
		pair<int, int> tree_center_gb = dynamics.grid->get_gridbased_position(tree.position);
		Cell cell = dynamics.grid->distribution[dynamics.grid->pos_2_idx(tree_center_gb)];
		cell.LAI = get_rand_float(0, 5);
		float prev_LAI = cell.LAI;

		// Test
		cell.add_tree_if_not_exists(&tree, dynamics.grid->cell_area);

		// Check
		if (cell.trees.size() == 0) {
			if (verbosity > 0) printf("Tree not added to grid cell. \n");
			success = false;
		}
		float LAI_averaged_over_cell = (tree.LAI * tree.crown_area) / dynamics.grid->cell_area;
		if (!is_float_equal((cell.LAI - prev_LAI), LAI_averaged_over_cell)) {
			if (verbosity > 0) printf("LAI not updated correctly (expected %f, got %f). \n", LAI_averaged_over_cell, cell.LAI - prev_LAI);
			success = false;
		}

		if (!success) failed_tests.push_back("add_tree_if_not_exists");
		return success;
	}
	void run_all() {
		printf("Beginning tests...\n");
		vector<string> failed_tests = {};
		int successes = 0;

		// Run tests
		successes += test_flammability(failed_tests);
		successes += test_add_tree_if_not_exists(failed_tests);

		printf("Completed all tests. ");
		if (failed_tests.size() > 0) {
			printf("Tests failed(%i / %i) :\n - %s \n", failed_tests.size(), successes + failed_tests.size(),
				help::join(&failed_tests, "\n - ").c_str());
		}
		else printf("All tests (%i) successful.", successes);
	}
	Dynamics dynamics;
	int verbosity = 0;
};

