#pragma once
#include "dynamics.h"

class Tests {
public:
	Tests() {};
	Tests(int _timestep, float _cell_width, float _self_ignition_factor, float _rainfall, float _seed_bearing_threshold,
		float _growth_rate_multiplier, float _unsuppressed_flammability, float _min_suppressed_flammability, float _max_suppressed_flammability,
		float _radius_suppr_flamm_min, float radius_range_suppr_flamm, float _max_dbh, float _saturation_threshold, float _fire_resistance_argmin,
		float _fire_resistance_argmax, float _fire_resistance_stretch, float _background_mortality, map<string, map<string, float>> _strategy_distribution_params,
		float _resource_grid_relative_size, float _mutation_rate, float _STR, int _verbosity, int grid_size, float dbh_q1, float dbh_q2, float growth_multiplier_stdev,
		float growth_multiplier_min, float growth_multiplier_max, int random_seed, int firefreq_random_seed, float enforce_no_recruits
	) {
		dynamics = Dynamics(_timestep, _cell_width, _self_ignition_factor, _rainfall, _seed_bearing_threshold,
			_growth_rate_multiplier, _unsuppressed_flammability, _min_suppressed_flammability, _max_suppressed_flammability,
			_radius_suppr_flamm_min, radius_range_suppr_flamm, _max_dbh, _saturation_threshold, _fire_resistance_argmin,
			_fire_resistance_argmax, _fire_resistance_stretch, _background_mortality, _strategy_distribution_params,
			_resource_grid_relative_size, _mutation_rate, _STR, _verbosity, random_seed, firefreq_random_seed, enforce_no_recruits
		);
		dynamics.init_state(grid_size, dbh_q1, dbh_q2, growth_multiplier_stdev, growth_multiplier_min, growth_multiplier_max);
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
	void create_tree_sapling(Tree& tree) {
		pair<float, float> position(help::get_rand_float(0, dynamics.grid->width_r), help::get_rand_float(0, dynamics.grid->width_r));
		float cell_width = dynamics.grid->cell_width;
		float max_radius = 0.95f * (cell_width / 2.0);
		float dbh = 5;
		Timer t; t.start();
		while (true) {
			tree = Tree(1, position, dbh, dynamics.seed_bearing_threshold, dynamics.pop->resprout_growthcurve, 1);
			if (tree.radius < max_radius || tree.dbh < 0) break;
			else dbh *= 0.9f;
			if (t.elapsedSeconds() > 1) {
				printf("Error: Tree creation taking too long. Terminating program. \n"); throw;
			}
		}
		if (verbosity > 0) printf("Created tree with dbh %f and radius %f (cell width = %f) \n", dbh, tree.radius, cell_width);
	}
	bool test_add_tree_if_not_present(vector<string>& failed_tests) {
		bool success = true;

		// Setup
		Tree tree;
		create_tree_sapling(tree);
		pair<int, int> tree_center_gb = dynamics.grid->get_gridbased_position(tree.position);
		Cell cell = dynamics.grid->distribution[dynamics.grid->pos_2_idx(tree_center_gb)];
		cell.set_LAI(3);
		float prev_LAI = cell.get_LAI();

		// Test
		cell.add_tree_if_not_present(&tree, dynamics.grid->cell_area, dynamics.grid->cell_halfdiagonal_sqrt);

		// Check
		if (cell.trees.size() == 0) {
			if (verbosity > 0) printf("Tree not added to grid cell. \n");
			success = false;
		}
		float true_LAI_averaged_over_cell = (tree.LAI * tree.crown_area) / dynamics.grid->cell_area;
		float LAI_diff = cell.get_LAI() - prev_LAI;
		if (!approx(LAI_diff, true_LAI_averaged_over_cell, 0.001f)) {
			if (verbosity > 0) printf("LAI not updated correctly (expected %f, got %f). \n", true_LAI_averaged_over_cell, LAI_diff);
			success = false;
		}

		if (!success) failed_tests.push_back("add_tree_if_not_present");
		return success;
	}
	bool test_remove_tree_if_sapling(vector<string>& failed_tests) {
		bool success = true;

		// Setup
		Tree tree;
		create_tree_sapling(tree);
		pair<int, int> tree_center_gb = dynamics.grid->get_gridbased_position(tree.position);
		Cell cell = dynamics.grid->distribution[dynamics.grid->pos_2_idx(tree_center_gb)];
		cell.set_LAI(3);
		float prev_LAI = cell.get_LAI();

		// Test
		cell.remove_tree_if_sapling(&tree, dynamics.grid->cell_area, dynamics.grid->cell_halfdiagonal_sqrt);

		// Check
		if (cell.tree_is_present(&tree)) {
			if (verbosity > 0) printf("Tree not removed from cell. \n");
			success = false;
		}
		float true_LAI_averaged_over_cell = (tree.LAI * tree.crown_area) / dynamics.grid->cell_area;
		float LAI_diff = prev_LAI - cell.get_LAI();
		if (!approx(LAI_diff, true_LAI_averaged_over_cell, 0.001f)) {
			if (verbosity > 0) printf("LAI not updated correctly (expected %f, got %f). \n", true_LAI_averaged_over_cell, LAI_diff);
			success = false;
		}

		if (!success) failed_tests.push_back("remove_tree_if_sapling");
		return success;
	}
	bool test_is_float_equal(vector<string>& failed_tests) {
		bool success = true;

		// Cases
		vector<float> c1 = { -4.5551, -4.5551, -4.5551, 0, 0, 0, 1, 1, 1 };
		vector<float> c2 = { -4.5550, -4.5551, -4.5552, -0.0001, 0, 0.0001, 0.9999, 1, 1.0001 };
		vector<bool> answers = { false, true, false, false, true, false, false, true, false };

		// Test and check
		for (int i = 0; i < c1.size(); i++) {
			if (is_float_equal(c1[i], c2[i]) != answers[i]) {
				if (verbosity > 0) printf(
					"Case %f, %f failed (Gave answer (%i) that should be (%i)) \n", c1[i], c2[i], is_float_equal(c1[i], c2[i]), (int)answers[i]
				);
				success = false;
			}
		}

		if (!success) failed_tests.push_back("is_float_equal");
		return success;
	}
	bool test_approx(vector<string>& failed_tests) {
		bool success = true;

		// Cases
		vector<float> c1 = { -4.5551, -4.5551, -4.565, 0, 0, 0, 1, 1, 1 };
		vector<float> c2 = { -4.5550, -4.5551, -4.55, -0.0001, 0, 0.02, 0.9999, 1, 1.02 };
		vector<bool> answers = { true, true, false, true, true, false, true, true, false };

		// Test and check
		for (int i = 0; i < c1.size(); i++) {
			if (approx(c1[i], c2[i], 0.01f) != answers[i]) {
				if (verbosity > 0) printf(
					"Case %f, %f failed (Gave answer (%i) that should be (%i)) \n", c1[i], c2[i], approx(c1[i], c2[i], 0.01f), (int)answers[i]
				);
				success = false;
			}
		}

		if (!success) failed_tests.push_back("approx");
		return success;
	}
	bool test_readable_number(vector<string>& failed_tests) {
		bool success = true;

		// Cases
		vector<int> nums = { 5005030, 1031066456 };
		vector<string> answers = { "5 005 030", "1 031 066 456"};

		// Test and check
		for (int i = 0; i < nums.size(); i++) {
			if (readable_number(nums[i]) != answers[i]) {
				if (verbosity > 0) printf(
					"Case %i failed (Gave answer %s that should be %s) \n", nums[i], readable_number(nums[i]).c_str(), answers[i].c_str()
				);
				success = false;
			}
		}

		if (!success) failed_tests.push_back("readable_number");
		return success;
	}
	void run_all() {
		printf("Beginning tests...\n");
		vector<string> failed_tests = {};
		int successes = 0;

		// Run tests
		successes += test_flammability(failed_tests);
		successes += test_add_tree_if_not_present(failed_tests);
		successes += test_remove_tree_if_sapling(failed_tests);
		successes += test_is_float_equal(failed_tests);
		successes += test_approx(failed_tests);
		successes += test_readable_number(failed_tests);

		printf("Completed all tests. ");
		if (failed_tests.size() > 0) {
			printf("\nTests failed(%i / %i) :\n - %s \n", failed_tests.size(), successes + failed_tests.size(),
				help::join(&failed_tests, "\n - ").c_str());
		}
		else printf("\nAll tests (%i) successful.", successes);
	}
	Dynamics dynamics;
	int verbosity = 0;
};

