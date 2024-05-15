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
		verbosity = _verbosity;
	}
		bool success = false;

		// Cases
		vector<float> c1 = { 2.0, 0.5, 1.5, 0.5 };
		map<vector<float>, float> cases;
		cases.insert(pair<vector<float>, float>(c1, 2.0));

		for (auto& _case : cases) {
			float result = help::get_lowest_solution_for_quadratic(_case.first[0], _case.first[1], _case.first[2], _case.first[3]);
			if (result != _case.second) {
				printf("Case %s failed (result obtained: %s, correct answer: %s) \n", help::join_as_string(_case.first, ", "), to_string(result), to_string(_case.second));
			}
		}

		if (!success) failed_tests.push_back("Linear PDF sampler");
		return success;
	}
	void run_all(int verbosity) {
		printf("Beginning tests...\n");
		vector<string> failed_tests = {};
		int successes = 0;

		// Run tests
		successes += test_quadratic_solver(failed_tests, verbosity);

		printf("Completed all tests. ");
		if (failed_tests.size() > 0) {
			printf("Tests failed(%i / %i) : \n %s \n", failed_tests.size(), successes,
				help::join(&failed_tests, "\n"));
		}
		else printf("All tests (%i) successful.", successes);
	}
	Dynamics dynamics;
	int verbosity = 0;
};

