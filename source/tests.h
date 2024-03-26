#pragma once
#include "dynamics.h"

class Tests {
public:
	Tests() {};
	bool test_quadratic_solver(vector<string>& failed_tests, int verbosity) {
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
};

