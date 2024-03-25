#pragma once
#include "dynamics.h"

class Tests {
public:
	Tests() {};
	bool test_linear_pdf_sampler(vector<string>& failed_tests, bool verbose = false) {
		printf("Testing..\n");
		failed_tests.push_back("Linear PDF sampler");
		return false;
	}
	void run_all(bool verbose = false) {
		printf("Beginning tests...\n");
		vector<string> failed_tests = {};
		int successes = 0;

		// Run tests
		successes += test_linear_pdf_sampler(failed_tests, verbose);

		printf("Completed all tests. ");
		if (failed_tests.size() > 0) {
			printf("Tests failed(%i / %i) : \n %s \n", failed_tests.size(), successes,
				help::join(&failed_tests, "\n"));
		}
		else printf("All tests (%i) successful.", successes);
	}
};

