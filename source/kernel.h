#pragma once
#include "helpers.h"


using namespace help;

class Kernel {
public:
	Kernel() = default;
	Kernel(int tree_id, float _q1, float _q2, float _min, float _max) {
		q1 = _q1;
		q2 = _q2;
		min = _min;
		max = _max;
		id = tree_id;
		prob_model = help::LinearProbabilityModel(q1, q2, min, max);
	}
	float sample() {
		return prob_model.sample();
	}
	float q1 = 0;
	float q2 = 0;
	float min = 0;
	float max = 0;
	float a = 0;
	float b = 0;
	help::LinearProbabilityModel prob_model;
	int id = -1;
};
