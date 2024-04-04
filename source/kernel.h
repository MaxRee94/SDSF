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
		init_linear_probability_model();
	}
	void init_linear_probability_model() {
		float xrange = (max - min);
		a = (q2 - q1) / xrange;
		b = q1;
		float cdf_of_xrange = 0.5 * a * xrange * xrange + b * xrange; // CDF(xrange)

		// Scale PDF such that it integrates to 1, i.e., so that CDF(xrange) = 1
		a /= cdf_of_xrange;
		b /= cdf_of_xrange;
	}
	float sample() {
		float cdf_of_dx = help::get_rand_float(0, 1); // Obtain CDF(dx)
		float dx = get_lowest_solution_for_quadratic(cdf_of_dx, a / 2.0, b, 0); // Get corresponding value dx
		return min + dx;
	}
	float q1 = 0;
	float q2 = 0;
	float min = 0;
	float max = 0;
	float a = 0;
	float b = 0;
	int id = -1;
};
