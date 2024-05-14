#pragma once
#include "helpers.h"


using namespace help;


class LinearDiffusionKernel : public LinearProbabilityModel {
public:
	LinearDiffusionKernel() = default;
	LinearDiffusionKernel(float _q1, float _q2, float _min, float _max) : 
		LinearProbabilityModel(_q1, _q2, _min, _max) {
	}
	float get_ld_dist() {
		return LinearProbabilityModel::linear_sample();
	}
};


class WindKernel : public PieceWiseLinearProbModel {
public:
	WindKernel() = default;
	WindKernel(float _dist_max, float _wspeed_gmean, float _wspeed_stdev, float _wind_direction, float _wind_direction_stdev, float _seed_tspeed = 0, float _abs_height = 0) :
		PieceWiseLinearProbModel(_dist_max)
	{
		dist_max = _dist_max;
		wspeed_gmean = _wspeed_gmean;
		wspeed_stdev = _wspeed_stdev;
		wind_direction = _wind_direction;
		wind_direction_stdev = _wind_direction_stdev;
		seed_tspeed = _seed_tspeed;
		abs_height = _abs_height;
	};
	float get_wind_dispersed_dist() {
		return sample();
	}
	float pdf(float x) override {
		x = max(0.1f, x);
		float natural_log = log((seed_tspeed * x) / (wspeed_gmean * abs_height));
		float second_fraction = (natural_log * natural_log) / (2.0f * wspeed_stdev * wspeed_stdev);
		float first_fraction = 1.0f / (x * sqrtf(2.0f * M_PI) * wspeed_stdev);
		return first_fraction * exp(-second_fraction);
	}
	void update(float tree_height) {
		abs_height = tree_height * 0.8f; // Constant factor for now; might change this to be a function of tree height.
		if ((abs_height - prev_build_height) > 2) { // We rebuild the kernel if the tree height has increased by >2 meters since the last build.
			build();
			prev_build_height = abs_height;
		}
	}
	float wspeed_gmean = 0;
	float wspeed_stdev = 0;
	float wind_direction = 0;
	float wind_direction_stdev = 0;
	float seed_tspeed = 0;
	float abs_height = 0;
	float prev_build_height = -100;
	float domain_size = 0;
	float dist_max = 0;
};


class Kernel : public LinearDiffusionKernel, public WindKernel {
public:
	Kernel() = default;
	Kernel(int _tree_id, float _q1, float _q2, float _min, float _max) :
		LinearDiffusionKernel(_q1, _q2, _min, _max)
	{
		id = _tree_id;
		type = "linear";
	}
	Kernel(int _tree_id, float _dist_max, float _wspeed_gmean, float _wspeed_stdev, float _wind_direction, float _wind_direction_stdev, float _seed_tspeed = 0, float _abs_height = 0) :
		WindKernel(_dist_max, _wspeed_gmean, _wspeed_stdev, _wind_direction, _wind_direction_stdev, _seed_tspeed, _abs_height)
	{
		id = _tree_id;
		type = "wind";
	}
	Kernel(int _tree_id, string _type) : type(_type), id(_tree_id) {
		// Arbitrary kernel without built-in functionality, currently only used for animal dispersal
	}
	float get_ld_dist() {
		return LinearDiffusionKernel::get_ld_dist();
	}
	float get_wind_dispersed_dist() {
		return WindKernel::get_wind_dispersed_dist();
	}
	float get_dist() {
		if (type == "linear") return get_ld_dist();
		else if (type == "wind") return get_wind_dispersed_dist();
	}
	int id = -1;
	string type = "none";
};
