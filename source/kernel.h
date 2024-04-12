#pragma once
#include "helpers.h"


using namespace help;


class LinearDiffusionKernel : public LinearProbabilityModel {
public:
	LinearDiffusionKernel() = default;
	LinearDiffusionKernel(int tree_id, float _q1, float _q2, float _min, float _max) : 
		LinearProbabilityModel(_q1, _q2, _min, _max) {}
	float get_ld_dist() {
		return LinearProbabilityModel::linear_sample();
	}
};


class WindKernel : public PieceWiseLinearProbModel {
public:
	WindKernel() = default;
	WindKernel(int _tree_id, float _dist_max, float _wspeed_gmean, float _wspeed_stdev, float _seed_tspeed, float _abs_height) :
		PieceWiseLinearProbModel(_dist_max)
	{
		dist_max = _dist_max;
		wspeed_gmean = _wspeed_gmean;
		wspeed_stdev = _wspeed_stdev;
		seed_tspeed = _seed_tspeed;
		abs_height = _abs_height;
		build();
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
	float wspeed_gmean = 0;
	float wspeed_stdev = 0;
	float seed_tspeed = 0;
	float abs_height = 0;
	float domain_size = 0;
	float dist_max = 0;
};


class Kernel : public LinearDiffusionKernel, public WindKernel {
public:
	Kernel() = default;
	Kernel(int _tree_id, float _q1, float _q2, float _min, float _max) :
		LinearDiffusionKernel(_tree_id, _q1, _q2, _min, _max)
	{
		type = "linear";
	}
	Kernel(int _tree_id, float _dist_max, float _wspeed_gmean, float _wspeed_stdev, float _seed_tspeed, float _abs_height) :
		WindKernel(_tree_id, _dist_max, _wspeed_gmean, _wspeed_stdev, _seed_tspeed, _abs_height)
	{
		id = _tree_id;
		type = "wind";
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
