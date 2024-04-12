#pragma once
#include "helpers.h"


using namespace help;


class LinearDiffusionKernel : LinearProbabilityModel {
public:
	LinearDiffusionKernel() = default;
	LinearDiffusionKernel(int tree_id, float _q1, float _q2, float _min, float _max) : 
		LinearProbabilityModel(_q1, _q2, _min, _max) {}
	float get_ld_dist() {
		return LinearProbabilityModel::linear_sample();
	}
};


class WindKernel : PieceWiseLinearProbModel {
public:
	WindKernel() = default;
	WindKernel(int _tree_id, float _dist_max, float _wspeed_gmean, float _wspeed_stdev, float _seed_tspeed, float _abscission_height) :
		PieceWiseLinearProbModel(_dist_max, resolution)
	{
		dist_max = _dist_max;
		wspeed_gmean = _wspeed_gmean;
		wspeed_stdev = _wspeed_stdev;
		seed_tspeed = _seed_tspeed;
		abscission_height = _abscission_height;
	};
	float get_wind_dispersed_dist() {
		return PieceWiseLinearProbModel::sample();
	}
	float wspeed_gmean = 0;
	float wspeed_stdev = 0;
	float seed_tspeed = 0;
	float abscission_height = 0;
	float domain_size = 0;
	float dist_max = 0;
};


class Kernel : LinearDiffusionKernel, WindKernel {
public:
	Kernel() = default;
	Kernel(int _tree_id, float _q1, float _q2, float _min, float _max) :
		LinearDiffusionKernel(_tree_id, _q1, _q2, _min, _max)
	{
		id = _tree_id;
	}
	Kernel(int _tree_id, float _dist_max, float _wspeed_gmean, float _wspeed_stdev, float _seed_tspeed, float _abscission_height) :
		WindKernel(_tree_id, _dist_max, _wspeed_gmean, _wspeed_stdev, _seed_tspeed, _abscission_height)
	{
		id = _tree_id;
	}
	float get_ld_dist() {
		return LinearDiffusionKernel::get_ld_dist();
	}
	float get_wind_dispersed_dist() {
		return WindKernel::get_wind_dispersed_dist();
	}
	int id = -1;
};
