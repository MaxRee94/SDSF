#pragma once
#include "helpers.h"
#include "grid_agent.forward.h"
#include "kernel.h"


using namespace help;

class Strategy {
public:
	Strategy() = default;
	Strategy(float _seed_mass, float _seed_density) : seed_mass(_seed_mass), seed_density(_seed_density) {}
	int id = -1;
	float seed_mass = 0;
	float seed_density = 0;
};


class Tree {
public:
	Tree() = default;
	Tree(pair<float, float> _position) : position(_position) {};
	Tree(int _id, pair<float, float> _position, float _radius, float _radius_tmin1) : position(_position), radius(_radius), radius_tmin1(_radius_tmin1) {
		id = _id;
	};
	Tree(int _id, pair<float, float> _position, float _radius, float _radius_tmin1, int _life_phase) :
		position(_position), radius(_radius), radius_tmin1(_radius_tmin1), life_phase(_life_phase)
	{
		id = _id;
	};
	bool operator==(const Tree& tree) const
	{
		return id == tree.id;
	}
	bool is_within_radius(pair<float, float> pos2) {
		float dist = help::get_dist(position, pos2);
		return dist < radius;
	}
	float radius = -1;
	float radius_tmin1 = -1;
	int life_phase = 0;
	int last_mortality_check = 0;
	pair<float, float> position = pair(0, 0);
	uint id = -1;
};


class Crop {
public:
	Crop() = default;
	Crop(Strategy &strategy, Tree& tree, float _mass_budget_factor) {
		seed_mass = strategy.seed_mass;
		mass_budget_factor = _mass_budget_factor;
		origin = tree.position;
		id = tree.id;
	}
	void compute_mass(Tree& tree) {
		float radius_avg = (tree.radius + tree.radius_tmin1) * 0.5;
		//printf("radius avg: %f, cur radius: %f, old radius: %f, cubed radius: %f \n", radius_avg, tree_radius, tree_radius_tmin1, cubed(radius_avg));
		mass = cubed(radius_avg) - (cubed(tree.radius) - cubed(tree.radius_tmin1));
		mass *= mass_budget_factor;
	}
	void compute_no_seeds() {
		no_seeds = mass / seed_mass;
	}
	void update(Tree& tree) {
		compute_mass(tree);
		compute_no_seeds();
	}
	float mass = 0;
	int no_seeds = 0;
	float seed_mass = 0;
	float mass_budget_factor = 0;
	Kernel* kernel = 0;
	pair<float, float> origin = pair<float, float>(0, 0);
	int id = -1;
};


class Population {
public:
	Population() = default;
	Population(float _max_radius, float _cellsize, float _radius_q1, float _radius_q2, float _mass_budget_factor) : max_radius(_max_radius),
		cellsize(_cellsize), radius_q1(_radius_q1), radius_q2(_radius_q2), mass_budget_factor(_mass_budget_factor)
	{
		radius_probability_model = help::LinearProbabilityModel(radius_q1, radius_q2, 0, max_radius);
	}
	Tree* add(pair<float, float> position, Strategy &strategy, float radius = -2) {
		// Create tree
		if (radius == -1) radius = max_radius;
		else if (radius == -2) radius = radius_probability_model.sample();
		float radius_tmin1 = radius * 0.9; // TEMP. TODO: Make dependent on growth curve.
		Tree tree(no_created_trees + 1, position, radius, radius_tmin1, 1);
		members[tree.id] = tree;
		no_created_trees++;

		// Create crop
		strategy.id = tree.id;
		strategies[tree.id] = strategy;
		Crop crop(strategy, tree, mass_budget_factor);
		crop.update(tree);
		crops[tree.id] = crop;

		return &members[tree.id];
	}
	Kernel* add_kernel_linear_diffusion(int tree_id, float q1, float q2, float min, float max) {
		Kernel kernel(tree_id, q1, q2, min, max);
		kernels[tree_id] = kernel;
		return &kernel;
	}
	Kernel* add_kernel(int tree_id, Kernel &kernel) {
		kernels[tree_id] = kernel;
		return &kernel;
	}
	Tree* get(int id) {
		return &members[id];
	}
	void get(vector<int>& ids, vector<Tree*> &trees) {
		for (int id : ids) trees.push_back(get(id));
	}
	Crop* get_crop(int id) {
		return &crops[id];
	}
	Strategy* get_strat(int id) {
		return &strategies[id];
	}
	int size() {
		return members.size();
	}
	void remove(Tree* tree) {
		members.erase(tree->id);
	}
	void remove(int id) {
		members.erase(id);
	}
	bool is_population_member(Tree* tree) {
		auto it = members.find(tree->id);
		return (it != members.end());
	}
	bool is_population_member(int tree_id) {
		return members.find(tree_id) != members.end();
	}
	unordered_map<int, Tree> members;
	unordered_map<int, Crop> crops;
	unordered_map<int, Kernel> kernels;
	unordered_map<int, Strategy> strategies;
	help::LinearProbabilityModel radius_probability_model;
	float max_radius = 0;
	float cellsize = 0;
	float radius_q1 = 0;
	float radius_q2 = 0;
	float seed_mass = 0;
	Tree removed_tree;
	int no_created_trees = 0;
	float mass_budget_factor = 0;
};
