#pragma once
#include "helpers.h"
#include "grid_agent.forward.h"


class Tree {
public:
	Tree() = default;
	Tree(pair<float, float> _position) : position(_position) {};
	Tree(int _id, pair<float, float> _position, float _radius) : position(_position), radius(_radius) {
		id = _id;
	};
	Tree(int _id, pair<float, float> _position, float _radius, vector<float> _strategy, int _life_phase):
		radius(_radius), strategy(_strategy), life_phase(_life_phase), position(_position)
	{
		id = _id;
	};
	bool operator==(const Tree& tree) const
	{
		return id == tree.id;
	}
	void grow(float max_radius) {
		// TEMP: constant growth rate. TODO: Make dependent on radius and life phase (resprout or seedling)
		radius = min(radius + 0.05 * max_radius, max_radius);
	}
	bool is_within_radius(pair<float, float> pos2) {
		float dist = help::get_dist(position, pos2);
		return dist < radius;
	}
	float radius = -1;
	vector<float> strategy = {};
	int life_phase = 0;
	pair<float, float> position = pair(0, 0);
	uint id = 0;
};


class Population {
public:
	Population() = default;
	Population(float _max_radius, float _cellsize, float _radius_q1, float _radius_q2) : max_radius(_max_radius),
		cellsize(_cellsize), radius_q1(_radius_q1), radius_q2(_radius_q2) 
	{
		//help::init_RNG();
	}
	Tree* add(pair<float, float> position, string radius="sampled") {
		float _radius = max_radius;
		if (radius == "sampled") _radius = help::sample_linear_distribution(radius_q1, radius_q2, 0, max_radius);
		Tree tree(no_created_trees + 1, position, _radius);
		members[tree.id] = tree;
		no_created_trees++;
		return &members[tree.id];
	}
	void store_positions() {
		for (auto& [id, tree] : members) {
			positions[tree.id] = tree.position;
		}
	}
	Tree* get(int id) {
		return &members[id];
	}
	void check_positions() {
		for (auto& [id, tree] : members) {
			if (tree.position != positions[tree.id]) {
				printf("id %i position (%f, %f) differs from stored position (%f, %f), id %i \n", tree.id, tree.position.first, tree.position.second, positions[tree.id].first, positions[tree.id].second, tree.id);
			}
		}
	}
	int size() {
		return members.size();
	}
	void grow() {
		for (auto& [id, tree] : members) {
			tree.grow(max_radius);
		}
	}
	void remove(Tree* tree) {
		members.erase(tree->id);
	}
	bool is_population_member(Tree* tree) {
		auto it = members.find(tree->id);
		return (it != members.end());
	}
	bool is_population_member(int tree_id) {
		return members.find(tree_id) != members.end();
	}
	unordered_map<int, Tree> members;
	float max_radius = 0;
	float cellsize = 0;
	float radius_q1 = 0;
	float radius_q2 = 0;
	Tree removed_tree;
	map<uint, pair<float, float>> positions;
	int no_created_trees = 0;
};
