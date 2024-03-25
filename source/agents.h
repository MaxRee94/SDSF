#pragma once
#include "helpers.h"


class Tree {
public:
	Tree() = default;
	Tree(pair<float, float> _position) : position(_position) {};
	Tree(pair<float, float> _position, float _radius) : position(_position), radius(_radius) {
		id = help::get_rand_uint(0, 10e20);
	};
	Tree(float _radius, vector<float> _strategy, int _life_phase, pair<float, float> _position):
		radius(_radius), strategy(_strategy), life_phase(_life_phase), position(_position)
	{
		id = help::get_rand_uint(0, 10e20);
	};
	bool operator==(const Tree& tree) const
	{
		return id == tree.id;
	}
	float radius = 0;
	vector<float> strategy = {};
	int life_phase = 0;
	pair<float, float> position = pair(0, 0);
	int id = 0;

};


class Population {
public:
	Population() = default;
	Population(float _max_radius, float _cellsize) : max_radius(_max_radius), cellsize(_cellsize) {
		//help::init_RNG();
	}
	Tree* add(pair<float, float> position) {
		Tree tree(position, max_radius);
		members.push_back(tree);
		return &members.back();
	}
	int size() {
		return members.size();
	}
	vector<Tree*> get_neighbors(Tree* base) {
		vector<Tree*> neighbors;
		float dist = (base->radius + max_radius * 1.1);
		for (auto& candidate : members) {
			if (candidate == *base) continue;
			if (help::get_dist(candidate.position, base->position) < dist) {
				neighbors.push_back(&candidate);
			}
		}
		return neighbors;
	}
	void remove(Tree* tree) {
		auto it = std::find(members.begin(), members.end(), *tree);
		if (it != members.end()) { members.erase(it); }
	}
	vector<Tree> members = {};
	float max_radius = 0;
	float cellsize = 0;
};
