#pragma once
#include "helpers.h"


class Tree {
public:
	Tree() = default;
	Tree(pair<float, float> _position) : position(_position) {};
	Tree(pair<float, float> _position, float _radius) : position(_position), radius(_radius) {};
	Tree(float _radius, vector<float> _strategy, int _life_phase, pair<float, float> _position):
		radius(_radius), strategy(_strategy), life_phase(_life_phase), position(_position)
	{
		id = help::get_rand_uint(0, 1e28);
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
	Population(float _mean_radius) : mean_radius(_mean_radius) {}
	Tree* add(pair<float, float> position) {
		Tree tree(position, mean_radius);
		members.push_back(tree);
		return &members.back();
	}
	int size() {
		return members.size();
	}

	void remove(Tree* tree) {
		auto it = std::find(members.begin(), members.end(), *tree);
		if (it != members.end()) { members.erase(it); }
	}
	vector<Tree> members = {};
	float mean_radius = 0;
};
