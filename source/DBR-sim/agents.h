#pragma once
#include "helpers.h"


class Tree {
public:
	Tree() = default;
	Tree(pair<float, float> _position) : position(_position) {};
	Tree(float _radius, vector<float> _strategy, int _life_phase, pair<float, float> _position):
		radius(_radius), strategy(_strategy), life_phase(_life_phase), position(_position) {};
	float radius = 0;
	vector<float> strategy = {};
	int life_phase = 0;
	pair<float, float> position = pair(0, 0);
};


class Population {
public:
	Population() = default;
	void add(Tree tree) {
		members.push_back(tree);
	}
	vector<Tree> members = {};
};
