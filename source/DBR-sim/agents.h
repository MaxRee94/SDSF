#pragma once
#include "helpers.h"

/**
 * @class Tree
 * @brief Represents a tree in the DBR simulation.
 * 
 * Tree objects store position, radius, life phase, and strategy information.
 * Used to model individual trees in the forest-savanna simulation.
 */
class Tree {
public:
	/**
	 * @brief Default constructor.
	 */
	Tree() = default;
	
	/**
	 * @brief Constructor with position.
	 * @param _position The (x, y) position coordinates of the tree.
	 */
	Tree(pair<float, float> _position) : position(_position) {};
	
	/**
	 * @brief Full constructor with all tree properties.
	 * @param _radius The radius of the tree.
	 * @param _strategy Vector of strategy values for the tree.
	 * @param _life_phase The life phase of the tree (0=sapling, 1=mature, etc.).
	 * @param _position The (x, y) position coordinates of the tree.
	 */
	Tree(float _radius, vector<float> _strategy, int _life_phase, pair<float, float> _position):
		radius(_radius), strategy(_strategy), life_phase(_life_phase), position(_position) {};
	
	float radius = 0;                ///< Radius of the tree
	vector<float> strategy = {};     ///< Strategy parameters for the tree
	int life_phase = 0;              ///< Current life phase of the tree
	pair<float, float> position = pair(0, 0); ///< (x, y) position coordinates
};

/**
 * @class Population
 * @brief Represents a population of trees in the DBR simulation.
 * 
 * Manages a collection of Tree objects and provides methods for population operations.
 */
class Population {
public:
	/**
	 * @brief Default constructor.
	 */
	Population() = default;
	
	/**
	 * @brief Add a tree to the population.
	 * @param tree The Tree object to add to the population.
	 */
	void add(Tree tree) {
		members.push_back(tree);
	}
	
	vector<Tree> members = {}; ///< Collection of Tree objects in the population
};