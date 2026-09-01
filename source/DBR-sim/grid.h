#pragma once
#include "agents.h"

/**
 * @class Cell
 * @brief Represents a single cell in the Grid.
 * 
 * Each cell can be in a savanna or forest state and can contain a tree.
 */
class Cell {
public:
	/**
	 * @brief Default constructor.
	 */
	Cell() = default;
	
	int state = 0;        ///< Cell state: 0 = savanna, 1 = forest
	Tree* tree = 0;      ///< Pointer to the tree in this cell (nullptr if no tree)
};

/**
 * @class Grid
 * @brief Represents the spatial grid for the DBR simulation.
 * 
 * Manages a 2D array of Cell objects and provides methods for grid operations,
 * tree cover calculation, and cell state management.
 */
class Grid {
public:
	/**
	 * @brief Default constructor.
	 */
	Grid() = default;
	
	/**
	 * @brief Constructor with grid size.
	 * @param _size The size of the grid (size x size cells).
	 */
	Grid(int _size) {
		size = _size;
		distribution = new Cell[size * size];
		no_savanna_cells = size;
		size_r = (float)size * cell_width;
	}
	
	/**
	 * @brief Reset the grid to initial state.
	 * 
	 * Sets all cells to savanna state and resets counters.
	 */
	void reset() {
		for (int i = 0; i < size * size; i++) {
			distribution[i].state = 0;
		}
		no_forest_cells = 0;
		no_savanna_cells = size * size;
	}
	
	/**
	 * @brief Set a cell to savanna state.
	 * @param idx The index of the cell to set to savanna.
	 */
	void set_to_savanna(int idx) {
		if (distribution[idx].state == 1) {
			no_forest_cells -= 1;
			no_savanna_cells += 1;
		}
		distribution[idx].state = 0;
	}
	
	/**
	 * @brief Recalculate cell state counts.
	 * 
	 * Iterates through all cells and recounts forest vs savanna cells.
	 */
	void redo_count() {
		no_savanna_cells = 0;
		no_forest_cells = 0;
		for (int i = 0; i < size * size; i++) {
			no_savanna_cells += !distribution[i].state;
			no_forest_cells += distribution[i].state;
		}
	}
	
	/**
	 * @brief Get the current tree cover as a fraction.
	 * @return float The fraction of cells that are forest (tree cover).
	 */
	float get_tree_cover() {
		return (float)no_forest_cells / (float)no_savanna_cells;
	}
	
	/**
	 * @brief Get the cell at a specific position.
	 * @param pos The (x, y) position coordinates.
	 * @return Cell* Pointer to the Cell at the given position.
	 */
	Cell* get_cell_at_position(pair<float, float> pos) {
		return &distribution[(int)(pos.second * size + pos.first)];
	}
	
	/**
	 * @brief Populate the tree domain in the grid.
	 * @param tree Pointer to the Tree object to populate.
	 */
	void populate_tree_domain(Tree* tree) {
		pair<float, float> tree_center = get_gridbased_position(tree);
		cout << "tree center: " << tree_center.first << ", " << tree_center.second << endl;
		for (int x = tree_center.first - (tree->radius / 2); x < tree_center.first + (tree->radius / 2); x++) {
			for (int y = tree_center.second - (tree->radius / 2); y < tree_center.second + (tree->radius / 2); y++) {
				if (help::get_dist(pair<float, float>(x, y), tree_center) < tree->radius) {
					set_to_forest(pair<int, int>(x, y), tree);
				}
			}
		}
	}
	
	/**
	 * @brief Get the current state distribution as an array.
	 * @return int* Array of cell states.
	 */
	int* get_state_distribution() {
		if (state_distribution = 0) state_distribution = new int[size * size];
		for (int i = 0; i < size * size; i++) {
			state_distribution[i] = distribution[i].state;
		}
		return state_distribution;
	}
	
	/**
	 * @brief Set a cell to forest state.
	 * @param idx The index of the cell to set to forest.
	 * @param tree Pointer to the Tree object to place in the cell.
	 */
	void set_to_forest(int idx, Tree* tree) {
		if (distribution[idx].state == 0) {
			no_savanna_cells -= 1;
			no_forest_cells += 1;
		}
		distribution[idx].state = 1;
		distribution[idx].tree = tree;
	}
	
	/**
	 * @brief Set a cell to forest state using grid coordinates.
	 * @param position_grid The (x, y) grid coordinates of the cell.
	 * @param tree Pointer to the Tree object to place in the cell.
	 */
	void set_to_forest(pair<int, int> position_grid, Tree* tree) {
		set_to_forest(position_grid.second * size + position_grid.first, tree);
	}
	
	/**
	 * @brief Get the grid-based position of a tree.
	 * @param tree Pointer to the Tree object.
	 * @return pair<int, int> The (x, y) grid coordinates corresponding to the tree's position.
	 */
	pair<int, int> get_gridbased_position(Tree* tree) {
		return pair<int, int>(tree->position.first / cell_width, tree->position.second / cell_width);
	}
	
	int size = 1000;               ///< Grid size (size x size cells)
	float size_r = 0;              ///< Real-world size of the grid
	float cell_width = 1.5;       ///< Width of each cell in real units
	Cell* distribution = 0;       ///< 2D array of Cell objects
	int* state_distribution = 0;  ///< Cached state distribution array
	int no_savanna_cells = 0;      ///< Number of savanna cells
	int no_forest_cells = 0;      ///< Number of forest cells
};