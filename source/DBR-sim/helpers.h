#pragma once
#include <vector>
#include <map>
#include <cstdlib>
#include <set>
#include <Windows.h>
#include <string>
#include <iostream>

#define VERBOSE false

using namespace std;

typedef unsigned int uint;

/**
 * @struct _comparator
 * @brief Comparator for sorting pairs by their second value.
 * 
 * Used to sort pairs in increasing order of their second value,
 * with the first value as a tiebreaker.
 */
struct _comparator {
	template <typename T>

	// Comparator function
	bool operator()(const T& l, const T& r) const
	{
		if (l.second != r.second) {
			return l.second < r.second;
		}
		return l.first < r.first;
	}
};

typedef std::set < std::pair<int, double>, _comparator> PairSet;

namespace help {

	/**
	 * @brief Initialize the random number generator.
	 * 
	 * Seeds the random number generator using the current time.
	 */
	void init_RNG();

	/**
	 * @brief Sort a map into a set of pairs ordered by value.
	 * @param _map The map to sort.
	 * @param _set The set to populate with sorted pairs.
	 */
	void sort(std::map<int, double>& _map, PairSet& _set);

	/**
	 * @brief Get the key for a given value in a map.
	 * @param _map Pointer to the map to search.
	 * @param value The value to find the key for.
	 * @return int The key corresponding to the value, or -1 if not found.
	 */
	int get_key(std::map<int, int>* _map, int value);

	/**
	 * @brief Fill a 2D double array with zeros.
	 * @param _array Pointer to the array to populate.
	 * @param dim_x The x dimension size.
	 * @param dim_y The y dimension size.
	 */
	void populate_with_zeroes(double* _array, int dim_x, int dim_y);
	
	/**
	 * @brief Fill a 2D unsigned int array with zeros.
	 * @param _array Pointer to the array to populate.
	 * @param dim_x The x dimension size.
	 * @param dim_y The y dimension size.
	 */
	void populate_with_zeroes(uint* _array, int dim_x, int dim_y);

	/**
	 * @brief Split a string into substrings using a separator.
	 * @param basestring The string to split.
	 * @param separator The separator string to split on.
	 * @param substrings Reference to vector to store the resulting substrings.
	 */
	void split(std::string basestring, std::string separator, vector<std::string>& substrings);

	/**
	 * @brief Check if a vector contains a specific integer.
	 * @param vec Pointer to the vector to search.
	 * @param item The integer to search for.
	 * @return bool true if the item is found, false otherwise.
	 */
	bool is_in(std::vector<int>* vec, int item);

	/**
	 * @brief Print the contents of a map to standard output.
	 * @param map Pointer to the map to print.
	 */
	void print_map(std::map<int, int>* map);

	/**
	 * @brief Print the contents of a vector to standard output.
	 * @param vec Pointer to the vector to print.
	 */
	void print_vector(std::vector<int>* vec);

	/**
	 * @brief Print the contents of a vector of pairs to standard output.
	 * @param vec Pointer to the vector of pairs to print.
	 */
	void print_pairs(std::vector<pair<int, int>>* vec);

	/**
	 * @brief Print a string (only if VERBOSE is true).
	 * @param str The string to print.
	 */
	void print(std::string);

	/**
	 * @brief Find all occurrences of a substring in a string.
	 * @param basestring The string to search in.
	 * @param target The substring to find.
	 * @return vector<size_t> Vector of positions where the target is found.
	 */
	vector<size_t> FindAll(std::string basestring, std::string target);

	/**
	 * @brief Check if a string contains a substring.
	 * @param basestring The string to search in.
	 * @param target The substring to find.
	 * @return bool true if the substring is found, false otherwise.
	 */
	bool is_in(std::string basestring, std::string target);

	/**
	 * @brief Replace all occurrences of a substring with another string.
	 * @param basestring The original string.
	 * @param toReplace The substring to replace.
	 * @param replaceWith The string to replace with.
	 * @return string The modified string with all replacements made.
	 */
	std::string replace_occurrences(std::string basestring, std::string toReplace, std::string replaceWith);

	/**
	 * @brief Compute fast inverse square root.
	 * @param n The number to compute the inverse square root of.
	 * @return float The approximate inverse square root of n.
	 */
	double fisqrt(float n);

	/**
	 * @brief Increment a key's value in a map.
	 * @param _map Pointer to the map to modify.
	 * @param key The key to increment.
	 */
	void increment_key(std::map<std::string, int>* _map, std::string key);

	/**
	 * @brief Get a float value from a string-keyed map.
	 * @param _map Pointer to the map to search.
	 * @param key The key to look up.
	 * @return float The value for the key, or 0 if not found.
	 */
	float get_value(std::map<std::string, float>* _map, std::string key);

	/**
	 * @brief Get an int value from a string-keyed map.
	 * @param _map Pointer to the map to search.
	 * @param key The key to look up.
	 * @return int The value for the key, or 0 if not found.
	 */
	int get_value(std::map<std::string, int>* _map, std::string key);

	/**
	 * @brief Get an int value from an int-keyed map.
	 * @param _map Pointer to the map to search.
	 * @param key The key to look up.
	 * @return int The value for the key, or -1 if not found.
	 */
	int get_value(std::map<int, int>* _map, int key);

	/**
	 * @brief Get a double value from an int-keyed map.
	 * @param _map Pointer to the map to search.
	 * @param key The key to look up.
	 * @return int The value for the key, or -1 if not found.
	 */
	int get_value(std::map<int, double>* _map, int key);

	/**
	 * @brief Get an int value from an uint32_t-keyed map.
	 * @param _map Pointer to the map to search.
	 * @param key The key to look up.
	 * @return int The value for the key, or -1 if not found.
	 */
	int get_value(std::map<uint32_t, uint32_t>* _map, uint32_t key);

	/**
	 * @brief Get free RAM memory information.
	 * @return vector<float> Vector containing RAM, VM, Pagefile sizes in GB, and memory percentage.
	 */
	vector<float> get_free_memory();

	/**
	 * @brief Calculate standard deviation of a distribution.
	 * @param distribution Pointer to the vector of data values.
	 * @param mean Optional pre-calculated mean value.
	 * @return double The standard deviation of the distribution.
	 */
	double get_stdev(vector<double>* distribution, double mean = -999999);

	/**
	 * @brief Calculate mean of a distribution.
	 * @param distribution Pointer to the vector of data values.
	 * @return double The mean of the distribution.
	 */
	double get_mean(vector<double>* distribution);

	/**
	 * @brief Get maximum value from a distribution.
	 * @param distribution Pointer to the vector of data values.
	 * @return double The maximum value in the distribution.
	 */
	double get_max(vector<double>* distribution);

	/**
	 * @brief Get minimum value from a distribution.
	 * @param distribution Pointer to the vector of data values.
	 * @return double The minimum value in the distribution.
	 */
	double get_min(vector<double>* distribution);

	/**
	 * @brief Generate a random float in a range.
	 * @param min The minimum value of the range.
	 * @param max The maximum value of the range.
	 * @return float A random float between min and max.
	 */
	float get_rand_float(float min, float max);

	/**
	 * @brief Generate a random unsigned integer in a range.
	 * @param min The minimum value of the range.
	 * @param max The maximum value of the range.
	 * @return uint A random unsigned integer between min and max.
	 */
	uint get_rand_uint(float min, float max);

	/**
	 * @brief Remove an item from a vector.
	 * @param vec Pointer to the vector to modify.
	 * @param item The item to remove.
	 * @throws Exception if the item is not found in the vector.
	 */
	void remove(vector<int>* vec, int item);

	/**
	 * @brief Calculate Euclidean distance between two points.
	 * @param p1 The first point (x, y).
	 * @param p2 The second point (x, y).
	 * @return float The Euclidean distance between p1 and p2.
	 */
	float get_dist(pair<float, float> p1, pair<float, float> p2);

	/**
	 * @brief Add padding to a string based on version number.
	 * @param basestring The base string to add padding to.
	 * @param version The version number to use for padding calculation.
	 * @return string The padded string.
	 */
	std::string add_padding(std::string basestring, int version);

	/**
	 * @brief Join a vector of integers into a string with separators.
	 * @param numbers The vector of integers to join.
	 * @param separator The separator string to use between numbers.
	 * @return string The joined string.
	 */
	std::string join_as_string(vector<int> numbers, std::string separator);

	/**
	 * @brief Join a vector of floats into a string with separators.
	 * @param numbers The vector of floats to join.
	 * @param separator The separator string to use between numbers.
	 * @return string The joined string.
	 */
	std::string join_as_string(vector<float> numbers, std::string separator);

	/**
	 * @brief Join a vector of pairs into a string with separators.
	 * @param numbers The vector of pairs to join.
	 * @param separator The separator string to use between pairs.
	 * @return string The joined string.
	 */
	std::string join_as_string(vector<pair<int, int>> numbers, std::string separator);

	/**
	 * @brief Join a vector of strings into a string with separators.
	 * @param strings Pointer to the vector of strings to join.
	 * @param separator The separator string to use between strings.
	 * @return string The joined string.
	 */
	std::string join(vector<std::string>* strings, std::string separator);

	/**
	 * @brief Remove the largest vector from a collection of vectors.
	 * @param vectors Pointer to the vector of vectors to modify.
	 * @param max_size Reference to store the size of the removed vector.
	 */
	void remove_largest_vector(vector<vector<int>>* vectors, int& max_size);

	/**
	 * @brief Check if a string ends with a specific ending.
	 * @param full_string The string to check.
	 * @param ending The ending to look for.
	 * @return bool true if the string ends with the ending, false otherwise.
	 */
	bool ends_with(std::string full_string, std::string ending);

	/**
	 * @brief Check if two vectors have any overlapping elements.
	 * @param larger_vector Pointer to the larger vector to check.
	 * @param smaller_vector Pointer to the smaller vector to check.
	 * @return bool true if there are overlapping elements, false otherwise.
	 */
	bool have_overlap(vector<int>* larger_vector, vector<int>* smaller_vector);

	/**
	 * @brief Append elements from one vector to another (pointer version).
	 * @param result The target vector to append to.
	 * @param vec2 Pointer to the source vector to append from.
	 */
	void append_vector(vector<int>& result, vector<int>* vec2);

	/**
	 * @brief Append elements from one vector to another (value version).
	 * @param result The target vector to append to.
	 * @param vec2 The source vector to append from.
	 */
	void append_vector(vector<int>& result, vector<int> vec2);

	/**
	 * @brief Append elements from one pair vector to another (pointer version).
	 * @param result The target vector to append to.
	 * @param vec2 Pointer to the source vector to append from.
	 */
	void append_vector(vector<pair<int, int>>& result, vector<pair<int, int>>* vec2);

	/**
	 * @brief Append elements from one pair vector to another (value version).
	 * @param result The target vector to append to.
	 * @param vec2 The source vector to append from.
	 */
	void append_vector(vector<pair<int, int>>& result, vector<pair<int, int>> vec2);

	/**
	 * @brief Append elements from one string vector to another.
	 * @param result The target vector to append to.
	 * @param vec2 Pointer to the source vector to append from.
	 */
	void append_vector(vector<string>& result, vector<string>* vec2);
