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

// Comparison function for sorting the
// set by increasing order of its pair's
// second value
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

	void init_RNG();

	void sort(std::map<int, double>& _map, PairSet& _set);

	int get_key(std::map<int, int>* _map, int value);

	void populate_with_zeroes(double* _array, int dim_x, int dim_y);
		
	void populate_with_zeroes(uint* _array, int dim_x, int dim_y);

	void split(std::string basestring, std::string separator, vector<std::string>& substrings);

	//Return whether the given vector <vec> contains the integer <item>
	bool is_in(std::vector<int>* vec, int item);

	void print_map(std::map<int, int>* map);

	void print_vector(std::vector<int>* vec);

	void print_pairs(std::vector<pair<int, int>>* vec);

	void print(std::string);

	vector<size_t> FindAll(std::string basestring, std::string target);

	bool is_in(std::string basestring, std::string target);

	std::string replace_occurrences(std::string basestring, std::string toReplace, std::string replaceWith);

	double fisqrt(float n);

	void increment_key(std::map<std::string, int>* map, std::string key);

	float get_value(std::map<std::string, float>* map, std::string key);

	int get_value(std::map<std::string, int>* map, std::string key);

	int get_value(std::map<int, int>* map, int key);

	int get_value(std::map<int, double>* map, int key);

	int get_value(std::map<uint32_t, uint32_t>* map, uint32_t key);

	float get_rand_float(float min, float max);

	uint get_rand_uint(float min, float max);

	void remove(vector<int>* vec, int item);

	std::string add_padding(std::string basestring, int version);

	std::string join_as_string(vector<int> numbers, std::string separator);

	std::string join_as_string(vector<float> numbers, std::string separator);

	std::string join_as_string(vector<pair<int, int>> numbers, std::string separator);

	std::string join(vector<std::string>* strings, std::string separator);

	void remove_largest_vector(vector<vector<int>>* vectors, int& max_size);

	bool ends_with(std::string full_string, std::string ending);

	bool have_overlap(vector<int>* larger_vector, vector<int>* smaller_vector);

	// Push back the items in vec2 to the vector <result>
	void append_vector(vector<int>& result, vector<int>* vec2);

	// Push back the items in vec2 to the vector <result>
	void append_vector(vector<int>& result, vector<int> vec2);

	// Push back the items in vec2 to the vector <result>
	void append_vector(vector<pair<int, int>>& result, vector<pair<int, int>>* vec2);

	// Push back the items in vec2 to the vector <result>
	void append_vector(vector<pair<int, int>>& result, vector<pair<int, int>> vec2);

	// Push back the items in vec2 to the vector <result>
	void append_vector(vector<string>& result, vector<string>* vec2);

	// Get free RAM memory
	vector<float> get_free_memory();

	// Get stdev
	double get_stdev(vector<double>* distribution, double mean = -999999);

	// Get mean
	double get_mean(vector<double>* distribution);

	// Get maximum
	double get_max(vector<double>* distribution);

	// Get maximum
	double get_min(vector<double>* distribution);
};
