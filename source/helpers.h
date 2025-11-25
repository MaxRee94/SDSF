#pragma once
#include <vector>
#include <map>
#include <cstdlib>
#include <set>
#include <Windows.h>
#include <string>
#include <iostream>
#include <queue>
#include <assert.h>
#include <unordered_map>
#include <random>
#include <chrono>
#include <numeric>

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <limits>
#include <concepts>

#include "timer.h"

using namespace std;
#define INV_RAND_MAX  1.0 / RAND_MAX


using namespace std::chrono;
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

typedef std::set < std::pair<int, int>, _comparator> PairIntSet;

template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U>& l, const std::pair<T, U>& r) {
	return { l.first + r.first,l.second + r.second };
}

template <typename T, typename U, typename V, typename W>
std::pair<T, U> operator+(const std::pair<T, U>& l, const std::pair<V, W>& r) {
	return { l.first + r.first,l.second + r.second };
}

template <typename T, typename U>
std::pair<T, U> operator-(const std::pair<T, U>& l, const std::pair<T, U>& r) {
	return { l.first - r.first,l.second - r.second };
}

template <typename T, typename U>
std::pair<T, U> operator*(const float& s, const std::pair<T, U>& p) {
	return { s * p.first, s * p.second };
}

template <typename T, typename U>
std::pair<T, U> operator*(const std::pair<T, U>& p, const float& s) {
	return s * p;
} 

template<std::floating_point value_t> [[nodiscard]]
constexpr bool is_float_equal(value_t l, value_t r)
{
	constexpr auto infinity = std::numeric_limits<value_t>::infinity();

	auto const min = std::nextafter(l, -infinity);
	auto const max = std::nextafter(l, infinity);
	return (min <= r && r <= max);
}

template<std::floating_point value_t> [[nodiscard]]
	constexpr bool approx(value_t l, value_t r, value_t precision)
	{
		auto const min = l - precision;
		auto const max = l + precision;
		return (min <= r && r <= max);
	}

namespace help {


	void init_RNG(int seed);

	void sort(std::map<int, double>& _map, PairSet& _set);

	void save_image(string name, shared_ptr<float[]> image);

	void sort(std::map<int, int>& _map, PairIntSet& _set);

	int get_key(std::map<int, int>* _map, int value);

	void populate_with_zeroes(double* _array, int dim_x, int dim_y);
		
	void populate_with_zeroes(uint* _array, int dim_x, int dim_y);

	void split(std::string basestring, std::string separator, vector<std::string>& substrings);

	float dot(pair<float, float> &p1, pair<float, float> &p2);

	vector<pair<int, int>> get_bbox(vector<pair<int, int>> positions2d);

	void normalize(pair<float, float>& vec, float length);

	void get_random_unit_vector(pair<float, float> &direction);

	int get_random_key(std::map<int, float>& map);

	void remove_from_vec(vector<int>* vec, int item);

	string readable_number(int number);

	void get_normal_distributed_direction(pair<float, float>& direction, float mean_direction, float direction_stdev);

	//Return whether the given vector <vec> contains the integer <item>
	bool is_in(std::vector<int>* vec, int item);

	bool is_in(std::map<int, int>& map, int item);

	bool is_in(std::map<int, float>& map, int item);

	bool is_in(std::string basestring, std::string target);

	template <typename T>
	vector<int> get_keys(map<int, T>& map);

	void print_map(std::map<int, int>* map);
	
	void print_map(std::map<int, float>* map);

	void print_vector(std::vector<int>* vec);

	void print_pairs(std::vector<pair<int, int>>* vec);

	void print(std::string);

	vector<size_t> FindAll(std::string basestring, std::string target);

	std::string replace_occurrences(std::string basestring, std::string toReplace, std::string replaceWith);

	double fisqrt(float n);

	void increment_key(std::map<std::string, int>* map, std::string key);

	float get_value(std::map<std::string, float>* map, std::string key);

	int get_value(std::map<std::string, int>* map, std::string key);

	int get_value(std::map<int, int>* map, int key);

	int get_value(std::map<int, double>* map, int key);

	float cubed(float val);

	int get_value(std::map<uint32_t, uint32_t>* map, uint32_t key);

	float get_rand_float(float min, float max);

	double _get_rand_double(double min, double max);

	double get_rand_double(double min, double max);

	uint get_rand_uint(int min, int max);

	int get_rand_int(int min, int max);

	void remove(vector<int>* vec, int item);

	float get_sigmoid(float x, float x_min, float x_max, float x_stretch);

	template <typename T>
	T pop(vector<T>* vec, int idx);

	template <typename T, typename U>
	pair<T, U> pop(map<T, U>* map, int idx);

	float get_dist(pair<float, float> p1, pair<float, float> p2);

	float get_manhattan_dist(pair<float, float> p1, pair<float, float> p2);

	std::string add_padding(std::string basestring, int version);
	
	std::string zfill(std::string basestring, int no_digits);

	std::string join_as_string(vector<int> numbers, std::string separator);

	std::string join_as_string(vector<float> numbers, std::string separator);

	std::string join_as_string(vector<pair<int, int>> numbers, std::string separator);

	std::string join(vector<std::string>* strings, std::string separator);

	std::string join(vector<std::string> strings, std::string separator);

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
	
	template <typename T>
	T get_max(vector<pair<T, T>>& distribution, int index);

	// Do binary search in array of floats or doubles
	template <typename T>
	int binary_search(T* arr, int size, T target);

	// Measure the number of microseconds elapsed since given time point
	int microseconds_elapsed_since(high_resolution_clock::time_point start_time);

	// Measure the number of milliseconds elapsed since given time point
	int milliseconds_elapsed_since(high_resolution_clock::time_point start_time);

	// Measure the number of seconds elapsed since given time point
	int seconds_elapsed_since(high_resolution_clock::time_point start_time);

	// Do linear search in array of floats or doubles
	template <typename T>
	int do_linear_search(T* arr, int size, T target);

	// Get minimum
	template <typename T>
	double get_min(vector<T>* distribution);

	template <typename T>
	T get_min(vector<pair<T, T>>& distribution, int index);

	// Get exponential function value
	float exponential_function(float x, float a, float b, float c);

	// Solve for x in quadratic function f(x), given coefficients (a and b), constant c, and value y = f(x)
	float get_lowest_solution_for_quadratic(float y, float a, float b, float c);

	class LinearProbabilityModel {
	public:
		LinearProbabilityModel();
		LinearProbabilityModel(float _q1, float _q2, float _min, float _max);
		float linear_sample();
		float q1 = 0;
		float q2 = 0;
		float min = 0;
		float max = 0;
		float a = 0;
		float b = 0;
	};

	class ProbModelPiece {
	public:
		ProbModelPiece();
		ProbModelPiece(float _xmin, float _xmax, float _ymin, float _ymax);
		float intersect(float cdf_y);
		void rescale(float factor);
		float xmin = 0;
		float xmax = 0;
		float ymin = 0;
		float ymax = 0;
		float ysize = 0;
		float xsize = 0;
	};

	class DiscreteProbModelPiece {
	public:
		DiscreteProbModelPiece();
		DiscreteProbModelPiece(int idx, float _ymin, float _ymax);
		int intersect(float cdf_y);
		void rescale(float factor);
		int idx = 0;
		float ymin = 0;
		float ymax = 0;
		float ysize = 0;
	};

	class PieceWiseLinearProbModel {
	public:
		PieceWiseLinearProbModel();
		PieceWiseLinearProbModel(float _xmax);
		virtual float pdf(float x);
		virtual float sample();
		virtual void build();
		float xmax = 0;
		float resolution = 1000.0f;
		float piece_width = 0.0f;
		shared_ptr<double[]> cdf = 0;
		bool built = false;
	};

	class NormalProbModel {
	public:
		NormalProbModel() = default;
		NormalProbModel(float mean, float stdev) {
			distribution = normal_distribution<float>(mean, stdev);
			generator = default_random_engine(rand());
		};
		virtual float get_normal_distr_sample() {
			return distribution(generator);
		}
		normal_distribution<float> distribution;
		default_random_engine generator;
	};

	class UniformProbModel {
	public:
		UniformProbModel() = default;
		UniformProbModel(float _min, float _max) {
			min = _min; max = _max;
		};
		virtual float get_uniform_sample() {
			return help::get_rand_float(min, max);
		}
		float min = 0;
		float max = 0;
	};

	class ConstantModel {
	public:
		ConstantModel() = default;
		ConstantModel(float _constant) {
			constant = _constant;
		};
		virtual float get_constant() {
			return constant;
		}
		float constant = 0;
	};

	class DiscreteProbabilityModel{
	public:
		DiscreteProbabilityModel() = default;
		DiscreteProbabilityModel(int _size) {
			size = _size;
			probabilities = std::make_shared<double[]>(_size);
			cdf = std::make_shared<double[]>(_size);
			id = rand();
		};
		/*~DiscreteProbabilityModel() {
			free();
		}*/
		/*void free() {
			printf("Freeing prob model %i \n", id);
			if (probabilities != nullptr) {
				printf("probs: %i\n", probabilities[0]);
				printf("probs: %i (ptr)\n", probabilities);
				delete[] probabilities;
				probabilities = nullptr;

				printf("probs after delete: %i\n", probabilities);
			}
			if (cdf != nullptr) {
				delete[] cdf;
				cdf = nullptr;
			}
		}*/
		void build_cdf() {
			double height = 0.0f;
			for (int i = 0; i < size; i++) {
				cdf[i] = height;
				height += probabilities[i];
			}
		}
		int sample() {
			double cdf_sample = help::get_rand_double(0.0f, 1.0);
			int idx = binary_search(cdf.get(), size, cdf_sample);
			if (idx != -1) return idx;
			else return uniform_rand_idx();
		}
		int uniform_rand_idx() {
			return help::get_rand_int(0, size - 1);
		}
		void set_probabilities(double* probs, double &integral) {
			integral = 0;
			for (int i = 0; i < size; i++) {
				probabilities[i] = probs[i];
				integral += probs[i];
			}
		}
		void set_probabilities(shared_ptr<float[]> probs, float& integral) {
			integral = 0;
			for (int i = 0; i < size; i++) {
				probabilities[i] = probs[i];
				integral += probs[i];
			}
		}
		void normalize(double integral) {
			double recipr = 1.0 / integral;
			for (int i = 0; i < size; i++) {
				probabilities[i] *= recipr;
			}
		}
		shared_ptr<double[]> probabilities = 0;
		shared_ptr<double[]> cdf = 0;
		int size = 0;
		int id = 0;
	};

	class SmallDiscreteProbabilityModel {
	public:
		SmallDiscreteProbabilityModel() = default;
		SmallDiscreteProbabilityModel(float _prob0, float _prob1, float _prob2) {
			prob0 = _prob0;
			prob1 = _prob1;
			prob2 = _prob2;
		};
		virtual float pick_outcome() {
			float uniform_sample = help::get_rand_float(0, 1);
			if (uniform_sample < prob0) {
				return 0;
			}
			else if (uniform_sample < (prob0 + prob1)) {
				return 1;
			}
			else {
				return 2;
			}
		}
		float prob0 = 0;
		float prob1 = 0;
		float prob2 = 0;
	};

	class GammaProbModel {
	public:
		GammaProbModel() = default;
		GammaProbModel(float shape, float scale) {
			generator = default_random_engine(rand());
			distribution = std::gamma_distribution<float>(shape, scale);
		};
		virtual float get_gamma_sample() {
			return distribution(generator);
		}
		gamma_distribution<float> distribution;
		default_random_engine generator;
	};

	class ProbModel : public SmallDiscreteProbabilityModel, public LinearProbabilityModel, public NormalProbModel, public UniformProbModel {
	public:
		ProbModel() = default;
		ProbModel(float _prob0, float _prob1, float _prob2) : SmallDiscreteProbabilityModel(_prob0, _prob1, _prob2) { type = "discrete"; };
		ProbModel(float _min, float _max) : UniformProbModel(_min, _max) { type = "uniform"; };
		ProbModel(float _mean, float _stdev, float _min, float _max, int dummy) : NormalProbModel(_mean, _stdev) {
			type = "normal";
			min_value = _min;
			max_value = _max;
		};
		ProbModel(float _q1, float _q2, float _min, float _max) : LinearProbabilityModel(_q1, _q2, _min, _max) { type = "linear"; };
		ProbModel(float _constant) { type = "constant"; constant = _constant; };
		float sample() {
			if (type == "uniform") {
				return UniformProbModel::get_uniform_sample();
			}
			else if (type == "normal") {
				float _sample = NormalProbModel::get_normal_distr_sample();
				int i = 0;
				while (true) {
					if (_sample >= min_value && _sample <= max_value) {
						return _sample;
					}
					if (i > 100) printf("Warning: Normal distribution sampling is requiring a lot of resamples (%i) to get a value within the given range\n", i);
					_sample = NormalProbModel::get_normal_distr_sample();
					i++;
				}
			}
			else if (type == "linear") {
				return LinearProbabilityModel::linear_sample();
			}
			else if (type == "discrete") {
				return SmallDiscreteProbabilityModel::pick_outcome();
			}
			else if (type == "constant") {
				return constant;
			}
			else {
				return -999;
			}
		}
		string type = "none";
		float min_value = 0;
		float max_value = 0;
		float constant = 0;
	};
};
