#include "helpers.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <filesystem.>
#include <time.h>
#include <string.h>


using namespace std;

#define RAND_DOUBLE_PRECISION 0.0001
#define INV_RAND_DOUBLE_PRECISION_PLUSONE 1.0 / (1.0 + RAND_DOUBLE_PRECISION)

void help::init_RNG(int seed) {
    if (seed == -999) {
        // Seed the random number generator with the current time
		srand(time(NULL));
	}
	else {
        // Seed the random number generator with the given seed
		srand(seed);
	}
    int x = rand();
    int z = x * 3;
}

float help::get_rand_float(float min, float max) {
    return min + (float)rand() * INV_RAND_MAX * (max - min);
}

double help::_get_rand_double(double min, double max) {
    return min + ((double)rand() * INV_RAND_MAX) * (max - min);
}

double help::get_rand_double(double min, double max) {
    float range = max - min;
    double val = help::_get_rand_double(0, range) + help::_get_rand_double(0, range * RAND_DOUBLE_PRECISION);
    val = min + INV_RAND_DOUBLE_PRECISION_PLUSONE * val;
    return val;
}

uint help::get_rand_uint(int min, int max) {
    return round(help::get_rand_float(min, max));
}

int help::get_rand_int(int min, int max) {
    return round(help::get_rand_float(min, max));
}

bool help::is_in(std::vector<int>* vec, int item) {
    return find(vec->begin(), vec->end(), item) != vec->end();
}

bool is_in(std::map<int, int>& map, int item) {
	return map.find(item) != map.end();
}

// Add padding as suffix to given basestring
std::string help::add_padding(std::string basestring, int version) {
    int padding = 4;
    if (version > 9) {
        if (version > 99) {
            if (version > 999) {
                if (version > 9999) {
                    padding = 0;
                }
                else padding = 1;
            }
            else padding = 2;
        }
        else padding = 3;
    }
    string pad = "";
    for (int i = 0; i < padding; i++) {
        pad += "0";
    }
    return basestring + pad;
}

std::string help::zfill(std::string digits, int no_digits) {
	while (digits.size() < no_digits) {
        digits = "0" + digits;
	}
    return digits;
}

float help::dot(pair<float, float> &p1, pair<float, float> &p2) {
	return p1.first * p2.first + p1.second * p2.second;
}

void help::normalize(pair<float, float> &vec, float length) {
    vec = vec * (1.0f / length);
}

void help::get_random_unit_vector(pair<float, float>& vec) {
    vec = pair<float, float>(help::get_rand_float(-1, 1), help::get_rand_float(-1, 1));
    float squared_length = dot(vec, vec);
    while (squared_length > 1) {
        vec = pair<float, float>(help::get_rand_float(-1, 1), help::get_rand_float(-1, 1));
        squared_length = dot(vec, vec);
    }
    normalize(vec, sqrtf(squared_length));
}

void help::remove_from_vec(vector<int>* vec, int item) {
	auto position = find(vec->begin(), vec->end(), item);
	if (position != vec->end()) {
		vec->erase(position);
	}
}

void help::print_map(std::map<int, int>* map) {
    int i = 0;
    for (auto const& [key, val] : (*map))
    {
        if (i > 0) {
            std::cout << ", ";
        }
        std::cout << key        // string (key)
            << ':'
            << val;        // string's value
        i++;
    }
    if (i == 0) {
        std::cout << "<empty map>" << endl;
    }
    else {
        std::cout << endl;
    }
}

string help::join_as_string(vector<int> numbers, string separator) {
    string result = "";
    for (auto number : numbers) {
        result += to_string(number) + separator;
    }
    return result;
}

string help::join_as_string(vector<float> numbers, string separator) {
    string result = "";
    for (auto number : numbers) {
        result += to_string(number) + separator;
    }
    return result;
}

string help::join_as_string(vector<pair<int, int>> numbers, string separator) {
    string result = "";
    for (auto pair : numbers) {
        result += "(" + to_string(pair.first) + ", " + to_string(pair.second) + ")" + separator;
    }
    return result;
}

string help::join(vector<string>* strings, string separator) {
    string result = "";
    for (int i = 0; i < strings->size(); i++) {
        result += strings->at(i);
        if (i < strings->size() - 1) result += separator;
    }
    return result;
}

string help::join(vector<string> strings, string separator) {
    string result = "";
    for (int i = 0; i < strings.size(); i++) {
        result += strings.at(i);
        if (i < strings.size() - 1) result += separator;
    }
    return result;
}

void help::print_vector(std::vector<int>* vec) {
    for (int i = 0; i < vec->size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << vec->at(i);
    }
    std::cout << endl;
}

void help::print(std::string str) {
#if verbosity
    std::cout << str;
#endif
}

void help::print_pairs(std::vector<pair<int, int>>* pairs) {
    for (int i = 0; i < pairs->size(); i++) {
        if (i > 0) std::cout << " ";
        std::cout << "(";
        std::cout << to_string(pairs->at(i).first / 6) + ", " + to_string(pairs->at(i).first % 6);
        //cout << to_string(pairs->at(i).first) + ", " + to_string(pairs->at(i).second);
        std::cout << ")";
    }
    std::cout << endl;
}

string help::replace_occurrences(string basestring, string toReplace, string replaceWith) {
    string newstring = basestring;
    int pos = newstring.find(toReplace);
    while (pos != string::npos) {
        newstring.replace(pos, toReplace.size(), replaceWith);
        pos = newstring.find(toReplace);
    }

    return newstring;
}


vector<size_t> help::FindAll(string basestring, string target) {
    vector<size_t> occurrences;
    size_t found = 0;
    while (true) {
        found = basestring.find(target, found);
        if (found != string::npos) {
            occurrences.push_back(found);
        }
        else {
            break;
        }
        found++;
    }

    return occurrences;
}

bool help::is_in(string basestring, string target) {
    size_t found = 0;
    found = basestring.find(target, found);
    if (found == string::npos) {
        return false;
    }
    else return true;
}


double help::fisqrt(float n)
{
    float y = n;
    long i = *(long*)&y;
    i = 0x5f3759df - (i >> 1);
    y = *(float*)&i;

    return y * (1.5f - ((n * 0.5f) * y * y));
}


void help::increment_key(std::map<std::string, int>* _map, std::string key) {
    if (_map->find(key) == _map->end()) {
        (*_map)[key] = 1;
    }
    else {
        (*_map)[key]++;
    }
}

float help::get_value(std::map<std::string, float>* _map, std::string key) {
    if (_map == 0) return 0;
    if (_map->find(key) == _map->end()) {
        return 0;
    }
    else {
        return _map->at(key);
    }
}

int help::get_value(std::map<std::string, int>* _map, std::string key) {
    if (_map == 0) return 0;
    if (_map->find(key) == _map->end()) {
        return 0;
    }
    else {
        return _map->at(key);
    }
}

int help::get_value(std::map<int, int>* _map, int key) {
    if (_map == 0) return -1;
    map<int, int>::iterator it = _map->find(key);
    if (it == _map->end()) {
        return -1;
    }
    else {
        return it->second;
    }
}

int help::get_value(std::map<int, double>* _map, int key) {
    if (_map == 0) return -1;
    map<int, double>::iterator it = _map->find(key);
    if (it == _map->end()) {
        return -1;
    }
    else {
        return it->second;
    }
}

int help::get_value(std::map<uint32_t, uint32_t>* _map, uint32_t key) {
    if (_map == 0) return -1;
    std::map<uint32_t, uint32_t>::iterator it = _map->find(key);
    if (it == _map->end()) {
        return -1;
    }
    else {
        return it->second;
    }
}

int help::get_key(std::map<int, int>* _map, int value) {
    if (_map == 0) return -1;
    for (auto const& [key, val] : (*_map))
    {
        if (val == value) {
            return key;
        }
    }
    return -1;
}

void help::populate_with_zeroes(double* _array, int dim_x, int dim_y) {
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            _array[x * dim_y + y] = 0.0;
        }
    }
}

void help::populate_with_zeroes(uint* _array, int dim_x, int dim_y) {
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            _array[x * dim_y + y] = 0.0;
        }
    }
}

// Function to sort the map according
// to value in a (key-value) pairs
void help::sort(std::map<int, double>& _map, PairSet& _set)
{
    _set = PairSet(_map.begin(), _map.end());
}


void help::split(string basestring, string separator, vector<string>& substrings) {
    vector<size_t> occurrences = help::FindAll(basestring, separator);
    if (occurrences.size() == 0) {
        substrings = { basestring };
        return;
    }
    substrings.push_back(basestring.substr(0, occurrences[0]));
    for (int i = 0; i < occurrences.size() - 1; i++) {
        int name_length = occurrences[i + 1] - (occurrences[i] + separator.size());
        string name = basestring.substr(occurrences[i] + separator.size(), name_length);
        substrings.push_back(name);
    }
    string end = basestring.substr(occurrences[occurrences.size() - 1] + separator.size(), string::npos);
    if (end != "") {
        substrings.push_back(end);
    }
}

// Remove the largest item from the given vector
void help::remove_largest_vector(vector<vector<int>>* vectors, int& max_size) {
    max_size = 0;
    int largest_item_idx = -1;
    for (int i = 0; i < vectors->size(); i++) {
        if (vectors->at(i).size() > max_size) {
            max_size = vectors->at(i).size();
            largest_item_idx = i;
        }
    }
    vectors->erase(vectors->begin() + largest_item_idx);
}

void help::remove(vector<int>* vec, int item) {
    auto position = find(vec->begin(), vec->end(), item);
    if (position != vec->end()) {
        int idx = position - vec->begin();
        vec->erase(vec->begin() + idx);
    }
    else throw("Error: Cannot remove item from vector because it is not present\n");
}

template <typename T>
T help::pop(vector<T>* vec, int idx) {
    T item = vec->at(idx);
    vec->erase(vec->begin() + idx);
    return item;
}
template pair<int, int> help::pop<pair<int, int>>(vector<pair<int, int>>* vec, int idx);
template float help::pop<float>(vector<float>* vec, int idx);
template int help::pop<int>(vector<int>* vec, int idx);

//template <typename T, typename U>
//pair<T, U> pop(map<T, U>* map, int idx);

float help::cubed(float val) {
    return val * val * val;
}

bool help::ends_with(string full_string, string ending) {
    if (full_string.length() >= ending.length()) {
        return (0 == full_string.compare(full_string.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}

string help::readable_number(int number) {
    string billions = number >= (int)1e9 ? to_string(number / (int)1e9) + " " : "";

    string millions = number >= (int)1e6 > 0 ? to_string((number % (int)1e9) / (int)1e6) : "";
    if (billions != "") millions = zfill(millions, 3) + " ";
    else if (millions != "") millions += " ";

    string thousands = number >= (int)1e3 > 0 ? to_string((number % (int)1e6) / (int)1e3) : "";
    if (millions != "") thousands = zfill(thousands, 3) + " ";
    else if (thousands != "") thousands += " ";

    string rest = to_string(number % (int)1e3);
    if (thousands != "") rest = zfill(rest, 3);
    else rest;

    string reformatted_num = billions + millions + thousands + rest;
    return reformatted_num;
}

int help::microseconds_elapsed_since(high_resolution_clock::time_point start_time) {
    high_resolution_clock::time_point stop_time = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop_time - start_time);
    return duration.count();
}

// Measure the number of milliseconds elapsed since given time point
int help::milliseconds_elapsed_since(high_resolution_clock::time_point start_time) {
    high_resolution_clock::time_point stop_time = high_resolution_clock::now();
    return (duration_cast<milliseconds>(stop_time - start_time)).count();
}

// Measure the number of seconds elapsed since given time point
int help::seconds_elapsed_since(high_resolution_clock::time_point start_time) {
    high_resolution_clock::time_point stop_time = high_resolution_clock::now();
    return (duration_cast<seconds>(stop_time - start_time)).count();
}

float help::get_sigmoid(float x, float x_min, float x_max, float x_stretch) {
    float x2 = 2.0f * log10(x * 1e6f) - 5 - x_min;
    return 1.0 / (1 + exp(- x2 * x_stretch));
}

bool help::have_overlap(vector<int>* larger_vector, vector<int>* smaller_vector) {
    for (auto& item : *smaller_vector) {
        if (help::is_in(larger_vector, item)) return true;
    }
    return false;
}

void help::append_vector(vector<int>& result, vector<int>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

void help::append_vector(vector<int>& result, vector<int> vec2) {
    help::append_vector(result, &vec2);
}

void help::append_vector(vector<pair<int, int>>& result, vector<pair<int, int>>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

void help::append_vector(vector<pair<int, int>>& result, vector<pair<int, int>> vec2) {
    help::append_vector(result, &vec2);
}

void help::append_vector(vector<string>& result, vector<string>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

vector<float> help::get_free_memory() {
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    unsigned long long RAM_bytes = status.ullAvailPhys;
    unsigned long long VM_bytes = status.ullAvailVirtual;
    unsigned long long Pagefile_bytes = status.ullAvailPageFile;
    auto percent_memory = status.dwMemoryLoad;

    float RAM_gigabytes = (float)RAM_bytes / (float)(1 << 30);
    float VM_gigabytes = (float)VM_bytes / (float)(1 << 30);
    float Pagefile_gigabytes = (float)Pagefile_bytes / (float)(1 << 30);
    return { RAM_gigabytes, VM_gigabytes, Pagefile_gigabytes, (float)percent_memory };
}

double help::get_mean(vector<double>* distribution) {
    double sum = 0;
    for (double& sample : *distribution) {
        sum += sample;
    }
    return sum / distribution->size();
}

double help::get_stdev(vector<double>* distribution, double mean) {
    if (mean == -999999) {
        mean = help::get_mean(distribution);
    }
    double sos = 0;
    for (auto& sample : *distribution) sos += (sample - mean) * (sample - mean);
    double variance = sos / (distribution->size() - 1);
    return sqrt(variance);
}


double help::get_max(vector<double>* distribution) {
    double max = -INFINITY;
    for (double value : *distribution) {
        if (value > max) max = value;
    }
    return max;
}

template <typename T>
double help::get_min(vector<T>* distribution) {
    double min = INFINITY;
    for (auto value : *distribution) {
        if (value < min) min = value;
    }
    return min;
}
template double help::get_min<double>(vector<double>* vec);
template double help::get_min<float>(vector<float>* vec);

float help::get_dist(pair<float, float> p1, pair<float, float> p2) {
    float xdif = (p1.first - p2.first);
    float ydif = (p1.second - p2.second);
    return sqrtf(xdif * xdif + ydif * ydif);
}

float help::get_manhattan_dist(pair<float, float> p1, pair<float, float> p2) {
	return abs(p1.first - p2.first) + abs(p1.second - p2.second);
}

float help::exponential_function(float x, float a, float b, float c) {
    return a * (x * x) + b * x + c;
}

float help::get_lowest_solution_for_quadratic(float y, float a, float b, float c) {
    c -= y; // Make rhs equal to 0
    float result1 = (-b - sqrtf(b * b - 4.0 * a * c)) / (2.0 * a);
    float result2 = (-b + sqrtf(b * b - 4.0 * a * c)) / (2.0 * a);
    return (result1 < result2 ? result1 : result2);
}

help::LinearProbabilityModel::LinearProbabilityModel() = default;
help::LinearProbabilityModel::LinearProbabilityModel(float _q1, float _q2, float _min, float _max) {
    q1 = _q1;
    q2 = _q2;
    min = _min;
    max = _max;

    float xrange = (max - min);
    a = (q2 - q1) / xrange;
    b = q1;
    float cdf_of_xrange = 0.5 * a * xrange * xrange + b * xrange; // CDF(xrange)

    // Scale PDF such that it integrates to 1, i.e., so that CDF(xrange) = 1
    a /= cdf_of_xrange;
    b /= cdf_of_xrange;
}
float help::LinearProbabilityModel::linear_sample() {
    float cdf_of_dx = help::get_rand_float(0, 1); // Obtain CDF(dx)
    float dx = get_lowest_solution_for_quadratic(cdf_of_dx, a / 2.0, b, 0); // Get corresponding value dx
    return min + dx;
}

help::ProbModelPiece::ProbModelPiece() = default;
help::ProbModelPiece::ProbModelPiece(float _xmin, float _xmax, float _ymin, float _ymax) {
    xmin = _xmin; xmax = _xmax;
    ymin = _ymin; ymax = _ymax;
    ysize = ymax - ymin;
    xsize = xmax - xmin;
};
float help::ProbModelPiece::intersect(float cdf_y) {
    if (cdf_y >= ymin && cdf_y <= ymax) {
        float y_diff = cdf_y - ymin;
        return xmin + (y_diff / ysize) * xsize;
    }
    return -1;
};
void help::ProbModelPiece::rescale(float factor) {
    ymin *= factor; ymax *= factor;
    ysize = ymax - ymin;
}


help::DiscreteProbModelPiece::DiscreteProbModelPiece() = default;
help::DiscreteProbModelPiece::DiscreteProbModelPiece(int _idx, float _ymin, float _ymax) {
    idx = _idx;
    ymin = _ymin; ymax = _ymax;
    ysize = ymax - ymin;
};
int help::DiscreteProbModelPiece::intersect(float cdf_y) {
    if (cdf_y >= ymin && cdf_y <= ymax) {
        return idx;
    }
    return -1;
};
void help::DiscreteProbModelPiece::rescale(float factor) {
    ymin *= factor; ymax *= factor;
    ysize = ymax - ymin;
}

template <typename T>
int help::binary_search(T* arr, int size, T target) {
    int l = 0;
    int r = size - 1;
    while (l <= r) {
        int m = l + (r - l) / 2;

        // Check if target is present at mid
        float next_element;
        if (m < size - 1) {
            next_element = arr[m + 1];
        }
        else {
            next_element = INFINITY;
        }
        if (target >= arr[m] && target <= next_element) {
            return m;
        }   

        // If target greater, ignore left half
        if (arr[m] < target)
            l = m + 1;

        // If target is smaller, ignore right half
        else
            r = m - 1;
    }

    // If we reach here, then element was not present
    return -1;
}
template int help::binary_search<double>(double* arr, int size, double target);
template int help::binary_search<float>(float* arr, int size, float target);

template <typename T>
int help::do_linear_search(T* arr, int size, T target) {
	for (int i = 0; i < size; i++) {
		if (arr[i] >= target) {
			return i;
		}
        /*else if (target < 0.48f && target > 0.47f && i >= 470000) {
            printf("Arr %i (%f) is not greater than target %f\n", i, arr[i], target);
        }*/
	}
	return -1;
}
template int help::do_linear_search<double>(double* arr, int size, double target);
template int help::do_linear_search<float>(float* arr, int size, float target);


help::PieceWiseLinearProbModel::PieceWiseLinearProbModel() = default;
help::PieceWiseLinearProbModel::PieceWiseLinearProbModel(float _xmax) {
    xmax = _xmax;
    piece_width = xmax / (float)resolution;
}
float help::PieceWiseLinearProbModel::pdf(float x) {
    return 0.01 * x * x; // Dummy return value; OVERRIDE THIS FUNCTION TO USE YOUR CUSTOM PDF.
}
void help::PieceWiseLinearProbModel::build() {
    cdf = make_shared<double[]>(resolution);
    float _xmin = 0.0f;
    float _xmax = piece_width;
    float cdf_maximum = 0.0f;
    float pdf_begin = pdf(0);
    float pdf_end = pdf(_xmax);
    int i = 0;
    while (true) {
        float piece_area = min(pdf_begin, pdf_end) * piece_width + 0.5f * abs(pdf_begin - pdf_end) * piece_width;
        float _cdf_min = cdf_maximum;
        cdf_maximum += piece_area;
        cdf[i] = _cdf_min;
        i++;
        if (i >= resolution) break;
        _xmin += piece_width;
        _xmax += piece_width;
        pdf_begin = pdf_end;
        pdf_end = pdf(_xmax);
    }

    // Rescale CDF pieces to ensure CDF maximum is exactly 1
    float scale_factor = 1.0f / cdf_maximum;
    for (int i = 0; i < resolution; i++) {
        cdf[i] *= scale_factor;
    }

    built = 1;
}
float help::PieceWiseLinearProbModel::sample() {
    double cdf_sample = help::get_rand_float(0.0f, 1.0f);
    float x = binary_search(cdf.get(), resolution, cdf_sample);
    x *= piece_width;
    return x;
}

void help::get_normal_distributed_direction(pair<float, float>& direction, float mean_direction, float direction_stdev) {
    NormalProbModel prob_model(mean_direction, direction_stdev);
    float phi = prob_model.get_normal_distr_sample();
    direction = pair<float, float>(cos(phi), sin(phi));
}


