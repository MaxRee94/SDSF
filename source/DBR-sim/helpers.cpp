#include "helpers.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <filesystem.>
#include <time.h>
#include <string.h>

using namespace std;

float INV_RAND_MAX = 1.0 / (float)RAND_MAX;

/**
 * @brief Initialize the random number generator with current time seed.
 * 
 * Uses srand() with the current time to seed the random number generator.
 * Also performs an initial random number generation to improve randomness.
 */
void help::init_RNG() {
    srand(time(NULL));
    int x = rand();
    int z = x * 3;
}

/**
 * @brief Generate a random float in the range [min, max).
 * 
 * @param min The minimum value (inclusive).
 * @param max The maximum value (exclusive).
 * @return float A random float in the specified range.
 */
float help::get_rand_float(float min, float max) {
    return min + (float)rand() * INV_RAND_MAX * (max - min);
}

/**
 * @brief Generate a random unsigned integer in the range [min, max].
 * 
 * @param min The minimum value (inclusive).
 * @param max The maximum value (inclusive).
 * @return uint A random unsigned integer in the specified range.
 */
uint help::get_rand_uint(float min, float max) {
    float float_rand_range = (float)rand() * INV_RAND_MAX * (max - min);
    return round(min + float_rand_range);
}

/**
 * @brief Check if a vector contains a specific integer.
 * 
 * @param vec Pointer to the vector to search.
 * @param item The integer to search for.
 * @return bool true if the item is found, false otherwise.
 */
bool help::is_in(std::vector<int>* vec, int item) {
    return find(vec->begin(), vec->end(), item) != vec->end();
}

/**
 * @brief Add padding to a string based on version number.
 * 
 * Adds zero-padding as a suffix to the base string. The amount of padding
 * depends on the version number: 4 digits for versions < 10, 3 digits for
 * versions 10-99, 2 digits for 100-999, 1 digit for 1000-9999, and 0 digits for >= 10000.
 * 
 * @param basestring The base string to add padding to.
 * @param version The version number to determine padding length.
 * @return string The padded string.
 */
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

/**
 * @brief Print the contents of a map to standard output.
 * 
 * @param map Pointer to the map to print.
 */
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

/**
 * @brief Join a vector of integers into a string with separators.
 * 
 * @param numbers The vector of integers to join.
 * @param separator The separator string to use between numbers.
 * @return string The joined string with all numbers separated by the separator.
 */
string help::join_as_string(vector<int> numbers, string separator) {
    string result = "";
    for (auto number : numbers) {
        result += to_string(number) + separator;
    }
    return result;
}

/**
 * @brief Join a vector of floats into a string with separators.
 * 
 * @param numbers The vector of floats to join.
 * @param separator The separator string to use between numbers.
 * @return string The joined string with all numbers separated by the separator.
 */
string help::join_as_string(vector<float> numbers, string separator) {
    string result = "";
    for (auto number : numbers) {
        result += to_string(number) + separator;
    }
    return result;
}

/**
 * @brief Join a vector of pairs into a string with separators.
 * 
 * @param numbers The vector of pairs to join.
 * @param separator The separator string to use between pairs.
 * @return string The joined string with all pairs formatted as (x,y) and separated by the separator.
 */
string help::join_as_string(vector<pair<int, int>> numbers, string separator) {
    string result = "";
    for (auto pair : numbers) {
        result += "(" + to_string(pair.first) + ", " + to_string(pair.second) + ")" + separator;
    }
    return result;
}

/**
 * @brief Join a vector of strings into a string with separators.
 * 
 * @param strings Pointer to the vector of strings to join.
 * @param separator The separator string to use between strings.
 * @return string The joined string with all strings separated by the separator.
 */
string help::join(vector<string>* strings, string separator) {
    string result = "";
    for (int i = 0; i < strings->size(); i++) {
        result += strings->at(i);
        if (i < strings->size() - 1) result += separator;
    }
    return result;
}

/**
 * @brief Print the contents of a vector to standard output.
 * 
 * @param vec Pointer to the vector to print.
 */
void help::print_vector(std::vector<int>* vec) {
    for (int i = 0; i < vec->size(); i++) {
        if (i > 0) std::cout << ", ";
        std::cout << vec->at(i);
    }
    std::cout << endl;
}

/**
 * @brief Print a string (only if VERBOSE is true).
 * 
 * @param str The string to print.
 */
void help::print(std::string str) {
#if VERBOSE
    std::cout << str;
#endif
}

/**
 * @brief Print the contents of a vector of pairs to standard output.
 * 
 * @param pairs Pointer to the vector of pairs to print.
 */
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

/**
 * @brief Replace all occurrences of a substring with another string.
 * 
 * @param basestring The original string.
 * @param toReplace The substring to replace.
 * @param replaceWith The string to replace with.
 * @return string The modified string with all replacements made.
 */
string help::replace_occurrences(string basestring, string toReplace, string replaceWith) {
    string newstring = basestring;
    int pos = newstring.find(toReplace);
    while (pos != string::npos) {
        newstring.replace(pos, toReplace.size(), replaceWith);
        pos = newstring.find(toReplace);
    }

    return newstring;
}


/**
 * @brief Find all occurrences of a substring in a string.
 * 
 * @param basestring The string to search in.
 * @param target The substring to find.
 * @return vector<size_t> Vector of positions where the target is found.
 */
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

/**
 * @brief Check if a string contains a substring.
 * 
 * @param basestring The string to search in.
 * @param target The substring to find.
 * @return bool true if the substring is found, false otherwise.
 */
bool help::is_in(string basestring, string target) {
    size_t found = 0;
    found = basestring.find(target, found);
    if (found == string::npos) {
        return false;
    }
    else return true;
}


/**
 * @brief Compute fast inverse square root using magic number approximation.
 * 
 * Uses a well-known fast inverse square root algorithm for performance-critical
 * applications. This is an approximate method that provides good performance
 * at the cost of some precision.
 * 
 * @param n The number to compute the inverse square root of.
 * @return float The approximate inverse square root of n.
 */
double help::fisqrt(float n)
{
    float y = n;
    long i = *(long*)&y;
    i = 0x5f3759df - (i >> 1);
    y = *(float*)&i;

    return y * (1.5f - ((n * 0.5f) * y * y));
}


/**
 * @brief Increment a key's value in a string-keyed integer map.
 * 
 * If the key doesn't exist, it's created with value 1.
 * If the key exists, its value is incremented by 1.
 * 
 * @param _map Pointer to the map to modify.
 * @param key The key to increment.
 */
void help::increment_key(std::map<std::string, int>* _map, std::string key) {
    if (_map->find(key) == _map->end()) {
        (*_map)[key] = 1;
    }
    else {
        (*_map)[key]++;
    }
}

/**
 * @brief Get a float value from a string-keyed map.
 * 
 * @param _map Pointer to the map to search.
 * @param key The key to look up.
 * @return float The value for the key, or 0 if the map is null or key not found.
 */
float help::get_value(std::map<std::string, float>* _map, std::string key) {
    if (_map == 0) return 0;
    if (_map->find(key) == _map->end()) {
        return 0;
    }
    else {
        return _map->at(key);
    }
}

/**
 * @brief Get an integer value from a string-keyed map.
 * 
 * @param _map Pointer to the map to search.
 * @param key The key to look up.
 * @return int The value for the key, or 0 if the map is null or key not found.
 */
int help::get_value(std::map<std::string, int>* _map, std::string key) {
    if (_map == 0) return 0;
    if (_map->find(key) == _map->end()) {
        return 0;
    }
    else {
        return _map->at(key);
    }
}

/**
 * @brief Get an integer value from an int-keyed map.
 * 
 * @param _map Pointer to the map to search.
 * @param key The key to look up.
 * @return int The value for the key, or -1 if the map is null or key not found.
 */
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

/**
 * @brief Get an integer value from an int-keyed double map.
 * 
 * Note: Returns int value even for double-keyed map.
 * 
 * @param _map Pointer to the map to search.
 * @param key The key to look up.
 * @return int The value for the key, or -1 if the map is null or key not found.
 */
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

/**
 * @brief Get an integer value from a uint32_t-keyed map.
 * 
 * Note: Returns int value even for uint32_t values.
 * 
 * @param _map Pointer to the map to search.
 * @param key The key to look up.
 * @return int The value for the key, or -1 if the map is null or key not found.
 */
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

/**
 * @brief Get the key for a given value in an int-keyed map.
 * 
 * @param _map Pointer to the map to search.
 * @param value The value to find the key for.
 * @return int The key corresponding to the value, or -1 if not found or map is null.
 */
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

/**
 * @brief Fill a 2D double array with zeros.
 * 
 * @param _array Pointer to the array to populate.
 * @param dim_x The x dimension size.
 * @param dim_y The y dimension size.
 */
void help::populate_with_zeroes(double* _array, int dim_x, int dim_y) {
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            _array[x * dim_y + y] = 0.0;
        }
    }
}

/**
 * @brief Fill a 2D unsigned int array with zeros.
 * 
 * @param _array Pointer to the array to populate.
 * @param dim_x The x dimension size.
 * @param dim_y The y dimension size.
 */
void help::populate_with_zeroes(uint* _array, int dim_x, int dim_y) {
    for (int x = 0; x < dim_x; x++) {
        for (int y = 0; y < dim_y; y++) {
            _array[x * dim_y + y] = 0.0;
        }
    }
}

// Function to sort the map according
// to value in a (key-value) pairs
/**
 * @brief Sort a map into a set of pairs ordered by value.
 * 
 * Creates a PairSet from a map, sorted by the values in ascending order.
 * 
 * @param _map The map to sort.
 * @param _set The set to populate with sorted pairs.
 */
void sort(std::map<int, double>& _map, PairSet& _set)
{
    _set = PairSet(_map.begin(), _map.end());
}


/**
 * @brief Split a string into substrings using a separator.
 * 
 * @param basestring The string to split.
 * @param separator The separator string to split on.
 * @param substrings Reference to vector to store the resulting substrings.
 */
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
/**
 * @brief Remove the largest vector from a collection of vectors.
 * 
 * Finds the vector with the most elements, removes it from the collection,
 * and sets the max_size parameter to its size.
 * 
 * @param vectors Pointer to the vector of vectors to modify.
 * @param max_size Reference to store the size of the removed vector.
 */
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

/**
 * @brief Remove an item from a vector.
 * 
 * @param vec Pointer to the vector to modify.
 * @param item The item to remove.
 * @throws Exception if the item is not found in the vector.
 */
void help::remove(vector<int>* vec, int item) {
    auto position = find(vec->begin(), vec->end(), item);
    if (position != vec->end()) {
        int idx = position - vec->begin();
        vec->erase(vec->begin() + idx);
    }
    else throw("Error: Cannot remove item from vector because it is not present\n");
}

/**
 * @brief Check if a string ends with a specific ending.
 * 
 * @param full_string The string to check.
 * @param ending The ending to look for.
 * @return bool true if the string ends with the ending, false otherwise.
 */
bool help::ends_with(string full_string, string ending) {
    if (full_string.length() >= ending.length()) {
        return (0 == full_string.compare(full_string.length() - ending.length(), ending.length(), ending));
    }
    else {
        return false;
    }
}

/**
 * @brief Check if two vectors have any overlapping elements.
 * 
 * @param larger_vector Pointer to the larger vector to check.
 * @param smaller_vector Pointer to the smaller vector to check.
 * @return bool true if there are overlapping elements, false otherwise.
 */
bool help::have_overlap(vector<int>* larger_vector, vector<int>* smaller_vector) {
    for (auto& item : *smaller_vector) {
        if (help::is_in(larger_vector, item)) return true;
    }
    return false;
}

/**
 * @brief Append elements from one vector to another (pointer version).
 * 
 * @param result The target vector to append to.
 * @param vec2 Pointer to the source vector to append from.
 */
void help::append_vector(vector<int>& result, vector<int>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

/**
 * @brief Append elements from one vector to another (value version).
 * 
 * @param result The target vector to append to.
 * @param vec2 The source vector to append from.
 */
void help::append_vector(vector<int>& result, vector<int> vec2) {
    help::append_vector(result, &vec2);
}

/**
 * @brief Append elements from one pair vector to another (pointer version).
 * 
 * @param result The target vector to append to.
 * @param vec2 Pointer to the source vector to append from.
 */
void help::append_vector(vector<pair<int, int>>& result, vector<pair<int, int>>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

/**
 * @brief Append elements from one pair vector to another (value version).
 * 
 * @param result The target vector to append to.
 * @param vec2 The source vector to append from.
 */
void help::append_vector(vector<pair<int, int>>& result, vector<pair<int, int>> vec2) {
    help::append_vector(result, &vec2);
}

/**
 * @brief Append elements from one string vector to another.
 * 
 * @param result The target vector to append to.
 * @param vec2 Pointer to the source vector to append from.
 */
void help::append_vector(vector<string>& result, vector<string>* vec2) {
    for (auto& item : *vec2) result.push_back(item);
}

/**
 * @brief Get free RAM memory information.
 * 
 * Uses Windows API to get memory information including available RAM,
 * virtual memory, and page file sizes.
 * 
 * @return vector<float> Vector containing:
 *         - Available RAM in GB
 *         - Available virtual memory in GB  
 *         - Available page file in GB
 *         - Memory load percentage
 */
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

/**
 * @brief Calculate standard deviation of a distribution.
 * 
 * @param distribution Pointer to the vector of data values.
 * @param mean Optional pre-calculated mean value. If not provided, it will be calculated.
 * @return double The standard deviation of the distribution.
 */
double help::get_stdev(vector<double>* distribution, double mean) {
    if (mean == -999999) {
        mean = help::get_mean(distribution);
    }
    double sos = 0;
    for (auto& sample : *distribution) sos += (sample - mean) * (sample - mean);
    double variance = sos / (distribution->size() - 1);
    return sqrt(variance);
}

/**
 * @brief Calculate mean of a distribution.
 * 
 * @param distribution Pointer to the vector of data values.
 * @return double The arithmetic mean of the distribution.
 */
double help::get_mean(vector<double>* distribution) {
    double sum = 0;
    for (double& sample : *distribution) {
        sum += sample;
    }
    return sum / distribution->size();
}

/**
 * @brief Get maximum value from a distribution.
 * 
 * @param distribution Pointer to the vector of data values.
 * @return double The maximum value in the distribution.
 */
double help::get_max(vector<double>* distribution) {
    double max = -INFINITY;
    for (double value : *distribution) {
        if (value > max) max = value;
    }
    return max;
}

/**
 * @brief Get minimum value from a distribution.
 * 
 * @param distribution Pointer to the vector of data values.
 * @return double The minimum value in the distribution.
 */
double help::get_min(vector<double>* distribution) {
    double min = INFINITY;
    for (double value : *distribution) {
        if (value < min) min = value;
    }
    return min;
}

/**
 * @brief Calculate Euclidean distance between two points.
 * 
 * @param p1 The first point (x, y).
 * @param p2 The second point (x, y).
 * @return float The Euclidean distance between p1 and p2.
 */
float help::get_dist(pair<float, float> p1, pair<float, float> p2) {
    float xdif = (p1.first - p2.first);
    float ydif = (p1.second - p2.second);
    return sqrtf(xdif * xdif + ydif * ydif);
}