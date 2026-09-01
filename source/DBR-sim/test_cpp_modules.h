#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <cassert>
#include "timer.h"
#include "agents.h"
#include "grid.h"
#include "state.h"
#include "helpers.h"

// Platform-specific includes for sleep functions
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

using namespace std;

/**
 * Comprehensive unit test suite for C++ modules in DBR-sim.
 * 
 * Tests are organized by module and cover:
 * - Timer class functionality
 * - Helper functions (math, string, vector operations)
 * - Grid and Cell classes
 * - Tree and Population classes
 * - State class
 */

class CppModuleTests {
public:
    CppModuleTests() = default;
    CppModuleTests(int verbosity) : verbosity(verbosity) {}
    
    // Run all tests and return list of failed tests
    vector<string> run_all() {
        vector<string> failed_tests;
        
        // Timer tests
        run_timer_tests(failed_tests);
        
        // Helper function tests
        run_helper_tests(failed_tests);
        
        // Grid tests
        run_grid_tests(failed_tests);
        
        // Agents tests
        run_agents_tests(failed_tests);
        
        // State tests
        run_state_tests(failed_tests);
        
        return failed_tests;
    }
    
    int verbosity = 0;

private:
    // ==================== TIMER TESTS ====================
    
    void run_timer_tests(vector<string>& failed_tests) {
        if (verbosity > 0) printf("Running Timer tests...\n");
        
        if (!test_timer_basic_operations()) {
            failed_tests.push_back("Timer::basic_operations");
        }
        if (!test_timer_elapsed_accuracy()) {
            failed_tests.push_back("Timer::elapsed_accuracy");
        }
        if (!test_timer_seconds_conversion()) {
            failed_tests.push_back("Timer::seconds_conversion");
        }
        if (!test_timer_start_stop()) {
            failed_tests.push_back("Timer::start_stop");
        }
        
        if (verbosity > 0) printf("Timer tests completed.\n");
    }

    bool test_timer_basic_operations() {
        Timer timer;
        timer.start();
        
        // Timer should be running and have non-negative elapsed time
        double elapsed = timer.elapsedMilliseconds();
        if (elapsed < 0) return false;
        
        return true;
    }

    bool test_timer_elapsed_accuracy() {
        Timer timer;
        timer.start();
        
        // Wait a small amount to ensure measurable elapsed time
        #ifdef _WIN32
        Sleep(50); // 50ms delay on Windows
        #else
        usleep(50000); // 50ms delay on Unix
        #endif
        
        double elapsed = timer.elapsedMilliseconds();
        // Should be at least 50ms (allow some tolerance)
        if (elapsed < 40) return false;
        
        return true;
    }

    bool test_timer_seconds_conversion() {
        Timer timer;
        timer.start();
        
        #ifdef _WIN32
        Sleep(100); // 100ms delay
        #else
        usleep(100000); // 100ms delay
        #endif
        
        timer.stop();
        double elapsed_ms = timer.elapsedMilliseconds();
        double elapsed_sec = timer.elapsedSeconds();
        
        // Check that conversion is correct (within tolerance)
        double expected_sec = elapsed_ms / 1000.0;
        if (fabs(elapsed_sec - expected_sec) > 0.01) return false;
        
        return true;
    }

    bool test_timer_start_stop() {
        Timer timer;
        timer.start();
        
        #ifdef _WIN32
        Sleep(20);
        #else
        usleep(20000);
        #endif
        
        double running_elapsed = timer.elapsedMilliseconds();
        timer.stop();
        double stopped_elapsed = timer.elapsedMilliseconds();
        
        // After stopping, elapsed time should still be accessible and >= running time
        if (stopped_elapsed < running_elapsed) return false;
        
        // Start again and check it resets properly
        timer.start();
        
        #ifdef _WIN32
        Sleep(10);
        #else
        usleep(10000);
        #endif
        
        double new_elapsed = timer.elapsedMilliseconds();
        // New elapsed should be less than previous (since we waited less time)
        // Note: this might not always be true due to timer precision, so be lenient
        if (new_elapsed > stopped_elapsed + 50) return false;
        
        return true;
    }

    // ==================== HELPER FUNCTION TESTS ====================
    
    void run_helper_tests(vector<string>& failed_tests) {
        if (verbosity > 0) printf("Running Helper function tests...\n");
        
        if (!test_rng_functions()) {
            failed_tests.push_back("help::random_number_generation");
        }
        if (!test_math_functions()) {
            failed_tests.push_back("help::math_functions");
        }
        if (!test_string_functions()) {
            failed_tests.push_back("help::string_functions");
        }
        if (!test_collection_functions()) {
            failed_tests.push_back("help::collection_functions");
        }
        if (!test_array_functions()) {
            failed_tests.push_back("help::array_functions");
        }
        if (!test_statistical_functions()) {
            failed_tests.push_back("help::statistical_functions");
        }
        
        if (verbosity > 0) printf("Helper function tests completed.\n");
    }

    bool test_rng_functions() {
        help::init_RNG();
        
        // Test get_rand_float produces values in range
        float val = help::get_rand_float(0.0f, 1.0f);
        if (val < 0.0f || val > 1.0f) return false;
        
        // Test get_rand_float with different ranges
        val = help::get_rand_float(5.0f, 10.0f);
        if (val < 5.0f || val > 10.0f) return false;
        
        // Test get_rand_uint
        unsigned int uint_val = help::get_rand_uint(0.0f, 100.0f);
        if (uint_val < 0 || uint_val > 100) return false;
        
        return true;
    }

    bool test_math_functions() {
        // Test get_dist
        pair<float, float> p1 = {0.0f, 0.0f};
        pair<float, float> p2 = {3.0f, 4.0f};
        float dist = help::get_dist(p1, p2);
        if (fabs(dist - 5.0f) > 0.001f) return false; // 3-4-5 triangle
        
        // Test fisqrt (fast inverse square root)
        float result = help::fisqrt(25.0f);
        if (fabs(result - 5.0f) > 0.1f) return false; // Approximate
        
        return true;
    }

    bool test_string_functions() {
        // Test replace_occurrences
        string test = "hello world";
        string result = help::replace_occurrences(test, "world", "cpp");
        if (result != "hello cpp") return false;
        
        // Test replace_occurrences with multiple occurrences
        string multi = "test test test";
        result = help::replace_occurrences(multi, "test", "pass");
        if (result != "pass pass pass") return false;
        
        // Test is_in
        if (!help::is_in("hello world", "hello")) return false;
        if (help::is_in("hello world", "cpp")) return false;
        
        // Test ends_with
        if (!help::ends_with("test.cpp", ".cpp")) return false;
        if (help::ends_with("test.cpp", ".py")) return false;
        
        // Test add_padding
        string padded = help::add_padding("test", 42);
        if (padded != "test0042") return false;
        
        padded = help::add_padding("file", 5);
        if (padded != "file0005") return false;
        
        return true;
    }

    bool test_collection_functions() {
        // Test is_in with vectors
        vector<int> vec = {1, 2, 3, 4, 5};
        if (!help::is_in(&vec, 3)) return false;
        if (help::is_in(&vec, 6)) return false;
        
        // Test remove
        vector<int> test_vec = {1, 2, 3, 4, 5};
        help::remove(&test_vec, 3);
        if (test_vec.size() != 4) return false;
        if (help::is_in(&test_vec, 3)) return false;
        
        // Test have_overlap
        vector<int> a = {1, 2, 3};
        vector<int> b = {3, 4, 5};
        vector<int> c = {6, 7, 8};
        if (!help::have_overlap(&a, &b)) return false;
        if (help::have_overlap(&a, &c)) return false;
        
        // Test append_vector
        vector<int> result;
        vector<int> vec1 = {1, 2};
        vector<int> vec2 = {3, 4};
        help::append_vector(result, &vec1);
        help::append_vector(result, &vec2);
        if (result.size() != 4) return false;
        if (result[0] != 1 || result[1] != 2 || result[2] != 3 || result[3] != 4) return false;
        
        return true;
    }

    bool test_array_functions() {
        // Test populate_with_zeroes for double
        double* double_arr = new double[10];
        help::populate_with_zeroes(double_arr, 2, 5);
        for (int i = 0; i < 10; i++) {
            if (double_arr[i] != 0.0) {
                delete[] double_arr;
                return false;
            }
        }
        delete[] double_arr;
        
        // Test populate_with_zeroes for uint
        unsigned int* uint_arr = new unsigned int[6];
        help::populate_with_zeroes(uint_arr, 2, 3);
        for (int i = 0; i < 6; i++) {
            if (uint_arr[i] != 0) {
                delete[] uint_arr;
                return false;
            }
        }
        delete[] uint_arr;
        
        return true;
    }

    bool test_statistical_functions() {
        vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
        
        // Test get_mean
        double mean = help::get_mean(&data);
        if (fabs(mean - 3.0) > 0.001) return false;
        
        // Test get_stdev
        double stdev = help::get_stdev(&data, mean);
        // For [1,2,3,4,5], stdev should be sqrt(2.5) ≈ 1.5811
        if (fabs(stdev - 1.5811) > 0.01) return false;
        
        // Test get_max and get_min
        if (help::get_max(&data) != 5.0) return false;
        if (help::get_min(&data) != 1.0) return false;
        
        // Test with single element
        vector<double> single = {42.0};
        if (help::get_mean(&single) != 42.0) return false;
        if (help::get_max(&single) != 42.0) return false;
        if (help::get_min(&single) != 42.0) return false;
        
        return true;
    }

    // ==================== GRID TESTS ====================
    
    void run_grid_tests(vector<string>& failed_tests) {
        if (verbosity > 0) printf("Running Grid tests...\n");
        
        if (!test_Grid_construction()) {
            failed_tests.push_back("Grid::construction");
        }
        if (!test_Grid_cell_management()) {
            failed_tests.push_back("Grid::cell_management");
        }
        if (!test_Grid_tree_cover()) {
            failed_tests.push_back("Grid::tree_cover");
        }
        if (!test_Grid_position_mapping()) {
            failed_tests.push_back("Grid::position_mapping");
        }
        if (!test_Grid_state_management()) {
            failed_tests.push_back("Grid::state_management");
        }
        
        if (verbosity > 0) printf("Grid tests completed.\n");
    }

    bool test_Grid_construction() {
        // Test default constructor
        Grid grid1;
        if (grid1.size != 1000) return false; // Default size should be 1000
        
        // Test constructor with size
        Grid grid2(100);
        if (grid2.size != 100) return false;
        if (grid2.distribution == nullptr) return false;
        if (grid2.no_savanna_cells != 100 * 100) return false;
        if (grid2.no_forest_cells != 0) return false;
        
        return true;
    }

    bool test_Grid_cell_management() {
        Grid grid(10); // 10x10 grid
        
        // Test reset
        grid.reset();
        if (grid.no_savanna_cells != 100 || grid.no_forest_cells != 0) return false;
        
        // Test set_to_forest with index
        Tree tree({1.0f, 1.0f});
        grid.set_to_forest(0, &tree);
        if (grid.no_forest_cells != 1 || grid.no_savanna_cells != 99) return false;
        
        // Test set_to_savanna
        grid.set_to_savanna(0);
        if (grid.no_forest_cells != 0 || grid.no_savanna_cells != 100) return false;
        
        // Test set_to_forest with position
        pair<int, int> pos = {1, 1};
        grid.set_to_forest(pos, &tree);
        if (grid.no_forest_cells != 1 || grid.no_savanna_cells != 99) return false;
        
        return true;
    }

    bool test_Grid_tree_cover() {
        Grid grid(10); // 10x10 = 100 cells
        
        // Initially should be 0 tree cover
        if (grid.get_tree_cover() != 0.0f) return false;
        
        Tree tree({0.0f, 0.0f});
        
        // Set exactly 50 cells to forest
        for (int i = 0; i < 50; i++) {
            grid.set_to_forest(i, &tree);
        }
        
        // Should have 50% tree cover
        float cover = grid.get_tree_cover();
        if (fabs(cover - 0.5f) > 0.01f) return false;
        
        // Test edge cases
        Grid empty_grid(10);
        if (empty_grid.get_tree_cover() != 0.0f) return false;
        
        return true;
    }

    bool test_Grid_position_mapping() {
        Grid grid(100);
        
        // Test get_cell_at_position
        pair<float, float> pos = {1.5f, 2.7f};
        Cell* cell = grid.get_cell_at_position(pos);
        if (cell == nullptr) return false;
        
        // Test get_gridbased_position for Tree
        Tree tree({1.5f, 2.7f});
        pair<int, int> grid_pos = grid.get_gridbased_position(&tree);
        // With cell_width = 1.5, position (1.5, 2.7) should map to grid (1, 1)
        if (grid_pos.first != 1 || grid_pos.second != 1) return false;
        
        return true;
    }

    bool test_Grid_state_management() {
        Grid grid(5); // 5x5 grid
        
        // Test get_state_distribution
        int* state_dist = grid.get_state_distribution();
        if (state_dist == nullptr) return false;
        
        // Initially all should be 0 (savanna)
        for (int i = 0; i < 25; i++) {
            if (state_dist[i] != 0) {
                delete[] state_dist;
                return false;
            }
        }
        delete[] state_dist;
        
        // Test redo_count
        Tree tree({0.0f, 0.0f});
        grid.set_to_forest(0, &tree);
        grid.set_to_forest(1, &tree);
        grid.set_to_forest(2, &tree);
        
        if (grid.no_forest_cells != 3) return false;
        
        // Manually change a cell state to make counts wrong
        grid.distribution[3].state = 1; // Set to forest without updating counts
        grid.no_forest_cells = 0; // Reset counter
        
        grid.redo_count();
        if (grid.no_forest_cells != 4) return false; // Should recount to 4
        
        return true;
    }

    // ==================== AGENTS TESTS ====================
    
    void run_agents_tests(vector<string>& failed_tests) {
        if (verbosity > 0) printf("Running Agents tests...\n");
        
        if (!test_Tree_construction()) {
            failed_tests.push_back("Tree::construction");
        }
        if (!test_Tree_defaults()) {
            failed_tests.push_back("Tree::defaults");
        }
        if (!test_Population_operations()) {
            failed_tests.push_back("Population::operations");
        }
        
        if (verbosity > 0) printf("Agents tests completed.\n");
    }

    bool test_Tree_construction() {
        pair<float, float> pos = {1.0f, 2.0f};
        
        // Test default constructor
        Tree tree1;
        if (tree1.radius != 0) return false;
        if (tree1.life_phase != 0) return false;
        if (tree1.position.first != 0 || tree1.position.second != 0) return false;
        
        // Test constructor with position
        Tree tree2(pos);
        if (tree2.position.first != 1.0f || tree2.position.second != 2.0f) return false;
        
        // Test full constructor
        vector<float> strategy = {0.1f, 0.2f, 0.7f};
        Tree tree3(5.0f, strategy, 1, pos);
        if (tree3.radius != 5.0f) return false;
        if (tree3.life_phase != 1) return false;
        if (tree3.strategy.size() != 3) return false;
        if (tree3.strategy[0] != 0.1f) return false;
        
        return true;
    }

    bool test_Tree_defaults() {
        Tree tree;
        
        // Check default values
        if (tree.radius != 0.0f) return false;
        if (tree.life_phase != 0) return false;
        if (tree.position.first != 0.0f || tree.position.second != 0.0f) return false;
        if (!tree.strategy.empty()) return false;
        
        return true;
    }

    bool test_Population_operations() {
        Population pop;
        
        // Test initial state
        if (pop.members.size() != 0) return false;
        
        // Test adding trees
        Tree tree1({1.0f, 2.0f});
        Tree tree2({3.0f, 4.0f});
        
        pop.add(tree1);
        if (pop.members.size() != 1) return false;
        
        pop.add(tree2);
        if (pop.members.size() != 2) return false;
        
        // Check that trees are stored correctly
        if (pop.members[0].position.first != 1.0f) return false;
        if (pop.members[1].position.second != 4.0f) return false;
        
        return true;
    }

    // ==================== STATE TESTS ====================
    
    void run_state_tests(vector<string>& failed_tests) {
        if (verbosity > 0) printf("Running State tests...\n");
        
        if (!test_State_construction()) {
            failed_tests.push_back("State::construction");
        }
        if (!test_State_grid_management()) {
            failed_tests.push_back("State::grid_management");
        }
        
        if (verbosity > 0) printf("State tests completed.\n");
    }

    bool test_State_construction() {
        // Test default constructor
        State state1;
        if (state1.grid.size != 1000) return false; // Default grid size
        if (state1.population.members.size() != 0) return false;
        
        // Test constructor with grid size
        State state2(100);
        if (state2.grid.size != 100) return false;
        if (state2.population.members.size() != 0) return false;
        
        return true;
    }

    bool test_State_grid_management() {
        State state(10); // 10x10 grid
        
        // Test that grid is properly initialized
        if (state.grid.size != 10) return false;
        if (state.grid.no_savanna_cells != 100) return false;
        if (state.grid.no_forest_cells != 0) return false;
        
        return true;
    }

public:
    // Utility function to run tests and print results
    void run_and_report() {
        printf("Beginning C++ module unit tests...\n");
        vector<string> failed_tests = run_all();
        
        int total_tests = 0;
        // Count total tests by counting test methods
        total_tests = 4 +  // Timer tests
                      6 +  // Helper tests
                      5 +  // Grid tests
                      3 +  // Agents tests
                      2;  // State tests
        
        int successful_tests = total_tests - failed_tests.size();
        
        printf("Completed C++ module tests. ");
        if (failed_tests.size() > 0) {
            printf("\nTests failed (%i / %i):\n - %s\n", 
                   (int)failed_tests.size(), total_tests, 
                   help::join(&failed_tests, "\n - ").c_str());
        } else {
            printf("All tests (%i) successful.\n", successful_tests);
        }
    }
};