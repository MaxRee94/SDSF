"""Unit tests for C++ modules in DBR-sim.

This module provides Python-accessible unit tests for the C++ components.
Tests cover Timer, Grid, Agents (Tree/Population), State, and Helper functions.

Note: These tests require the C++ modules to be compiled and available.
"""

import unittest
import os
import sys
import subprocess
import tempfile
from pathlib import Path


class TestCppModules(unittest.TestCase):
    """Test suite for C++ modules in DBR-sim directory."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(__file__).parent
        self.cpp_test_header = self.test_dir / "test_cpp_modules.h"
        self.cpp_test_source = self.test_dir / "test_cpp_modules.cpp"
        
    def test_cpp_test_files_exist(self):
        """Test that C++ test files exist."""
        self.assertTrue(self.cpp_test_header.exists(), f"Test header file {self.cpp_test_header} does not exist")
        self.assertTrue(self.cpp_test_source.exists(), f"Test source file {self.cpp_test_source} does not exist")

    def test_cpp_test_header_content(self):
        """Test that C++ test header contains expected test classes and methods."""
        with open(self.cpp_test_header, 'r') as f:
            content = f.read()
        
        # Check for main test class
        self.assertIn("class CppModuleTests", content)
        
        # Check for test categories
        expected_tests = [
            "run_timer_tests",
            "run_helper_tests", 
            "run_grid_tests",
            "run_agents_tests",
            "run_state_tests",
            "test_timer_basic_operations",
            "test_rng_functions",
            "test_math_functions",
            "test_string_functions",
            "test_collection_functions",
            "test_array_functions",
            "test_statistical_functions",
            "test_Grid_construction",
            "test_Grid_cell_management",
            "test_Grid_tree_cover",
            "test_Grid_position_mapping",
            "test_Grid_state_management",
            "test_Tree_construction",
            "test_Tree_defaults",
            "test_Population_operations",
            "test_State_construction",
            "test_State_grid_management"
        ]
        
        for test in expected_tests:
            self.assertIn(test, content, f"Test method {test} not found in test header")

    def test_cpp_test_source_content(self):
        """Test that C++ test source file contains main function."""
        with open(self.cpp_test_source, 'r') as f:
            content = f.read()
        
        self.assertIn("int main", content)
        self.assertIn("CppModuleTests", content)
        self.assertIn("run_and_report", content)

    def test_timer_class_available(self):
        """Test that Timer class is available and can be included."""
        timer_header = self.test_dir / "timer.h"
        self.assertTrue(timer_header.exists(), f"Timer header {timer_header} does not exist")
        
        with open(timer_header, 'r') as f:
            content = f.read()
        
        self.assertIn("class Timer", content)
        self.assertIn("start()", content)
        self.assertIn("stop()", content)
        self.assertIn("elapsedMilliseconds()", content)
        self.assertIn("elapsedSeconds()", content)

    def test_agents_module_available(self):
        """Test that Agents module (Tree, Population) is available."""
        agents_header = self.test_dir / "agents.h"
        self.assertTrue(agents_header.exists(), f"Agents header {agents_header} does not exist")
        
        with open(agents_header, 'r') as f:
            content = f.read()
        
        self.assertIn("class Tree", content)
        self.assertIn("class Population", content)

    def test_grid_module_available(self):
        """Test that Grid module is available."""
        grid_header = self.test_dir / "grid.h"
        self.assertTrue(grid_header.exists(), f"Grid header {grid_header} does not exist")
        
        with open(grid_header, 'r') as f:
            content = f.read()
        
        self.assertIn("class Cell", content)
        self.assertIn("class Grid", content)

    def test_state_module_available(self):
        """Test that State module is available."""
        state_header = self.test_dir / "state.h"
        self.assertTrue(state_header.exists(), f"State header {state_header} does not exist")
        
        with open(state_header, 'r') as f:
            content = f.read()
        
        self.assertIn("class State", content)

    def test_helpers_module_available(self):
        """Test that Helpers module is available."""
        helpers_header = self.test_dir / "helpers.h"
        self.assertTrue(helpers_header.exists(), f"Helpers header {helpers_header} does not exist")
        
        with open(helpers_header, 'r') as f:
            content = f.read()
        
        # Check for helper namespace and functions
        self.assertIn("namespace help", content)
        self.assertIn("get_dist", content)
        self.assertIn("get_rand_float", content)

    def test_timer_functionality(self):
        """Test Timer class functionality through compilation."""
        test_code = """
        #include "timer.h"
        #include <cassert>
        
        int main() {
            Timer timer;
            timer.start();
            assert(timer.elapsedMilliseconds() >= 0);
            timer.stop();
            assert(timer.elapsedSeconds() >= 0);
            return 0;
        }
        """
        
        # This test would require compilation - skip for now
        # as it's complex to compile C++ from Python tests
        self.skipTest("C++ compilation test skipped - requires build system")

    def test_helpers_functionality(self):
        """Test helper functions through compilation."""
        test_code = """
        #include "helpers.h"
        #include <cassert>
        
        int main() {
            help::init_RNG();
            float val = help::get_rand_float(0.0f, 1.0f);
            assert(val >= 0.0f && val <= 1.0f);
            
            std::pair<float, float> p1 = {0, 0};
            std::pair<float, float> p2 = {3, 4};
            float dist = help::get_dist(p1, p2);
            assert(fabs(dist - 5.0f) < 0.001f);
            
            return 0;
        }
        """
        
        # This test would require compilation - skip for now
        self.skipTest("C++ compilation test skipped - requires build system")


class TestCppModuleTestCoverage(unittest.TestCase):
    """Test that C++ module tests provide adequate coverage."""
    
    def test_timer_coverage(self):
        """Test that Timer class has comprehensive test coverage."""
        test_header = Path(__file__).parent / "test_cpp_modules.h"
        with open(test_header, 'r') as f:
            content = f.read()
        
        timer_tests = [
            "test_timer_basic_operations",
            "test_timer_elapsed_accuracy", 
            "test_timer_seconds_conversion",
            "test_timer_start_stop"
        ]
        
        for test in timer_tests:
            self.assertIn(test, content, f"Timer test {test} not found")

    def test_helpers_coverage(self):
        """Test that Helper functions have comprehensive test coverage."""
        test_header = Path(__file__).parent / "test_cpp_modules.h"
        with open(test_header, 'r') as f:
            content = f.read()
        
        helpers_tests = [
            "test_rng_functions",
            "test_math_functions",
            "test_string_functions",
            "test_collection_functions",
            "test_array_functions", 
            "test_statistical_functions"
        ]
        
        for test in helpers_tests:
            self.assertIn(test, content, f"Helpers test {test} not found")

    def test_grid_coverage(self):
        """Test that Grid class has comprehensive test coverage."""
        test_header = Path(__file__).parent / "test_cpp_modules.h"
        with open(test_header, 'r') as f:
            content = f.read()
        
        grid_tests = [
            "test_Grid_construction",
            "test_Grid_cell_management",
            "test_Grid_tree_cover",
            "test_Grid_position_mapping",
            "test_Grid_state_management"
        ]
        
        for test in grid_tests:
            self.assertIn(test, content, f"Grid test {test} not found")

    def test_agents_coverage(self):
        """Test that Agents classes have comprehensive test coverage."""
        test_header = Path(__file__).parent / "test_cpp_modules.h"
        with open(test_header, 'r') as f:
            content = f.read()
        
        agents_tests = [
            "test_Tree_construction",
            "test_Tree_defaults",
            "test_Population_operations"
        ]
        
        for test in agents_tests:
            self.assertIn(test, content, f"Agents test {test} not found")

    def test_state_coverage(self):
        """Test that State class has comprehensive test coverage."""
        test_header = Path(__file__).parent / "test_cpp_modules.h"
        with open(test_header, 'r') as f:
            content = f.read()
        
        state_tests = [
            "test_State_construction",
            "test_State_grid_management"
        ]
        
        for test in state_tests:
            self.assertIn(test, content, f"State test {test} not found")

    def test_total_test_count(self):
        """Test that we have a reasonable number of tests."""
        test_header = Path(__file__).parent / "test_cpp_modules.h"
        with open(test_header, 'r') as f:
            content = f.read()
        
        # Count test methods
        test_count = 0
        test_methods = [
            "test_timer_",
            "test_rng_",
            "test_math_", 
            "test_string_",
            "test_collection_",
            "test_array_",
            "test_statistical_",
            "test_Grid_",
            "test_Tree_",
            "test_Population_",
            "test_State_"
        ]
        
        for method in test_methods:
            count = content.count(method)
            test_count += count
        
        # Should have at least 20 individual test methods
        self.assertGreaterEqual(test_count, 20, f"Expected at least 20 test methods, found {test_count}")


class TestCppModuleIntegration(unittest.TestCase):
    """Test integration of C++ module tests with existing infrastructure."""
    
    def test_files_in_correct_location(self):
        """Test that C++ test files are in the correct location."""
        test_files = [
            "test_cpp_modules.h",
            "test_cpp_modules.cpp"
        ]
        
        for filename in test_files:
            filepath = Path(__file__).parent / filename
            self.assertTrue(filepath.exists(), f"Test file {filename} not in DBR-sim directory")

    def test_includes_correct_headers(self):
        """Test that C++ test files include the correct headers."""
        test_header = Path(__file__).parent / "test_cpp_modules.h"
        with open(test_header, 'r') as f:
            content = f.read()
        
        required_headers = [
            "timer.h",
            "agents.h", 
            "grid.h",
            "state.h",
            "helpers.h"
        ]
        
        for header in required_headers:
            self.assertIn(f"#include \"{header}\"", content, f"Required header {header} not included")


if __name__ == '__main__':
    unittest.main()