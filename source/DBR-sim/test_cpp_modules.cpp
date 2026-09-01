#include "test_cpp_modules.h"

/**
 * @file test_cpp_modules.cpp
 * @brief Main entry point for C++ module unit tests.
 * 
 * This file provides a standalone main function for running the comprehensive
 * C++ module test suite defined in test_cpp_modules.h.
 */

/**
 * @brief Main function for running C++ module tests.
 * 
 * @param argc Number of command line arguments.
 * @param argv Array of command line argument strings.
 * @return int Exit code (0 for success, non-zero for failure).
 */
int main(int argc, char* argv[]) {
    int verbosity = 0;
    if (argc > 1 && string(argv[1]) == "-v") {
        verbosity = 1;
    }
    
    CppModuleTests tests(verbosity);
    tests.run_and_report();
    
    return 0;
}