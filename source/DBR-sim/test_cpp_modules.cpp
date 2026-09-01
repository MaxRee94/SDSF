#include "test_cpp_modules.h"

// This file provides a simple main function for running C++ module tests
// It can be used standalone or integrated with the Python test infrastructure

int main(int argc, char* argv[]) {
    int verbosity = 0;
    if (argc > 1 && string(argv[1]) == "-v") {
        verbosity = 1;
    }
    
    CppModuleTests tests(verbosity);
    tests.run_and_report();
    
    return 0;
}