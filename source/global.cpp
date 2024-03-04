#include <pybind11/pybind11.h>
#include <iostream>


using namespace std;
namespace py = pybind11;

void say_hello() {
    printf("Hi..\n");
}

int add(int a, int b) {
    return a + b;
}

PYBIND11_MODULE(dbr_cpp, module) {
    module.doc() = "DBR-cpp module (contains python extensions written in c++)"; // optional module docstring

    module.def("say_hello", &say_hello);
    module.def("add", &add);
}