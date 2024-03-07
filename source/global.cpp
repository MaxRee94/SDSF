#include <pybind11/pybind11.h>
#include "state.h"

using namespace std;
namespace py = pybind11;


void check_communication() {
    printf("Verified python->cpp communication.\n");
}

PYBIND11_MODULE(dbr_cpp, module) {
    module.doc() = "DBR-cpp module (contains python extensions written in c++)"; // optional module docstring

    module.def("check_communication", &check_communication);

    py::class_<State>(module, "State")
        .def(py::init<>())
        .def("populate_grid", &State::populate_grid);
}