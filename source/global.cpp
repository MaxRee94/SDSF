#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "state.h"


using namespace std;
namespace py = pybind11;


void check_communication() {
    printf("Verified python->cpp communication.\n");
}

py::array_t<int> as_numpy_array(int* distribution, int width) {
    constexpr size_t element_size = sizeof(int);
    size_t shape[2]{ width, width };
    size_t strides[2]{ width * element_size, element_size };
    auto numpy_array = py::array_t<int>(shape, strides);
    auto setter = numpy_array.mutable_unchecked<2>();

    for (size_t i = 0; i < numpy_array.shape(0); i++)
    {
        for (size_t j = 0; j < numpy_array.shape(1); j++)
        {
            setter(j, i) = distribution[i * width + j] * ((float)j / (1000.0 / 255.0));
            //setter(i, j) = distribution[i * width + j];
        }
    }
    return numpy_array;
}


PYBIND11_MODULE(dbr_cpp, module) {
    module.doc() = "DBR-cpp module (contains python extensions written in c++)";

    module.def("check_communication", &check_communication);

    py::class_<State>(module, "State")
        .def(py::init<>())
        .def(py::init<const int&, const float&, const float&>())
        .def("populate_grid", &State::populate_grid)
        .def("set_tree_cover", &State::set_tree_cover)
        .def_readwrite("grid", &State::grid);

    py::class_<Grid>(module, "Grid")
        .def(py::init<>())
        .def(py::init<const int&, const float&>())
        .def_readwrite("size", &Grid::size)
        .def("get_distribution", [](Grid &grid) {
            int* state_distribution = grid.get_state_distribution();
            return as_numpy_array(state_distribution, grid.size);
        });
}