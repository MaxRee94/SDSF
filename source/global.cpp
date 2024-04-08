#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "tests.h"


using namespace std;
namespace py = pybind11;


void check_communication() {
    printf("Verified python->cpp communication.\n");
}

py::array_t<int> as_2d_numpy_array(int* distribution, int width) {
    constexpr size_t element_size = sizeof(int);
    size_t shape[2]{ width, width };
    size_t strides[2]{ width * element_size, element_size };
    auto numpy_array = py::array_t<int>(shape, strides);
    auto setter = numpy_array.mutable_unchecked<2>();

    for (size_t i = 0; i < numpy_array.shape(0); i++)
    {
        for (size_t j = 0; j < numpy_array.shape(1); j++)
        {
            setter(i, j) = distribution[i * width + j];
        }
    }
    return numpy_array;
}

py::array_t<float> as_1d_numpy_array(float* distribution, int size) {
    constexpr size_t element_size = sizeof(float);
    size_t shape[1]{ size };
    size_t strides[1]{ element_size };
    auto numpy_array = py::array_t<float>(shape, strides);
    auto setter = numpy_array.mutable_unchecked<1>();

    for (size_t i = 0; i < numpy_array.shape(0); i++)
    {
        setter(i) = distribution[i];
    }
    return numpy_array;
}


PYBIND11_MODULE(dbr_cpp, module) {
    module.doc() = "DBR-cpp module (contains python extensions written in c++)";

    module.def("check_communication", &check_communication);

    py::class_<State>(module, "State")
        .def(py::init<>())
        .def(py::init<const int&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const float&>())
        .def("repopulate_grid", &State::repopulate_grid)
        .def("set_tree_cover", &State::set_tree_cover)
        .def_readwrite("grid", &State::grid)
        .def_readwrite("population", &State::population);

    py::class_<Grid>(module, "Grid")
        .def(py::init<>())
        .def(py::init<const int&, const float&>())
        .def("get_tree_cover", &Grid::get_tree_cover)
        .def_readwrite("width", &Grid::width)
        .def("get_distribution", [](Grid &grid, bool &collect_states) {
            int* state_distribution = grid.get_state_distribution(collect_states);
            return as_2d_numpy_array(state_distribution, grid.width);
        });

    py::class_<Dynamics>(module, "Dynamics")
        .def(py::init<>())
        .def(py::init<const int&, const float&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const int&>())
        .def_readwrite("time", &Dynamics::time)
        .def_readwrite("state", &Dynamics::state)
        .def_readwrite("timestep", &Dynamics::timestep)
        .def_readwrite("seeds_dispersed", &Dynamics::seeds_dispersed)
        .def_readwrite("fire_spatial_extent", &Dynamics::fire_spatial_extent)
        .def("init_state", &Dynamics::init_state)
        .def("set_global_kernel", &Dynamics::set_global_kernel)
        .def("update", &Dynamics::update)
        .def("simulate_fires", &Dynamics::burn)
        .def("get_firefree_intervals", [](Dynamics& dynamics) {
            float* intervals = dynamics.get_firefree_intervals();
            return as_1d_numpy_array(intervals, dynamics.grid->no_cells);
        });

    py::class_<Population>(module, "Population")
        .def(py::init<>())
        .def(py::init<const float&, const float&, const float&, const float&, const float&>())
        .def("size", &Population::size);
    
    py::class_<Tests>(module, "Tests")
        .def(py::init<>())
        .def("run_all", &Tests::run_all);
}