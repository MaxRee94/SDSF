#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "state.h"


using namespace std;
namespace py = pybind11;


void check_communication() {
    printf("Verified python->cpp communication.\n");
}

//py::array_t<double> get_distribution(const Grid &grid) {
//    double* testarray = new double[grid.gridsize * grid.gridsize];
//
//    constexpr size_t element_size = sizeof(double);
//    size_t shape[2]{ grid.gridsize, grid.gridsize };
//    size_t strides[2]{ grid.gridsize * element_size, element_size };
//    auto numpy_array = py::array_t<double>(shape, strides);
//    auto view = numpy_array.mutable_unchecked<2>();
//
//    for (size_t i = 0; i < numpy_array.shape(0); i++)
//    {
//        for (size_t j = 0; j < numpy_array.shape(1); j++)
//        {
//            testarray[i * grid.gridsize + j] = 0;
//            if (i > 4 * grid.gridsize) testarray[i * grid.gridsize + j] = 1;
//            view(i, j) = testarray[i * grid.gridsize + j];
//        }
//    }
//
//    delete[] testarray;
//    return numpy_array;
//}

PYBIND11_MODULE(dbr_cpp, module) {
    module.doc() = "DBR-cpp module (contains python extensions written in c++)"; // optional module docstring

    module.def("check_communication", &check_communication);

    py::class_<State>(module, "State")
        .def(py::init<>())
        .def("populate_grid", &State::populate_grid)
        .def_readwrite("grid", &State::grid);

    py::class_<Grid>(module, "Grid")
        .def(py::init<>())
        .def("get_distribution", [](const Grid &grid) {
            int gridsize = grid.gridsize;
            double* testarray = new double[gridsize * gridsize];

            constexpr size_t element_size = sizeof(double);
            size_t shape[2]{ gridsize, gridsize };
            size_t strides[2]{ gridsize * element_size, element_size };
            auto numpy_array = py::array_t<double>(shape, strides);
            auto view = numpy_array.mutable_unchecked<2>();

            for (size_t i = 0; i < numpy_array.shape(0); i++)
            {
                for (size_t j = 0; j < numpy_array.shape(1); j++)
                {
                    testarray[i * gridsize + j] = 0;
                    if (i > gridsize/2) testarray[i * gridsize + j] = 1;
                    view(i, j) = testarray[i * gridsize + j];
                }
            }

            delete[] testarray;
            return numpy_array;
        }
    );

}