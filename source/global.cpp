#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "tests.h"


using namespace std;
namespace py = pybind11;


void check_communication() {
    printf("Verified python->cpp communication.\n");
}

void convert_from_numpy_array(py::array_t<float>& img, float*& cover_image, int& width, int& height) {
	auto buf1 = img.request();
	float* ptr = (float*)buf1.ptr;
	width = buf1.shape[0];
	height = buf1.shape[1];
	cover_image = new float[width * height];
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			cover_image[x * height + y] = ptr[x * height + y];
		}
	}
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

py::array_t<float> as_2d_pairwise_numpy_array(float* distribution1, float* distribution2, int size) {
    constexpr size_t element_size = sizeof(float);
    size_t shape[2]{ size, 2 };
    size_t strides[2]{ 2 * element_size, element_size };
    auto numpy_array = py::array_t<float>(shape, strides);
    auto setter = numpy_array.mutable_unchecked<2>();

    for (size_t i = 0; i < numpy_array.shape(0); i++)
    {
        setter(i, 0) = distribution1[i];
        setter(i, 1) = distribution2[i + 1];
    }
    return numpy_array;
}

py::array_t<float> as_2d_numpy_array(float* distribution, int width) {
    constexpr size_t element_size = sizeof(float);
    size_t shape[2]{ width, width };
    size_t strides[2]{ width * element_size, element_size };
    auto numpy_array = py::array_t<float>(shape, strides);
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
    module.def("init_RNG", &help::init_RNG);

    py::class_<State>(module, "State")
        .def(py::init<>())
        .def(py::init<const int&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const float&, map<string, map<string, float>>&, const float&>())
        .def("repopulate_grid", &State::repopulate_grid)
        .def("set_tree_cover", &State::set_tree_cover)
        .def("set_cover_from_image", [](State& state, py::array_t<float>& img) {
            float* cover_image;
            int width, height;
            convert_from_numpy_array(img, cover_image, width, height);
            state.set_cover_from_image(cover_image, width, height);
        })
        .def("get_tree_positions", [](State& state) {
            float* tree_positions_x = new float[state.population.size()];
            float* tree_positions_y = new float[state.population.size()];
            state.get_tree_positions(tree_positions_x, tree_positions_y);
			return as_2d_pairwise_numpy_array(tree_positions_x, tree_positions_y, state.population.size());
		})
        .def_readwrite("grid", &State::grid)
        .def_readwrite("population", &State::population)
        .def_readwrite("initial_tree_cover", &State::initial_tree_cover);

    py::class_<Grid>(module, "Grid")
        .def(py::init<>())
        .def(py::init<const int&, const float&>())
        .def("get_tree_cover", &Grid::get_tree_cover)
        .def_readwrite("width", &Grid::width)
        .def_readwrite("width_r", &Grid::width_r)
        .def_readwrite("tree_cover", &Grid::tree_cover)
        .def("get_distribution", [](Grid &grid, bool &collect_states) {
            int* state_distribution = grid.get_state_distribution(collect_states);
            return as_2d_numpy_array(state_distribution, grid.width);
        });

    py::class_<Dynamics>(module, "Dynamics")
        .def(py::init<>())
        .def(py::init<const int&, const float&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const map<string, map<string, float>>&, const float&, const float&, const int& >())
        .def_readwrite("time", &Dynamics::time)
        .def_readwrite("state", &Dynamics::state)
        .def_readwrite("timestep", &Dynamics::timestep)
        .def_readwrite("seeds_dispersed", &Dynamics::seeds_dispersed)
        .def_readwrite("fire_spatial_extent", &Dynamics::fire_spatial_extent)
        .def("init_state", &Dynamics::init_state)
        .def("update", &Dynamics::update)
        .def("simulate_fires", &Dynamics::burn)
        .def("get_firefree_intervals", [](Dynamics& dynamics) {
            float* intervals = dynamics.get_firefree_intervals();
            return as_1d_numpy_array(intervals, dynamics.grid->no_cells);
        })
        .def("set_global_linear_kernel", &Dynamics::set_global_linear_kernel)
        .def("set_global_wind_kernel", &Dynamics::set_global_wind_kernel)
        .def("set_global_animal_kernel", [](Dynamics& dynamics, const py::dict& _animal_dispersal_map) {
            std::map<string, std::map<string, float>> animal_dispersal_map = py::cast<std::map<string, std::map<string, float>>>(_animal_dispersal_map);
            dynamics.set_global_animal_kernel(animal_dispersal_map);
        })
        .def("set_global_kernels", [](Dynamics& dynamics, const py::dict& _nonanimal_kernel_params, const py::dict& _animal_kernel_params) {
            std::map<string, std::map<string, float>> nonanimal_kernel_params = py::cast<std::map<string, std::map<string, float>>>(_nonanimal_kernel_params);
            std::map<string, std::map<string, float>> animal_kernel_params = py::cast<std::map<string, std::map<string, float>>>(_animal_kernel_params);
            dynamics.set_global_kernels(nonanimal_kernel_params, animal_kernel_params);
        })
        .def("get_resource_grid_colors", [](Dynamics& dynamics, string& species, string& type, int& verbosity) {
            int* color_distribution = dynamics.resource_grid.get_color_distribution(species, type, verbosity);
            return as_2d_numpy_array(color_distribution, dynamics.resource_grid.width);
        });

    py::class_<DiscreteProbabilityModel>(module, "DiscreteProbabilityModel")
        .def(py::init<>())
        .def(py::init<const int&>())
        .def("set_probabilities", [](DiscreteProbabilityModel& model, py::array_t<float>& arr) {
            auto buf1 = arr.request();
            float* ptr = (float*)buf1.ptr;
            size_t size = buf1.shape[0];
			for (int x = 0; x < size; x++) {
                model.probabilities[x] = ptr[x];
            }
            model.build_cdf();
		})
        .def("sample", [](DiscreteProbabilityModel& model) {
            return model.sample();
        });

    py::class_<Population>(module, "Population")
        .def(py::init<>())
        .def(py::init<const float&, const float&, const float&, const float&, const float&, const map<string, map<string, float>>&, const float& >())
        .def("size", &Population::size);
    
    py::class_<Tests>(module, "Tests")
        .def(py::init<>())
        .def("run_all", &Tests::run_all);

    py::class_<PieceWiseLinearProbModel>(module, "PieceWiseLinearProbModel")
        .def(py::init<>())
        .def(py::init<const float&>())
        .def("sample", &PieceWiseLinearProbModel::sample);

    py::class_<Kernel>(module, "Kernel")
        .def(py::init<>())
        .def(py::init < const int&, const float&, const float&, const float&, const float&>())
        .def(py::init < const int&, const float&, const float&, const float&, const float&, const float&, const float&, const float&>())
        .def("get_dist", &Kernel::get_dist)
        .def("build", &Kernel::build);
}