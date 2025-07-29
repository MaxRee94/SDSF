#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "tests.h"


using namespace std;
namespace py = pybind11;


void check_communication() {
    printf("Verified python->cpp communication.\n");
}

void convert_from_numpy_array(py::array_t<float>& img, shared_ptr<float[]>& cover_image, int& width, int& height) {
	auto buf1 = img.request();
	float* ptr = (float*)buf1.ptr;
	width = buf1.shape[0];
    height = width;
    cover_image = make_shared<float[]>(width * height);
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			cover_image[y * width + x] = ptr[y * width + x];
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

py::array_t<int> as_2d_numpy_array(shared_ptr<int[]> distribution, int width) {
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
        setter(i, 1) = distribution2[i];
    }
    return numpy_array;
}

py::array_t<float> as_state_report_numpy_array(float* distribution, int rows) {
    constexpr size_t element_size = sizeof(float);
    const int cols = 4;
    size_t shape[2]{ rows, cols };
    size_t strides[2]{ cols * element_size, element_size };
    auto numpy_array = py::array_t<float>(shape, strides);
    auto setter = numpy_array.mutable_unchecked<2>();

    for (size_t i = 0; i < numpy_array.shape(0); i++)
    {
        for (int j = 0; j < cols; j++) {
            setter(i, j) = distribution[i * cols + j];
        }
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

py::array_t<float> as_2d_numpy_array(shared_ptr<float[]> distribution, int width) {
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

py::array_t<float> as_1d_numpy_array(shared_ptr<float[]> distribution, int size) {
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
        .def(py::init<const int&, const float&, const float&, const float&, const float&,
            const float&, const float&, map<string, map<string, float>>&, const float&,
            const float&, const float&, const float&>())
        .def("repopulate_grid", &State::repopulate_grid)
        .def("set_tree_cover", &State::set_tree_cover)
        .def("set_cover_from_image", [](State& state, py::array_t<float>& img, float& target_cover) {
            shared_ptr<float[]> cover_image;
            int width, height;
            convert_from_numpy_array(img, cover_image, width, height);
            if (target_cover > 0) state.set_cover_from_image(cover_image, width, height, target_cover);
            else state.set_cover_from_image(cover_image, width, height);
        })
        .def("get_state_table", [](State& state) {
            int number_of_values_per_tree = 4;
            float* state_table = new float[state.population.size() * number_of_values_per_tree];
            state.get_state_table(state_table);
            py::array_t<float> np_arr = as_state_report_numpy_array(state_table, state.population.size());
            delete[] state_table;
            return np_arr;
		})
        .def("get_tree_sizes", [](State& state) {
            float* tree_sizes = new float[state.population.size()];
            state.get_tree_sizes(tree_sizes);
            py::array_t<float> np_arr = as_1d_numpy_array(tree_sizes, state.population.size());
            delete[] tree_sizes;
            return np_arr;
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
        .def("get_distribution", [](Grid &grid, int &collect_states) {
            shared_ptr<int[]> state_distribution = grid.get_state_distribution(collect_states);
            return as_2d_numpy_array(state_distribution, grid.width);
        });

    py::class_<Dynamics>(module, "Dynamics")
        .def(py::init<>())
        .def(py::init<const int&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const map<string, map<string, float>>&, const int&, const float&, const float&, const int&, const int&, const int&,
            const float&, const int& >())
        .def_readwrite("time", &Dynamics::time)
        .def_readwrite("state", &Dynamics::state)
        .def_readwrite("timestep", &Dynamics::timestep)
        .def_readwrite("seeds_produced", &Dynamics::seeds_produced)
        .def_readwrite("fire_spatial_extent", &Dynamics::fire_spatial_extent)
        .def_readwrite("max_dbh", &Dynamics::max_dbh)
        .def("init_state", &Dynamics::init_state)
        .def("update", &Dynamics::update)
        .def("simulate_fires", &Dynamics::burn)
        .def("get_firefree_intervals", [](Dynamics& dynamics, string& type) {
            shared_ptr<float[]> intervals = dynamics.get_firefree_intervals(type);
            py::array_t<float> np_arr = as_1d_numpy_array(intervals, dynamics.grid->no_cells);
            return np_arr;
        })
        .def("get_no_recruits", [](Dynamics& dynamics, string& type) {
			float no_recruits = dynamics.get_no_recruits(type);
			return no_recruits;
		})
        .def("free", &Dynamics::free)
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
            shared_ptr<int[]> color_distribution = dynamics.resource_grid.get_color_distribution(species, type, verbosity);
            return as_2d_numpy_array(color_distribution, dynamics.resource_grid.width);
        })
        .def("get_resource_grid_lookup_table", [](Dynamics& dynamics, string& species) {
            shared_ptr<float[]> lookup_table = dynamics.resource_grid.get_lookup_table(species);
            return as_2d_numpy_array(lookup_table, dynamics.resource_grid.width * dynamics.resource_grid.width);
        })
        .def("set_resource_grid_lookup_table", [](Dynamics& dynamics, py::array_t<float>& lookup_table, string& species) {
            shared_ptr<float[]> float_array;
            int width, height;
			convert_from_numpy_array(lookup_table, float_array, width, height);
			dynamics.resource_grid.set_dist_lookup_table(float_array, species);
		})
        .def("get_fraction_time_spent_moving", [](Dynamics& dynamics) {
			return dynamics.fraction_time_spent_moving;
		})
        .def("get_fraction_seedlings_dead_due_to_shade", [](Dynamics& dynamics) {
            return (float)dynamics.no_seedlings_dead_due_to_shade / (float)dynamics.no_germination_attempts;
        })
        .def("get_fraction_seedlings_outcompeted", [](Dynamics& dynamics) {
            return (float)dynamics.no_seedling_competitions / (float)dynamics.no_germination_attempts;
        })
        .def("get_fraction_seedlings_outcompeted_by_older_trees", [](Dynamics& dynamics) {
            return (float)dynamics.no_competitions_with_older_trees / (float)dynamics.no_germination_attempts;
        })
        .def("get_fraction_cases_seedling_competition_and_shading", [](Dynamics& dynamics) {
            return (float)dynamics.no_cases_seedling_competition_and_shading / (float)dynamics.no_germination_attempts;
        })
        .def("get_fraction_cases_oldstem_competition_and_shading", [](Dynamics& dynamics) {
            return (float)dynamics.no_cases_oldstem_competition_and_shading / (float)dynamics.no_germination_attempts;
        })
        .def("get_no_germination_attempts", [](Dynamics& dynamics) {
			return dynamics.no_germination_attempts;
		})
        .def("get_no_fire_induced_deaths", [](Dynamics& dynamics) {
            return dynamics.no_fire_induced_deaths;
		})
        .def("get_no_fire_induced_topkills", [](Dynamics& dynamics) {
            return dynamics.no_fire_induced_topkills;
        })
        .def("get_no_fire_induced_nonseedling_topkills", [](Dynamics& dynamics) {
            return dynamics.no_fire_induced_nonseedling_topkills;
        })
        .def("get_initial_no_dispersals", [](Dynamics& dynamics) {
            return dynamics.initial_no_effective_dispersals;
        })
        .def("precompute_resourcegrid_lookup_table", [](Dynamics& dynamics, string& species) {
            dynamics.resource_grid.precompute_dist_lookup_table(species);
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
        .def(py::init<const float&, const float&, const float&, const float&, const map<string, map<string, float>>&, const float&, const float&, const float&,
            const float&, const float& >())
        .def("size", &Population::size);
    
    py::class_<Tests>(module, "Tests")
        .def(py::init<>())
        .def(py::init<const int&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&,
            const float&, const map<string, map<string, float>>&, const float&, const float&, const float&, const int&, const int&, const float&, const float&,
            const float&, const float&, const float&, const int&, const int&, const float&, const int& >()
        )
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