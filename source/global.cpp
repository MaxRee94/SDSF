#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "tests.h"


using namespace std;
namespace py = pybind11;


void check_communication() {
    printf("Verified python->cpp communication.\n");
}

void convert_from_numpy_array(py::array_t<float>& img, shared_ptr<float[]>& output_image, int& width, int& height) {
	auto buf1 = img.request();
	float* ptr = (float*)buf1.ptr;
	width = buf1.shape[0];
    height = width;
    output_image = make_shared<float[]>(width * height);
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			output_image[y * width + x] = ptr[y * width + x];
		}
	}
}

py::list convert_from_idx_pair_vector_to_position_pylist(const vector<pair<int, int>>& vec, Grid& grid) {
    py::list pylist;
    for (size_t j = 0; j < vec.size(); j++) {
        pair<int, int> _pair = vec[j];
        pair<int, int> pos_1_cpp = grid.idx_2_pos(_pair.first);
        pair<int, int> pos_2_cpp = grid.idx_2_pos(_pair.second);
        py::tuple pos_1 = py::make_tuple(pos_1_cpp.first, pos_1_cpp.second);
        py::tuple pos_2 = py::make_tuple(pos_2_cpp.first, pos_2_cpp.second);
        py::tuple neighbors = py::make_tuple(pos_1, pos_2);
        pylist.append(neighbors);
    }
	return pylist;
}

py::list convert_from_idx_vector_to_position_pylist(const vector<int>& vec, Grid& grid) {
    py::list pylist;
    for (size_t j = 0; j < vec.size(); j++) {
        pair<int, int> pos_cpp = grid.idx_2_pos(vec[j]);
        py::tuple pos = py::make_tuple(pos_cpp.first, pos_cpp.second);
        pylist.append(pos);
    }
    return pylist;
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

py::array_t<float> as_1d_numpy_array(vector<float> distribution) {
    constexpr size_t element_size = sizeof(float);
    size_t shape[1]{ distribution.size() };
    size_t strides[1]{ element_size };
    auto numpy_array = py::array_t<float>(shape, strides);
    auto setter = numpy_array.mutable_unchecked<1>();

    for (size_t i = 0; i < numpy_array.shape(0); i++)
    {
        setter(i) = distribution[i];
    }
    return numpy_array;
}

Dynamics create_dynamics(py::dict dict) {
    auto get = [&](const std::string& key) {
        if (!dict.contains(key)) throw std::runtime_error("Missing key: " + key);
        return dict[key.c_str()];
    };

    return Dynamics(
        get("timestep").cast<int>(),
        get("cell_width").cast<float>(),
        get("self_ignition_factor").cast<float>(),
        get("rainfall").cast<float>(),
        get("seed_bearing_threshold").cast<float>(),
        get("growth_rate_multiplier").cast<float>(),
        get("unsuppressed_flammability").cast<float>(),
        get("max_dbh").cast<float>(),
        get("saturation_threshold").cast<float>(),
        get("fire_resistance_params").cast<map<string, float>>(),
        get("background_mortality").cast<float>(),
        get("strategy_distribution_params").cast<map<string, map<string, float>>>(),
        get("resource_grid_width").cast<int>(),
        get("mutation_rate").cast<float>(),
        get("STR").cast<float>(),
        get("verbosity").cast<int>(),
        get("random_seed").cast<int>(),
        get("firefreq_random_seed").cast<int>(),
        get("enforce_no_recruits").cast<float>(),
        get("animal_group_size").cast<int>(),
		get("display_fire_effects").cast<bool>()
    );
}


PYBIND11_MODULE(dbr_cpp, module) {
    module.doc() = "DBR-cpp module (contains python extensions written in c++)";

    module.def("check_communication", &check_communication);
    module.def("init_RNG", &help::init_RNG);
    module.def("create_dynamics", &create_dynamics, "Create a Dynamics object from a Python dictionary");

    py::class_<State>(module, "State")
        .def(py::init<>())
        .def(py::init<const int&, const float&, const float&, const float&, const float&,
            const float&, const float&, map<string, map<string, float>>&, const float&,
            const float&, const float&, const float&, const float&, const float&>())
        .def("repopulate_grid", [](State& state, int verbosity) {state.repopulate_grid(verbosity); })
        .def("set_tree_cover", &State::set_tree_cover)
        .def("set_cover_from_image", [](State& state, py::array_t<float>& img, float& override_image_treecover) {
            shared_ptr<float[]> cover_image;
            int width, height;
            convert_from_numpy_array(img, cover_image, width, height);
            if (override_image_treecover >= 0) state.set_cover_from_image(cover_image, width, height, override_image_treecover);
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
        .def("get_tree_ages", [](State& state) {
            float* tree_ages = new float[state.population.size()];
            state.get_tree_ages(tree_ages);
            py::array_t<float> np_arr = as_1d_numpy_array(tree_ages, state.population.size());
            delete[] tree_ages;
            return np_arr;
        })
        .def_readwrite("grid", &State::grid)
        .def_readwrite("population", &State::population)
        .def_readwrite("initial_tree_cover", &State::initial_tree_cover);

    py::class_<Grid>(module, "Grid")
        .def(py::init<>())
        .def(py::init<const int&, const float&>())
        .def("get_tree_cover", &Grid::get_tree_cover)
		.def("reset_state_distr", &Grid::reset_state_distr)
		.def("redo_count", &Grid::redo_count)
        .def_readwrite("width", &Grid::width)
        .def_readwrite("width_r", &Grid::width_r)
        .def_readwrite("tree_cover", &Grid::tree_cover)
        .def("get_distribution", [](Grid& grid, int& collect_states) {
            shared_ptr<int[]> state_distribution = grid.get_state_distribution(collect_states);
            return as_2d_numpy_array(state_distribution, grid.width);
        })
        .def("set_grass_carrying_capacity", [](Grid& grid, py::array_t<float>& py_image) {
            shared_ptr<float[]> image;
            int _width, _height;
            convert_from_numpy_array(py_image, image, _width, _height);
            grid.set_grass_carrying_capacity(image);
        })
        .def("set_local_growth_multipliers", [](Grid& grid, py::array_t<float>& py_image) {
            shared_ptr<float[]> image;
            int _width, _height;
            convert_from_numpy_array(py_image, image, _width, _height);
            grid.set_local_growth_multipliers(image);
        })
        .def("get_aggr_tree_LAI_distribution", [](Grid& grid) {
            shared_ptr<float[]> aggr_tree_LAI_distribution = grid.get_aggr_tree_LAI_distribution();
            return as_2d_numpy_array(aggr_tree_LAI_distribution, grid.width);
        })
        .def("get_fuel_distribution", [](Grid& grid) {
            shared_ptr<float[]> fuel_load_distribution = grid.get_fuel_load_distribution();
            return as_2d_numpy_array(fuel_load_distribution, grid.width);
        });

    py::class_<Dynamics>(module, "Dynamics")
        .def(py::init<>())
        //.def(py::init<const int&, const float&, const float&, const float&, const float&, const float&, const float&,
        //    const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&, const float&,
        //    const float&, const map<string, map<string, float>>&, const int&, const float&, const float&, const int&, const int&, const int&,
        //    const float&, const int& >())
        .def(py::init<>())
        .def_readwrite("time", &Dynamics::time)
        .def_readwrite("state", &Dynamics::state)
        .def_readwrite("timestep", &Dynamics::timestep)
        .def_readwrite("seeds_produced", &Dynamics::seeds_produced)
        .def_readwrite("max_dbh", &Dynamics::max_dbh)
        .def("disperse_within_forest", [](Dynamics& dynamics, py::array_t<float>& img) {
            shared_ptr<float[]> mask;
            int width, height;
            convert_from_numpy_array(img, mask, width, height);
            dynamics.disperse_within_forest(mask);
        })
        .def("prune", [](Dynamics& dynamics, py::array_t<float>& img) {
            shared_ptr<float[]> mask;
            int width, height;
            convert_from_numpy_array(img, mask, width, height);
            dynamics.prune(mask);
        })
        .def("init_state", &Dynamics::init_state)
        .def("get_fires", [](Dynamics& dynamics) {
            py::array_t<float> np_arr = as_1d_numpy_array(dynamics.fires);
            return np_arr;
        })
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
        .def("get_basal_area", [](Dynamics& dynamics) {
            float basal_area = dynamics.get_basal_area();
            return basal_area;
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
        .def("get_patches", [](Dynamics& dynamics) {
            map<int, Patch> _forest_patches = dynamics.state.grid.forest_patches;
            map<int, Patch> _savanna_patches = dynamics.state.grid.savanna_patches;
            map<int, Patch> _patches = {};
			_patches.merge(_forest_patches);
            _patches.merge(_savanna_patches);
            Grid& grid = dynamics.state.grid;
            py::list patches;
            for (auto const& [id, _patch] : _patches) {
                py::dict patch;
                py::list perimeter = convert_from_idx_pair_vector_to_position_pylist(_patch.perimeter, dynamics.state.grid);
                py::list cells = convert_from_idx_vector_to_position_pylist(_patch.cells, dynamics.state.grid);

                // Assemble patch info into a dictionary
                patch["type"] = py::str(_patch.type);
                patch["id"] = py::int_(_patch.id);
                patch["area"] = py::float_(_patch.area);
                patch["cells"] = cells;
                patch["perimeter"] = perimeter;
                patch["perimeter_length"] = _patch.perimeter_length;
                patch["centroid"] = py::make_tuple(_patch.centroid_x, _patch.centroid_y);

                patches.append(patch);
			}
            return patches;
        })
        .def("get_forest_perimeter_length", [](Dynamics& dynamics) {
		    return dynamics.state.grid.get_forest_perimeter_length();
        })
        .def("get_perimeter_area_ratio", [](Dynamics& dynamics) {
            return dynamics.state.grid.get_perimeter_area_ratio();
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
		})
        .def("disperse", &Dynamics::disperse)
        .def("burn", &Dynamics::burn)
        .def("grow", &Dynamics::grow)
        .def("induce_background_mortality", &Dynamics::induce_background_mortality)
        .def("report_state", &Dynamics::report_state)
        .def("update_firefree_interval_averages", &Dynamics::update_firefree_interval_averages)
        .def("update_forest_patch_detection", &Dynamics::update_forest_patch_detection)
        ;

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