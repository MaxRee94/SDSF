"""DBR-sim defaults and constants."""
import argparse
import string
import json
import os


# global constants
cwd = os.getcwd()
if cwd.endswith("source"):
    cwd = cwd + "/DBR-sim"
REPOSITORY_BASEDIR = os.path.dirname(os.path.dirname(cwd))
DATA_IN_DIR = REPOSITORY_BASEDIR + "/data_in"
DATA_OUT_DIR = REPOSITORY_BASEDIR + "/data_out"
BUILD_DIR = REPOSITORY_BASEDIR + "/build"


defaults = {
    "grid_width": 1000,
    "treecover": 0.5,
    "cellsize": 1,
    "max_radius": 6.0,
    "image_width": 1000,
    "timestep": 1,
    "timelimit": 1e32,
    "self_ignition_factor": 140,
    "unsuppressed_flammability": 0.4,
    "rainfall": 0.1,
    "test": "none",
    "radius_q1": 1,
    "radius_q2": 0,
    "verbosity": 0,
    "seed_bearing_threshold": 0.4,
    "mass_budget_factor": 1,
    "dispersal_mode": "all",
    "linear_diffusion_q1": 1,
    "linear_diffusion_q2": 0,
    "wind_dispersal_params": [30, 5, 0.1, 15],
    "animal_dispersal_params": '{"bird": {"popsize": 20, "fruitsize_pref": 0.02}}',
    "multi_disperser_params": f"{DATA_IN_DIR}/multi_disperser_params.json",
    "dispersal_min": 0,
    "dispersal_max": 300,
    "growth_rate_multiplier": 0.1,
    "seed_mass": 0.01,
    "flammability_coefficients_and_constants": [0.1, 0.35, 0, 1],
    "saturation_threshold": 3,
    "fire_resistance_params": [8.5, 50, 2.857],
    "constant_mortality": 0.03,
    "csv_path": "",
    "headless": False,
    "max_timesteps": 1e9,
    "resource_grid_relative_size": f"{DATA_IN_DIR}/strategy_distribution_params.json",
    "resource_grid_relative_size": 0.1,
}

gui_defaults = {
    "grid_width": defaults["grid_width"],
    "treecover": defaults["treecover"],
    "setting2": [
        "default1",
        "default2",
        "default3"
    ]
}

_parameter_config = {
    "grid_width": {
        "keys": {
            "cli": ["--grid_width", "-gw"]
        },
        "settings": {
            "type": int,
            "help": (
                "The width (measured by the number of cells along the horizontal- or vertical axis) of the spatial domain."
            ),
            "default": defaults["grid_width"],
        },
    },
    "treecover": {
        "keys": {
            "cli": ["--treecover", "-tc"]
        },
        "settings": {
            "type": float,
            "help": (
                "The minimal fraction of the spatial domain occupied by tree cells."
            ),
            "default": defaults["treecover"],
        },
    },
    "cellsize": {
        "keys": {
            "cli": ["--cellsize", "-cs"]
        },
        "settings": {
            "type": float,
            "help": (
                "The width (in meters) of each grid cell."
            ),
            "default": defaults["cellsize"],
        },
    },
    "verbosity": {
        "keys": {
            "cli": ["--verbosity", "-vb"]
        },
        "settings": {
            "type": int,
            "help": (
                "Controls the amount of informational progress updates that are printed to the terminal."
            ),
        },
    },
    "max_radius": {
        "keys": {
            "cli": ["--max_radius", "-mr"]
        },
        "settings": {
            "type": float,
            "help": (
                "The mean radius (in meters) of each tree in the initial timestep."
            ),
            "default": defaults["max_radius"],
        },
    },
    "image_width": {
        "keys": {
            "cli": ["--image_width", "-iw"]
        },
        "settings": {
            "type": int,
            "help": (
                "The width (in pixels) of the viewer image (which displays the simulated spatial domain)."
            ),
            "default": defaults["image_width"],
        },
    },
    "timestep": {
        "keys": {
            "cli": ["--timestep", "-ts"]
        },
        "settings": {
            "type": int,
            "help": (
                "The duration of a single timestep, in whole years."
            ),
            "default": defaults["timestep"],
        },
    },
    "timelimit": {
        "keys": {
            "cli": ["--timelimit", "-tl"]
        },
        "settings": {
            "type": int,
            "help": (
                "The maxmimum number of seconds that the simulation is allowed to run before termination."
            ),
            "default": defaults["timelimit"],
        },
    },
    "self_ignition_factor": {
        "keys": {
            "cli": ["--self_ignition_factor", "-sif"]
        },
        "settings": {
            "type": float,
            "help": (
                "Determines the number of cells ignited per year by multiplication with rainfall levels."
            ),
            "default": defaults["self_ignition_factor"],
        },
    },
    "unsuppressed_flammability": {
        "keys": {
            "cli": ["--unsuppressed_flammability", "-uf"]
        },
        "settings": {
            "type": float,
            "help": (
                "Flammability of savanna cells and forest cells that do not suppress fire (WARNING: ideally, this should not be a constant value but should be contingent on rainfall levels. Interpret results with care.)."
            ),
            "default": defaults["unsuppressed_flammability"],
        },
    },
    "rainfall": {
        "keys": {
            "cli": ["--rainfall", "-rf"]
        },
        "settings": {
            "type": float,
            "help": (
                "Rainfall levels, kept constant across all timesteps."
            ),
            "default": defaults["rainfall"],
        },
    },
    "test": {
        "keys": {
            "cli": ["--test", "-test"]
        },
        "settings": {
            "type": str,
            "help": (
                "Whether or not to run tests. Possible options: 'none', 'all'."
            ),
            "default": defaults["test"],
        },
    },
    "radius_q1": {
        "keys": {
            "cli": ["--radius_q1", "-rq1"]
        },
        "settings": {
            "type": float,
            "help": (
                "The relative probability of occurrence in the initial state of the smallest tree radius."
            ),
            "default": defaults["radius_q1"],
        },
    },
    "radius_q2": {
        "keys": {
            "cli": ["--radius_q2", "-rq2"]
        },
        "settings": {
            "type": float,
            "help": (
                "The relative probability of occurrence in the initial state of the largest tree radius."
            ),
            "default": defaults["radius_q2"],
        },
    },
    "seed_bearing_threshold": {
        "keys": {
            "cli": ["--seed_bearing_threshold", "-sbt"]
        },
        "settings": {
            "type": float,
            "help": (
                "The fraction of the tree maximum radius above which a tree becomes seed-bearing."
            ),
            "default": defaults["seed_bearing_threshold"],
        },
    },
    "mass_budget_factor": {
        "keys": {
            "cli": ["--mass_budget_factor", "-mbf"]
        },
        "settings": {
            "type": float,
            "help": (
                "Factor determining (along with radius) the total annual tree seed biomass."
            ),
            "default": defaults["mass_budget_factor"],
        },
    },
    "dispersal_mode": {
        "keys": {
            "cli": ["--dispersal_mode", "-dm"]
        },
        "settings": {
            "type": str,
            "help": (
                "Dispersal mode. Currently only 'linear diffusion' is possible, which uses a linear probability model to sample dispersal distances, in uniformly sampled directions."
            ),
            "default": defaults["dispersal_mode"],
        },
    },
    "linear_diffusion_q1": {
        "keys": {
            "cli": ["--linear_diffusion_q1", "-ldq1"]
        },
        "settings": {
            "type": float,
            "help": (
                "The relative probability of sampling short distances when using a linear diffusion model for seed dispersal."
            ),
            "default": defaults["linear_diffusion_q1"],
        },
    },
    "linear_diffusion_q2": {
        "keys": {
            "cli": ["--linear_diffusion_q2", "-ldq2"]
        },
        "settings": {
            "type": float,
            "help": (
                "The relative probability of sampling short distances when using a linear diffusion model for seed dispersal."
            ),
            "default": defaults["linear_diffusion_q2"],
        },
    },
    "wind_dispersal_params": {
        "keys": {
            "cli": ["--wind_dispersal_params", "-wdp"]
        },
        "settings": {
            "nargs": "*",
            "type": float,
            "help": (
                ("Parameters for the wind dispersal model, in the order [windspeed geometric mean, windspeed stdev, seed terminal speed, abscission height].")
            ),
            "default": defaults["wind_dispersal_params"],
        },
    },
    "dispersal_min": {
        "keys": {
            "cli": ["--dispersal_min", "-dmin"]
        },
        "settings": {
            "type": float,
            "help": (
                "The minimum dispersal distance."
            ),
            "default": defaults["dispersal_min"],
        },
    },
    "dispersal_max": {
        "keys": {
            "cli": ["--dispersal_max", "-dmax"]
        },
        "settings": {
            "type": float,
            "help": (
                "The maximum dispersal distance."
            ),
            "default": defaults["dispersal_max"],
        },
    },
    "growth_rate_multiplier": {
        "keys": {
            "cli": ["--growth_rate_multiplier", "-grm"]
        },
        "settings": {
            "type": float,
            "help": (
                "Constant multiplied with the growth rates of all trees."
            ),
            "default": defaults["growth_rate_multiplier"],
        },
    },
    "seed_mass": {
        "keys": {
            "cli": ["--seed_mass", "-sm"]
        },
        "settings": {
            "type": float,
            "help": (
                "Seed mass."
            ),
            "default": defaults["seed_mass"],
        },
    },
    "flammability_coefficients_and_constants": {
        "keys": {
            "cli": ["--flammability_coefficients_and_constants", "-fcac"]
        },
        "settings": {
            "nargs": "*",
            "type": float,
            "help": (
                ("Coefficients and constants [minimum_flammability, maximum_flammability, minimum_radius, radius_range] which determine the flammability" +
                 " function's output given the tree radius.")
            ),
            "default": defaults["flammability_coefficients_and_constants"],
        },
    },
    "saturation_threshold": {
        "keys": {
            "cli": ["--saturation_threshold", "-st"]
        },
        "settings": {
            "type": float,
            "help": (
                "The minimum number of trees in a cell to prohibit seed germination due to competition for locally available resources."
            ),
            "default": defaults["saturation_threshold"],
        },
    },
    "fire_resistance_params": {
        "keys": {
            "cli": ["--fire_resistance_params", "-frp"]
        },
        "settings": {
            "nargs": "*",
            "type": float,
            "help": (
                ("Parameters to tree mortality sigmoid function, which takes bark thickness as input and returns survival probability.")
            ),
            "default": defaults["fire_resistance_params"],
        },
    },
    "constant_mortality": {
        "keys": {
            "cli": ["--constant_mortality", "-cm"]
        },
        "settings": {
            "type": float,
            "help": (
                "Constant background mortality rate, independent of fire risk."
            ),
            "default": defaults["constant_mortality"],
        },
    },
    "headless": {
        "keys": {
            "cli": ["--headless", "-hl"]
        },
        "settings": {
            "action": "store_true",
            "help": (
                "Run DBR-sim in headless mode."
            ),
        },
    },
    "csv_path": {
        "keys": {
            "cli": ["--csv_path", "-csv"]
        },
        "settings": {
            "type": str,
            "help": (
                "CSV path to export data to."
            ),
            "default": defaults["csv_path"],
        },
    },
    "max_timesteps": {
        "keys": {
            "cli": ["--max_timesteps", "-mt"]
        },
        "settings": {
            "type": str,
            "help": (
                "Maximum number of timesteps to run simulation for."
            ),
            "default": defaults["max_timesteps"],
        },
    },
    "animal_dispersal_params": {
        "keys": {
            "cli": ["--animal_dispersal_params", "-adp"]
        },
        "settings": {
            "nargs": "*",
            "type": json.loads,
            "help": (
                ('Parameters for animal, to be entered as a json string in the format {"animal1": {"popsize": <your int>, "fruitsize_pref": <your float>}}.')
            ),
            "default": defaults["animal_dispersal_params"],
        },
    },
    "multi_disperser_params": {
        "keys": {
            "cli": ["--multi_disperser_params", "-mdp"]
        },
        "settings": {
            "type": str,
            "help": (
                "Path to a json file containing parameters for multiple dispersers."
            ),
            "default": defaults["multi_disperser_params"],
        },
    },
    "strategy_distribution_params": {
        "keys": {
            "cli": ["--strategy_distribution_params", "-sdp"]
        },
        "settings": {
            "type": str,
            "help": (
                "Path to a json file containing parameters for the distribution of strategies."
            ),
            "default": defaults["strategy_distribution_params"],
        },
    },
    "resource_grid_relative_size": {
        "keys": {
            "cli": ["--resource_grid_relative_size", "-rgrs"]
        },
        "settings": {
            "type": float,
            "help": (
                "The relative size (in number of cells along the vertical- or horizontal axis) of the resource grid versus the regular grid."
            ),
            "default": defaults["resource_grid_relative_size"],
        },
    }
}


class ParameterConfig():
    """Wrapper class for _parameter_config."""

    def __init__(self):
        """Initialize and load data from _parameter_config."""
        self.data = {}
        self.load()

    def load(self):
        """Update `self.data` by loading parameter config data."""
        for name, cfg in _parameter_config.copy().items():
            _cfg = cfg.copy()

            self.data[name] = _cfg

    def __iter__(self):
        """Yield all parameter names and -configurations.

        Yields:
            name(str): Name of the parameter.
            cfg(dict): Contains the configuration of the parameter.
        """
        for name, cfg in self.data.items():
            yield name, cfg

    def get(self, name):
        """Return the configuration associated with the given name of the parameter.

        Args:
            name(str): The name of the parameter to return the configuration for.
        Returns:
            (dict): The configuration associated with the name of the parameter.
        """
        return self.data[name]
