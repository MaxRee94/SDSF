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
DATA_INTERNAL_DIR = REPOSITORY_BASEDIR + "/data_internal"  
BUILD_DIR = REPOSITORY_BASEDIR + "/build"
PERLIN_NOISE_DIR = DATA_IN_DIR + "/state patterns/perlin_noise"
CONTROLLED_PATTERN_DIR = DATA_IN_DIR + "/state patterns/controlled_patterns"
LEGEND_PATH = DATA_OUT_DIR + "/legends"


defaults = {
    "grid_width": 960,
    "treecover": 0.5,
    "cell_width": 1,
    "max_dbh": 44.3, # (Close to) Theoretical maximum due to density-dependent constraint on LAI (see 'Tree Allometric Relationships v03.xlsx' for details)
    "image_width": 1000,
    "timestep": 1,
    "patch_width": 170,
    "enforce_no_recruits": -1,
    "noise_octaves": 5,
    "timelimit": 1e32,
    "self_ignition_factor": 3,
    "unsuppressed_flammability": 0.5,
    "rainfall": 0.1,
    "tests": "none",
    "random_seed": -999,
    "firefreq_random_seed": 0,
    "termination_conditions": "all",
    "STR": 10000, # The number of seeds produced by a tree with a dbh of 30 cm
    "dbh_q1": 1,
    "dbh_q2": 0,
    "verbosity": 0,
    "growth_rate_multiplier_params": [0.5, 0.5, 2.13],
    "seed_bearing_threshold": 0.25, # dbh fraction of theoretical maximum height. We assume all trees are seed bearing beyond this height, based on Minor and Kobe (2018), Figure 5.
    "dispersal_mode": "all",
    "multi_disperser_params": f"multi_disperser_params.json",
    "dispersal_min": 0,
    "dispersal_max": 300,
    "growth_rate_multiplier": 0.4,
    "saturation_threshold": 3,
    "fire_resistance_params": "{\"argmin\": 8.5, \"argmax\": 50, \"stretch\": 2.857}", # See 'notes/fire_resistance_threshold_curve.xlsx' for details
    "background_mortality": 0.01,
    "csv_path": "",
    "headless": False,
    "max_timesteps": 100,
    "strategy_distribution_params": f"windkernel.json",
    "resource_grid_width": 48,
    "initial_pattern_image": "none",
    "mutation_rate": 0, # We do not incorporate mutation in this study.
    #"batch_parameters": "{\"control_variable\": \"growth_rate_multiplier_params-><idx>0\", \"control_value\": 0.99}"
    "batch_parameters": "",
    "report_state": False,
    "animal_group_size": 10,
    "image_size":(1000, 1000),
    "mean_radius":30,
    "cv_radius":0,
    "mean_distance":100,
    "cv_distance":0,
    "sine_amp1":5,
    "sine_wave1":6,
    "sine_amp2":2,
    "sine_wave2":12,
    "draw_stripes":False,
    "stripe_angle_deg":0,
    "stripe_mean_length":200,
    "enforce_distance_uniformity":True,
    "stripe_std_length":20,
    "sin_stripe":False,
    "sin_stripe_amp":10,
    "sin_stripe_wavelength":100,
    "grid_type":"square",
    "rotate_randomly":False,
    "suppress_distance_warning": False,  # If True, suppresses the warning about distance uniformity in the cpg script.
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
    "report_state": {
        "keys": {
            "cli": ["--report_state", "-repst"]
        },
        "settings": {
            "action": "store_true",
            "help": (
                "Export state report in each time step. Default false."
            ),
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
    "cell_width": {
        "keys": {
            "cli": ["--cell_width", "-cs"]
        },
        "settings": {
            "type": float,
            "help": (
                "The width (in meters) of each grid cell."
            ),
            "default": defaults["cell_width"],
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
            "default": 0,
        },
    },
    "max_dbh": {
        "keys": {
            "cli": ["--max_dbh", "-md"]
        },
        "settings": {
            "type": float,
            "help": (
                "The maximum dbh (in cm) of each tree in the initial timestep."
            ),
            "default": defaults["max_dbh"],
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
    "tests": {
        "keys": {
            "cli": ["--tests", "-tests"]
        },
        "settings": {
            "type": str,
            "help": (
                "Whether or not to run unit or output tests. Possible options: 'none', 'unit_tests', 'output_tests'."
            ),
            "default": defaults["tests"],
        },
    },
    "dbh_q1": {
        "keys": {
            "cli": ["--dbh_q1", "-rq1"]
        },
        "settings": {
            "type": float,
            "help": (
                "The relative probability of occurrence in the initial state of the smallest tree radius."
            ),
            "default": defaults["dbh_q1"],
        },
    },
    "dbh_q2": {
        "keys": {
            "cli": ["--dbh_q2", "-rq2"]
        },
        "settings": {
            "type": float,
            "help": (
                "The relative probability of occurrence in the initial state of the largest tree radius."
            ),
            "default": defaults["dbh_q2"],
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
            "type": str,
            "help": (
                ("Parameters to tree mortality sigmoid function, which takes bark thickness as input and returns survival probability.")
            ),
            "default": defaults["fire_resistance_params"],
        },
    },
    "background_mortality": {
        "keys": {
            "cli": ["--background_mortality", "-cm"]
        },
        "settings": {
            "type": float,
            "help": (
                "Constant background mortality rate, independent of fire risk."
            ),
            "default": defaults["background_mortality"],
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
    "resource_grid_width": {
        "keys": {
            "cli": ["--resource_grid_width", "-rgw"]
        },
        "settings": {
            "type": int,
            "help": (
                "The relative size (in number of cells along the vertical- or horizontal axis) of the resource grid versus the regular grid."
            ),
            "default": defaults["resource_grid_width"],
        },
    },
    "initial_pattern_image": {
        "keys": {
            "cli": ["--initial_pattern_image", "-ipi"]
        },
        "settings": {
            "type": str,
            "help": (
                "Path to an image file with an initial tree cover pattern."
            ),
            "default": defaults["initial_pattern_image"],
        },
    },
    "mutation_rate": {
        "keys": {
            "cli": ["--mutation_rate", "-murate"]
        },
        "settings": {
            "type": float,
            "help": (
                "The relative size (in number of cells along the vertical- or horizontal axis) of the resource grid versus the regular grid."
            ),
            "default": defaults["mutation_rate"],
        },
    },
    "STR": {
        "keys": {
            "cli": ["--STR", "-fcm"]
        },
        "settings": {
            "type": float,
            "help": (
                "Multiplier to artifically increase- or decrease the number seeds produced for all trees."
            ),
            "default": defaults["STR"],
        },
    },
    "termination_conditions": {
        "keys": {
            "cli": ["--termination_conditions", "-tcond"]
        },
        "settings": {
            "type": str,
            "help": (
                "Condition for terminating the simulation. Possible options: 'all', 'timelimit'."
            ),
            "default": defaults["termination_conditions"],
        },
    },
    "patch_width": {
        "keys": {
            "cli": ["--patch_width", "-pw"]
        },
        "settings": {
            "type": float,
            "help": (
                "Parameter to the perlin noise function."
            ),
            "default": defaults["patch_width"],
        },
    },
    "noise_octaves": {
        "keys": {
            "cli": ["--noise_octaves", "-noict"]
        },
        "settings": {
            "type": int,
            "help": (
                "Parameter to the perlin noise function."
            ),
            "default": defaults["noise_octaves"],
        },
    },
    "random_seed": {
        "keys": {
            "cli": ["--random_seed", "-rseed"]
        },
        "settings": {
            "type": int,
            "help": (
                "Random seed used by the c++ component of this program. A value of -999 (default) indicates each run will use a different unique seed."
            ),
            "default": defaults["random_seed"],
        },
    },
    "firefreq_random_seed": {
        "keys": {
            "cli": ["--firefreq_random_seed", "-rseedff"]
        },
        "settings": {
            "type": int,
            "help": (
                "Random seed used for the fire frequency distribution. If -999 is given, a random seed will be generated on the fly. Default: {}.".format(defaults["firefreq_random_seed"])
            ),
            "default": defaults["firefreq_random_seed"],
        },
    },
    "enforce_no_recruits": {
        "keys": {
            "cli": ["--enforce_no_recruits", "-enr"]
        },
        "settings": {
            "type": int,
            "help": (
                "Enforce a fraction of the total number of produced seeds to become recruits. -1 (default) means no specific fraction is enforced."
            ),
            "default": defaults["enforce_no_recruits"],
        },
    },
    "animal_group_size": {
        "keys": {
            "cli": ["--animal_group_size", "-ags"]
        },
        "settings": {
            "type": int,
            "help": (
                "The number of animals per simulated animal unit (a cluster/flock of animals), which moves and disperses seeds in unison."
            ),
            "default": defaults["animal_group_size"],
        },
    },
    "growth_rate_multiplier_params": {
        "keys": {
            "cli": ["--growth_rate_multiplier_params", "-grmp"]
        },
        "settings": {
            "nargs": "*",
            "type": float,
            "help": (
                ("Parameters to normal distribution (stdev and minimum) used to sample tree dbh growth rates.")
            ),
            "default": defaults["growth_rate_multiplier_params"],
        },
    },
    "batch_parameters": {
        "keys": {
            "cli": ["--batch_parameters", "-bp"]
        },
        "settings": {
            "type": str,
            "help": (
                ("Parameter settings set through a batch script.")
            ),
            "default": defaults["batch_parameters"],
        },
    },
    "image_size": {
        "keys": {
            "cli": ["--image_size", "-isz"]
        },
        "settings": {
            "type": tuple,
            "help": "Size of the generated image (width, height).",
            "default": defaults["image_size"],
        }
    },
    "mean_radius": {
        "keys": {
            "cli": ["--mean_radius", "-mr"]
        },
        "settings": {
            "type": float,
            "help": "Mean radius of the disks.",
            "default": defaults["mean_radius"],
        },
    },
    "cv_radius": {
        "keys": {
            "cli": ["--cv_radius", "-cvr"]
        },
        "settings": {
            "type": float,
            "help": "Coefficient of variation (CV) for the disk radius.",
            "default": defaults["cv_radius"],
        },
    },
    "mean_distance": {
        "keys": {
            "cli": ["--mean_distance", "-mdist"]
        },
        "settings": {
            "type": float,
            "help": "Mean distance between disks.",
            "default": defaults["mean_distance"],
        },
    },
    "cv_distance": {
        "keys": {
            "cli": ["--cv_distance", "-cvd"]
        },
        "settings": {
            "type": float,
            "help": "Coefficient of variation (CV) for the inter-disk distance.",
            "default": defaults["cv_distance"],
        },
    },
    "sine_amp1": {
        "keys": {
            "cli": ["--sine_amp1", "-sa1"]
        },
        "settings": {
            "type": float,
            "help": "Amplitude of the first sine wave used to modulate disk perimeters.",
            "default": defaults["sine_amp1"],
        },
    },
    "sine_wave1": {
        "keys": {
            "cli": ["--sine_wave1", "-sw1"]
        },
        "settings": {
            "type": float,
            "help": "Number of sine wave cycles (wave number) around disk perimeter for sine wave 1.",
            "default": defaults["sine_wave1"],
        },
    },
    "sine_amp2": {
        "keys": {
            "cli": ["--sine_amp2", "-sa2"]
        },
        "settings": {
            "type": float,
            "help": "Amplitude of the second sine wave used to modulate disk perimeters.",
            "default": defaults["sine_amp2"],
        },
    },
    "sine_wave2": {
        "keys": {
            "cli": ["--sine_wave2", "-sw2"]
        },
        "settings": {
            "type": float,
            "help": "Number of sine wave cycles (wave number) around disk perimeter for sine wave 2.",
            "default": defaults["sine_wave2"],
        },
    },
    "draw_stripes": {
        "keys": {
            "cli": ["--draw_stripes", "-ds"]
        },
        "settings": {
            "type": bool,
            "help": "Whether to draw convex-hull-based stripes on the disks.",
            "default": defaults["draw_stripes"],
        },
    },
    "stripe_angle_deg": {
        "keys": {
            "cli": ["--stripe_angle_deg", "-sad"]
        },
        "settings": {
            "type": float,
            "help": "Angle in degrees of the stripes relative to horizontal.",
            "default": defaults["stripe_angle_deg"],
        },
    },
    "stripe_mean_length": {
        "keys": {
            "cli": ["--stripe_mean_length", "-sml"]
        },
        "settings": {
            "type": float,
            "help": "Mean length of stripes drawn within each disk.",
            "default": defaults["stripe_mean_length"],
        },
    },
    "enforce_distance_uniformity": {
        "keys": {
            "cli": ["--enforce_distance_uniformity", "-edu"]
        },
        "settings": {
            "type": bool,
            "help": "Whether to enforce uniform spacing between disk centers.",
            "default": defaults["enforce_distance_uniformity"],
        },
    },
    "stripe_std_length": {
        "keys": {
            "cli": ["--stripe_std_length", "-ssl"]
        },
        "settings": {
            "type": float,
            "help": "Standard deviation of stripe lengths.",
            "default": defaults["stripe_std_length"],
        },
    },
    "sin_stripe": {
        "keys": {
            "cli": ["--sin_stripe", "-ss"]
        },
        "settings": {
            "type": bool,
            "help": "If True, apply sinusoidal curvature to the stripes.",
            "default": defaults["sin_stripe"],
        },
    },
    "sin_stripe_amp": {
        "keys": {
            "cli": ["--sin_stripe_amp", "-ssa"]
        },
        "settings": {
            "type": float,
            "help": "Amplitude of sinusoidal modulation of stripes.",
            "default": defaults["sin_stripe_amp"],
        },
    },
    "sin_stripe_wavelength": {
        "keys": {
            "cli": ["--sin_stripe_wavelength", "-ssw"]
        },
        "settings": {
            "type": float,
            "help": "Wavelength of sinusoidal modulation of stripes.",
            "default": defaults["sin_stripe_wavelength"],
        },
    },
    "grid_type": {
        "keys": {
            "cli": ["--grid_type", "-gt"]
        },
        "settings": {
            "type": str,
            "help": "Type of grid layout to use for disk placement ('hex' or 'square').",
            "default": defaults["grid_type"],
        },
    },
    "rotate_randomly": {
        "keys": {
            "cli": ["--rotate_randomly", "-rr"]
        },
        "settings": {
            "type": bool,
            "help": "If True, randomly rotate disks to avoid orientation artifacts.",
            "default": defaults["rotate_randomly"],
        }
    },
    "suppress_distance_warning": {
        "keys": {
            "cli": ["--suppress_distance_warning", "-sdw"]
        },
        "settings": {
            "type": bool,
            "help": "If True, no error is raised if the given mean_distance (used for generating disk patterns) is automatically changed to generate a perfect grid.",
            "default": defaults["suppress_distance_warning"],
        }
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


def get_all_defaults():
    """Return all default parameters and their values.

    Returns:
        dict: A dictionary containing all default parameters and their values.
    """
    return defaults.copy()