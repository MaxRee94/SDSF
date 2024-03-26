"""DBR-sim defaults and constants."""
import argparse
import string

defaults = {
    "gridsize": 1000,
    "treecover": 0.5,
    "cellsize": 1.5,
    "max_radius": 6.0,
    "image_width": 1000,
    "timestep": 1,
    "timelimit": 1e8,
    "self_ignition_factor": 20.0,
    "unsuppressed_flammability": -999.0,
    "suppressed_flammability": -999.0,
    "rainfall": 0.1,
    "test": "none",
    "radius_q1": 1,
    "radius_q2": 0,
    "verbosity": 0,
}

gui_defaults = {
    "gridsize": defaults["gridsize"],
    "treecover": defaults["treecover"],
    "setting2": [
        "default1",
        "default2",
        "default3"
    ]
}

_parameter_config = {
    "gridsize": {
        "keys": {
            "cli": ["--gridsize", "-gs"]
        },
        "settings": {
            "type": int,
            "help": (
                "The number of grid cells to be used along one axis of the spatial domain."
            ),
            "default": defaults["gridsize"],
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
    "suppressed_flammability": {
        "keys": {
            "cli": ["--suppressed_flammability", "-sf"]
        },
        "settings": {
            "type": float,
            "help": (
                "Flammability of forest cells that suppress fire (WARNING: ideally, this should not be a constant value but should be contingent on rainfall levels. Interpret results with care.)."
            ),
            "default": defaults["suppressed_flammability"],
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
                "The probability of occurrence in the initial state of the smallest tree radius."
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
                "The probability of occurrence in the initial state of the largest tree radius."
            ),
            "default": defaults["radius_q2"],
        },
    },
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
