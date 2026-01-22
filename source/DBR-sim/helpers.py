import numpy as np
from rdfpy import rdf
import sys
import json
import argparse
from config import ParameterConfig


def add_kwargs(parser):
    """Dynamically add arguments using `ParameterConfig` defined in `config.py`.

    Arguments are dynamically loaded from a single config to reduce code duplication
    between cli, api and gui. Please view `config.py` > `ParameterConfig`
    for reference of all arguments and their configurations.

    Args:
        parser(argparse.ArgumentParser instance): The argument parser to add the arguments to.
    """
    for _, arg_config in ParameterConfig():
        keys = arg_config["keys"]["cli"]
        settings = arg_config["settings"]
        parser.add_argument(*keys, **settings)


def overwrite_from_global_arguments(args, global_args):
    """Overwrite arguments from global arguments.

    Args:
        args (dict): Dictionary of arguments to overwrite.
        global_args (dict): Dictionary of  global arguments.
    """
    for key, value in args.items():
        if global_args.get(key):
            args[key] = global_args[key]
    return args


def inject_unknown_args(args, unknown_args):
    # Process unknown arguments: assume form --key value
    it = iter(unknown_args)
    for token in it:
        if token.startswith("--"):
            key = token.lstrip("-")
            try:
                value = next(it)
                if value.startswith("--"):  # standalone flag
                    raise ValueError(f"Standalone flag detected ({value}, for key {key}) where value expected.")
                else:
                    if value.isnumeric():
                        value = int(value)
                    else:
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                    setattr(args, key, value)
            except StopIteration:
                setattr(args, key, True)
        else:
            raise ValueError(f"Unexpected token {token} in unknown arguments.")

    return args


def parse_args(parser=None):
    """Parse commandline arguments.

    Returns:
        kwargs(dict): Dictionary containing the cli args and their values.
    """
    if parser is None:
        # Create a new ArgumentParser instance if none is provided
        parser = argparse.ArgumentParser()
    
    add_kwargs(parser)
    args, unknown_args = parser.parse_known_args()

    # Handle unknown arguments (if any)
    if args.allow_unknown_args and unknown_args:
        kwargs = inject_unknown_args(args, unknown_args)
    else:
        # Re-start parser without unknown args allowed to raise error
        _ = parser.parse_args()
    
    kwargs = vars(args)
        
    return kwargs


def load_json_strings_if_any(kwargs):
    # Convert JSON strings to dicts
    for k, v in kwargs.items():
        if type(v) == str and "{" in v:
            kwargs[k] = json.loads(v) 
    return kwargs


class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)


def get_max(val1, val2):
    if val1 > val2:
        return val1
    else:
        return val2
    

def get_2d_dist(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def apply_user_args_to_configuration(args, cfg):
    for key, value in vars(args).items():
        setattr(cfg, key, value)
    return cfg


def compute_radial_distribution_function(dynamics, stepsize=0.02):
    tree_positions = dynamics.state.get_tree_positions()
    
    print("Computing radial distribution function...")
    g_r, radii = rdf(tree_positions, stepsize * dynamics.state.grid.width_r)
    print("Finished computing radial distribution function.")
    
    return g_r, radii


def is_number(inputString):
    return all(char.isdigit() for char in inputString)
