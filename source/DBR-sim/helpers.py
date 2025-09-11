import numpy as np
from rdfpy import rdf
import sys
import json
from argparse import ArgumentParser
from config import ParameterConfig



def add_kwargs(parser):
    """Dynamically add arguments using `ParameterConfig` defined in `config.py`.

    Arguments are dynamically loaded from a single config to reduce code duplication
    between cli, api and gui. Please view `config.py` > `ParameterConfig`
    for reference of all arguments and their configurations.

    Args:
        parser(ArgumentParser instance): The argument parser to add the arguments to.
    """
    for _, arg_config in ParameterConfig():
        keys = arg_config["keys"]["cli"]
        settings = arg_config["settings"]
        parser.add_argument(*keys, **settings)


def parse_args(parser=None):
    """Parse commandline arguments.

    Returns:
        kwargs(dict): Dictionary containing the cli args and their values.
    """
    if parser is None:
        # Create a new ArgumentParser instance if none is provided
        parser = ArgumentParser()
    
    add_kwargs(parser)
    kwargs = vars(parser.parse_args())
    
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


def compute_radial_distribution_function(dynamics, stepsize=0.02):
    tree_positions = dynamics.state.get_tree_positions()
    
    print("Computing radial distribution function...")
    g_r, radii = rdf(tree_positions, stepsize * dynamics.state.grid.width_r)
    print("Finished computing radial distribution function.")
    
    return g_r, radii

def is_number(inputString):
    return all(char.isdigit() for char in inputString)
