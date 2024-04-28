import numpy as np
from rdfpy import rdf



def get_max(val1, val2):
    if val1 > val2:
        return val1
    else:
        return val2
    

def compute_radial_distribution_function(dynamics, stepsize=0.02):
    tree_positions = dynamics.state.get_tree_positions()
    
    print("Computing radial distribution function...")
    g_r, radii = rdf(tree_positions, stepsize * dynamics.state.grid.width_r)
    print("Finished computing radial distribution function.")
    
    return g_r, radii
