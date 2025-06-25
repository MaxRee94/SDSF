import numpy as np
from rdfpy import rdf
import sys


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
    

def compute_radial_distribution_function(dynamics, stepsize=0.02):
    tree_positions = dynamics.state.get_tree_positions()
    
    print("Computing radial distribution function...")
    g_r, radii = rdf(tree_positions, stepsize * dynamics.state.grid.width_r)
    print("Finished computing radial distribution function.")
    
    return g_r, radii

def is_number(inputString):
    return all(char.isdigit() for char in inputString)
