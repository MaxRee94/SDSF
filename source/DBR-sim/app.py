import sys

from cv2 import textureFlattening
import visualization as vis
import numpy as np

sys.path.append(r"F:\Development\DBR-sim\build")
from x64.Release import dbr_cpp as cpp



def init(timestep=None, gridsize=None, cellsize=None, mean_radius=None, treecover=None, **_):
    print("timestep:", timestep)
    dynamics = cpp.Dynamics(timestep)
    dynamics.init_state(gridsize, cellsize, mean_radius)
    dynamics.state.set_tree_cover(treecover)
    
    return dynamics


def main(**kwargs):
    dynamics = init(**kwargs)
    cpp.check_communication()
    
    vis.visualize(dynamics.state.grid, image_width=kwargs["image_width"])

