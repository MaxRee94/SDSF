import sys

from cv2 import textureFlattening
import visualization as vis
import numpy as np


def main(gridsize=None, treecover=None, cellsize=None, mean_radius=None):
    print("Launching app...")
    sys.path.append(r"F:\Development\DBR-sim\build")
    from x64.Release import dbr_cpp as cpp
    #from x64.Debug import dbr_cpp as cpp

    cpp.check_communication()

    state = cpp.State(gridsize[0], cellsize[0], mean_radius[0])
    state.set_tree_cover(treecover[0])
    distr = np.ndarray(shape=(1000, 1000), dtype=np.double)

    #help(cpp.get_distribution)
    print("Obtained grid distribution.")
    #print(distribution[600])
    vis.visualize(state.grid)

