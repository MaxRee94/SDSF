import sys
import visualization as vis
import numpy as np


def main(gridsize=None):
    print("Launching app...")
    sys.path.append(r"F:\Development\DBR-sim\build")
    from x64.Release import dbr_cpp as cpp
    cpp.check_communication()

    state = cpp.State(gridsize[0])
    #print("Grid size:", state.grid.gridsize)
    #state.set_tree_cover(0.5)
    state.populate_grid()
    distr = np.ndarray(shape=(1000, 1000), dtype=np.double)

    #help(cpp.get_distribution)
    distribution = state.grid.get_distribution()
    print(distribution[600])
    #vis.visualize(state.grid)

