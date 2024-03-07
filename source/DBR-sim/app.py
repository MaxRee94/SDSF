import sys


def main(gridsize=None):
    print("Launching app...")
    sys.path.append(r"F:\Development\DBR-sim\build")
    from x64.Release import dbr_cpp as cpp
    cpp.check_communication()

    state = cpp.State()
    #print("Grid size:", state.grid.gridsize)
    state.populate_grid()

