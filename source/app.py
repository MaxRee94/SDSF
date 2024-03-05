import sys


def main(gridsize=None):
    print("Launching app...")
    print("Grid size:", gridsize)
    sys.path.append("../build")
    from x64.Release import dbr_cpp
    dbr_cpp.say_hello()


