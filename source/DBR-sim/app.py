import sys


def main(gridsize=None):
    print("Launching app...")
    print("Grid size:", gridsize)
    sys.path.append(r"F:\Development\DBR-sim\build")
    #sys.path.append(".../build")
    print(sys.path)
    from x64.Release import dbr_cpp
    dbr_cpp.say_hello()


