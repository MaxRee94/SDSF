import sys

import cv2
import visualization as vis
import numpy as np
import time

sys.path.append(r"F:\Development\DBR-sim\build")
from x64.Release import dbr_cpp as cpp



def init(timestep=None, gridsize=None, cellsize=None, max_radius=None,
         treecover=None, self_ignition_factor=None, flammability=None,
         rainfall=None, unsuppressed_flammability=None, suppressed_flammability=None,
         verbosity=None, radius_q1=None, radius_q2=None,
         **_):
    cpp.check_communication()
    dynamics = cpp.Dynamics(
        timestep, unsuppressed_flammability, suppressed_flammability,
        self_ignition_factor, rainfall, verbosity
    )
    dynamics.init_state(gridsize, cellsize, max_radius, radius_q1, radius_q2)
    dynamics.state.set_tree_cover(treecover)
    
    return dynamics


def termination_condition_satisfied(dynamics, start_time, user_args):
    satisfied = False
    cause = ""
    if (time.time() - start_time > user_args["timelimit"]):
        cause = f"Exceeded timelimit ({user_args['timelimit']})"
    if (dynamics.state.population.size() == 0):
        cause = f"Population collapse."
    #print("Pop size: ", dynamics.state.population.size())
    satisfied = len(cause) > 0
    if satisfied:
        print("Termination condition reached. Cause:", cause)

    return satisfied


def updateloop(dynamics, **user_args):
    start = time.time()
    is_within_timelimit = True;
    while not termination_condition_satisfied(dynamics, start, user_args):
        dynamics.update()
        vis.visualize(dynamics.state.grid, user_args["image_width"])

    cv2.destroyAllWindows()


def do_tests(**user_args):
    tests = cpp.Tests()
    tests.run_all(True)


def main(**user_args):
    if user_args["test"] == "all":
        do_tests(**user_args)
        return
        
    dynamics = init(**user_args)
    updateloop(dynamics, **user_args)


    

