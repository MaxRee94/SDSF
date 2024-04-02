import sys

import cv2
import visualization as vis
import numpy as np
import time
import os
from pathlib import Path

from helpers import *

sys.path.append(r"F:\Development\DBR-sim\build")
#from x64.Debug import dbr_cpp as cpp
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
    condition = ""
    if (time.time() - start_time > user_args["timelimit"]):
        condition = f"Exceeded timelimit ({user_args['timelimit']} seconds)."
    if (dynamics.state.population.size() == 0):
        condition = f"Population collapse."
    #print("Pop size: ", dynamics.state.population.size())
    satisfied = len(condition) > 0
    if satisfied:
        print("Simulation terminated. Cause:", condition)

    return satisfied


def updateloop(dynamics, **user_args):
    start = time.time()
    is_within_timelimit = True
    print("Creating color dict...")
    color_dict = vis.get_color_dict(dynamics.state.population.size() + 1, begin=0.2, end=0.7)
    collect_states = True
    print("Visualizing state at t = 0")
    img = vis.visualize(
        dynamics.state.grid, user_args["image_width"], collect_states=collect_states,
        color_dict=color_dict
    )
    imagepath = os.path.join(str(Path(os.getcwd()).parent.parent), "data_out/image_timeseries/" + str(dynamics.time) + ".png")
    vis.save_image(img, imagepath)
    print("Beginning simulation...")
    while not termination_condition_satisfied(dynamics, start, user_args):
        dynamics.update()
        img = vis.visualize(
            dynamics.state.grid, user_args["image_width"], collect_states=collect_states,
            color_dict=color_dict
        )
        imagepath = os.path.join(str(Path(os.getcwd()).parent.parent), "data_out/image_timeseries/" + str(dynamics.time) + ".png")
        vis.save_image(img, imagepath)

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
 


    

