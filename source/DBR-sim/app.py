from dataclasses import dataclass
import sys

import cv2
import visualization as vis
import file_handling as io
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
         rainfall=None, unsuppressed_flammability=None, 
         verbosity=None, radius_q1=None, radius_q2=None, seed_bearing_threshold=None,
         mass_budget_factor=None, dispersal_mode=None, linear_diffusion_q1=None, linear_diffusion_q2=None,
         dispersal_min=None, dispersal_max=None, growth_rate_multiplier=None, seed_mass=None,
         flammability_coefficients_and_constants=None,
         **_):
    cpp.check_communication()
    dynamics = cpp.Dynamics(
        timestep, cellsize, self_ignition_factor, rainfall, seed_bearing_threshold,
        mass_budget_factor, growth_rate_multiplier, unsuppressed_flammability, flammability_coefficients_and_constants[0],
        flammability_coefficients_and_constants[1], flammability_coefficients_and_constants[2], 
        flammability_coefficients_and_constants[3], max_radius, verbosity
    )
    dynamics.init_state(gridsize, radius_q1, radius_q2, seed_mass)
    dynamics.state.set_tree_cover(treecover)
    print("disp mode: ", dispersal_mode)
    if (dispersal_mode == "linear_diffusion"):
        dynamics.set_global_kernel(linear_diffusion_q1, linear_diffusion_q2, dispersal_min, dispersal_max)
    
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
    # dynamics.state.population.size() + 1
    no_colors = 1000
    color_dict = vis.get_color_dict(no_colors, begin=0.2, end=0.7)
    collect_states = True
    print("Visualizing state at t = 0")
    img = vis.visualize(
        dynamics.state.grid, user_args["image_width"], collect_states=collect_states,
        color_dict=color_dict
    )
    imagepath = os.path.join(str(Path(os.getcwd()).parent.parent), "data_out/image_timeseries/" + str(dynamics.time) + ".png")
    vis.save_image(img, imagepath)
    print("Beginning simulation...")
    datapath = None
    treecover_graph = vis.Graphs(dynamics, "Tree cover")
    while not termination_condition_satisfied(dynamics, start, user_args):
        dynamics.update()
        img = vis.visualize(
            dynamics.state.grid, user_args["image_width"], collect_states=collect_states,
            color_dict=color_dict
        )
        imagepath = os.path.join(str(Path(os.getcwd()).parent.parent), "data_out/image_timeseries/" + str(dynamics.time) + ".png")
        vis.save_image(img, imagepath, get_max(1000, img.shape[0]))
        datapath = io.export_state(dynamics, datapath)
        treecover_graph.update()

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
 


    

