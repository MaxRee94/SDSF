from dataclasses import dataclass
import sys

import cv2
import numpy as np
import time
import os
from pathlib import Path
from matplotlib import pyplot as plt
import json

import visualization as vis
import file_handling as io
from helpers import *
from config import *

sys.path.append(BUILD_DIR)
#from x64.Debug import dbr_cpp as cpp
from x64.Release import dbr_cpp as cpp


def set_dispersal_kernel(
        dynamics, dispersal_mode, linear_diffusion_q1, linear_diffusion_q2, dispersal_min, dispersal_max,
        wind_dispersal_params, multi_disperser_params, animal_dispersal_params
    ):
    if (dispersal_mode == "linear_diffusion"):
        dynamics.set_global_linear_kernel(linear_diffusion_q1, linear_diffusion_q2, dispersal_min, dispersal_max)
    elif (dispersal_mode == "wind"):
        dynamics.set_global_wind_kernel(*wind_dispersal_params);
    elif (dispersal_mode == "animal"):
        dynamics.set_global_animal_kernel(animal_dispersal_params)
    elif (dispersal_mode == "all"):
        with open(multi_disperser_params, "r") as mdp_jsonfile:
            multi_disperser_params = json.load(mdp_jsonfile)
        animal_dispersal_params = multi_disperser_params["animal"]
        multi_disperser_params.pop("animal")
        dynamics.set_global_kernels(multi_disperser_params, animal_dispersal_params)

    return dynamics


def init(
    timestep=None, grid_width=None, cellsize=None, max_radius=None, image_width=None,
    treecover=None, self_ignition_factor=None, flammability=None,
    rainfall=None, unsuppressed_flammability=None, 
    verbosity=None, radius_q1=None, radius_q2=None, seed_bearing_threshold=None,
    mass_budget_factor=None, dispersal_mode=None, linear_diffusion_q1=None, linear_diffusion_q2=None,
    dispersal_min=None, dispersal_max=None, growth_rate_multiplier=None, seed_mass=None,
    flammability_coefficients_and_constants=None, saturation_threshold=None, fire_resistance_params=None,
    constant_mortality=None, headless=False, wind_dispersal_params=None, animal_dispersal_params=None,
    multi_disperser_params=None, strategy_distribution_params=None, resource_grid_relative_size=None,
    **user_args
    ):
    
    # Obtain strategy distribution parameters
    with open(strategy_distribution_params, "r") as sdp_jsonfile:
        strategy_distribution_params = json.load(sdp_jsonfile)

    # Initialize dynamics object and state
    dynamics = cpp.Dynamics(
        timestep, cellsize, self_ignition_factor, rainfall, seed_bearing_threshold,
        mass_budget_factor, growth_rate_multiplier, unsuppressed_flammability, flammability_coefficients_and_constants[0],
        flammability_coefficients_and_constants[1], flammability_coefficients_and_constants[2], 
        flammability_coefficients_and_constants[3], max_radius, saturation_threshold, fire_resistance_params[0],
        fire_resistance_params[1], fire_resistance_params[2], constant_mortality, strategy_distribution_params, 
        resource_grid_relative_size, verbosity
    )
    dynamics.init_state(grid_width, radius_q1, radius_q2, seed_mass)
    dynamics.state.set_tree_cover(treecover)
    dynamics = set_dispersal_kernel(
        dynamics, dispersal_mode, linear_diffusion_q1, linear_diffusion_q2, dispersal_min, dispersal_max,
        wind_dispersal_params, multi_disperser_params, animal_dispersal_params
    )
    
    # Create a color dictionary
    no_colors = 100
    color_dict = vis.get_color_dict(no_colors, begin=0.3, end=0.6)
    
    # Visualize the initial state
    collect_states = True
    print("Visualizing state at t = 0")
    if not headless:
        # Get a color image representation of the initial state and show it.
        img = vis.visualize(
        dynamics.state.grid, image_width, collect_states=collect_states,
        color_dict=color_dict
        )
    else:
        # Get a color image representation of the initial state
        img = vis.get_image(dynamics.state.grid, collect_states, color_dict)
   
    # Export image file
    imagepath = os.path.join(DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
    vis.save_image(img, imagepath)

    return dynamics, color_dict


def termination_condition_satisfied(dynamics, start_time, user_args):
    satisfied = False
    condition = ""
    if (time.time() - start_time > user_args["timelimit"]):
        condition = f"Exceeded timelimit ({user_args['timelimit']} seconds)."
    if (dynamics.state.population.size() == 0):
        condition = "Population collapse."
    if (dynamics.state.grid.tree_cover >= 0.95):
        condition = "Tree cover converged to >95%."
    if (dynamics.state.grid.tree_cover <= 0.05):
        condition = "Tree cover converged to <5%."
    if (dynamics.time >= user_args["max_timesteps"]):
        condition = f"Maximum number of timesteps ({user_args['max_timesteps']}) reached."
    #print("Pop size: ", dynamics.state.population.size())
    satisfied = len(condition) > 0
    if satisfied:
        print("Simulation terminated. Cause:", condition)

    return satisfied


def updateloop(dynamics, color_dict, **user_args):
    start = time.time()
    print("Beginning simulation...")
    csv_path = user_args["csv_path"]
    init_csv = True
    collect_states = True
    if not user_args["headless"]:
        graphs = vis.Graphs(dynamics)
    while not termination_condition_satisfied(dynamics, start, user_args):
        dynamics.update()
        if user_args["headless"]:
            # Get a color image representation of the initial state
            img = vis.get_image(dynamics.state.grid, collect_states, color_dict)
        else:
            # Get a color image representation of the initial state and show it.
            img = vis.visualize(
                dynamics.state.grid, user_args["image_width"], collect_states=collect_states,
                color_dict=color_dict
            )
        imagepath = os.path.join(DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
        vis.save_image(img, imagepath, get_max(1000, img.shape[0]))
        csv_path = io.export_state(dynamics, csv_path, init_csv)
        init_csv = False
        if not user_args["headless"]:
            graphs.update()

    if not user_args["headless"]:
        cv2.destroyAllWindows()
    
    return dynamics


def do_tests(**user_args):
    tests = cpp.Tests()
    tests.run_all(True)
 

def main(**user_args):
    # cpp.init_RNG()
    # dist_max = 200
    # windspeed_gmean = 20
    # windspeed_stdev = 6
    # seed_terminal_speed = 1.5
    # abscission_height = 15
    # wind_kernel = cpp.Kernel(1, dist_max, windspeed_gmean, windspeed_stdev, seed_terminal_speed, abscission_height)
    # vis.visualize_kernel(wind_kernel)
    # return

    if user_args["test"] == "all":
        do_tests(**user_args)
        return
        
    dynamics, color_dict = init(**user_args)
    return updateloop(dynamics, color_dict, **user_args)
 


    

