from dataclasses import dataclass
import sys
from turtle import color

import cv2
import numpy as np
import time
import os
from pathlib import Path
from matplotlib import pyplot as plt
import json
import scipy.stats as stats
import math

import visualization as vis
import file_handling as io
from helpers import *
from config import *

sys.path.append(BUILD_DIR)
#from x64.Debug import dbr_cpp as cpp
from x64.Release import dbr_cpp as cpp


def set_dispersal_kernel(
        dynamics, dispersal_mode, multi_disperser_params
    ):
    with open(multi_disperser_params, "r") as mdp_jsonfile:
        multi_disperser_params = json.load(mdp_jsonfile)
    if (dispersal_mode == "linear_diffusion"):
        dynamics.set_global_linear_kernel(
            multi_disperser_params["linear"]["q1"], multi_disperser_params["linear"]["q2"],
            multi_disperser_params["linear"]["min"], multi_disperser_params["linear"]["max"]
        )
    elif (dispersal_mode == "wind"):
        dynamics.set_global_wind_kernel(
            multi_disperser_params["wind"]["wspeed_gmean"], multi_disperser_params["wind"]["wspeed_stdev"],
            multi_disperser_params["wind"]["wind_direction"], multi_disperser_params["wind"]["wind_direction_stdev"]
            
        );
    elif (dispersal_mode == "animal"):
        dynamics.set_global_animal_kernel(multi_disperser_params["animal"])
    elif (dispersal_mode == "all"):
        animal_dispersal_params = multi_disperser_params["animal"]
        multi_disperser_params.pop("animal")
        dynamics.set_global_kernels(multi_disperser_params, animal_dispersal_params)

    return dynamics


def init_tests(
    timestep=None, grid_width=None, cellsize=None, max_dbh=None, image_width=None,
    treecover=None, self_ignition_factor=None, flammability=None,
    rainfall=None, unsuppressed_flammability=None, 
    verbosity=None, dbh_q1=None, dbh_q2=None, seed_bearing_threshold=None,
    dispersal_mode=None, linear_diffusion_q1=None, linear_diffusion_q2=None,
    dispersal_min=None, dispersal_max=None, growth_rate_multiplier=None, 
    flammability_coefficients_and_constants=None, saturation_threshold=None, fire_resistance_params=None,
    constant_mortality=None, headless=False, wind_dispersal_params=None, animal_dispersal_params=None,
    multi_disperser_params=None, strategy_distribution_params=None, resource_grid_relative_size=None,
    initial_pattern_image=None, mutation_rate=None, **user_args
):
    # Obtain strategy distribution parameters
    with open(strategy_distribution_params, "r") as sdp_jsonfile:
        strategy_distribution_params = json.load(sdp_jsonfile)    

    tests = cpp.Tests(timestep, cellsize, self_ignition_factor, rainfall, seed_bearing_threshold,
        growth_rate_multiplier, unsuppressed_flammability, flammability_coefficients_and_constants[0],
        flammability_coefficients_and_constants[1], flammability_coefficients_and_constants[2], 
        flammability_coefficients_and_constants[3], max_dbh, saturation_threshold, fire_resistance_params[0],
        fire_resistance_params[1], fire_resistance_params[2], constant_mortality, strategy_distribution_params, 
        resource_grid_relative_size, mutation_rate, verbosity, grid_width, dbh_q1, dbh_q2)
    
    return tests


def init(
    timestep=None, grid_width=None, cellsize=None, max_dbh=None, image_width=None,
    treecover=None, self_ignition_factor=None, flammability=None,
    rainfall=None, unsuppressed_flammability=None, 
    verbosity=None, dbh_q1=None, dbh_q2=None, seed_bearing_threshold=None,
    dispersal_mode=None, linear_diffusion_q1=None, linear_diffusion_q2=None,
    dispersal_min=None, dispersal_max=None, growth_rate_multiplier=None, 
    flammability_coefficients_and_constants=None, saturation_threshold=None, fire_resistance_params=None,
    constant_mortality=None, headless=False, wind_dispersal_params=None, animal_dispersal_params=None,
    multi_disperser_params=None, strategy_distribution_params=None, resource_grid_relative_size=None,
    initial_pattern_image=None, mutation_rate=None,
    **user_args
    ):
    
    # Obtain strategy distribution parameters
    with open(strategy_distribution_params, "r") as sdp_jsonfile:
        strategy_distribution_params = json.load(sdp_jsonfile)

    # Initialize dynamics object and state
    dynamics = cpp.Dynamics(
        timestep, cellsize, self_ignition_factor, rainfall, seed_bearing_threshold,
        growth_rate_multiplier, unsuppressed_flammability, flammability_coefficients_and_constants[0],
        flammability_coefficients_and_constants[1], flammability_coefficients_and_constants[2], 
        flammability_coefficients_and_constants[3], max_dbh, saturation_threshold, fire_resistance_params[0],
        fire_resistance_params[1], fire_resistance_params[2], constant_mortality, strategy_distribution_params, 
        resource_grid_relative_size, mutation_rate, verbosity
    )
    dynamics.init_state(grid_width, dbh_q1, dbh_q2)
    
    # Set dispersal kernel
    dynamics = set_dispersal_kernel(dynamics, dispersal_mode, multi_disperser_params)
    
    # Set initial tree cover
    if initial_pattern_image == "none":
        dynamics.state.set_tree_cover(treecover)
    else:
        img = cv2.imread(f"{DATA_IN_DIR}/state patterns/" + initial_pattern_image, cv2.IMREAD_GRAYSCALE)
        img = cv2.resize(img, (dynamics.state.grid.width, dynamics.state.grid.width), interpolation=cv2.INTER_LINEAR)
        #(thresh, img) = cv2.threshold(img, 128, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
        print("Setting tree cover from image...")
        print(img / 255)
        dynamics.state.set_cover_from_image(img / 255)
    dynamics.state.repopulate_grid(0)
    
    # Create a color dictionary
    no_colors = 100
    color_dict = vis.get_color_dict(no_colors, begin=0.2, end=0.5)
    color_dict_recruitment = vis.get_color_dict(no_colors, begin=0.2, end=0.5, distr_type="recruitment")
    color_dict_fire_freq = vis.get_color_dict(10, begin=0.2, end=0.5, distr_type="fire_freq")
    color_dict_blackwhite = vis.get_color_dict(no_colors, distr_type="blackwhite")
    color_dicts = {}
    color_dicts["normal"] = color_dict
    color_dicts["recruitment"] = color_dict_recruitment
    color_dicts["fire_freq"] = color_dict_fire_freq
    color_dicts["blackwhite"] = color_dict_blackwhite
    
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
        img = vis.get_image_from_grid(dynamics.state.grid, collect_states, color_dict)
   
    # Export image file
    imagepath = os.path.join(DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
    vis.save_image(img, imagepath)

    return dynamics, color_dicts


def termination_condition_satisfied(dynamics, start_time, user_args):
    satisfied = False
    condition = ""
    if (time.time() - start_time > user_args["timelimit"]):
        condition = f"Exceeded timelimit ({user_args['timelimit']} seconds)."
    if (dynamics.state.population.size() == 0):
        condition = "Population collapse."
    if (dynamics.state.grid.tree_cover >= 0.9):
        condition = "Tree cover converged to >90%."
    if (dynamics.state.grid.tree_cover <= 0.05):
        condition = "Tree cover converged to <5%."
    if (dynamics.time >= user_args["max_timesteps"]):
        condition = f"Maximum number of timesteps ({user_args['max_timesteps']}) reached."
    #print("Pop size: ", dynamics.state.population.size())
    satisfied = len(condition) > 0
    if satisfied:
        print("Simulation terminated. Cause:", condition)

    return satisfied


def do_visualizations(dynamics, fire_freq_arrays, fire_no_timesteps, verbose, color_dicts, collect_states, visualization_types, user_args):
    if ("recruitment" in visualization_types):
        print("Saving recruitment img...") if verbose else None
        recruitment_img = vis.get_image_from_grid(dynamics.state.grid, 0, color_dicts["recruitment"])
        imagepath_recruitment = os.path.join(DATA_OUT_DIR, "image_timeseries/recruitment/" + str(dynamics.time) + ".png")
        vis.save_image(recruitment_img, imagepath_recruitment, get_max(1000, recruitment_img.shape[0]))
    
    if ("fire_freq" in visualization_types):
        fire_freq_arrays.append(dynamics.state.grid.get_distribution(0) == -5)
        if dynamics.time > fire_no_timesteps:
            print("Saving fire frequency img...") if verbose else None
            fire_freq_img = vis.get_fire_freq_image(fire_freq_arrays[-fire_no_timesteps:], color_dicts["fire_freq"], dynamics.state.grid.width, fire_no_timesteps)
            imagepath_fire_freq = os.path.join(DATA_OUT_DIR, "image_timeseries/fire_frequencies/" + str(dynamics.time) + ".png")
            vis.save_image(fire_freq_img, imagepath_fire_freq, get_max(1000, fire_freq_img.shape[0]))

    print("-- Visualizing image...") if verbose else None
    if user_args["headless"]:
        # Get a color image representation of the initial state
        img = vis.get_image_from_grid(dynamics.state.grid, collect_states, color_dicts["normal"])
    else:
        # Get a color image representation of the initial state and show it.
        img = vis.visualize(
            dynamics.state.grid, user_args["image_width"], collect_states=collect_states,
            color_dict=color_dicts["normal"]
        )

    print("-- Saving image...") if verbose else None
    imagepath = os.path.join(DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
    vis.save_image(img, imagepath, get_max(1000, img.shape[0]))
    
    if ("fuel" in visualization_types):
        print("-- Saving fuel image...") if verbose else None
        fuel_img = vis.get_image_from_grid(dynamics.state.grid, 2, color_dicts["blackwhite"])
        imagepath_fuel = os.path.join(DATA_OUT_DIR, "image_timeseries/fuel/" + str(dynamics.time) + ".png")
        vis.save_image(fuel_img, imagepath_fuel, get_max(1000, fuel_img.shape[0]))

    if ("tree_LAI" in visualization_types):
        print("-- Saving tree LAI image...") if verbose else None
        tree_LAI_img = vis.get_image_from_grid(dynamics.state.grid, 1, color_dicts["blackwhite"], invert=True)
        imagepath_fuel = os.path.join(DATA_OUT_DIR, "image_timeseries/tree_LAI/" + str(dynamics.time) + ".png")
        vis.save_image(tree_LAI_img, imagepath_fuel, get_max(1000, tree_LAI_img.shape[0]))


def updateloop(dynamics, color_dicts, **user_args):
    start = time.time()
    print("Beginning simulation...")
    csv_path = user_args["csv_path"]
    visualization_types = [] # Options: "fire_freq", "recruitment", "fuel", "tree_LAI"
    #visualization_types = ["fire_freq", "recruitment", "fuel", "tree_LAI"] # Options: "fire_freq", "recruitment", "fuel", "tree_LAI"
    init_csv = True
    export_animal_resources = True
    collect_states = 1
    fire_no_timesteps = 1
    verbose = user_args["verbosity"]
    fire_freq_arrays = []
    if not user_args["headless"]:
        graphs = vis.Graphs(dynamics)
    while not termination_condition_satisfied(dynamics, start, user_args):
        print("-- Starting update (calling from python)") if verbose else None
        dynamics.update()
        print("-- Finished update") if verbose else None
        
        # Do visualizations
        do_visualizations(dynamics, fire_freq_arrays, fire_no_timesteps, verbose, color_dicts, collect_states, visualization_types, user_args)
        
        print("-- Exporting state data...") if verbose else None
        csv_path = io.export_state(dynamics, csv_path, init_csv)
        init_csv = False
        
        print("-- Saving tree positions...") if verbose else None
        io.save_tree_positions(dynamics)

        print("-- Showing graphs...") if verbose else None
        if not user_args["headless"]:
            graphs.update()
        
        break

        if export_animal_resources and (user_args["dispersal_mode"] == "all" or user_args["dispersal_mode"] == "animal"):
            # Get color image representations of the resource grid from the last iteration
            cover_path = os.path.join(DATA_OUT_DIR, "image_timeseries/cover/" + str(dynamics.time) + ".png")
            fruits_path = os.path.join(DATA_OUT_DIR, "image_timeseries/fruits/" + str(dynamics.time) + ".png")
            visits_path = os.path.join(DATA_OUT_DIR, "image_timeseries/visits/" + str(dynamics.time) + ".png")
            distance_path = os.path.join(DATA_OUT_DIR, "image_timeseries/distance/" + str(dynamics.time) + ".png")
            distance_coarse_path = os.path.join(DATA_OUT_DIR, "image_timeseries/distance_coarse/" + str(dynamics.time) + ".png")
            k_coarse_path = os.path.join(DATA_OUT_DIR, "image_timeseries/k_coarse/" + str(dynamics.time) + ".png")
            k_path = os.path.join(DATA_OUT_DIR, "image_timeseries/k/" + str(dynamics.time) + ".png")
            vis.save_resource_grid_colors(dynamics, "Turdus pilaris", "cover", cover_path, user_args["resource_grid_relative_size"])
            vis.save_resource_grid_colors(dynamics, "Turdus pilaris", "fruits", fruits_path, user_args["resource_grid_relative_size"])
            vis.save_resource_grid_colors(dynamics, "Turdus pilaris", "visits", visits_path, user_args["resource_grid_relative_size"])
            vis.save_resource_grid_colors(dynamics, "Turdus pilaris", "distance_single", distance_path, user_args["resource_grid_relative_size"])
            vis.save_resource_grid_colors(dynamics, "Turdus pilaris", "distance_single_coarse", distance_coarse_path, user_args["resource_grid_relative_size"])
            vis.save_resource_grid_colors(dynamics, "Turdus pilaris", "k", k_path, user_args["resource_grid_relative_size"])
            vis.save_resource_grid_colors(dynamics, "Turdus pilaris", "k_coarse", k_coarse_path, user_args["resource_grid_relative_size"])

    if not user_args["headless"]:
        cv2.destroyAllWindows()
        dynamics.free()
    
    return dynamics


def test_kernel():
    cpp.init_RNG()
    dist_max = 5000
    windspeed_gmean = 20
    windspeed_stdev = 3
    seed_terminal_speed = 0.65
    abscission_height = 30
    wind_kernel = cpp.Kernel(1, dist_max, windspeed_gmean, windspeed_stdev, 0, 3600, seed_terminal_speed, abscission_height)
    wind_kernel.build()
    vis.visualize_kernel(wind_kernel, "Wind kernel. d_max = {}, w_gmean = {}, \n w_stdev = {}, v_t = {}, h = {}".format(
        dist_max, windspeed_gmean, windspeed_stdev, seed_terminal_speed, abscission_height)
    )
    return


def test_discrete_probmodel():
    cpp.init_RNG()
    # Get a list of probabilities (normal distributed)
    x_interval = np.linspace(-3, 3, 100)
    probabilities = np.array(stats.norm.pdf(x_interval))
    probabilities = list(probabilities / np.sum(probabilities))
    #print("probs: ", probabilities)
    
    # Create the discrete probability model
    probmodel = cpp.DiscreteProbabilityModel(100)
    probmodel.set_probabilities(probabilities)
    
    # Plot samples from the model
    samples = [probmodel.sample() for i in range(100000)]
    plt.hist(samples, bins=50)
    plt.show()


def main(**user_args):
    #test_kernel()
    #test_discrete_probmodel()
    #return

    if user_args["test"] == "all":
        tests = init_tests(**user_args)
        tests.run_all()
    else:
        dynamics, color_dicts = init(**user_args)
        return updateloop(dynamics, color_dicts, **user_args)
 


    

