from dataclasses import dataclass
import sys

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


def unpack_control_keys(control_variable):
    control_keys = control_variable.split("->")
    return control_keys


def set_dispersal_kernel(
        dynamics, dispersal_mode, multi_disperser_params
    ):
    animal_species = []
    with open(os.path.join(DATA_IN_DIR, multi_disperser_params), "r") as mdp_jsonfile:
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
        animal_species = [species for species in animal_dispersal_params.keys() if species != "population"]
        multi_disperser_params.pop("animal")
        dynamics.set_global_kernels(multi_disperser_params, animal_dispersal_params)

    return dynamics, animal_species


def init_tests(
    timestep=None, grid_width=None, cellsize=None, max_dbh=None, image_width=None,
    treecover=None, self_ignition_factor=None, flammability=None,
    rainfall=None, unsuppressed_flammability=None, 
    verbosity=None, dbh_q1=None, dbh_q2=None, seed_bearing_threshold=None,
    dispersal_mode=None, linear_diffusion_q1=None, linear_diffusion_q2=None,
    dispersal_min=None, dispersal_max=None, growth_rate_multiplier=None, 
    flammability_coefficients_and_constants=None, saturation_threshold=None, fire_resistance_params=None,
    constant_mortality=None, headless=False, wind_dispersal_params=None, animal_dispersal_params=None,
    multi_disperser_params=None, strategy_distribution_params=None, resource_grid_width=None,
    initial_pattern_image=None, mutation_rate=None, growth_rate_multiplier_params=None, **user_args
):
    # Obtain strategy distribution parameters
    with open(os.path.join(DATA_IN_DIR, strategy_distribution_params), "r") as sdp_jsonfile:
        strategy_distribution_params = json.load(sdp_jsonfile)    

    tests = cpp.Tests(timestep, cellsize, self_ignition_factor, rainfall, seed_bearing_threshold,
        growth_rate_multiplier, unsuppressed_flammability, flammability_coefficients_and_constants[0],
        flammability_coefficients_and_constants[1], flammability_coefficients_and_constants[2], 
        flammability_coefficients_and_constants[3], max_dbh, saturation_threshold, fire_resistance_params[0],
        fire_resistance_params[1], fire_resistance_params[2], constant_mortality, strategy_distribution_params, 
        resource_grid_width, mutation_rate, verbosity, grid_width, dbh_q1, dbh_q2, growth_rate_multiplier_params[0], growth_rate_multiplier_params[1])
    
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
    multi_disperser_params=None, strategy_distribution_params=None, resource_grid_width=None,
    initial_pattern_image=None, mutation_rate=None, STR=None,
    ):
    
    # Obtain strategy distribution parameters
    with open(os.path.join(DATA_IN_DIR, strategy_distribution_params), "r") as sdp_jsonfile:
        strategy_distribution_params = json.load(sdp_jsonfile)
    if batch_parameters and "strategies->" in batch_parameters["control_variable"]:
        control_keys = unpack_control_keys(batch_parameters["control_variable"])
        control_keys.remove("strategies")
        strategy_distribution_params[control_keys[0]][control_keys[1]] = batch_parameters["control_value"]

    # Initialize dynamics object and state
    print("Verbosity type: ", type(verbosity))
    dynamics = cpp.Dynamics(
        timestep, cellsize, self_ignition_factor, rainfall, seed_bearing_threshold,
        growth_rate_multiplier, unsuppressed_flammability, flammability_coefficients_and_constants[0],
        flammability_coefficients_and_constants[1], flammability_coefficients_and_constants[2], 
        flammability_coefficients_and_constants[3], max_dbh, saturation_threshold, fire_resistance_params[0],
        fire_resistance_params[1], fire_resistance_params[2], constant_mortality, strategy_distribution_params, 
        resource_grid_width, mutation_rate, STR, verbosity
    )
    dynamics.init_state(grid_width, dbh_q1, dbh_q2, growth_rate_multiplier_params[0], growth_rate_multiplier_params[1])
    
    # Set dispersal kernel
    dynamics, animal_species = set_dispersal_kernel(dynamics, dispersal_mode, multi_disperser_params)
    print("species: ", animal_species)
    
    # Set initial tree cover
    if initial_pattern_image == "none":
        dynamics.state.set_tree_cover(treecover)
    else:
        path = f"{DATA_IN_DIR}/state patterns/" + initial_pattern_image
        if "perlin_noise" == initial_pattern_image:
            print("Setting tree cover using perlin noise function...")
            path = f"{PERLIN_NOISE_DIR}/" + initial_pattern_image + ".png"
            vis.generate_perlin_noise_image(path, frequency=user_args["noise_frequency"], octaves=user_args["noise_octaves"])
        else:
            print("Setting tree cover from image...")
        img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        img = vis.get_thresholded_image(img, treecover * img.shape[0] * img.shape[0] * 255 )
        cv2.imwrite(path.replace(".png", "_thresholded.png"), img)
        img = cv2.resize(img, (dynamics.state.grid.width, dynamics.state.grid.width), interpolation=cv2.INTER_LINEAR) 
        dynamics.state.set_cover_from_image(img / 255, -1)
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
    
    if dispersal_mode == "all" or dispersal_mode == "animal":

        # Export resource grid lookup table
        for species in animal_species:
            lookup_table, fpath = io.get_lookup_table(species, resource_grid_width * resource_grid_width)
            if lookup_table is None:
                print(f"Lookup table file {fpath} not found. Creating new one...")
                dynamics.precompute_resourcegrid_lookup_table(species)
                lookup_table = dynamics.get_resource_grid_lookup_table(species)
                io.export_lookup_table(lookup_table, species)
            else:
                print(f"Lookup table file {fpath} found. Loading...")
                dynamics.set_resource_grid_lookup_table(lookup_table, species)

    return dynamics, color_dicts


def termination_condition_satisfied(dynamics, start_time, user_args):
    satisfied = False
    condition = ""
    if (dynamics.time >= int(user_args["max_timesteps"])):
        condition = f"Maximum number of timesteps ({user_args['max_timesteps']}) reached."
    if (user_args["termination_conditions"] == "all" and dynamics.state.grid.get_tree_cover() > 0.98):
        condition = f"Tree cover exceeds 98%."
    satisfied = len(condition) > 0
    if satisfied:
        print("\nSimulation terminated. Cause:", condition)

    return satisfied


def do_visualizations(dynamics, fire_freq_arrays, fire_no_timesteps, verbose, color_dicts, collect_states, visualization_types, user_args):
    if ("recruitment" in visualization_types):
        print("Saving recruitment img...") if verbose else None
        recruitment_img = vis.get_image_from_grid(dynamics.state.grid, 0, color_dicts["recruitment"])
        imagepath_recruitment = os.path.join(DATA_OUT_DIR, "image_timeseries/recruitment/" + str(dynamics.time) + ".png")
        vis.save_image(recruitment_img, imagepath_recruitment, get_max(1000, recruitment_img.shape[0]), interpolation="none")
    
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
    visualization_types = ["fire_freq"] # Options: "fire_freq", "recruitment", "fuel", "tree_LAI"
    #visualization_types = ["fire_freq", "recruitment", "fuel", "tree_LAI"] # Options: "fire_freq", "recruitment", "fuel", "tree_LAI"
    #visualization_types = ["recruitment"] # Options: "fire_freq", "recruitment", "fuel", "tree_LAI"
    init_csv = True
    prev_tree_cover = [user_args["treecover"]] * 60
    slope = 0
    export_animal_resources = True
    collect_states = 1
    fire_no_timesteps = 200
    verbose = user_args["verbosity"]
    fire_freq_arrays = []
    color_dict_fire_freq = vis.get_color_dict(fire_no_timesteps, begin=0.2, end=0.5, distr_type="fire_freq")
    color_dicts["fire_freq"] = color_dict_fire_freq
    if not user_args["headless"]:
        graphs = vis.Graphs(dynamics)
    while True:
        print("-- Starting update (calling from python)") if verbose else None
        dynamics.update()
        print("-- Finished update") if verbose else None
        
        # Track tree cover trajectory and evaluate termination conditions
        prev_tree_cover.append(dynamics.state.grid.get_tree_cover())
        if dynamics.time > 10:
            slope = (prev_tree_cover[-1] - prev_tree_cover[-10]) / 10
        else:
            slope = 0
        do_terminate = termination_condition_satisfied(dynamics, start, user_args)
        
        # Do visualizations
        do_visualizations(dynamics, fire_freq_arrays, fire_no_timesteps, verbose, color_dicts, collect_states, visualization_types, user_args)
        
        print("-- Exporting state data...") if verbose else None
        csv_path = io.export_state(dynamics, csv_path, init_csv, tree_cover_slope=slope)
        init_csv = False
        
        print("-- Saving tree positions...") if verbose else None
        io.save_tree_positions(dynamics)

        print("-- Showing graphs...") if verbose else None
        if not user_args["headless"]:
            graphs.update()

        if export_animal_resources and (user_args["dispersal_mode"] == "all" or user_args["dispersal_mode"] == "animal"):
            # Get color image representations of the resource grid from the last iteration
            cover_path = os.path.join(DATA_OUT_DIR, "image_timeseries/cover/" + str(dynamics.time) + ".png")
            fruits_path = os.path.join(DATA_OUT_DIR, "image_timeseries/fruits/" + str(dynamics.time) + ".png")
            visits_path = os.path.join(DATA_OUT_DIR, "image_timeseries/visits/" + str(dynamics.time) + ".png")
            distance_path = os.path.join(DATA_OUT_DIR, "image_timeseries/distance/" + str(dynamics.time) + ".png")
            distance_coarse_path = os.path.join(DATA_OUT_DIR, "image_timeseries/distance_coarse/" + str(dynamics.time) + ".png")
            k_coarse_path = os.path.join(DATA_OUT_DIR, "image_timeseries/k_coarse/" + str(dynamics.time) + ".png")
            k_path = os.path.join(DATA_OUT_DIR, "image_timeseries/k/" + str(dynamics.time) + ".png")
            vis.save_resource_grid_colors(dynamics, "Turdus merula", "cover", cover_path)
            vis.save_resource_grid_colors(dynamics, "Turdus merula", "fruits", fruits_path)
            vis.save_resource_grid_colors(dynamics, "Turdus merula", "visits", visits_path)
            vis.save_resource_grid_colors(dynamics, "Turdus merula", "distance_single", distance_path)
            vis.save_resource_grid_colors(dynamics, "Turdus merula", "distance_single_coarse", distance_coarse_path)
            vis.save_resource_grid_colors(dynamics, "Turdus merula", "k", k_path)
            vis.save_resource_grid_colors(dynamics, "Turdus merula", "k_coarse", k_coarse_path)

        if do_terminate:
            break

    if not user_args["headless"]:
        cv2.destroyAllWindows()
        dynamics.free()
    
    return dynamics, slope


def test_kernel():
    cpp.init_RNG()
    dist_max = 200
    windspeed_gmean = 10
    windspeed_stdev = 5
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

    #vis.visualize_legend()
    #return

    #test_kernel()
    #test_discrete_probmodel()
    #return

    assert (user_args["grid_width"] % user_args["resource_grid_width"] == 0), "Resource grid width must be a divisor of grid width."

    # Set number of cells
    #user_args["grid_width"] = int( user_args["grid_width"] / user_args["cellsize"])

    if user_args["test"] == "all":
        tests = init_tests(**user_args)
        tests.run_all()
    else:
        dynamics, color_dicts = init(**user_args)
        return updateloop(dynamics, color_dicts, **user_args)
 


    

