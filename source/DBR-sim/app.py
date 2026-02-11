from ast import arg
from dataclasses import dataclass
import sys
from tkinter import E, N, Y

import cv2
import numpy as np
import copy
import time
import os
from pathlib import Path
from matplotlib import pyplot as plt
import json
import scipy.stats as stats
import random
import math
from types import SimpleNamespace

import visualization
import file_handling as io
from helpers import *
from config import *
import disk_pattern_generator as _dpg
import sine_pattern_generator as spg
import simple_noise_generator as sng

#from x64.Debug import dbr_cpp as cpp
from x64.Release import dbr_cpp as cpp


def set_dispersal_kernel(
        dynamics, dispersal_mode, multi_disperser_params
    ):
    """Set the parameters for the dispersal kernel associated with the given dispersal_mode.

    The parameters are loaded from the given JSON file into a dictionary. This JSON file can contain parameters
    for multiple dispersal modes, though only the parameters of the requested dispersal mode will actually be used.
    
    Examples:
    For wind dispersal, the parameters may be scalars defining wind speed and standard deviation, whereas
    for animal dispersal, it will be separate dictionaries per species, containing parameter values for flight speed 
    and such.

    Input:
        dynamics: The C++ object that contains the core logic of SDSF.
        dispersal_mode (string): the dispersal mode used, such as animal or wind.
        multi_disperser_params (str): Name of the JSON-file containing the parameters of the dispersal kernel. 
    Return:
        dynamics: The C++ object that contains the core logic of SDSF, now including or more initialized dispersal kernels.
        animal_species (list): List of animal species for which the dispersal kernel has been initialized. If you did not
            specify 'animal' or 'all' as your dispersal_mode, then this will be an empty list.
    """

    with open(os.path.join(cfg.DATA_IN_DIR, multi_disperser_params), "r") as mdp_jsonfile:
        multi_disperser_params = json.load(mdp_jsonfile)

    animal_species = []
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
        animal_dispersal_params = multi_disperser_params["animal"]
        animal_species = [species for species in animal_dispersal_params.keys() if species != "population"]
        dynamics.set_global_animal_kernel(multi_disperser_params["animal"])
    elif (dispersal_mode == "all"):
        animal_dispersal_params = multi_disperser_params["animal"]
        animal_species = [species for species in animal_dispersal_params.keys() if species != "population"]
        multi_disperser_params.pop("animal")
        dynamics.set_global_kernels(multi_disperser_params, animal_dispersal_params)

    return dynamics, animal_species


def set_initial_tree_cover(dynamics, cfg, color_dicts):
    print("Initial pattern image:", cfg.initial_pattern_image)
    if cfg.initial_pattern_image == "ctrl":
        _cfg = copy.copy(cfg)
        img, cfg.cover_img_path, benchmark_cover = cfg.vis.generate_disk_pattern(**vars(_cfg))
        img = cv2.resize(img, (cfg.grid_width, cfg.grid_width), interpolation=cv2.INTER_NEAREST)
        img = img / 255

        # If the user has set the override_image_treecover to 2, we will use the benchmark cover value from the controllable pattern image.
        # The benchmark cover is the cover for a version of the pattern produced using the given parameters, but with a sine amplitude set to 0 (i.e., with circular disks).
        if cfg.override_image_treecover == 2:
            cfg.override_image_treecover = benchmark_cover
    elif cfg.initial_pattern_image == "perlin_noise":
        cfg.cover_img_path = f"{cfg.PERLIN_NOISE_DIR}/" + cfg.initial_pattern_image + ".png"
        noise_frequency = 5.0 / cfg.patch_width # Convert patch width to noise frequency
        noise_frequency = round(noise_frequency, 2) # Conform noise frequency to 2 decimal places, to ensure periodicity of the noise pattern
        img = cfg.vis.generate_perlin_noise_image(cfg.cover_img_path, frequency=noise_frequency, octaves=cfg.noise_octaves)
    elif cfg.initial_pattern_image == "none":
        cfg.cover_img_path = f"{cfg.SIMPLE_PATTERNS_DIR}/uniform.png" # Pure white image
        img = cv2.imread(cfg.cover_img_path, cv2.IMREAD_GRAYSCALE)
    else:
        print("Setting tree cover from image...")

        cfg.cover_img_path = f"{cfg.DATA_IN_DIR}/state_patterns/" + cfg.initial_pattern_image
        img = cv2.imread(cfg.cover_img_path, cv2.IMREAD_GRAYSCALE)
        print(f"Path of the generated {cfg.initial_pattern_image} pattern image: ", cfg.cover_img_path) if cfg.initial_pattern_image in ["ctrl", "perlin_noise"] else None
        
    if ("perlin_noise" in cfg.cover_img_path) and (not "thresholded.png" in cfg.cover_img_path):
        img = cfg.vis.get_thresholded_image(img, cfg.treecover * img.shape[0] * img.shape[0] * 255 )
        cv2.imwrite(cfg.cover_img_path.replace(".png", "_thresholded.png"), img)
 
    img = cv2.resize(img, (dynamics.state.grid.width, dynamics.state.grid.width), interpolation=cv2.INTER_NEAREST)

    # Do burn-in procedure to stabilize initial forest density and demography    
    dynamics, cfg = do_burn_in(dynamics, cfg, img, color_dicts, target_treecover=cfg.treecover)

    return dynamics, cfg


def init(cfg):

    # Make any changes to cfg if needed
    cfg.patch_width = round(cfg.patch_width, 2)

    # Set random seed
    if (cfg.random_seed == -999):
        cfg.random_seed = random.randint(0, 100000000)
        cfg.rng = np.random.default_rng(cfg.random_seed)
        print(f"Generated new global random seed ({cfg.random_seed})")
    else:
        cfg.rng = np.random.default_rng(cfg.random_seed)
        print(f"Using given global random seed ({cfg.random_seed})")

    # Set random seed for fire frequency probability distribution. If -999 is given, a random seed will be generated. Otherwise, the given seed will be used.
    if cfg.firefreq_random_seed == -999:
        cfg.firefreq_random_seed = random.randint(0, 1000000)
        print(f"Generated new random seed ({cfg.firefreq_random_seed}) for fire frequency probability distribution: ", )
    else:
        print(f"Using given random seed ({cfg.firefreq_random_seed}) for fire frequency probability distribution.")

    # Initialize visualization tool
    cfg.vis = visualization.Visualiser(cfg)
    
    # Initialize disk pattern generator
    cfg.dpg = _dpg.DiskPatternGenerator(cfg)

    # Initialize dynamics object and state
    dynamics = cpp.create_dynamics(vars(cfg))
    dynamics.init_state(
        cfg.grid_width,
        cfg.dbh_q1,
        cfg.dbh_q2,
        cfg.growth_rate_multiplier_params[0],
        cfg.growth_rate_multiplier_params[1],
        cfg.growth_rate_multiplier_params[2],
        cfg.minimum_patch_size,
        cfg.LAI_aggregation_radius
    )
    
    # Set dispersal kernel
    dynamics, animal_species = set_dispersal_kernel(dynamics, cfg.dispersal_mode, cfg.multi_disperser_params)
    
    # Create color dictionaries for visualizations
    no_colors = 100
    if cfg.display_fire_effects == 1:
        color_dict = cfg.vis.get_color_dict(no_colors, begin=0.2, end=0.5, distr_type="normal_with_fire_effects")
    else:
        color_dict = cfg.vis.get_color_dict(no_colors, begin=0.2, end=0.5, distr_type="normal")
    color_dict_recruitment = cfg.vis.get_color_dict(no_colors, begin=0.2, end=0.5, distr_type="recruitment")
    color_dict_fire_freq = cfg.vis.get_color_dict(10, begin=0.2, end=0.5, distr_type="fire_freq")
    color_dict_blackwhite = cfg.vis.get_color_dict(no_colors, distr_type="blackwhite")
    color_dict_colored_patches = cfg.vis.get_color_dict(no_colors, distr_type="colored_patches")
    color_dicts = {}
    color_dicts["normal"] = color_dict
    color_dicts["recruitment"] = color_dict_recruitment
    color_dicts["fire_freq"] = color_dict_fire_freq
    color_dicts["blackwhite"] = color_dict_blackwhite
    color_dicts["colored_patches"] = color_dict_colored_patches
    
    # Set input maps
    print("Setting heterogeneity maps...") if cfg.verbosity > 0 else None
    dynamics, cfg = io.set_heterogeneity_maps(dynamics, cfg)
    
    # Set initial tree cover
    print("Setting initial tree cover...") if cfg.verbosity > 0 else None
    dynamics, cfg = set_initial_tree_cover(dynamics, cfg, color_dicts)
    
    # Visualize the initial state
    collect_states = True
    print("Visualizing state at t = 0")
    if not cfg.headless:
        # Get a color image representation of the initial state and show it.
        img = cfg.vis.visualize(
            dynamics.state.grid, cfg.image_width, collect_states=collect_states,
            color_dict=color_dict
        )
    else:
        # Get a color image representation of the initial state
        img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dict, collect_states=collect_states)
   
    # Export image file
    imagepath = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
    cfg.vis.save_image(img, imagepath)
    
    if cfg.dispersal_mode == "all" or cfg.dispersal_mode == "animal":

        # Export resource grid lookup table
        print("animal species: ", animal_species)
        for species in animal_species:
            lookup_table, fpath = io.get_lookup_table(species, cfg.grid_width, cfg.resource_grid_width * cfg.resource_grid_width)
            if lookup_table is None:
                print(f"Lookup table file {fpath} not found. Creating new one... (might take some time).")
                dynamics.precompute_resourcegrid_lookup_table(species)
                lookup_table = dynamics.get_resource_grid_lookup_table(species)
                io.export_lookup_table(lookup_table, cfg.grid_width, species)
            else:
                print(f"Lookup table file {fpath} found. Loading...")
                dynamics.set_resource_grid_lookup_table(lookup_table, species)

    return dynamics, color_dicts, cfg


def termination_condition_satisfied(dynamics, start_time, cfg):
    satisfied = False
    condition = ""
    if (dynamics.time >= int(cfg.max_timesteps)):
        condition = f"Maximum number of timesteps ({cfg.max_timesteps}) reached."
    if (cfg.termination_conditions == "all" and dynamics.state.grid.get_tree_cover() > 0.98):
        condition = f"Tree cover exceeds 98%."
    satisfied = len(condition) > 0
    if satisfied:
        print("\nSimulation terminated. Cause:", condition)

    return satisfied


def do_burn_in(dynamics, cfg, forest_mask, color_dicts, target_treecover=1):
    """Perform burn-in procedure to stabilize initial forest density and demography."""
    
    print(f"Beginning burn-in for {cfg.burnin_duration} timesteps...")
    init_csv = True
    cfg.treesize_bins = "initialize"
    visualization_types = ["aggr_tree_LAI", "colored_patches"]
    fire_freq_arrays = []
    fire_no_timesteps = 1
    patch_colors = {}
    while dynamics.time < cfg.burnin_duration or dynamics.state.grid.get_tree_cover() > target_treecover:
        print(f"Burn-in timestep:    {dynamics.time})")
        dynamics.disperse_within_forest(forest_mask)
        dynamics.grow()
        dynamics.induce_background_mortality()
        dynamics.state.repopulate_grid(cfg.verbosity)
        dynamics.prune(forest_mask)
        dynamics.state.repopulate_grid(cfg.verbosity)
        dynamics.state.grid.reset_state_distr()
        dynamics.update_forest_patch_detection()
        dynamics.report_state()
        cfg = io.export_state(dynamics, path=cfg.csv_path, init_csv=init_csv, cfg=cfg)
        init_csv = False

        # Feedback control on tree cover during burn-in
        if dynamics.state.grid.get_tree_cover() > target_treecover:
            prune_fraction = dynamics.state.grid.get_tree_cover() - target_treecover
            dynamics.prune_randomly(prune_fraction)
        
        # Obtain patches from simulation
        patches = dynamics.get_patches() 

        # Do visualizations
        if not cfg.headless:
            img = cfg.vis.visualize(
                dynamics.state.grid, cfg.image_width, collect_states=1,
                color_dict=color_dicts["normal"]
            )
        dynamics.time += 1

    # Remove seedlings, then let trees disperse for one time step to create a realistic initial distribution of seedlings.
    dynamics.time -= 1
    dynamics.remove_trees_up_to_age(1)
    dynamics.state.repopulate_grid(0)

    dynamics.time = 0

    return dynamics, cfg


def do_iteration(dynamics, cfg):
    """Run a single simulation iteration."""

    # -- PREPARE -- #
    print(f"\nTime: {dynamics.time}")

    if cfg.verbosity > 0:
        print("Resetting state distr...")

    dynamics.state.grid.reset_state_distr()

    # -- DISPERSE -- #
    if cfg.disperse:
        t0 = time.time()
        if cfg.verbosity > 0:
            print("Beginning dispersal...")

        if dynamics.time > 0:
            dynamics.disperse()

        t1 = time.time()
        if cfg.verbosity > 0:
            print(f"Dispersal took {t1 - t0:.6f} seconds. Beginning burn...")

    # -- BURN -- #
    if cfg.burn:
        t0 = time.time()
        dynamics.burn()
        t1 = time.time()

        if cfg.verbosity > 0:
            print(f"Burns took {t1 - t0:.6f} seconds. Beginning growth...")

    # -- GROW -- #
    if cfg.grow:
        t0 = time.time()
        dynamics.grow()
        t1 = time.time()

        if cfg.verbosity > 0:
            print(f"Growth took {t1 - t0:.6f} seconds.")

    # -- BACKGROUND MORTALITY -- #
    dynamics.induce_background_mortality()

    if cfg.verbosity > 0:
        print("Induced background mortality. Repopulating grid...")

    # Post-simulation cleanup and reporting
    dynamics.state.repopulate_grid(cfg.verbosity)

    if cfg.verbosity > 1:
        print("Redoing grid count...")

    dynamics.state.grid.redo_count()
    dynamics.report_state()
    dynamics.update_firefree_interval_averages()
    dynamics.update_forest_patch_detection()

    return dynamics


def updateloop(dynamics, color_dicts, cfg):
    print("Beginning simulation...")
    start = time.time()
    treesize_bins = cfg.treesize_bins
    visualization_types = ["fire_freq", "recruitment", "fuel", "tree_LAI", "aggr_tree_LAI", "colored_patches", "fuel_penetration"]
    init_csv = False
    prev_tree_cover = [cfg.treecover] * 60
    slope = 0
    largest_absolute_slope = 0
    export_animal_resources = False
    patches = {}
    collect_states = 1
    fire_no_timesteps = 1
    patch_colors = {}
    verbose = cfg.verbosity
    fire_freq_arrays = []
    color_dict_fire_freq = cfg.vis.get_color_dict(fire_no_timesteps, begin=0.2, end=0.5, distr_type="fire_freq")
    color_dicts["fire_freq"] = color_dict_fire_freq

    if not cfg.headless:
        graphs = visualization.Graphs(dynamics)

    while True:
        print("-- Starting iteration...") if verbose else None
        dynamics = do_iteration(dynamics, cfg)
        print("-- Finished update") if verbose else None
    
        # Obtain initial number of recruits
        if dynamics.time == 0:
            initial_no_dispersals = dynamics.get_initial_no_dispersals()
    
        # Track tree cover trajectory
        prev_tree_cover.append(dynamics.state.grid.get_tree_cover())
        if dynamics.time > 10:
            slope = (prev_tree_cover[-1] - cfg.treecover) / dynamics.time
            single_tstep_slope = prev_tree_cover[-1] - prev_tree_cover[-2]
            if abs(single_tstep_slope) > largest_absolute_slope:
                largest_absolute_slope = abs(single_tstep_slope)
        else:
            slope = 0

        # Track treecover trajectory
        do_terminate = termination_condition_satisfied(dynamics, start, cfg)
    
        # Obtain patches from simulation
        old_no_patches = len(patches)
        patches = dynamics.get_patches() 
        patch_count_change =  len(patches) - old_no_patches
    
        # Do visualizations
        if cfg.export_visualizations:
            visualization.do_visualizations(
                dynamics, fire_freq_arrays, fire_no_timesteps, verbose,
                color_dicts, collect_states, visualization_types,
                patches, patch_count_change, patch_colors, cfg
            )
    
        print("-- Exporting state_data...") if verbose else None
        cfg = io.export_state(
            dynamics, path=cfg.csv_path, init_csv=init_csv,
            tree_cover_slope=slope, cfg=cfg
        )

        print("-- Saving tree positions...") if verbose else None
        # io.save_state(dynamics)
    
        if cfg.report_state == "True" or cfg.report_state is True:
            io.update_state_report(dynamics)

        print("-- Showing graphs...") if verbose else None
        if not cfg.headless:
            graphs.update()

        if export_animal_resources and (cfg.dispersal_mode == "all" or cfg.dispersal_mode == "animal"):
            # Get color image representations of the resource grid from the last iteration
            cover_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/cover/" + str(dynamics.time) + ".png")
            fruits_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fruits/" + str(dynamics.time) + ".png")
            visits_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/visits/" + str(dynamics.time) + ".png")
            distance_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/distance/" + str(dynamics.time) + ".png")
            distance_coarse_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/distance_coarse/" + str(dynamics.time) + ".png")
            k_coarse_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/k_coarse/" + str(dynamics.time) + ".png")
            k_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/k/" + str(dynamics.time) + ".png")

            cfg.vis.save_resource_grid_colors(dynamics, "Turdus merula", "cover", cover_path)
            cfg.vis.save_resource_grid_colors(dynamics, "Turdus merula", "fruits", fruits_path)
            cfg.vis.save_resource_grid_colors(dynamics, "Turdus merula", "visits", visits_path)
            cfg.vis.save_resource_grid_colors(dynamics, "Turdus merula", "distance_single", distance_path)
            cfg.vis.save_resource_grid_colors(dynamics, "Turdus merula", "distance_single_coarse", distance_coarse_path)
            cfg.vis.save_resource_grid_colors(dynamics, "Turdus merula", "k", k_path)
            cfg.vis.save_resource_grid_colors(dynamics, "Turdus merula", "k_coarse", k_coarse_path)

        if do_terminate:
            break
        
        dynamics.time += 1

    if not cfg.headless:
        cv2.destroyAllWindows()
        dynamics.free()

    if not hasattr(cfg, "dependent_vars"):
        cfg.dependent_vars = {
            "tree_cover_slope": slope, 
            "largest_absolute_slope": largest_absolute_slope, 
            "initial_no_dispersals": initial_no_dispersals
        }

    return dynamics, cfg


def strategy_distribution_params_are_loaded(strategy_distribution_params):
    return type(strategy_distribution_params) == dict


def main(**user_args): 
    global cfg

    if user_args["verbosity"] >= 0:
        print(f"Loading dbr_cpp from {cfg.BUILD_DIR} and preparing model run... (this may take some time)")

    # Set any given batch parameters
    if not strategy_distribution_params_are_loaded(user_args["strategy_distribution_params"]):
        with open(os.path.join(cfg.DATA_IN_DIR, user_args["strategy_distribution_params"]), "r") as sdp_jsonfile:
            user_args["strategy_distribution_params"] = json.load(sdp_jsonfile)

    args = SimpleNamespace(**user_args)
    cfg = apply_user_args_to_configuration(args, cfg)

    dynamics, color_dicts, cfg = init(cfg)
    return updateloop(dynamics, color_dicts, cfg)
 


    

