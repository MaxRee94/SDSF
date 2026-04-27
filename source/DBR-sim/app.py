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
import helpers as h
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


def create_controlled_patches_image(cfg, dynamics):
    """Create an image with patches whose perimeter-area, size, inter-patch distance
    is determined by params in cfg."""

    print("Generating controlled pattern image...")

    _cfg = copy.copy(cfg)
    img, cfg.cover_img_path, benchmark_cover = cfg.vis.generate_disk_pattern(**vars(_cfg))
    img = cv2.resize(img, (cfg.grid_width, cfg.grid_width), interpolation=cv2.INTER_NEAREST)

    # If the user has set the override_image_treecover to 2, we will use the benchmark cover value from
    # the controllable pattern image. The benchmark cover is the tree cover you would get if you would
    # use circular patches.
    if cfg.override_image_treecover == 2:
        cfg.override_image_treecover = benchmark_cover

    return img, dynamics, cfg


def create_perlin_noise_image(cfg, dynamics):
    """Create perlin noise image according to the parameters set in the cfg."""

    print("Generating perlin noise image...")

    cfg.cover_img_path = f"{cfg.PERLIN_NOISE_DIR}/" + cfg.initial_pattern_image + ".png"
    noise_frequency = 5.0 / round(cfg.patch_width, 2) # Convert patch width to noise frequency
    noise_frequency = round(noise_frequency, 2) # Conform noise frequency to 2 decimal places, to ensure periodicity of the noise pattern
    img = cfg.vis.generate_perlin_noise_image(cfg.cover_img_path, frequency=noise_frequency, octaves=cfg.noise_octaves)

    # Ensure the image is thresholded
    img = cfg.vis.get_thresholded_image(img * 255, cfg.treecover * img.shape[0] * img.shape[0] * 255 )
    cv2.imwrite(cfg.cover_img_path.replace(".png", "_thresholded.png"), img)

    return img, dynamics, cfg


def load_homogeneous_image(cfg, dynamics):
    """Load a pure white image so that no specific tree cover pattern is applied at the start of the simulation (i.e.,
    trees are placed randomly)."""
    
    print("Loading pure white image...")

    cfg.cover_img_path = f"{cfg.SIMPLE_PATTERNS_DIR}/uniform.png" # Pure white image
    img = cv2.imread(cfg.cover_img_path, cv2.IMREAD_GRAYSCALE)

    return img, dynamics, cfg


def load_specified_existing_image(cfg, dynamics):
    """Load the image whose name was specified in cfg.initial_pattern_image."""

    cfg.cover_img_path = f"{cfg.DATA_IN_DIR}/state_patterns/" + cfg.initial_pattern_image
    img = cv2.imread(cfg.cover_img_path, cv2.IMREAD_GRAYSCALE)
    print(f"Path of the generated {cfg.initial_pattern_image} pattern image: ", cfg.cover_img_path) if cfg.initial_pattern_image in ["ctrl", "perlin_noise"] else None

    return img, dynamics, cfg


def set_initial_tree_cover(dynamics, cfg, color_dicts):
    """Generate- or load a black-and-white image and use it to set an initial pattern of forest cover."""
    
    if cfg.initial_pattern_image == "ctrl":
        img, dynamics, cfg = create_controlled_patches_image(cfg, dynamics)
    elif cfg.initial_pattern_image == "perlin_noise":
        img, dynamics, cfg = create_perlin_noise_image(cfg, dynamics)
    elif cfg.initial_pattern_image == "none":
        img, dynamics, cfg = load_homogeneous_image(cfg, dynamics)
    else:
        img, dynamics, cfg = load_specified_existing_image(cfg, dynamics)
 
    # Resize the image to the dimensions of the grid
    img = cv2.resize(img, (dynamics.state.grid.width, dynamics.state.grid.width), interpolation=cv2.INTER_NEAREST)

    # Do burn-in procedure to allow forest demography to organically reach an equilibrium.   
    dynamics, cfg = do_burn_in(dynamics, cfg, img, color_dicts, target_treecover=cfg.treecover)

    return dynamics


def set_random_seeds(cfg):
    """Set random seeds for reproducibility. If -999 is given as the seed, a random seed will be generated and used.
    
    Args:
        cfg (namespace): 
            The configuration object containing the random seed settings. cfg.random_seed is the global
            random seed, which is used for all stochastic processes in the model except for the fire frequency
            probability distribution. cfg.firefreq_random_seed is the random seed specifically for the fire frequency
            probability distribution, which can be set separately to allow for independent reproducibility of temporal
            patterns of fire frequency across different runs.
    """

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

    return cfg


def generate_and_or_load_animal_dist_lookup_tables(dynamics, cfg, animal_species):
    # Export resource grid lookup table
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


def init(cfg):
    """Initialize the model.

    This includes setting random seeds, initializing the dynamics object, setting the parameters for the dispersal kernel,
    creating an initial pattern of tree cover, and visualizing the initial state of the model. The function returns the
    initialized dynamics object."""

    # Set random seeds
    cfg = set_random_seeds(cfg)

    # Initialize visualization tool
    cfg.vis = visualization.Visualiser(cfg)
    
    # Initialize disk pattern generator
    cfg.dpg = _dpg.DiskPatternGenerator(cfg)

    # Initialize dynamics object (containing process-related functionality) 
    # and state (containing state variables and data-management functionality).
    dynamics = cpp.create_dynamics(vars(cfg))
    dynamics.init_state(
        cfg.grid_width,
        cfg.growth_rate_multiplier_params[0],
        cfg.growth_rate_multiplier_params[1],
        cfg.growth_rate_multiplier_params[2],
        cfg.minimum_patch_size,
        cfg.local_neighborhood_radius
    )
    
    # Set the parameters for the dispersal kernel associated with the given dispersal_mode.
    dynamics, animal_species = set_dispersal_kernel(dynamics, cfg.dispersal_mode, cfg.multi_disperser_params)
    
    # Create color dictionaries to be used in visualizations
    visualization.create_color_dict(cfg)

    # Set input maps
    print("Setting heterogeneity maps...") if cfg.verbosity > 0 else None
    dynamics = io.set_heterogeneity_maps(dynamics, cfg)

    # Set initial tree cover
    print("Setting initial tree cover...") if cfg.verbosity > 0 else None
    dynamics = set_initial_tree_cover(dynamics, cfg, cfg.color_dicts)
    
    # Visualize the initial state
    visualization.visualize_initial_state(dynamics, cfg)
    
    if cfg.dispersal_mode in ["animal", "all"]:
        # Export resource grid lookup table
        generate_and_or_load_animal_dist_lookup_tables(dynamics, cfg, animal_species)

    # Initialize variables for tracking compute time, tree cover trajectory, and patch dynamics
    cfg.visualization_types = [
        "fire_freq", "recruitment", "fuel", "tree_LAI", "aggr_tree_LAI", "colored_patches",
        "fuel_penetration", "stand_density"
    ]
    cfg.computer_start_time = time.time()
    cfg.init_csv = False
    cfg.prev_tree_cover = [cfg.treecover] * 60
    cfg.treecover_slope = 0
    cfg.largest_absolute_slope = 0
    cfg.export_animal_resources = False
    cfg.patches = {}
    cfg.collect_states = 1
    cfg.fire_no_timesteps = 1
    cfg.patch_colors = {}
    cfg.fire_freq_arrays = []

    return dynamics


def termination_condition_satisfied(dynamics, start_time, cfg):
    """Check if any of the termination conditions specified in cfg have been satisfied.
    return 'True' if the simulation should be terminated, and 'False' otherwise."""

    satisfied = False
    satisifed_condition = ""
    if (dynamics.time >= int(cfg.max_timesteps)):
        satisifed_condition = f"Maximum number of timesteps ({cfg.max_timesteps}) reached."
    if (cfg.termination_conditions == "any" and dynamics.state.grid.get_tree_cover() > 0.99):
        satisifed_condition = f"Tree cover exceeds 99%."

    satisfied = len(satisifed_condition) > 0 # If 'satisfied_contion' is not an empty string, then a
                                             # termination condition was satisfied.
    if satisfied:
        print("\nSimulation terminated. Cause:", satisifed_condition)

    return satisfied


def do_burn_in(dynamics, cfg, forest_mask, color_dicts, target_treecover=1):
    """Perform burn-in procedure to stabilize initial forest density and demography.
    
    Forest is dispersed and allowed to grow for a number of timesteps, while feedback control is applied to keep
    tree cover close to the target tree cover. This allows the age distribution of trees to stabilize and prevents
    the model from starting with an unrealistically high tree cover (if the initial pattern image is very dense)
    or low tree cover (if the initial pattern image is very sparse). The burn-in procedure ends when either the
    specified burn-in duration has been reached, or when the tree cover exceeds the target tree cover, whichever
    comes later.
    """

    print(f"Beginning burn-in for {cfg.burnin_duration} timesteps...")
    cfg.init_csv = True
    cfg.treesize_bins = "initialize"
    #visualization_types = ["aggr_tree_LAI", "colored_patches"]
    visualization_types = []
    fire_freq_arrays = []
    fire_no_timesteps = 1
    patch_colors = {}
    while dynamics.time < cfg.burnin_duration or (cfg.initial_pattern_image == "none" and dynamics.state.grid.get_tree_cover() > target_treecover):
        print(f"Burn-in timestep:    {dynamics.time})")
        dynamics.disperse_within_forest(forest_mask)
        dynamics.grow()
        dynamics.induce_background_mortality()
        print("repopulating grid..")
        dynamics.state.repopulate_grid(cfg.verbosity)
        #dynamics.state.repopulate_grid(cfg.verbosity)
        if dynamics.time > (cfg.burnin_duration - 2):
            print("pruning..")
            dynamics.prune(forest_mask)
            print("repopulating grid..")
            dynamics.state.repopulate_grid(cfg.verbosity)
            print("resetting state distr..")
            dynamics.state.grid.reset_state_distr()
        
        dynamics.report_state()
        cfg = io.export_state(dynamics, path=cfg.csv_path, init_csv=cfg.init_csv, cfg=cfg)
        cfg.init_csv = False

        # Feedback control on tree cover during burn-in
        if cfg.initial_pattern_image == "none" and dynamics.state.grid.get_tree_cover() > target_treecover:
            prune_fraction = dynamics.state.grid.get_tree_cover() - target_treecover
            dynamics.prune_randomly(prune_fraction)

        # Do visualizations
        if ((dynamics.time % cfg.visualization_interval) == 0) and not cfg.headless:
            img = cfg.vis.visualize(
                dynamics.state, cfg.image_width, collect_states=1,
                color_dict=color_dicts.normal, cheap_visualization=True
            )
        dynamics.time += 1

    # Remove seedlings, then let trees disperse for one time step to create a realistic
    # initial distribution of seedlings.
    dynamics.time -= 1
    dynamics.remove_trees_up_to_age(1)
    dynamics.state.repopulate_grid(0)
    dynamics.update_forest_patch_detection()

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

    # -- MORTALITY TEMPLATE -- #
    dynamics.invoke_mortality_template()

    if cfg.verbosity > 0:
        print("Induced background mortality. Repopulating grid...")

    # Post-simulation cleanup and reporting
    dynamics.state.repopulate_grid(cfg.verbosity)
    dynamics.state.grid.redo_count()
    dynamics.report_state()
    dynamics.update_firefree_interval_averages()
    dynamics.update_forest_patch_detection()

    return dynamics


def export_animal_resources(dynamics):
    """Get color image representations of the resource grid."""

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


def do_update(dynamics, cfg):
    """Perform a single update of the model.

    This includes running the processes of dispersal, growth, and burning, as well as updating the
    state variables, creating visualizations, and exporting data to files. The function also checks if
    any of the termination conditions have been satisfied, and returns 'True' if the simulation should
    be terminated, and 'False' otherwise."""

    print("-- Starting iteration...") if cfg.verbosity else None
    dynamics = do_iteration(dynamics, cfg)
    print("-- Finished update") if cfg.verbosity else None
    
    # Obtain initial number of recruits
    if dynamics.time == 0:
        cfg.initial_no_dispersals = dynamics.get_initial_no_dispersals()
    
    # Track tree cover trajectory
    cfg.prev_tree_cover.append(dynamics.state.grid.get_tree_cover())
    if dynamics.time > 10:
        cfg.treecover_slope = (cfg.prev_tree_cover[-1] - cfg.treecover) / dynamics.time
        single_tstep_slope = cfg.prev_tree_cover[-1] - cfg.prev_tree_cover[-2]
        if abs(single_tstep_slope) > cfg.largest_absolute_slope:
            cfg.largest_absolute_slope = abs(single_tstep_slope)
    else:
        cfg.treecover_slope = 0

    # Track treecover trajectory
    do_terminate = termination_condition_satisfied(dynamics, cfg.computer_start_time, cfg)
    
    # Obtain patches from simulation
    cfg.old_no_patches = len(cfg.patches)
    cfg.patches = dynamics.get_patches() 
    cfg.patch_count_change =  len(cfg.patches) - cfg.old_no_patches
    
    # Create visualizations and export to image files
    if (dynamics.time % cfg.visualization_interval) == 0:
        visualization.do_visualizations(
            dynamics, cfg.fire_freq_arrays, cfg.fire_no_timesteps, cfg.verbosity,
            cfg.color_dicts, cfg.collect_states, cfg.visualization_types,
            cfg.patches, cfg.patch_count_change, cfg.patch_colors, cfg
        )
    
    # Export state variables to csv file
    print("-- Exporting state_data...") if cfg.verbosity else None
    cfg = io.export_state(
        dynamics, path=cfg.csv_path, init_csv=cfg.init_csv,
        tree_cover_slope=cfg.treecover_slope, cfg=cfg
    )

    # Update graphs 
    print("-- Showing graphs...") if cfg.verbosity else None
    if not cfg.headless:
        cfg.graphs.update()

    if export_animal_resources and (cfg.dispersal_mode == "all" or cfg.dispersal_mode == "animal"):
        export_animal_resources(dynamics)

    if do_terminate:
        return dynamics, True
        
    dynamics.time += 1

    return dynamics, False


def updateloop(dynamics, cfg):
    """Run the main update loop of the model, which repeatedly calls
    'do_update' until a termination condition is satisfied."""

    if not cfg.headless:
        cfg.graphs = visualization.Graphs(dynamics)

    print("Beginning simulation...")
    terminate_simulation = False
    while not terminate_simulation:
        dynamics, terminate_simulation = do_update(dynamics, cfg)

    if not cfg.headless:
        cv2.destroyAllWindows()
        dynamics.free()

    if not hasattr(cfg, "dependent_vars"):
        cfg.dependent_vars = {
            "tree_cover_slope": cfg.treecover_slope, 
            "largest_absolute_slope": cfg.largest_absolute_slope,
            "initial_no_dispersals": cfg.initial_no_dispersals
        }

    return dynamics, cfg


def main(**user_args):
    """Initialize the model, run the simulation, and return the final state of the model and other tracked variables."""

    global cfg

    if user_args["verbosity"] >= 0:
        print(f"Loading dbr_cpp from {cfg.BUILD_DIR} and preparing model run... (this may take some time)")

    # Set any given batch parameters
    if not h.strategy_distribution_params_are_loaded(user_args["strategy_distribution_params"]):
        with open(os.path.join(cfg.DATA_IN_DIR, user_args["strategy_distribution_params"]), "r") as sdp_jsonfile:
            user_args["strategy_distribution_params"] = json.load(sdp_jsonfile)

    args = SimpleNamespace(**user_args)
    cfg = h.apply_user_args_to_configuration(args, cfg)

    dynamics = init(cfg)
    return updateloop(dynamics, cfg)
 


    

