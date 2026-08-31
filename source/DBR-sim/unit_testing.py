"""Unit testing module for DBR-sim.

Contains Python-based unit tests for individual components (kernel, probability models)
and integration with C++-based test suite.

Run via: python __main__.py --tests unit_tests [--additional_args]
"""

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
import random
import math

import visualization as vis
import file_handling as io
from helpers import *
from config import *

#from x64.Debug import dbr_cpp as cpp
from x64.Release import dbr_cpp as cpp


def test_kernel():
    """Test the wind dispersal kernel visualization."""
    cpp.init_RNG(42)  # Initialize with seed 42
    dist_max = 200
    windspeed_gmean = 10
    windspeed_stdev = 5
    seed_terminal_speed = 0.65
    abscission_height = 30
    wind_kernel = cpp.Kernel(1, dist_max, windspeed_gmean, windspeed_stdev, 0, 3600, seed_terminal_speed, abscission_height)
    wind_kernel.build()
    # Initialize visualization for kernel visualization
    import visualization
    vis = visualization.Visualiser(cfg)
    vis.visualize_kernel(wind_kernel, "Wind kernel. d_max = {}, w_gmean = {}, \n w_stdev = {}, v_t = {}, h = {}".format(
        dist_max, windspeed_gmean, windspeed_stdev, seed_terminal_speed, abscission_height)
    )
    return


def test_discrete_probmodel():
    """Test the discrete probability model sampling."""
    cpp.init_RNG(42)  # Initialize with seed 42
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


def init_tests(
    timestep=None, grid_width=None, cell_width=None, max_dbh=None, image_width=None,
    treecover=None, self_ignition_factor=None, flammability=None,
    rainfall=None, unsuppressed_flammability=None, 
    verbosity=None, seed_bearing_threshold=None,
    dispersal_mode=None, linear_diffusion_q1=None, linear_diffusion_q2=None,
    dispersal_min=None, dispersal_max=None, growth_rate_multiplier=None, 
    flammability_coefficients_and_constants=None, saturation_threshold=None, fire_resistance_params=None,
    background_mortality=None, headless=False, wind_dispersal_params=None, animal_dispersal_params=None,
    multi_disperser_params=None, strategy_distribution_params=None, resource_grid_width=None,
    initial_pattern_image=None, mutation_rate=None, growth_rate_multiplier_params=None, 
    random_seed=None, firefreq_random_seed=None, enforce_no_recruits=None, **user_args
):
    """
    Initialize unit tests with the given parameters.

    Parameters
    ----------
    timestep : float, optional
        Duration of a single timestep in years.
    grid_width : int, optional
        Width of the spatial domain in cells.
    cell_width : float, optional
        Width of each grid cell in meters.
    max_dbh : float, optional
        Maximum dbh (diameter at breast height) in cm.
    image_width : int, optional
        Width of the visualization image in pixels.
    treecover : float, optional
        Fraction of domain occupied by tree cells.
    self_ignition_factor : float, optional
        Expected number of ignitions per year per km^2.
    flammability : float, optional
        Flammability parameter.
    rainfall : float, optional
        Rainfall parameter.
    unsuppressed_flammability : float, optional
        Flammability of savanna cells.
    verbosity : int, optional
        Verbosity level.
    seed_bearing_threshold : float, optional
        Fraction of max radius above which a tree is seed-bearing.
    dispersal_mode : str, optional
        Dispersal mode.
    linear_diffusion_q1 : float, optional
        Linear diffusion parameter 1.
    linear_diffusion_q2 : float, optional
        Linear diffusion parameter 2.
    dispersal_min : float, optional
        Minimum dispersal distance.
    dispersal_max : float, optional
        Maximum dispersal distance.
    growth_rate_multiplier : float, optional
        Growth rate multiplier.
    flammability_coefficients_and_constants : list, optional
        Flammability model coefficients.
    saturation_threshold : float, optional
        Tree density saturation threshold.
    fire_resistance_params : list, optional
        Fire resistance sigmoid parameters.
    background_mortality : float, optional
        Background mortality rate.
    headless : bool, optional
        Run in headless mode.
    wind_dispersal_params : dict, optional
        Wind dispersal parameters.
    animal_dispersal_params : dict, optional
        Animal dispersal parameters.
    multi_disperser_params : str, optional
        Path to multi-disperser parameters JSON file.
    strategy_distribution_params : str, optional
        Path to strategy distribution parameters JSON file.
    resource_grid_width : int, optional
        Width of the resource grid.
    initial_pattern_image : str, optional
        Path to initial pattern image.
    mutation_rate : float, optional
        Mutation rate.
    growth_rate_multiplier_params : list, optional
        Growth rate multiplier distribution parameters.
    random_seed : int, optional
        Random seed.
    firefreq_random_seed : int, optional
        Fire frequency random seed.
    enforce_no_recruits : int, optional
        Recruit enforcement parameter.
    **user_args : dict
        Additional user arguments.

    Returns
    -------
    cpp.Tests
        Initialized test suite object.
    """
    print("Starting unit tests...")

    # Set defaults for optional parameters
    if flammability_coefficients_and_constants is None:
        flammability_coefficients_and_constants = [0, 0, 0, 0]
    if fire_resistance_params is None:
        fire_resistance_params = [8.5, 50, 2.857]
    if growth_rate_multiplier_params is None:
        growth_rate_multiplier_params = [0, 1.0, 1.0]

    # Obtain strategy distribution parameters
    with open(os.path.join(cfg.DATA_IN_DIR, strategy_distribution_params), "r") as sdp_jsonfile:
        strategy_distribution_params = json.load(sdp_jsonfile)    

    tests = cpp.Tests(timestep, cell_width, self_ignition_factor, rainfall, seed_bearing_threshold,
        growth_rate_multiplier, unsuppressed_flammability, flammability_coefficients_and_constants[0],
        flammability_coefficients_and_constants[1], flammability_coefficients_and_constants[2], 
        flammability_coefficients_and_constants[3], max_dbh, saturation_threshold, fire_resistance_params[0],
        fire_resistance_params[1], fire_resistance_params[2], background_mortality, strategy_distribution_params, 
        resource_grid_width, mutation_rate, verbosity, grid_width, growth_rate_multiplier_params[0],
        growth_rate_multiplier_params[1], growth_rate_multiplier_params[2], random_seed, firefreq_random_seed,
        enforce_no_recruits
    )
    
    return tests


def main(batch_parameters=None, **user_args): 
    """
    Run unit tests.

    Runs both Python-based unit tests (test_kernel, test_discrete_probmodel)
    and C++-based tests via the Tests suite.

    Parameters
    ----------
    batch_parameters : dict, optional
        Batch parameters (unused).
    **user_args : dict
        User arguments passed to test initialization.
    """
    # Run Python-based unit tests
    print("Running Python unit tests...")
    try:
        test_kernel()
        print("[PASS] test_kernel")
    except Exception as e:
        print(f"[FAIL] test_kernel: {e}")
    
    try:
        test_discrete_probmodel()
        print("[PASS] test_discrete_probmodel")
    except Exception as e:
        print(f"[FAIL] test_discrete_probmodel: {e}")

    # Run C++-based tests
    print("\nRunning C++ unit tests...")
    tests = init_tests(**user_args)
    tests.run_all()
    print("[COMPLETE] C++ unit tests")
