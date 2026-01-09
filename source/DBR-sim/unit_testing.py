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


def init_tests(
    timestep=None, grid_width=None, cell_width=None, max_dbh=None, image_width=None,
    treecover=None, self_ignition_factor=None, flammability=None,
    rainfall=None, unsuppressed_flammability=None, 
    verbosity=None, dbh_q1=None, dbh_q2=None, seed_bearing_threshold=None,
    dispersal_mode=None, linear_diffusion_q1=None, linear_diffusion_q2=None,
    dispersal_min=None, dispersal_max=None, growth_rate_multiplier=None, 
    flammability_coefficients_and_constants=None, saturation_threshold=None, fire_resistance_params=None,
    background_mortality=None, headless=False, wind_dispersal_params=None, animal_dispersal_params=None,
    multi_disperser_params=None, strategy_distribution_params=None, resource_grid_width=None,
    initial_pattern_image=None, mutation_rate=None, growth_rate_multiplier_params=None, 
    random_seed=None, firefreq_random_seed=None, enforce_no_recruits=None, **user_args
):
    print("Starting unit tests...")

    # Obtain strategy distribution parameters
    with open(os.path.join(cfg.DATA_IN_DIR, strategy_distribution_params), "r") as sdp_jsonfile:
        strategy_distribution_params = json.load(sdp_jsonfile)    

    tests = cpp.Tests(timestep, cell_width, self_ignition_factor, rainfall, seed_bearing_threshold,
        growth_rate_multiplier, unsuppressed_flammability, flammability_coefficients_and_constants[0],
        flammability_coefficients_and_constants[1], flammability_coefficients_and_constants[2], 
        flammability_coefficients_and_constants[3], max_dbh, saturation_threshold, fire_resistance_params[0],
        fire_resistance_params[1], fire_resistance_params[2], background_mortality, strategy_distribution_params, 
        resource_grid_width, mutation_rate, verbosity, grid_width, dbh_q1, dbh_q2, growth_rate_multiplier_params[0],
        growth_rate_multiplier_params[1], growth_rate_multiplier_params[2], random_seed, firefreq_random_seed,
        enforce_no_recruits
    )
    
    return tests


def main(batch_parameters=None, **user_args): 

    #test_kernel()
    #test_discrete_probmodel()
    #return

    tests = init_tests(**user_args)
    tests.run_all()
