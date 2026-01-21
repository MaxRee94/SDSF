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


def unpack_control_keys(control_variable):
    control_keys = control_variable.split("->")
    return control_keys


def set_dispersal_kernel(
        dynamics, dispersal_mode, multi_disperser_params
    ):
    animal_species = []
    with open(os.path.join(cfg.DATA_IN_DIR, multi_disperser_params), "r") as mdp_jsonfile:
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
        animal_dispersal_params = multi_disperser_params["animal"]
        animal_species = [species for species in animal_dispersal_params.keys() if species != "population"]
        dynamics.set_global_animal_kernel(multi_disperser_params["animal"])
    elif (dispersal_mode == "all"):
        animal_dispersal_params = multi_disperser_params["animal"]
        animal_species = [species for species in animal_dispersal_params.keys() if species != "population"]
        multi_disperser_params.pop("animal")
        dynamics.set_global_kernels(multi_disperser_params, animal_dispersal_params)

    print("finished setting dispersal kernel.")

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
        img = sng.generate(scale=2, grid_width=cfg.grid_width, binary_connectivity=cfg.treecover, cfg=cfg)
        cfg.cover_img_path = f"{cfg.SIMPLE_PATTERNS_DIR}/whitenoise.png"
        cv2.imwrite(cfg.cover_img_path, img)
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
                print(f"Lookup table file {fpath} not found. Creating new one...")
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
    if (cfg.termination_conditions == "all" and dynamics.state.grid.get_tree_cover() > 0.93):
        condition = f"Tree cover exceeds 93%."
    satisfied = len(condition) > 0
    if satisfied:
        print("\nSimulation terminated. Cause:", condition)

    return satisfied


def do_visualizations(dynamics, fire_freq_arrays, fire_no_timesteps, verbose, color_dicts, collect_states, visualization_types, patches, patch_count_change, patch_color_ids, cfg):
    if ("recruitment" in visualization_types):
        print("Saving recruitment img...") if verbose else None
        recruitment_img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["recruitment"], collect_states=0)
        imagepath_recruitment = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/recruitment/" + str(dynamics.time) + ".png")
        cfg.vis.save_image(recruitment_img, imagepath_recruitment, get_max(1000, recruitment_img.shape[0]), interpolation="none")
    
    if ("fire_freq" in visualization_types):
        fire_freq_arrays.append(dynamics.state.grid.get_distribution(0) == -5)
        if dynamics.time > fire_no_timesteps:
            print("Saving fire frequency img...") if verbose else None
            fire_freq_img = cfg.vis.get_fire_freq_image(fire_freq_arrays[-fire_no_timesteps:], color_dicts["fire_freq"], dynamics.state.grid.width, fire_no_timesteps)
            imagepath_fire_freq = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fire_frequencies/" + str(dynamics.time) + ".png")
            cfg.vis.save_image(fire_freq_img, imagepath_fire_freq, get_max(1000, fire_freq_img.shape[0]))

    def get_bbox(positions2d):
        min_x = min(positions2d, key=lambda item: item[0])[0]
        max_x = max(positions2d, key=lambda item: item[0])[0]
        min_y = min(positions2d, key=lambda item: item[1])[1]
        max_y = max(positions2d, key=lambda item: item[1])[1]
        return (min_x, min_y), (max_x, max_y)

    if ("colored_patches" in visualization_types):
        print("Creating colored patches image...") if verbose else None

        # Initialize color index array
        patch_colors_indices = np.zeros((dynamics.state.grid.width, dynamics.state.grid.width), dtype=int) -1
        
        # Sort patches so that larger patches are assigned colors first
        patches = sorted(patches, key=lambda patch: patch["area"], reverse=True)

        minimum_considered_patch_size = 0 # m^2
        forest_area = 0
        central_patch = -1
        max_no_colors = 1000
        offset = -10

        for i, patch in enumerate(patches):
            if (len(patches) > 20) and dynamics.time == 0:
                if i >= 100:
                    print("Breaking patch coloring loop at patch no", i, f"for performance reasons (number of patches = {len(patches)}).")
                    break
                if i % 15 == 0:
                    print(f"Coloring patch no {i}/{len(patches)}... Will be cut off at 100.")
            if patch["area"] < minimum_considered_patch_size: # Only consider patches of a certain size
                continue
            patch_id = patch["id"]

            if not patch_color_ids.get(str(patch_id)):
                if (len(patches) > 30) and dynamics.time > 0:
                    color_idx = cfg.vis.get_random_color_index(patch_color_ids.values(), max_no_colors, offset)
                else:
                    color_idx = cfg.vis.get_most_distinct_index(patch_color_ids.values(), max_no_colors, offset)
                patch_color_ids[str(patch_id)] = color_idx # Assign a new color index to the patch
            patch_color_id = patch_color_ids[str(patch_id)]
            col=color_dicts["colored_patches"][patch_color_id]
            forest_area += patch["area"] if patch["type"] == "forest" else 0
            for cell in patch["cells"]:
                patch_colors_indices[cell[1]][cell[0]] = patch_color_id

        colored_patches_img = cfg.vis.get_image(patch_colors_indices, color_dicts["colored_patches"], dynamics.state.grid.width)
        cfg.show_edges = False # Hardcoded for now
        if cfg.show_edges:
            for patch in patches:
                if patch["area"] < minimum_considered_patch_size: # Only consider patches of a certain size
                    continue
                for edge in patch["perimeter"]:
                    point1 = (edge[0][0], edge[0][1])
                    point2 = (edge[1][0], edge[1][1])
                    if get_2d_dist(point1, point2) > 0.5*dynamics.state.grid.width:
                        continue # Skip drawing wrap-around edges for now.
                    cv2.line(colored_patches_img, point1, point2, (0, 0, 0), 1)  # black line, thickness=2
        imagepath_colored_patches = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/colored_patches/" + str(dynamics.time) + ".png")
        cfg.vis.save_image(colored_patches_img, imagepath_colored_patches, get_max(1000, colored_patches_img.shape[0]), interpolation="none")

    print("-- Visualizing image...") if verbose else None
    if cfg.export_visualizations:
        # Get a color image representation of the state
        img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["normal"], collect_states=1)
    else:
        # Get a color image representation of the state and show it.
        img = cfg.vis.visualize(
            dynamics.state.grid, cfg.image_width, collect_states=collect_states,
            color_dict=color_dicts["normal"]
        )

    print("-- Saving image...") if verbose else None
    imagepath = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
    cfg.vis.save_image(img, imagepath, get_max(1000, img.shape[0]), interpolation="none")
    
    # if ("fuel" in visualization_types):
    #     print("-- Saving fuel image...") if verbose else None
    #     fuel_img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="fuel")
    #     imagepath_fuel = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fuel/" + str(dynamics.time) + ".png")
    #     cfg.vis.save_image(fuel_img, imagepath_fuel, get_max(1000, fuel_img.shape[0]))

    if ("aggr_tree_LAI" in visualization_types):
        print("-- Saving aggregated tree LAI image...") if verbose else None
        aggr_tree_LAI_img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="aggr_tree_LAI", invert=False)
        imagepath_aggr_tree_LAI = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/aggr_tree_LAI/" + str(dynamics.time) + ".png")
        cfg.vis.save_image(aggr_tree_LAI_img, imagepath_aggr_tree_LAI, get_max(1000, aggr_tree_LAI_img.shape[0]))

    # if ("fuel_penetration" in visualization_types):
    #     print("-- Saving fuel penetration image...") if verbose else None
    #     fuel_penetration_img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="fuel_penetration")
    #     imagepath_fuel_penetration = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fuel_penetration/" + str(dynamics.time) + ".png")
    #     cfg.vis.save_image(fuel_penetration_img, imagepath_fuel_penetration, get_max(1000, fuel_penetration_img.shape[0]))


def do_burn_in(dynamics, cfg, forest_mask, color_dicts, target_treecover=1):
    """Perform burn-in procedure to stabilize initial forest density and demography."""
    
    print(f"Beginning burn-in for {cfg.burnin_duration} timesteps...")
    init_csv = True
    cfg.treesize_bins = "initialize"
    cur_treecover = 0
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
        cur_treecover = dynamics.state.grid.get_tree_cover()

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
            do_visualizations(
                dynamics, fire_freq_arrays, fire_no_timesteps, verbose,
                color_dicts, collect_states, visualization_types,
                patches, patch_count_change, patch_colors, cfg
            )
    
        print("-- Exporting state_data...") if verbose else None
        cfg = io.export_state(
            dynamics, cfg.csv_path, init_csv,
            tree_cover_slope=slope,
            cfg=SimpleNamespace(**vars(cfg))
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

    return dynamics, slope, largest_absolute_slope, initial_no_dispersals, cfg


def test_kernel():
    cpp.init_RNG()
    dist_max = 200
    windspeed_gmean = 10
    windspeed_stdev = 5
    seed_terminal_speed = 0.65
    abscission_height = 30
    wind_kernel = cpp.Kernel(1, dist_max, windspeed_gmean, windspeed_stdev, 0, 3600, seed_terminal_speed, abscission_height)
    wind_kernel.build()
    cfg.vis.visualize_kernel(wind_kernel, "Wind kernel. d_max = {}, w_gmean = {}, \n w_stdev = {}, v_t = {}, h = {}".format(
        dist_max, windspeed_gmean, windspeed_stdev, seed_terminal_speed, abscission_height)
    )
    return


def strategy_distribution_params_are_loaded(strategy_distribution_params):
    return type(strategy_distribution_params) == dict


def main(batch_parameters=None, **user_args): 
    
    global cfg

    # Set any given batch parameters
    if not strategy_distribution_params_are_loaded(user_args["strategy_distribution_params"]):
        with open(os.path.join(cfg.DATA_IN_DIR, user_args["strategy_distribution_params"]), "r") as sdp_jsonfile:
            user_args["strategy_distribution_params"] = json.load(sdp_jsonfile)
    if batch_parameters:
        print("-- Setting batch parameters: ", batch_parameters)
        if (type(batch_parameters) == str):
            batch_parameters = json.loads(batch_parameters)
        if "strategies->" in batch_parameters["control_variable"]:
            control_keys = unpack_control_keys(batch_parameters["control_variable"])
            control_keys.remove("strategies")
            user_args["strategy_distribution_params"][control_keys[0]][control_keys[1]] = batch_parameters["control_value"]
        elif "<idx>" in batch_parameters["control_variable"]:
            control_keys = unpack_control_keys(batch_parameters["control_variable"])
            idx = int(control_keys[1].split("<idx>")[1])
            user_args[control_keys[0]][idx] = batch_parameters["control_value"]

    args = SimpleNamespace(**user_args)
    cfg = apply_user_args_to_configuration(args, cfg)

    dynamics, color_dicts, cfg = init(cfg)
    return updateloop(dynamics, color_dicts, cfg)
 


    

