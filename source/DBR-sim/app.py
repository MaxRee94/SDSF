from dataclasses import dataclass
import sys
from tkinter import E, N, Y

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
from types import SimpleNamespace

import visualization as vis
import file_handling as io
from helpers import *
from config import *
import controllable_pattern_generator as cpg

#from x64.Debug import dbr_cpp as cpp
from x64.Release import dbr_cpp as cpp

print(cpp.__file__)

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


def set_initial_tree_cover(dynamics, args):
    print("Initial pattern image:", args.initial_pattern_image)
    if args.initial_pattern_image == "none":
        dynamics.state.set_tree_cover(args.treecover)
    elif args.initial_pattern_image == "ctrl":
        img, path, benchmark_cover = vis.generate_controllable_pattern_image(**vars(args))
        # If the user has set the override_image_treecover to 2, we will use the benchmark cover value from the controllable pattern image.
        # The benchmark cover is the cover for a version of the pattern produced using the given parameters, but with a sine amplitude set to 0 (i.e., with circular disks).
        if args.override_image_treecover == 2:
            args.override_image_treecover = benchmark_cover
    elif args.initial_pattern_image == "perlin_noise":
        path = f"{cfg.DATA_IN_DIR}/state_patterns/" + args.initial_pattern_image
        if "perlin_noise" == args.initial_pattern_image:
            print("Setting tree cover using perlin noise function...")
            path = f"{cfg.PERLIN_NOISE_DIR}/" + args.initial_pattern_image + ".png"
            noise_frequency = 5.0 / args.patch_width # Convert patch width to noise frequency
            noise_frequency = round(noise_frequency, 2) # Conform noise frequency to 2 decimal places, to ensure periodicity of the noise pattern
            vis.generate_perlin_noise_image(path, frequency=noise_frequency, octaves=args.noise_octaves)
    else:
        print("Setting tree cover from image...")
    
    if args.initial_pattern_image != "none":
        img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        print(f"Path of the generated {args.initial_pattern_image} pattern image: ", path) if args.initial_pattern_image in ["ctrl", "perlin_noise"] else None
        if img is None and args.initial_pattern_image == "perlin_noise": # Hotfix; sometimes perlin noise-generated images fail to load properly
            print("Image not found. Generating new one...")
            vis.generate_perlin_noise_image(path, frequency=noise_frequency, octaves=args.noise_octaves)
            
            img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
            if img is None:
                time.sleep(1)
                img = cv2.imread(path, cv2.IMREAD_GRAYSCALE) # Wait for the image to be generated before trying to load it again
            print("Path perlin noise image: (attempt 2)", path)
            if img is None and args.initial_pattern_image == "perlin_noise": # Hotfix; sometimes perlin noise-generated images fail to load properly
                print("Image not found. Generating new one...")
                vis.generate_perlin_noise_image(path, frequency=noise_frequency, octaves=args.noise_octaves)
                time.sleep(1) # Wait for the image to be generated before trying to load it again
                img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        
        if ("perlin_noise" in path) and (not "thresholded.png" in path):
            img = vis.get_thresholded_image(img, args.treecover * img.shape[0] * img.shape[0] * 255 )
            cv2.imwrite(path.replace(".png", "_thresholded.png"), img)
 
        img = cv2.resize(img, (dynamics.state.grid.width, dynamics.state.grid.width), interpolation=cv2.INTER_LINEAR)
        print("Created image. Setting cover...")
        dynamics.state.set_cover_from_image(img / 255, args.override_image_treecover)
    dynamics.state.repopulate_grid(0)

    return dynamics, args


def init(user_args):
    args = SimpleNamespace(**user_args)

    # Set random seed for fire frequency probability distribution. If -999 is given, a random seed will be generated. Otherwise, the given seed will be used.
    if args.firefreq_random_seed == -999:
        args.firefreq_random_seed = random.randint(0, 1000000)
        print("Generated random seed for fire frequency probability distribution: ", args.firefreq_random_seed)
    print("Using fire frequency random seed: ", args.firefreq_random_seed)
    # Initialize dynamics object and state
    print("type:", type(user_args["fire_resistance_params"]))
    print("Initializing dynamics object with args:\n", user_args)
    print("LAI aggregation radius:", args.LAI_aggregation_radius)
    dynamics = cpp.create_dynamics(user_args)
    dynamics.init_state(
        args.grid_width,
        args.dbh_q1,
        args.dbh_q2,
        args.growth_rate_multiplier_params[0],
        args.growth_rate_multiplier_params[1],
        args.growth_rate_multiplier_params[2],
        args.minimum_patch_size,
        args.LAI_aggregation_radius
    )

    # Set 
    
    # Set dispersal kernel
    dynamics, animal_species = set_dispersal_kernel(dynamics, args.dispersal_mode, args.multi_disperser_params)
    
    # Set initial tree cover
    dynamics, args = set_initial_tree_cover(dynamics, args)
    
    # Create color dictionaries for visualizations
    no_colors = 100
    if args.display_fire_effects == 1:
        color_dict = vis.get_color_dict(no_colors, begin=0.2, end=0.5, distr_type="normal_with_fire_effects")
    else:
        color_dict = vis.get_color_dict(no_colors, begin=0.2, end=0.5, distr_type="normal")
    color_dict_recruitment = vis.get_color_dict(no_colors, begin=0.2, end=0.5, distr_type="recruitment")
    color_dict_fire_freq = vis.get_color_dict(10, begin=0.2, end=0.5, distr_type="fire_freq")
    color_dict_blackwhite = vis.get_color_dict(no_colors, distr_type="blackwhite")
    color_dict_colored_patches = vis.get_color_dict(no_colors, distr_type="colored_patches")
    color_dicts = {}
    color_dicts["normal"] = color_dict
    color_dicts["recruitment"] = color_dict_recruitment
    color_dicts["fire_freq"] = color_dict_fire_freq
    color_dicts["blackwhite"] = color_dict_blackwhite
    color_dicts["colored_patches"] = color_dict_colored_patches
    
    # Visualize the initial state
    collect_states = True
    print("Visualizing state at t = 0")
    if not args.headless:
        # Get a color image representation of the initial state and show it.
        img = vis.visualize(
            dynamics.state.grid, args.image_width, collect_states=collect_states,
            color_dict=color_dict
        )
    else:
        # Get a color image representation of the initial state
        img = vis.get_image_from_grid(dynamics.state.grid, color_dict, collect_states=collect_states)
   
    # Export image file
    imagepath = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
    vis.save_image(img, imagepath)
    
    if args.dispersal_mode == "all" or args.dispersal_mode == "animal":

        # Export resource grid lookup table
        print("animal species: ", animal_species)
        for species in animal_species:
            lookup_table, fpath = io.get_lookup_table(species, args.grid_width, args.resource_grid_width * args.resource_grid_width)
            #lookup_table = None # Hotfix; lookup table might lead to memory leaks (?)
            if lookup_table is None:
                print(f"Lookup table file {fpath} not found. Creating new one...")
                dynamics.precompute_resourcegrid_lookup_table(species)
                lookup_table = dynamics.get_resource_grid_lookup_table(species)
                io.export_lookup_table(lookup_table, args.grid_width, species)
            else:
                print(f"Lookup table file {fpath} found. Loading...")
                dynamics.set_resource_grid_lookup_table(lookup_table, species)

    return dynamics, color_dicts


def termination_condition_satisfied(dynamics, start_time, user_args):
    satisfied = False
    condition = ""
    if (dynamics.time >= int(user_args["max_timesteps"])):
        condition = f"Maximum number of timesteps ({user_args['max_timesteps']}) reached."
    if (user_args["termination_conditions"] == "all" and dynamics.state.grid.get_tree_cover() > 0.93):
        condition = f"Tree cover exceeds 93%."
    satisfied = len(condition) > 0
    if satisfied:
        print("\nSimulation terminated. Cause:", condition)

    return satisfied


def do_visualizations(dynamics, fire_freq_arrays, fire_no_timesteps, verbose, color_dicts, collect_states, visualization_types, patches, patch_color_ids, user_args):
    if ("recruitment" in visualization_types):
        print("Saving recruitment img...") if verbose else None
        recruitment_img = vis.get_image_from_grid(dynamics.state.grid, color_dicts["recruitment"], collect_states=0)
        imagepath_recruitment = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/recruitment/" + str(dynamics.time) + ".png")
        vis.save_image(recruitment_img, imagepath_recruitment, get_max(1000, recruitment_img.shape[0]), interpolation="none")
    
    if ("fire_freq" in visualization_types):
        fire_freq_arrays.append(dynamics.state.grid.get_distribution(0) == -5)
        if dynamics.time > fire_no_timesteps:
            print("Saving fire frequency img...") if verbose else None
            fire_freq_img = vis.get_fire_freq_image(fire_freq_arrays[-fire_no_timesteps:], color_dicts["fire_freq"], dynamics.state.grid.width, fire_no_timesteps)
            imagepath_fire_freq = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fire_frequencies/" + str(dynamics.time) + ".png")
            vis.save_image(fire_freq_img, imagepath_fire_freq, get_max(1000, fire_freq_img.shape[0]))

    def get_bbox(positions2d):
        min_x = min(positions2d, key=lambda item: item[0])[0]
        max_x = max(positions2d, key=lambda item: item[0])[0]
        min_y = min(positions2d, key=lambda item: item[1])[1]
        max_y = max(positions2d, key=lambda item: item[1])[1]
        return (min_x, min_y), (max_x, max_y)

    if ("colored_patches" in visualization_types):
        print("Creating colored patches image...") if verbose else None

        # Modify color array to give patches distinct colors
        patch_colors_indices = dynamics.state.grid.get_distribution(False)
        
        # Sort patches so that larger patches are assigned colors first
        patches = sorted(patches, key=lambda patch: patch["area"], reverse=True)

        minimum_considered_patch_size = 0 # m^2
        forest_area = 0
        central_patch = -1
        max_no_colors = 1000
        offset = -10
        for i, patch in enumerate(patches):
            if patch["area"] < minimum_considered_patch_size: # Only consider patches of a certain size
                continue
            patch_id = patch["id"]

            if not patch_color_ids.get(str(patch_id)):
                if (len(patches) > 30) and dynamics.time > 1:
                    color_idx = vis.get_random_color_index(patch_color_ids.values(), max_no_colors, offset)
                else:    
                    color_idx = vis.get_most_distinct_index(patch_color_ids.values(), max_no_colors, offset)
                patch_color_ids[str(patch_id)] = color_idx # Assign a new color index to the patch
            patch_color_id = patch_color_ids[str(patch_id)]
            col=color_dicts["colored_patches"][patch_color_id]
            forest_area += patch["area"] if patch["type"] == "forest" else 0
            for cell in patch["cells"]:
                patch_colors_indices[cell[1]][cell[0]] = patch_color_id
        
        colored_patches_img = vis.get_image(patch_colors_indices, color_dicts["colored_patches"], dynamics.state.grid.width)
        user_args["show_edges"] = False # Hardcoded for now
        if user_args["show_edges"]:
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
        vis.save_image(colored_patches_img, imagepath_colored_patches, get_max(1000, colored_patches_img.shape[0]), interpolation="none")


    print("-- Visualizing image...") if verbose else None
    if user_args["headless"]:
        # Get a color image representation of the initial state
        img = vis.get_image_from_grid(dynamics.state.grid, color_dicts["normal"], collect_states=1)
    else:
        # Get a color image representation of the initial state and show it.
        img = vis.visualize(
            dynamics.state.grid, user_args["image_width"], collect_states=collect_states,
            color_dict=color_dicts["normal"]
        )

    print("-- Saving image...") if verbose else None
    imagepath = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
    vis.save_image(img, imagepath, get_max(1000, img.shape[0]), interpolation="none")
    
    if ("fuel" in visualization_types):
        print("-- Saving fuel image...") if verbose else None
        fuel_img = vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="fuel")
        imagepath_fuel = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fuel/" + str(dynamics.time) + ".png")
        vis.save_image(fuel_img, imagepath_fuel, get_max(1000, fuel_img.shape[0]))
        print("fuel done.")

    if ("aggr_tree_LAI" in visualization_types):
        print("-- Saving aggregated tree LAI image...") if verbose else None
        aggr_tree_LAI_img = vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="aggr_tree_LAI", invert=False)
        imagepath_aggr_tree_LAI = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/aggr_tree_LAI/" + str(dynamics.time) + ".png")
        vis.save_image(aggr_tree_LAI_img, imagepath_aggr_tree_LAI, get_max(1000, aggr_tree_LAI_img.shape[0]))

    if ("fuel_penetration" in visualization_types):
        print("-- Saving fuel penetration image...") if verbose else None
        fuel_penetration_img = vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="fuel_penetration")
        imagepath_fuel_penetration = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fuel_penetration/" + str(dynamics.time) + ".png")
        vis.save_image(fuel_penetration_img, imagepath_fuel_penetration, get_max(1000, fuel_penetration_img.shape[0]))


def updateloop(dynamics, color_dicts, **user_args):
    start = time.time()
    print("Beginning simulation...")
    csv_path = user_args["csv_path"]
    visualization_types = ["fire_freq", "recruitment", "fuel", "tree_LAI", "aggr_tree_LAI", "colored_patches", "fuel_penetration"]
    init_csv = True
    prev_tree_cover = [user_args["treecover"]] * 60
    slope = 0
    largest_absolute_slope = 0
    export_animal_resources = False
    collect_states = 1
    fire_no_timesteps = 1
    patch_colors = {}
    verbose = user_args["verbosity"]
    fire_freq_arrays = []
    color_dict_fire_freq = vis.get_color_dict(fire_no_timesteps, begin=0.2, end=0.5, distr_type="fire_freq")
    color_dicts["fire_freq"] = color_dict_fire_freq
    if not user_args["headless"]:
        graphs = vis.Graphs(dynamics)
    while True:
        print("-- Starting update (calling from python)") if verbose else None
        try:
            dynamics.update()
        except Exception as e:
            print("Error during dynamics update:\n", e)
            exit(1)
        print("-- Finished update") if verbose else None
        
        # Obtain initial number of recruits
        if dynamics.time == 1:
            initial_no_dispersals = dynamics.get_initial_no_dispersals()
        
        # Track tree cover trajectory and evaluate termination conditions
        prev_tree_cover.append(dynamics.state.grid.get_tree_cover())
        if dynamics.time > 10:
            slope = (prev_tree_cover[-1] - user_args["treecover"]) / dynamics.time
            single_tstep_slope = prev_tree_cover[-1] - prev_tree_cover[-2]
            if abs(single_tstep_slope) > largest_absolute_slope:
                largest_absolute_slope = abs(single_tstep_slope)
        else:
            slope = 0
        do_terminate = termination_condition_satisfied(dynamics, start, user_args)
        
        # WIP: Obtain forest patch perimeters from simulation
        patches = dynamics.get_patches()
        print("----- Number of patches: ", len(patches)) 
        if verbose:
            for patch in patches:
                print(f"No cells in patch {patch['id']}: {len(patch['cells'])}, example cell position: ({patch['cells'][0]}), \n" +
                    f"centroid: ({patch['centroid'][0], patch['centroid'][1]}, area: {patch['area']} m^2, perimeter length: {patch['perimeter_length']} m.")
        
        # Do visualizations
        if not user_args["headless"]:
            do_visualizations(dynamics, fire_freq_arrays, fire_no_timesteps, verbose, color_dicts, collect_states, visualization_types, patches, patch_colors, user_args)
        
        print("-- Exporting state_data...") if verbose else None
        csv_path = io.export_state(dynamics, csv_path, init_csv, tree_cover_slope=slope, args=SimpleNamespace(**user_args))
        init_csv = False
        
        
        print("-- Saving tree positions...") if verbose else None
        #io.save_state(dynamics)
        
        if user_args["report_state"] == "True" or user_args["report_state"] == True:
            io.update_state_report(dynamics)

        print("-- Showing graphs...") if verbose else None
        if not user_args["headless"]:
            graphs.update()

        if export_animal_resources and (user_args["dispersal_mode"] == "all" or user_args["dispersal_mode"] == "animal"):
            # Get color image representations of the resource grid from the last iteration
            cover_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/cover/" + str(dynamics.time) + ".png")
            fruits_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fruits/" + str(dynamics.time) + ".png")
            visits_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/visits/" + str(dynamics.time) + ".png")
            distance_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/distance/" + str(dynamics.time) + ".png")
            distance_coarse_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/distance_coarse/" + str(dynamics.time) + ".png")
            k_coarse_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/k_coarse/" + str(dynamics.time) + ".png")
            k_path = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/k/" + str(dynamics.time) + ".png")
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
    
    return dynamics, slope, largest_absolute_slope, initial_no_dispersals


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




def main(batch_parameters=None, **user_args): 

    #vis.visualize_legend()
    #return

    # Set any given batch parameters
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

    user_args["patch_width"] = round(user_args["patch_width"], 2)

    dynamics, color_dicts = init(user_args)
    return updateloop(dynamics, color_dicts, **user_args)
 


    

