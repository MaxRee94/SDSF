from cgitb import lookup
import csv
import os
from pathlib import Path
import datetime
from shutil import ExecError
import cv2
import time
from config import *
import numpy as np

import sine_pattern_generator as spg



def get_tree_sizes(dynamics):
    _tree_sizes = dynamics.state.get_tree_sizes()
    counts, bins = np.histogram(_tree_sizes, bins=5)
    return counts


def get_firefree_interval_stats(dynamics, _type):
    intervals = dynamics.get_firefree_intervals(_type)
    stdev = np.std(intervals)
    mean = np.mean(intervals)
    return mean, stdev


def load_json_config(file_path):
    """Load a JSON configuration file.
    
    Params:
        file_path (str): Path to the JSON file.
        
    Returns:
        config (dict): Loaded configuration dictionary.
    """
    with open(file_path, 'r') as json_file:
        config = json.load(json_file)
    return config


def set_heterogeneity_maps(dynamics, args):
    """Set input maps from args if provided.
    
    Params:
        dynamics (Dynamics): Dynamics object
        args(SimpleNamespace): User arguments
    """
    input_map_dir = f"{cfg.DATA_IN_DIR}/heterogeneity/"
    json_path = os.path.join(input_map_dir, args.heterogeneity)
    heterogeneity_map_args = load_json_config(json_path)
    if heterogeneity_map_args.get("grass_carrying_capacity"):
        m_args = heterogeneity_map_args["grass_carrying_capacity"]
        grass_carcap_dir = os.path.join(input_map_dir, "grass_carrying_capacity")
        if m_args.get("filename"):
            path = os.path.join(grass_carcap_dir, m_args["filename"])
            image = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        else:
            image = spg.generate((dynamics.state.grid.width, dynamics.state.grid.width), **m_args)
            cv2.imwrite(os.path.join(grass_carcap_dir, "generated_grass_carrying_capacity.png"), image)

        image = cv2.resize(image, (dynamics.state.grid.width, dynamics.state.grid.width), interpolation=cv2.INTER_NEAREST)
        image = image / 255  # Normalize to 0-1
        dynamics.state.grid.set_grass_carrying_capacity(image)

    return dynamics, args


def export_state(
        dynamics, path="", init_csv=True, control_variable=None, control_value=None, tree_cover_slope=0,
        extra_parameters="", secondary_variable=None, secondary_value=None, dependent_var=None, dependent_val=None, initial_no_dispersals=None, dependent_result_range_stdev=None,
        args=None
    ):
    fieldnames = [
        "time", "treecover", "slope", "population_size", "#seeds_produced", "fires", "top_kills", "nonseedling_top_kills", "deaths",
        "trees[dbh_0-20%]", "trees[dbh_20-40%]", "trees[dbh_40-60%]", "trees[dbh_60-80%]", "trees[dbh_80-100%]", "extra_parameters",
        "firefree_interval_mean", "firefree_interval_stdev", "firefree_interval_full_sim_mean", "firefree_interval_full_sim_stdev", "time_spent_moving",
        "shaded_out", "outcompeted_by_seedlings", "outcompeted_by_oldstems", "initial_no_dispersals",
        "recruits", "germination_attempts", "oldstem_competition_and_shading", "seedling_competition_and_shading", "perimeter_length", "perimeter-area_ratio"
    ]
    if dependent_var:
        fieldnames.insert(0, dependent_var)
    if secondary_variable:
        if secondary_variable == "treecover":
            secondary_variable = "initial_tree_cover"
        fieldnames.insert(0, secondary_variable)
    if initial_no_dispersals:
        fieldnames.insert(-3, "initial_number_of_dispersals")
    if control_variable:
        if control_variable == "treecover":
            control_variable = "initial tree cover"
        fieldnames.insert(0, control_variable)
    if dependent_result_range_stdev:
        fieldnames.insert(0, "dependent_result_range_stdev")
    if not args.rotate_randomly:
        fieldnames.insert(3, "global_rotation_offset")
    if init_csv and not os.path.exists(path):
        if path == "":
            print("\n\nExport location:", cfg.EXPORT_DIR)
            path = os.path.join(cfg.EXPORT_DIR, "Simulation_" + str(datetime.datetime.now()).replace(":", "-") + ".csv")
        with open(path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
        

    with open(path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        tree_sizes = get_tree_sizes(dynamics)
        firefree_interval_mean, firefree_interval_stdev = get_firefree_interval_stats(dynamics, "current_iteration")
        firefree_interval_fullsim_mean, firefree_interval_fullsim_stdev = get_firefree_interval_stats(dynamics, "average")
        fires = dynamics.get_fires()
        fires = "|".join([str(fire) for fire in fires])
        result = {
            "time": str(dynamics.time),
            "treecover": str(dynamics.state.grid.get_tree_cover()), 
            "slope": str(tree_cover_slope),
            "population_size": str(dynamics.state.population.size()),
            "#seeds_produced": str(dynamics.seeds_produced),
            "fires": fires,
            "trees[dbh_0-20%]": tree_sizes[0],
            "trees[dbh_20-40%]": tree_sizes[1],
            "trees[dbh_40-60%]": tree_sizes[2],
            "trees[dbh_60-80%]": tree_sizes[3],
            "trees[dbh_80-100%]": tree_sizes[4],
            "top_kills": str(dynamics.get_no_fire_induced_topkills()),
            "deaths": str(dynamics.get_no_fire_induced_deaths()),
            "extra_parameters": extra_parameters,
            "firefree_interval_mean": firefree_interval_mean,
            "firefree_interval_stdev": firefree_interval_stdev,
            "firefree_interval_full_sim_mean": firefree_interval_fullsim_mean,
            "firefree_interval_full_sim_stdev": firefree_interval_fullsim_stdev,
            "recruits": str(dynamics.get_no_recruits("all")),
            "initial_no_dispersals": str(dynamics.get_initial_no_dispersals()),
            "time_spent_moving": str(dynamics.get_fraction_time_spent_moving()),
            "shaded_out": str(dynamics.get_fraction_seedlings_dead_due_to_shade()),
            "nonseedling_top_kills": str(dynamics.get_no_fire_induced_nonseedling_topkills()),
            "outcompeted_by_seedlings": str(dynamics.get_fraction_seedlings_outcompeted()),
            "outcompeted_by_oldstems": str(dynamics.get_fraction_seedlings_outcompeted_by_older_trees()),
            "germination_attempts": str(dynamics.get_no_germination_attempts()),
            "oldstem_competition_and_shading": str(dynamics.get_fraction_cases_oldstem_competition_and_shading()),
            "seedling_competition_and_shading": str(dynamics.get_fraction_cases_seedling_competition_and_shading()),
            "perimeter-area_ratio": str(dynamics.get_perimeter_area_ratio()),
            "perimeter_length": str(dynamics.get_forest_perimeter_length())
        }
        if not args.rotate_randomly:
            result["global_rotation_offset"] = str(args.global_rotation_offset)
        if control_variable:
            result[control_variable] = control_value
        if secondary_variable:
            result[secondary_variable] = secondary_value
        if dependent_var:
            result[dependent_var] = dependent_val
        if initial_no_dispersals:
            result["initial_number_of_dispersals"] = initial_no_dispersals
        if dependent_result_range_stdev:
            result["dependent_result_range_stdev"] = dependent_result_range_stdev
        writer.writerow(result)
        
    return path


def get_lookup_table(species, width, rcg_width):
    path = os.path.join(cfg.DATA_INTERNAL_DIR, f"lookup_table_{species}_width-{width}_rcg_width-{rcg_width}.npy")
    if os.path.exists(path):
        return np.load(path), path
    else:
        return None, path

def export_lookup_table(lookup_table, width, species):
    path = os.path.join(cfg.DATA_INTERNAL_DIR, f"lookup_table_{species}_width-{width}_rcg_width-{lookup_table.shape[0]}.npy")
    print("path: ", path)
    np.save(path, lookup_table)

def import_image(fpath):
    img = cv2.imread(fpath)
    return img

def save_numpy_array_to_file(array, path):
    np.save(path, array)

def save_state(dynamics):
    state_table = dynamics.state.get_state_table()
    time = str(dynamics.time).zfill(4)
    save_numpy_array_to_file(state_table, f"{cfg.DATA_OUT_DIR}/tree_positions/tree_positions_iter_{time}.npy")

def load_tree_dbh_values(_time):
    if os.path.exists(cfg.TREE_DBH_FILE):
        if _time == 1:
            os.remove(cfg.TREE_DBH_FILE)
            return {_time : {}}
    else:
        return {_time : {}}
    
     
    with open (cfg.TREE_DBH_FILE, "r") as dbh_json_file:
        tree_dbh_values = json.load(dbh_json_file)
            
    tree_dbh_values[_time] = {}
    
    return tree_dbh_values
    
def save_tree_dbh_values(tree_dbh_values):
    with open(cfg.TREE_DBH_FILE, "w") as dbh_json_file:
        json.dump(tree_dbh_values, dbh_json_file)

def update_state_report(dynamics):
    state_report_file = f"{cfg.DATA_OUT_DIR}/state_reports/state_report.npy"
    current_state = dynamics.state.get_state_table()
    germ_index = 3 
    death_index = germ_index + 1
    tree_dbh_values = load_tree_dbh_values(dynamics.time)

    if os.path.exists(state_report_file) and dynamics.time == 1:
        os.remove(state_report_file)
    elif os.path.exists(state_report_file):
        state_report = np.load(state_report_file)
        
        # Update the state report with newly germinated trees
        for i in range(current_state.shape[0]):
            _id = int(current_state[i, 0])
            if _id not in state_report[:, 0]:
                # Use tree format: id, x, y, time of germination, time of death
                new_tree = list(current_state[i])[:3] # Get the tree id, x and y coordinates
                new_tree = np.concatenate((new_tree, np.array([dynamics.time, 99999])), axis=0) # Add time of germination and placeholder for time of death
                state_report = np.vstack([state_report, new_tree])
            
        # Update dbh file
        for i in range(current_state.shape[0]):
            _id = int(current_state[i, 0])
            cur_dbh = float(current_state[i][3])
            tree_dbh_values[dynamics.time][_id] = cur_dbh
        
        # Update the state report with tree deaths
        for i in range(state_report.shape[0]):
            _id = int(state_report[i, 0])
            if (_id not in current_state[:, 0]) and (int(state_report[i][death_index]) == 99999):
                state_report[i][death_index] = dynamics.time - 1
    
    if dynamics.time == 1:
        state_report = np.zeros((current_state.shape[0], 5))
        state_report[:,:3] = current_state[:,:3]
        state_report[:,-1] = 99999

    save_numpy_array_to_file(state_report, state_report_file)
    save_tree_dbh_values(tree_dbh_values)




