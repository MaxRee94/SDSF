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


TREE_DBH_FILE = f"{DATA_OUT_DIR}/state_reports/tree_dbh_values.json"
EXPORT_LOCATION = r"F:/Development/DBR-sim/data_out/state data"


def get_tree_sizes(dynamics):
    _tree_sizes = dynamics.state.get_tree_sizes()
    counts, bins = np.histogram(_tree_sizes, bins=5)
    return counts


def get_firefree_interval_stats(dynamics, _type):
    intervals = dynamics.get_firefree_intervals(_type)
    stdev = np.std(intervals)
    mean = np.mean(intervals)
    return mean, stdev


def export_state(
        dynamics, path="", init_csv=True, control_variable=None, control_value=None, tree_cover_slope=0,
        extra_parameters="", secondary_variable=None, secondary_value=None, dependent_var=None, dependent_val=None, initial_no_dispersals=None, dependent_result_range_stdev=None
    ):
    fieldnames = [
        "time", "tree cover", "slope", "population size", "#seeds produced", "fire mean spatial extent", "#top kills", "#deaths",
        "#trees[dbh 0-20%]", "#trees[dbh 20-40%]", "#trees[dbh 40-60%]", "#trees[dbh 60-80%]", "#trees[dbh 80-100%]", "extra_parameters",
        "firefree interval mean", "firefree interval stdev", "firefree interval full sim mean", "firefree interval full sim stdev", "time_spent_moving",
        "shaded_out", "outcompeted_by_seedlings", "outcompeted_by_oldstems", "initial_no_dispersals",
        "recruits", "germination_attempts", "oldstem_competition_and_shading", "seedling_competition_and_shading"
    ]
    if dependent_var:
        fieldnames.insert(0, dependent_var)
    if secondary_variable:
        if secondary_variable == "treecover":
            secondary_variable = "initial tree cover"
        fieldnames.insert(0, secondary_variable)
    if initial_no_dispersals:
        fieldnames.insert(-3, "initial_no_dispersals")
    if control_variable:
        if control_variable == "treecover":
            control_variable = "initial tree cover"
        fieldnames.insert(0, control_variable)
    if dependent_result_range_stdev:
        fieldnames.insert(0, "dependent_result_range_stdev")
    if init_csv and not os.path.exists(path):
        if path == "":
            path = os.path.join(EXPORT_LOCATION, "Simulation_" + str(datetime.datetime.now()).replace(":", "-") + ".csv")
        with open(path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

    with open(path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        tree_sizes = get_tree_sizes(dynamics)
        firefree_interval_mean, firefree_interval_stdev = get_firefree_interval_stats(dynamics, "current_iteration")
        firefree_interval_fullsim_mean, firefree_interval_fullsim_stdev = get_firefree_interval_stats(dynamics, "average")
        result = {
            "time": str(dynamics.time),
            "tree cover": str(dynamics.state.grid.get_tree_cover()), 
            "slope": str(tree_cover_slope),
            "population size": str(dynamics.state.population.size()),
            "#seeds produced": str(dynamics.seeds_produced),
            "fire mean spatial extent": str(dynamics.fire_spatial_extent),
            "#trees[dbh 0-20%]": tree_sizes[0],
            "#trees[dbh 20-40%]": tree_sizes[1],
            "#trees[dbh 40-60%]": tree_sizes[2],
            "#trees[dbh 60-80%]": tree_sizes[3],
            "#trees[dbh 80-100%]": tree_sizes[4],
            "#top kills": str(dynamics.get_no_fire_induced_topkills()),
            "#deaths": str(dynamics.get_no_fire_induced_deaths()),
            "extra_parameters": extra_parameters,
            "firefree interval mean": firefree_interval_mean,
            "firefree interval stdev": firefree_interval_stdev,
            "firefree interval full sim mean": firefree_interval_fullsim_mean,
            "firefree interval full sim stdev": firefree_interval_fullsim_stdev,
            "recruits": str(dynamics.get_no_recruits("all")),
            "initial_no_dispersals": str(dynamics.get_initial_no_dispersals()),
            "time_spent_moving": str(dynamics.get_fraction_time_spent_moving()),
            "shaded_out": str(dynamics.get_fraction_seedlings_dead_due_to_shade()),
            "outcompeted_by_seedlings": str(dynamics.get_fraction_seedlings_outcompeted()),
            "outcompeted_by_oldstems": str(dynamics.get_fraction_seedlings_outcompeted_by_older_trees()),
            "germination_attempts": str(dynamics.get_no_germination_attempts()),
            "oldstem_competition_and_shading": str(dynamics.get_fraction_cases_oldstem_competition_and_shading()),
            "seedling_competition_and_shading": str(dynamics.get_fraction_cases_seedling_competition_and_shading())
        }
        if control_variable:
            result[control_variable] = control_value
        if secondary_variable:
            result[secondary_variable] = secondary_value
        if dependent_var:
            result[dependent_var] = dependent_val
        if initial_no_dispersals:
            result["initial number of dispersals"] = initial_no_dispersals
        if dependent_result_range_stdev:
            result["dependent_result_range_stdev"] = dependent_result_range_stdev
        writer.writerow(result)
        
    return path


def get_lookup_table(species, width, rcg_width):
    path = os.path.join(DATA_INTERNAL_DIR, f"lookup_table_{species}_width-{width}_rcg_width-{rcg_width}.npy")
    if os.path.exists(path):
        return np.load(path), path
    else:
        return None, path

def export_lookup_table(lookup_table, width, species):
    path = os.path.join(DATA_INTERNAL_DIR, f"lookup_table_{species}_width-{width}_rcg_width-{lookup_table.shape[0]}.npy")
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
    save_numpy_array_to_file(state_table, f"{DATA_OUT_DIR}/tree_positions/tree_positions_iter_{time}.npy")

def load_tree_dbh_values(_time):
    if os.path.exists(TREE_DBH_FILE):
        if _time == 1:
            os.remove(TREE_DBH_FILE)
            return {_time : {}}
    else:
        return {_time : {}}
    
     
    with open (TREE_DBH_FILE, "r") as dbh_json_file:
        tree_dbh_values = json.load(dbh_json_file)
            
    tree_dbh_values[_time] = {}
    
    return tree_dbh_values
    
def save_tree_dbh_values(tree_dbh_values):
    with open(TREE_DBH_FILE, "w") as dbh_json_file:
        json.dump(tree_dbh_values, dbh_json_file)

def update_state_report(dynamics):
    state_report_file = f"{DATA_OUT_DIR}/state_reports/state_report.npy"
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




