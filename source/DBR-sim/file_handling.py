from cgitb import lookup
import csv
import os
from pathlib import Path
import datetime
import cv2
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
        extra_parameters="", secondary_variable=None, secondary_value=None, dependent_var=None, dependent_val=None, initial_no_recruits=None
    ):
    fieldnames = [
        "time", "tree cover", "slope", "population size", "#seeds produced", "fire mean spatial extent",
        "#trees[dbh 0-20%]", "#trees[dbh 20-40%]", "#trees[dbh 40-60%]", "#trees[dbh 60-80%]", "#trees[dbh 80-100%]", "extra_parameters",
        "firefree interval mean", "firefree interval stdev", "firefree interval full sim mean", "firefree interval full sim stdev", "time_spent_moving",
        "shaded_out", "outcompeted_by_seedlings", "outcompeted_by_oldstems",
        "recruits", "germination_attempts", "oldstem_competition_and_shading", "seedling_competition_and_shading"
    ]
    if dependent_var:
        fieldnames.insert(0, dependent_var)
    if secondary_variable:
        if secondary_variable == "treecover":
            secondary_variable = "initial tree cover"
        fieldnames.insert(0, secondary_variable)
    if initial_no_recruits:
        fieldnames.insert(-3, "initial_no_recruits")
    if control_variable:
        if control_variable == "treecover":
            control_variable = "initial tree cover"
        fieldnames.insert(0, control_variable)
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
            "extra_parameters": extra_parameters,
            "firefree interval mean": firefree_interval_mean,
            "firefree interval stdev": firefree_interval_stdev,
            "firefree interval full sim mean": firefree_interval_fullsim_mean,
            "firefree interval full sim stdev": firefree_interval_fullsim_stdev,
            "recruits": str(dynamics.get_no_recruits("all")),
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
        if initial_no_recruits:
            result["initial_no_recruits"] = initial_no_recruits
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

def save_tree_positions(dynamics):
    tree_positions = dynamics.state.get_tree_positions()
    time = str(dynamics.time).zfill(4)
    save_numpy_array_to_file(tree_positions, f"{DATA_OUT_DIR}/tree_positions/tree_positions_iter_{time}.npy")




