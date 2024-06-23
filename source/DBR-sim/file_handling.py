import csv
import os
from pathlib import Path
import datetime
import cv2
from config import *
import numpy as np



EXPORT_LOCATION = r"F:/Development/DBR-sim/data_out/state data"


def get_tree_sizes(dynamics):
    _tree_sizes = dynamics.state.get_tree_sizes()
    counts, bins = np.histogram(_tree_sizes, bins=5)
    return counts


def export_state(dynamics, path="", init_csv=True, control_variable=None, control_value=None, tree_cover_slope=0, extra_parameters=""):
    fieldnames = [
        "time", "tree cover", "slope", "population size", "#seeds produced", "fire mean spatial extent",
        "#trees[dbh 0-20%]", "#trees[dbh 20-40%]", "#trees[dbh 40-60%]", "#trees[dbh 60-80%]", "#trees[dbh 80-100%]", "extra_parameters"
    ]
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
                "extra_parameters": extra_parameters
            }
        if control_variable:
            result[control_variable] = control_value
        writer.writerow(result)
        
    return path


def import_image(fpath):
    img = cv2.imread(fpath)
    return img

def save_numpy_array_to_file(array, path):
    np.save(path, array)

def save_tree_positions(dynamics):
    tree_positions = dynamics.state.get_tree_positions()
    time = str(dynamics.time).zfill(4)
    save_numpy_array_to_file(tree_positions, f"{DATA_OUT_DIR}/tree_positions/tree_positions_iter_{time}.npy")




