import config
import app
from argparse import ArgumentParser
import file_handling as _io
from time import time
import os
import visualization as vis
from helpers import *
import numpy as np
import random
import time
import json


def get_csv_parent_dir(run):
    with open("F:/Development/DBR-sim/data_out/state data/tmp/batchfolder_lookup_table.json", "r") as lookup_table_file:
        folder_lookup_table = json.load(lookup_table_file)

    csv_parent_dir = folder_lookup_table[run]
    if not os.path.exists(csv_parent_dir):
        os.makedirs(csv_parent_dir)
    
    return csv_parent_dir


def main(process_index=None, control_variable=None, control_range=None, extra_parameters=None, run=0, no_processes=7):
    csv_parent_dir = get_csv_parent_dir(run)
    print("csv parent dir: ", csv_parent_dir)
    
    no_colors = 100
    if "->" in control_variable:
        sim_name = control_variable.split("->")[-2]
    else:
        sim_name = control_variable
    color_dict = vis.get_color_dict(no_colors, begin=0.2, end=0.5)
    color_dict[0] = np.array((170, 255, 255), np.uint8)
    params = config.defaults
    params["headless"] = True
    if extra_parameters:
        print("Extra parameters:", extra_parameters) # Example of usage: {\"verbosity\":1}
        extra_parameters = json.loads(extra_parameters)
        for key, val in extra_parameters.items():
            params[key] = val
            print("extra parameter:", val, f"(key = {key})")
            
    control_value = control_range[0] + process_index * (control_range[2] / 7)
    total_results_csv = csv_parent_dir + "/{}_results.csv".format(csv_parent_dir.split("state data/")[1])
    init_csv = True
    i = 0
    time_budget_per_run = 60 * 60
    while control_value < control_range[1]:
        singlerun_name = f"/{sim_name}={str(control_value)}_process_{str(process_index)}"
        print(f"\n ------- Beginning simulation with {singlerun_name} ------- \n")
        singlerun_csv_path = csv_parent_dir + singlerun_name + ".csv"
        singlerun_image_path = singlerun_csv_path.replace(".csv", ".png")
        params["csv_path"] = singlerun_csv_path
        if "->" in control_variable:
            params["batch_parameters"] = {"control_variable": control_variable, "control_value": control_value}
        else:
            params[control_variable] = control_value
        
        # Run the simulation and append its results to the total results csv
        run_starttime = time.time()
        dynamics, color_dicts = app.init(**params) 
        dynamics, tree_cover_slope = app.updateloop(dynamics, color_dicts, **params)

        _io.export_state(
            dynamics, total_results_csv, init_csv, control_variable=control_variable, control_value=control_value,
            tree_cover_slope=tree_cover_slope, extra_parameters=str(extra_parameters)
        )
        init_csv = False
        
        # Get the next control value
        control_value = control_range[0] + control_range[2] * (i+1) + process_index * (control_range[2] / no_processes)
        i+=1
        
        # Get a color image representation of the final state
        img = vis.get_image_from_grid(dynamics.state.grid, True, color_dict)
        vis.save_image(img, singlerun_image_path, get_max(dynamics.state.grid.width, 1000))
        
        # Free memory
        dynamics.free()
        print(f"\n ------- Simulation with {control_variable}={str(control_value)} (process idx {str(process_index)}) complete. -------- \n")
        

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-pi', '--process_index', type=int)
    parser.add_argument('-cv', '--control_variable', type=str)
    parser.add_argument('-cr', '--control_range', type=float, nargs="*", help="Format: min max stepsize")
    parser.add_argument('-np', '--no_processes', type=int, help="Number of processes to run simultaneously.")
    parser.add_argument('-r', '--run', type=str, default=1, help="Unique run identifier in case of multiple runs per batch")
    parser.add_argument(
        '-ep', '--extra_parameters', type=str,
        help=r"Json string containing key-value pairs of custom parameters. Keys should be surrounded by double quotes preceded by a backslash. Example: {\"verbosity\":1}"
    )
    args = parser.parse_args()
    main(**vars(args))






