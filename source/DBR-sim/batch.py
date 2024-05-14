from msilib import Control
import config
import app
from argparse import ArgumentParser
import file_handling as _io
from time import time
import os
import visualization as vis
from helpers import *
import numpy as np
import time


def main(process_index=None, control_variable=None, control_range=None, extra_parameters=None):
    batch_no = 1
    csv_parent_dir = "F:/Development/DBR-sim/data_out/state data/batch_1"
    while os.path.exists(csv_parent_dir):
        batch_no += 1
        csv_parent_dir = csv_parent_dir.split("batch_")[0] + "batch_" + str(batch_no).zfill(6)
    if process_index == 0:
        os.makedirs(csv_parent_dir)
    else:
        csv_parent_dir = csv_parent_dir.split("batch_")[0] + "batch_" + str(batch_no - 1).zfill(6)
    
    no_colors = 100
    color_dict = vis.get_color_dict(no_colors, begin=0.3, end=0.6)
    color_dict[-5] = np.array((0,0,0), np.uint8)
    color_dict[-6] = np.array((0,0,0), np.uint8)
    params = config.defaults
    params["headless"] = True
    params["max_timesteps"] = 1000
    if extra_parameters:
        print("no params:", int(len(extra_parameters) / 2))
        for i in range(int(len(extra_parameters) / 2)):
            key = extra_parameters[i * 2]
            _type = key.split("*")[0]
            key = key.split("*")[1]
            val = extra_parameters[i * 2 + 1]
            if _type == "float":
                val = float(val)
            elif _type == "int":
                val = int(val)
            elif _type == "list":
                val = [float(v) for v in val.split("_")]
            params[key] = val
            print("treecover:", params[key], f"(key = {key})")
            
    control_value = control_range[0] + process_index * (control_range[2] / 7)
    total_results_csv = csv_parent_dir + "/{}_results.csv".format(csv_parent_dir.split("state data/")[1])
    init_csv = True
    i = 0
    time_budget_per_run = 60 * 60
    while control_value < control_range[1]:
        print(f"\n ------- Beginning simulation with {control_variable}={str(control_value)} (process idx {str(process_index)}) ------- \n")
        singlerun_csv_path = csv_parent_dir + f"/{control_variable}={str(control_value)}_process_{str(process_index)}.csv"
        singlerun_image_path = singlerun_csv_path.replace(".csv", ".png")
        params["csv_path"] = singlerun_csv_path
        params[control_variable] = control_value
        print("csv path: ", params["csv_path"])
        print("image path: ", singlerun_image_path)
        
        run_starttime = time.time()
        dynamics = app.main(**params)
            
        try:
            _io.export_state(dynamics, total_results_csv, init_csv, control_variable=control_variable, control_value=control_value)
        except:
            time.sleep(3)
            _io.export_state(dynamics, total_results_csv, init_csv, control_variable=control_variable, control_value=control_value)
        init_csv = False
        control_value = control_range[0] + control_range[2] * i + process_index * (control_range[2] / 7)
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
    parser.add_argument('-cr', '--control_range', type=float, nargs="*")
    parser.add_argument('-ep', '--extra_parameters', type=str, nargs="*")
    args = parser.parse_args()
    main(**vars(args))






