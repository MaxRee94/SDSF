import config
import app
from argparse import ArgumentParser
import file_handling as _io
from time import time
import os
import visualization as vis
from helpers import *
import numpy as np


def main(process_index):
    batch_no = 1
    csv_parent_dir = "F:/Development/DBR-sim/data_out/state data/batch_1"
    while os.path.exists(csv_parent_dir):
        batch_no += 1
        csv_parent_dir = csv_parent_dir.split("batch_")[0] + "batch_" + str(batch_no).zfill(6)
    if process_index == 1:
        os.makedirs(csv_parent_dir)
    else:
        csv_parent_dir = csv_parent_dir.split("batch_")[0] + "batch_" + str(batch_no - 1).zfill(6)
    
    no_colors = 100
    color_dict = vis.get_color_dict(no_colors, begin=0.3, end=0.6)
    color_dict[0] = np.array((0, 0, 0), np.uint8)
    control_range = [0.454, 0.47, 0.001]
    control_variable = "treecover"
    params = config.defaults
    params["headless"] = True
    params["max_timesteps"] = 1000
    control_value = control_range[0]
    total_results_csv = csv_parent_dir + "/{}_results.csv".format(csv_parent_dir.split("state data/")[1])
    init_csv = True
    i = 0
    while control_value < control_range[1]:
        print(f"\nBeginning simulation with {control_variable}={str(control_value)} (process idx {str(process_index)})...")
        singlerun_csv_path = csv_parent_dir + f"/{control_variable}={str(control_value)}_process_{str(process_index)}.csv"
        singlerun_image_path = singlerun_csv_path.replace(".csv", ".png")
        params["csv_path"] = singlerun_csv_path
        params[control_variable] = control_value
        print("csv path: ", params["csv_path"])
        print("image path: ", singlerun_image_path)
        dynamics = app.main(**params)
        try:
            _io.export_state(dynamics, total_results_csv, init_csv)
        except:
            time.sleep(3)
        init_csv = False
        control_value = control_range[0] + control_range[2] * i + process_index * (control_range[2] / 8)
        i+=1
        
        # Get a color image representation of the final state
        img = vis.get_image(dynamics.state.grid, True, color_dict)
        vis.save_image(img, singlerun_image_path, get_max(dynamics.state.grid.width, 1000))
        

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-pi', '--process_index', type=int)
    args = parser.parse_args()
    main(args.process_index)






