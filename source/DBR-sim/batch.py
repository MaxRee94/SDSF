
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
import statistics as stats


def get_csv_parent_dir(run):
    with open("F:/Development/DBR-sim/data_out/state data/tmp/batchfolder_lookup_table.json", "r") as lookup_table_file:
        folder_lookup_table = json.load(lookup_table_file)

    csv_parent_dir = folder_lookup_table[run]
    if not os.path.exists(csv_parent_dir):
        os.makedirs(csv_parent_dir)
    
    return csv_parent_dir


def get_next_control_value(i, control_variable, control_value, control_range, process_index, no_processes, dynamics, tree_cover_slope):
    control_value = control_range[0] + control_range[2] * (i+1) + process_index * (control_range[2] / no_processes)
    return control_value


def get_singlerun_name(batch_type, sim_name, control_value, csv_parent_dir, secondary_variable=None, secondary_value=None, process_index=None, no_runs_for_current_parameter_set=None):
    if batch_type == "range" or batch_type == "constant":
        singlerun_name = f"{sim_name}={str(control_value)}_process_{str(process_index)}_run_{str(no_runs_for_current_parameter_set)}"
        print(f"\n ------- Beginning simulation with {singlerun_name} ------- \n")
        singlerun_csv_path = csv_parent_dir + "/" + singlerun_name + ".csv"
        singlerun_image_path = singlerun_csv_path.replace(".csv", ".png")
    elif batch_type == "saddle_search":
        singlerun_name = f"{sim_name}={str(control_value)}, {secondary_variable}={str(secondary_value)}, process {str(process_index)}"
        print(f"\n ------- Beginning simulation with {singlerun_name} ------- \n")
        singlerun_csv_path = csv_parent_dir + "/" + singlerun_name + ".csv"
        singlerun_image_path = singlerun_csv_path.replace(".csv", ".png")
    
    return singlerun_name, singlerun_csv_path, singlerun_image_path


def execute_single_run(
        params, control_variable, control_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict,
        extra_parameters, dependent_var, opti_mode, statistic, batch_type, init_csv, secondary_variable=None, secondary_value=None, export_to_parent_csv=True,
        no_runs_for_current_parameter_set=None
    ):
    
    params[control_variable] = control_value
    singlerun_name, singlerun_csv_path, singlerun_image_path = get_singlerun_name(
        batch_type, sim_name, control_value, csv_parent_dir, secondary_variable=secondary_variable, secondary_value=secondary_value, process_index=process_index,
        no_runs_for_current_parameter_set=no_runs_for_current_parameter_set
    )
    params["csv_path"] = singlerun_csv_path
    if "->" in control_variable:
        print("----- batch parameters -----: ", control_variable, control_value, "-----------------")
        params["batch_parameters"] = {"control_variable": control_variable, "control_value": control_value}
    else:
        params[control_variable] = control_value

    # Run the simulation and append its results to the total results csv
    run_starttime = time.time()
    dynamics, tree_cover_slope, largest_absolute_slope, initial_no_recruits = app.main(**params) 
    print(f"\n ------- Simulation {singlerun_name} complete. -------- \n")
    
    return dynamics, tree_cover_slope, largest_absolute_slope, initial_no_recruits, singlerun_name, singlerun_csv_path, singlerun_image_path


            
def iterate_across_range(params, control_variable, control_range, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict,
                        extra_parameters, batch_type, dependent_var, opti_mode, statistic, secondary_variable, secondary_range, attempts=None):
    init_csv = True
    i = 0
    time_budget_per_run = 60 * 60
    no_runs_for_current_parameter_set = 0
    control_value = control_range[0] + process_index * (control_range[2] / no_processes)
    params_json_path = csv_parent_dir + "/parameters.json"
    with open(params_json_path, "w") as params_json_file:
        json.dump(params, params_json_file)
    while (control_value < control_range[1]) or (batch_type == "constant"):       
        run_starttime = time.time()
        if (batch_type == "range" or batch_type == "constant"):
            # Run the simulation and append its results to the total results csv
            dynamics, tree_cover_slope, largest_absolute_slope, initial_no_recruits, singlerun_name, singlerun_csv_path, singlerun_image_path = execute_single_run(
                params, control_variable, control_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict, extra_parameters,
                dependent_var, opti_mode, statistic, batch_type, init_csv, no_runs_for_current_parameter_set=no_runs_for_current_parameter_set
            )
            _io.export_state(
                dynamics, total_results_csv, init_csv, initial_no_recruits=initial_no_recruits, control_variable=control_variable, control_value=control_value,
                tree_cover_slope=tree_cover_slope, extra_parameters=str(extra_parameters)
            )
        elif (batch_type == "saddle_search"):
            secondary_value, dynamics, tree_cover_slope, largest_absolute_slope, guess_result, initial_no_recruits, singlerun_name, singlerun_csv_path, singlerun_image_path = execute_saddle_search(
                params, control_variable, control_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict, extra_parameters,
                dependent_var, opti_mode, statistic, secondary_variable, secondary_range, attempts
            )
            _io.export_state(
                dynamics, total_results_csv, init_csv, control_variable=control_variable, control_value=control_value,
                tree_cover_slope=tree_cover_slope, extra_parameters=str(extra_parameters), secondary_variable=secondary_variable, secondary_value=secondary_value,
                dependent_var=dependent_var, dependent_val=guess_result
            )
        else:
            raise RuntimeError("Invalid batch type '{}'. Exiting..".format(batch_type))
        
        init_csv = False

        # Get a color image representation of the final state
        img = vis.get_image_from_grid(dynamics.state.grid, True, color_dict)
        vis.save_image(img, singlerun_image_path, get_max(dynamics.state.grid.width, 1000))
        
        # Get the next control value
        no_runs_for_current_parameter_set += 1
        if batch_type != "constant" and (batch_type == "saddle_search" or no_runs_for_current_parameter_set >= no_reruns):
            no_runs_for_current_parameter_set = 0
            control_value = get_next_control_value(i, control_variable, control_value, control_range, process_index, no_processes, dynamics, largest_absolute_slope)
        i += no_runs_for_current_parameter_set == 0
 
        # Free memory
        dynamics.free()
        init_csv = False
        
        # If the batch type is constant, we terminate when we've performed a number of simulations equal to 'no_reruns'
        if batch_type == "constant" and no_runs_for_current_parameter_set > no_reruns:
            break
    
    print("Batch complete. Exiting..")


def main(
        process_index=None, control_variable=None, control_range=None, extra_parameters=None, run=0, no_processes=7, no_reruns=None, batch_type=None,
         dependent_var=None, opti_mode=None, statistic=None, secondary_variable=None, secondary_range=None, attempts=None
    ):
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
            
    total_results_csv = csv_parent_dir + "/{}_results.csv".format(csv_parent_dir.split("state data/")[1])
        
    iterate_across_range(params, control_variable, control_range, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, 
                         color_dict, extra_parameters, batch_type, dependent_var, opti_mode, statistic, secondary_variable, secondary_range, attempts=attempts)
    
        

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-pi', '--process_index', type=int)
    parser.add_argument('-cv', '--control_variable', default="constant_batch", type=str)
    parser.add_argument('-cr', '--control_range', type=float, nargs="*", default=[-1, -1, -1], help="Format: min max stepsize")
    parser.add_argument('-np', '--no_processes', type=int, help="Number of processes to run simultaneously.")
    parser.add_argument('-nrr', '--no_reruns', type=int, default=1, help="Number of reruns to perform per parameter set.")
    parser.add_argument('-r', '--run', type=str, default=1, help="Unique run identifier in case of multiple runs per batch")
    parser.add_argument('-bt', '--batch_type', type=str, default="range", help="Type of batch to run. Options: 'range', 'saddle_search', 'constant'")
    parser.add_argument('-dvar', '--dependent_var', type=str, default="tree_cover_slope", help="Dependent variable to optimize")
    parser.add_argument('-om', '--opti_mode', type=str, default="minimize", help="Mode of optimization. Options: 'maximize', 'minimize'")
    parser.add_argument('-st', '--statistic', type=str, default="mean", help="Statistic to optimize. Options: 'max', 'min'")
    parser.add_argument('-secvar', '--secondary_variable', type=str, default="self_ignition_factor", help="Secondary variable of saddle search")
    parser.add_argument('-svr', '--secondary_range', type=float, nargs="*", default=[0, 2], help="Value range of the secondary variable of saddle search")
    parser.add_argument('-na', '--attempts', type=str, default="5", help="Number of attempts (saddle search)")
    parser.add_argument(
        '-ep', '--extra_parameters', type=str,
        help=r"Json string containing key-value pairs of custom parameters. Keys should be surrounded by double quotes preceded by a backslash. Example: {\"verbosity\":1}"
    )
    args = parser.parse_args()
    main(**vars(args))






