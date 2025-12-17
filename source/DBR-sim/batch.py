import __main__
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
import traceback
import json
import statistics as stats
from types import SimpleNamespace


def get_csv_parent_dir(run):
    with open(os.path.join(config.cfg.DATA_OUT_DIR, "state_data/tmp/batchfolder_lookup_table.json"), "r") as lookup_table_file:
        folder_lookup_table = json.load(lookup_table_file)

    csv_parent_dir = folder_lookup_table[run]
    print("csv parent dir: ", csv_parent_dir)
    if not os.path.isdir(csv_parent_dir):
        os.makedirs(csv_parent_dir)
    
    return csv_parent_dir


def get_next_control_value(i, control_variable, control_value, control_range, process_index, no_processes, dynamics, tree_cover_slope):
    control_value = control_range[0] + control_range[2] * (i+1) + process_index * (control_range[2] / no_processes)
    return control_value


def get_singlerun_name(batch_type, sim_name, control_value, csv_parent_dir, secondary_variable=None, secondary_value=None, process_index=None, no_runs_for_current_parameter_set=None, rerun_idx=None):
    if batch_type == "range" or batch_type == "constant":
        if secondary_variable:
            singlerun_name = f"{sim_name}={str(control_value)}_process_{str(process_index)}_2nd_var({secondary_variable})={str(secondary_value)}_run_{str(rerun_idx+1)}" 
        else:
            singlerun_name = f"{sim_name}={str(control_value)}_process_{str(process_index)}_run_{str(rerun_idx+1)}"
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
        extra_parameters, dependent_var, opti_mode, statistic, batch_type, init_csv, export_to_parent_csv=True,
        secondary_variable=None, secondary_value=None, no_runs_for_current_parameter_set=None, rerun_idx=None
    ):
    
    _params = params.copy()
    _params[control_variable] = control_value
    if secondary_variable:
        _params[secondary_variable] = secondary_value
    print("rerun idx:", rerun_idx)
    singlerun_name, singlerun_csv_path, singlerun_image_path = get_singlerun_name(
        batch_type, sim_name, control_value, csv_parent_dir, secondary_variable=secondary_variable, secondary_value=secondary_value, process_index=process_index,
        no_runs_for_current_parameter_set=no_runs_for_current_parameter_set, rerun_idx=rerun_idx
    )
    _params["csv_path"] = singlerun_csv_path
    if "->" in control_variable:
        print("----- batch parameters -----: ", control_variable, control_value, "-----------------")
        _params["batch_parameters"] = {"control_variable": control_variable, "control_value": control_value}
    else:
        _params[control_variable] = control_value

    # Run the simulation and append its results to the total results csv
    run_starttime = time.time()
    dynamics, tree_cover_slope, largest_absolute_slope, initial_no_dispersals, params = app.main(**_params) 
    print(f"\n ------- Simulation {singlerun_name} complete. -------- \n")
    
    return dynamics, tree_cover_slope, largest_absolute_slope, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, params


def execute_multiple_runs(params, primary_variable, primary_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict, extra_parameters,
        dependent_var, opti_mode, statistic, batch_type, init_csv, secondary_value=None, secondary_variable=None, dependent_result_range_mean=1, dependent_result_range_stdev=1
    ):
    
    params[secondary_variable] = secondary_value
    params["report_state"] = "False"
    dependent_values = []
    for i in range(no_reruns):
        print(f"Beginning run {i+1} of {no_reruns}")
        dynamics, tree_cover_slope, largest_absolute_slope, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, params = execute_single_run(
            params, primary_variable, primary_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict, extra_parameters,
            dependent_var, opti_mode, statistic, batch_type, init_csv, secondary_value=secondary_value, secondary_variable=secondary_variable, export_to_parent_csv=False
        )
        if dependent_var == "tree_cover_slope":
            dependent_values.append(tree_cover_slope)
        else:
            dependent_values.append(getattr(dynamics, dependent_var))
        print("Current {}: ".format(dependent_var), tree_cover_slope)
        if i < (no_reruns - 1):
            dynamics.free()

    print("Dependent values: ", dependent_values)
    result = 0
    if (statistic == "mean"):
        result = abs(stats.mean(dependent_values))
    elif (statistic == "stdev"):
        result = stats.stdev(dependent_values)
    elif (statistic == "stdev_and_mean" and opti_mode == "minimize"):
        result = (abs(stats.mean(dependent_values)) / dependent_result_range_mean) - (stats.stdev(dependent_values) / dependent_result_range_stdev)
    else:
        raise RuntimeError("Invalid statistic '{}'".format(statistic))
    
    if (dependent_var == "tree_cover_slope"):
        tree_cover_slope = stats.mean(dependent_values)

    dependent_result_range_stdev = stats.stdev(dependent_values)

    return dynamics, tree_cover_slope, largest_absolute_slope, result, initial_no_dispersals, stats.stdev(dependent_values), singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, params


def interpolate_secondary_value(argmin, argmax, argmin_result, argmax_result, opti_mode):
    if opti_mode == "maximize":
        return (argmin_result + argmax_result) / 2
    elif opti_mode == "minimize":
        cumulative_deviation = argmin_result + argmax_result
        argmin_deviation = argmin_result / cumulative_deviation
        return argmin + argmin_deviation * (argmax - argmin) # Linear interpolation between argmin and argmax


def hillclimb(latest_arg, latest_result, previous_arg, previous_result, secondary_range, opti_mode, stepsize):
    stepsize = stepsize * (secondary_range[1] - secondary_range[0])
    print("Stepsize: ", stepsize)
    if opti_mode == "minimize":
        if latest_result > previous_result:
            # Result in direction of latest argument is worse, since latest result is larger.
            # Therefore, we should move back in the direction of the previous argument value
            print("Moving toward previous best argument value (latest result = {}, latest arg: {}, previous best result = {}, previous best arg: {})".format(
                    latest_result, latest_arg, previous_result, previous_arg
                )
            )
            diff = previous_arg - latest_arg
        else:
            # The latest argument value is better, since the result is smaller.
            # Therefore, we should continue moving in the same direction
            print("Moving toward latest argument value, which is better than previous (latest result = {}, latest arg: {}, previous best result = {}, previous best arg: {})".format(
                    latest_result, latest_arg, previous_result, previous_arg
                )
            )
            diff = latest_arg - previous_arg
            
    print("Diff: {}, abs(diff): {}".format(diff, abs(diff)))
    direction = diff / abs(diff)
    val = latest_arg + stepsize * direction
    if (val < secondary_range[0]):
        val = (latest_arg + secondary_range[0]) / 2
        print("-- latest arg: ", latest_arg, ", prev arg: ", previous_arg, ", lower bound: ", secondary_range[0], ", avg:", val)      
    elif (val > secondary_range[1]):
        val = (latest_arg + secondary_range[1]) / 2
        print("-- latest arg: ", latest_arg, ", prev arg: ", previous_arg, ", upper bound: ", secondary_range[1], ", avg:", val)
    
    return val

def execute_saddle_search(
        params, primary_variable, primary_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict,
        extra_parameters, dependent_var, opti_mode, statistic, secondary_variable, secondary_range, no_attempts
    ):
    no_attempts = int(no_attempts)
    print("\nExecuting saddle search with primary variable {} = {}".format(primary_variable, primary_value), "and secondary variable {} in {}".format(secondary_variable, secondary_range))
    init_csv = True
    
    # Perform the first two attempts
    argmin = secondary_range[0]
    argmax = secondary_range[1]
    dynamics, tree_cover_slope, largest_absolute_slope, argmin_result, initial_no_dispersals, argmin_results_stdev, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, params = execute_multiple_runs(
        params, primary_variable, primary_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv,
        color_dict, extra_parameters, dependent_var, opti_mode, statistic, "saddle_search", init_csv, secondary_value=argmin, secondary_variable=secondary_variable
    )
    dynamics, tree_cover_slope, largest_absolute_slope, argmax_result, initial_no_dispersals, argmax_results_stdev, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, params = execute_multiple_runs(
        params, primary_variable, primary_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv,
        color_dict, extra_parameters, dependent_var, opti_mode, statistic, "saddle_search", init_csv, secondary_value=argmax, secondary_variable=secondary_variable
    )
    secondary_value_trajectory = {argmin:argmin_result, argmax:argmax_result}
    negative_value_trajectory = {}
    if dependent_var == "tree_cover_slope":
        if argmin_result < 0:
            negative_value_trajectory = {argmin:argmin_result}
        if argmax_result < 0:
            negative_value_trajectory = {argmax:argmax_result}
    print("Initial argmin_result and argmax_result:", argmin_result, argmax_result)
    cur_secondary_value = interpolate_secondary_value(secondary_range[0], secondary_range[1], argmin_result, argmax_result, opti_mode)
    
    # Establish a secondary range-wide mean and stdev, for use in the stdev_and_mean statistic as a normalization factor
    dependent_result_range_mean = (argmin_result + argmax_result) / 2
    dependent_result_range_stdev = (argmin_results_stdev + argmax_results_stdev) / 2
    
    # Perform the remaining attempts using a hillclimbing algorithm
    # -- Initialize 'previous' variables; secondary value with the most optimal result is chosen as the initial 'previous' value
    stepsize = 0.5
    min_and_max_argvalues = np.array([argmin, argmax])
    min_and_max_results = np.array([argmin_result, argmax_result])
    if opti_mode == "minimize":
        prev_result = min_and_max_results[min_and_max_results.argmin()]
        prev_secondary_value = min_and_max_argvalues[min_and_max_results.argmin()]
    else:
        raise RuntimeError("Currently unsupported opti mode '{}'".format(opti_mode))
    
    # -- Run the hillclimbing algorithm
    best_result = prev_result
    best_secondary_value = prev_secondary_value
    best_negative_secondary_value = argmax
    best_negative_secondary_result = argmax_result
    best_positive_secondary_value = argmin
    best_positive_secondary_result = argmin_result
    for i in range(no_attempts - 2):
        print("Value trajectory so far: ", secondary_value_trajectory)
        dynamics, tree_cover_slope, largest_absolute_slope, cur_secondary_value_result, cur_initial_no_dispersals, _, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, params = execute_multiple_runs(
            params, primary_variable, primary_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv,
            color_dict, extra_parameters, dependent_var, opti_mode, statistic, "saddle_search", init_csv, secondary_value=cur_secondary_value, secondary_variable=secondary_variable,
            dependent_result_range_mean=dependent_result_range_mean, dependent_result_range_stdev=dependent_result_range_stdev
        )
        secondary_value_trajectory[cur_secondary_value] = cur_secondary_value_result
        
        # Get smallest and largest secondary values
        secondary_values = np.array([cur_secondary_value, prev_secondary_value])
        prev_and_current_results = np.array([cur_secondary_value_result, prev_result])
        argmax = secondary_values[np.argmax(secondary_values)]
        argmin = secondary_values[np.argmin(secondary_values)]
        
        # Obtain new secondary value
        if opti_mode == "minimize" and dependent_var == "tree_cover_slope":
            if tree_cover_slope > 0:
                # Get average of latest argument and argument corresponding to best negative result
                _secondary_value = (best_negative_secondary_value + cur_secondary_value) / 2
                print("Moving in negative direction of dependent variable")
            else:
                # Get average of latest argument and argument corresponding to best positive result
                _secondary_value = (best_positive_secondary_value + cur_secondary_value) / 2
                print("Moving in positive direction of dependent variable")
        else:
            _secondary_value = hillclimb(cur_secondary_value, cur_secondary_value_result, best_secondary_value, best_result, secondary_range, opti_mode, stepsize)
        prev_secondary_value = cur_secondary_value
        prev_result = cur_secondary_value_result
        stepsize *= 0.8 # Reduce stepsize asymptotically
        
        # Update best result and best secondary value
        if opti_mode == "minimize":
            if cur_secondary_value_result < best_result:
                best_result = cur_secondary_value_result
                best_secondary_value = cur_secondary_value
                initial_no_dispersals = cur_initial_no_dispersals
                if (dependent_var == "tree_cover_slope"):
                    largest_absolute_slope = best_result
                    if tree_cover_slope > 0:
                        best_positive_secondary_value = cur_secondary_value
                        best_positive_secondary_result = cur_secondary_value_result
                    else:
                        best_negative_secondary_value = cur_secondary_value
                        best_negative_secondary_result = cur_secondary_value_result
        else:
            raise RuntimeError("Currently unsupported opti mode '{}'".format(opti_mode))
        
        # Update secondary value
        cur_secondary_value = _secondary_value
        
    
    return best_secondary_value, dynamics, tree_cover_slope, largest_absolute_slope, cur_secondary_value_result, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, params


def get_secondary_value_if_applicable(secondary_variable, secondary_range):
    if secondary_variable is None:
        yield None
        return
    else:
        secondary_control_value = secondary_range[0]
        while secondary_control_value <= secondary_range[2]:
            yield secondary_control_value
            secondary_control_value += secondary_range[1] # New value is old value + step size


def yield_rerun_idx_if_not_saddle_search(batch_type, no_reruns):
    if batch_type == "saddle_search":
        yield None
        return
    else:
        if no_reruns == 0 or no_reruns == 1:
            yield 0
            return
        else:
            for rerun_idx in range(no_reruns):
                yield rerun_idx

            
def iterate_across_range(params, control_variable, control_range, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict,
                        extra_parameters, batch_type, dependent_var, opti_mode, statistic, secondary_variable, secondary_range, attempts=None):
    init_csv = True
    i = 0
    time_budget_per_run = 60 * 60
    no_runs_for_current_parameter_set = 0
    params["report_state"] = "False"
    control_value_minimum = control_range[0] + process_index * (control_range[2] / (no_processes) )
    control_value = control_value_minimum
    params_json_path = csv_parent_dir + "/parameters.json"
    with open(params_json_path, "w") as params_json_file:
        json.dump(params, params_json_file)
        
    for rerun_idx in yield_rerun_idx_if_not_saddle_search(batch_type, no_reruns):
        control_value = control_value_minimum
        print("\n\n\n\n ------- Received rerun index: ", rerun_idx, "for batch type: ", batch_type, "and no_reruns: ", no_reruns)
        while (control_value < control_range[1]) or (batch_type == "constant"):    
            print("\n\n\n\n---new control val:", control_value)
            run_starttime = time.time()
            if (batch_type == "range" or batch_type == "constant"):
                for secondary_value in get_secondary_value_if_applicable(secondary_variable, secondary_range):
                    # Run the simulation and append its results to the total results csv
                    dynamics, tree_cover_slope, largest_absolute_slope, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, _params = execute_single_run(
                        params, control_variable, control_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict, 
                        extra_parameters, dependent_var, opti_mode, statistic, batch_type, init_csv, secondary_variable=secondary_variable, secondary_value=secondary_value, 
                        no_runs_for_current_parameter_set=1, rerun_idx=rerun_idx
                    )
                    params = vars(_params) # Update params in case they were modified during the run

                    _io.export_state(
                        dynamics, total_results_csv, init_csv, initial_no_dispersals=initial_no_dispersals, control_variable=control_variable, control_value=control_value,
                        tree_cover_slope=tree_cover_slope, secondary_variable=secondary_variable, secondary_value=secondary_value, extra_parameters=str(extra_parameters),
                        args=SimpleNamespace(**params)
                    )
            elif (batch_type == "saddle_search"):
                secondary_value, dynamics, tree_cover_slope, largest_absolute_slope, guess_result, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, params = execute_saddle_search(
                    params, control_variable, control_value, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, color_dict, extra_parameters,
                    dependent_var, opti_mode, statistic, secondary_variable, secondary_range, attempts
                )
                _io.export_state(
                    dynamics, total_results_csv, init_csv, control_variable=control_variable, control_value=control_value,
                    tree_cover_slope=tree_cover_slope, extra_parameters=str(extra_parameters), secondary_variable=secondary_variable, secondary_value=secondary_value,
                    dependent_var=dependent_var, dependent_val=guess_result, dependent_result_range_stdev=dependent_result_range_stdev
                )
            else:
                raise RuntimeError("Invalid batch type '{}'. Exiting..".format(batch_type))
        
            init_csv = False

            # Get a color image representation of the final state
            img = vis.get_image_from_grid(dynamics.state.grid, color_dict, collect_states=1)
            vis.save_image(img, singlerun_image_path, get_max(dynamics.state.grid.width, 1000))
        
            # Get the next control value
            no_runs_for_current_parameter_set += 1
            if batch_type != "constant" and (batch_type in ["saddle_search", "range"]):
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

def _ensure_correct_arg_datatype(element):
    if not type(element) == str:
        return element
    element = element.strip()
    if element == "True":
        result = True
    elif element == "False":
        element = False
    elif element == "None":
        element = None
    elif is_number(element):
        if "." in element:
            element = float(element)
        else:
            element = int(element)
    elif "[" in element:
        stringlist = [x.strip() for x in element.split(",")]
        element = [_ensure_correct_arg_datatype(x) for x in stringlist]
    
    return element


def ensure_correct_arg_datatype(args):
    """Ensure that the arguments are of the correct datatype.

    Args:
        args (dict): Dictionary containing the arguments and their values.

    Returns:
        dict: Dictionary with the arguments and their values, with the correct datatypes.
    """
    print("args: ", args)
    for key in args:
        args[key] = _ensure_correct_arg_datatype(args[key])
    
    return args
 

def main(
        process_index=None, control_variable=None, control_range=None, extra_parameters=None, run=0, no_processes=7, no_reruns=None, batch_type=None,
         dependent_var=None, opti_mode=None, statistic=None, secondary_variable=None, secondary_range=None, attempts=None, **params
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
    params["headless"] = True
    if extra_parameters:
        print("Extra parameters:", extra_parameters) # Example of usage: {\"verbosity\":1}
        extra_parameters = ensure_correct_arg_datatype(extra_parameters)
        for key, val in extra_parameters.items():
            params[key] = val
            print("extra parameter:", val, f"(key = {key})")
            
    total_results_csv = csv_parent_dir + "/{}_results.csv".format(csv_parent_dir.split("state_data/")[1])
        
    iterate_across_range(params, control_variable, control_range, csv_parent_dir, process_index, no_processes, no_reruns, sim_name, total_results_csv, 
                         color_dict, extra_parameters, batch_type, dependent_var, opti_mode, statistic, secondary_variable, secondary_range, 
                         attempts=attempts)
    
        
if __name__ == "__main__":
    print(sys.argv)


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
    parser.add_argument('-secvar', '--secondary_variable', type=str, default=None, help="Secondary variable")
    parser.add_argument('-svr', '--secondary_range', type=float, nargs="*", default=[0, 2, 2], help="Value range of the secondary variable")
    parser.add_argument('-na', '--attempts', type=str, default="5", help="Number of attempts (saddle search)")
    parser.add_argument(
        '-ep', '--extra_parameters', type=str,
        help=r"Json string containing key-value pairs of custom parameters. Keys should be surrounded by double quotes preceded by a backslash. Example: {\"verbosity\":1}"
    )
    args = SimpleNamespace(**parse_args(parser))
    try:
        main(**vars(args))
    except Exception as e:
        print("======== ERROR =========")
        parent_dir = get_csv_parent_dir(args.run)
        csv_file = parent_dir + "/Error_log.txt";
        with open(csv_file, 'a') as f:
            f.write(str(e))
            f.write(traceback.format_exc())
            
        traceback.print_exc()
        time.sleep(10)

