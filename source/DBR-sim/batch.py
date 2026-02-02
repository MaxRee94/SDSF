import __main__
from operator import ge
import config
import app
from argparse import ArgumentParser
from multiprocessing import Process, Value
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
import logging
import statistics as stats
from types import SimpleNamespace
from copy import deepcopy

logger = logging.getLogger(__name__)


class Jobs:
    def __init__(self, reruns=None, arguments=None, rng=None, **args):
        self.defaults = config.get_all_defaults()
        self.default_values = list(self.defaults.values())
        self.n_reruns = reruns
        self.no_finished_rerun_cycles = 0
        self.rng = rng
        self.sampling_mode = "regular" # alternative: "random"
        self.jobs, self.job_indices = self.parse(arguments)

    def get_key_idx(self, key):
        return list(self.defaults.keys()).index(key)
    
    @staticmethod
    def get_vec(arg_cfg):
        if arg_cfg["interpolation"] == "linear":
            minim, maxim, stepsize = arg_cfg["range"] + [arg_cfg["stepsize"]]
            no_steps = int((maxim - minim) / stepsize)
            vec = list(np.arange(minim, maxim + stepsize/2, stepsize))
            no_digits_after_comma = max(digits_after_decimal(minim), digits_after_decimal(maxim), digits_after_decimal(stepsize))
            vec = [float(round(v, no_digits_after_comma)) for v in vec]
        elif arg_cfg["interpolation"] == "categorical":
            vec = arg_cfg["range"]

        return vec

    def convert_value_set_to_namespace(self, job_values):
        job = self.convert_job_to_dict(job_values)

        return SimpleNamespace(**job)

    def convert_job_to_dict(self, job_values):
        job = {}
        for i, (key, _) in enumerate(self.defaults.items()):
            job[key] = job_values[i]

        return job
    
    def derive_unique_job_count(self, arg_changes):
        job_count = 0
        for key, arg_cfg in arg_changes.items():
            vec = self.get_vec(arg_cfg)
            if job_count == 0:
                job_count = len(vec)
            else:
                job_count *= len(vec)

        return job_count

    def switch_to_random_sampling(self, job_count, generator_cfg):
        logger.warning(
            f"Number of jobs ({job_count}) exceeds 1 million. \n" +
            "I will therefore not precompute all job arguments. " +
            "Instead, I will switch to random parameter sampling."
        )
        self.sampling_mode = "random"
        jobs = []
        generator_cfg.rng = self.rng
        generator_cfg.low = 0
        generator_cfg.high = job_count - 1
        job_idx_generator = get_random_int_generator
        
        return jobs, job_idx_generator

    def parse(self, arg_changes):
        job_count = self.derive_unique_job_count(arg_changes)
        logger.info("Expecting {} unique jobs.".format(job_count))
        generator_cfg = SimpleNamespace()
        if job_count < 1e6:
            jobs = self.parse_arg_values(arg_changes)
            generator_cfg.n = len(jobs)
            job_idx_generator = midpoint_gap_indices
        else:
            jobs, job_idx_generator = self.switch_to_random_sampling(job_count, generator_cfg)

        job_indices = list(job_idx_generator(generator_cfg))

        return jobs, job_indices

    def parse_arg_values(self, arg_changes):
        value_sets = self.default_values.copy()
        for key, arg_cfg in arg_changes.items():
            vec = self.get_vec(arg_cfg)
            idx = self.get_key_idx(key)
            value_sets = self.add_range(value_sets, idx, vec)

        # If we're doing a multi-dimensional sensitivity analysis, we shuffle the order of the jobs
        if not type(value_sets[0]) == list:
            self.rng.shuffle(value_sets)

        return value_sets

    def get_specific_job(self, idx, return_dict=False):
        if return_dict:
            return self.convert_job_to_dict(self.jobs[idx])
        else:
            return self.convert_value_set_to_namespace(self.jobs[idx])

    def count(self):
        return self.n_reruns * len(self.jobs)

    def get(self, n_finished_simulations):
        if n_finished_simulations >= self.count():
            return None

        idx = self.job_indices[n_finished_simulations % len(self.jobs)]
        job = self.get_specific_job(idx)

        return job

    @staticmethod
    def apply_single_arg_change(argsets, is_single_argset, idx, value):
        if is_single_argset:
            argsets[idx] = value
        else:
            for argset in argsets:
                argset[idx] = value

        return argsets

    def add_range(self, argsets, idx, vec):
        """Add a range of values (i.e., a vector) to the argument sets at position idx.
        
        vec (list): Range of values to add.
        idx (int): Index into each argument list (='argset') of the argument to modify.
        argsets (list): list of argument lists (='argset'), or a single, flat argument list if no arg changes have yet been made.
        """

        is_single_argset = type(argsets[0]) != list
        if is_single_argset:
            expanded_argsets = [argsets.copy() for i in range(len(vec))]
        else:
            expanded_argsets = []
            for i in range(len(vec)):
                expanded_argsets += deepcopy(argsets)

        stepsize = len(argsets)
        for i, value in enumerate(vec):
            argsets = self.apply_single_arg_change(deepcopy(argsets), is_single_argset, idx, value)
            begin = i * stepsize
            end = begin + stepsize
            if is_single_argset:
                _argsets = [argsets]
            else:
                _argsets = argsets
                
            expanded_argsets[begin:end] = _argsets

        return expanded_argsets


def get_new_batch_folder():
    batch_no = 1
    batch_folder = config.cfg.DATA_OUT_DIR + "/state_data/batch_000001"
    new_batch_folders = []
    while os.path.exists(batch_folder):
        batch_no += 1
        batch_folder = batch_folder.split("batch_")[0] + "batch_" + str(batch_no).zfill(6)

    return batch_folder


def determine_csv_parent_dir(batch_cfg):
    batch_folder = get_new_batch_folder()
    if not os.path.isdir(batch_folder):
        os.makedirs(batch_folder)
    batch_cfg.csv_parent_dir = batch_folder    

    return batch_cfg


def get_next_control_value(i, control_variable, control_value, control_range, proc_id, no_processes, dynamics, tree_cover_slope):
    control_value = control_range[0] + control_range[2] * (i+1) + proc_id * (control_range[2] / no_processes)
    return control_value


def get_singlerun_name(type, sim_name, control_value, csv_parent_dir, secondary_variable=None, secondary_value=None, proc_id=None, no_runs_for_current_parameter_set=None, rerun_idx=None):
    if type == "range" or type == "constant":
        if secondary_variable:
            singlerun_name = f"{sim_name}={str(control_value)}_process_{str(proc_id)}_2nd_var({secondary_variable})={str(secondary_value)}_run_{str(rerun_idx+1)}" 
        else:
            singlerun_name = f"{sim_name}={str(control_value)}_process_{str(proc_id)}_run_{str(rerun_idx+1)}"
        print(f"\n ------- Beginning simulation with {singlerun_name} ------- \n")
        singlerun_csv_path = csv_parent_dir + "/" + singlerun_name + ".csv"
        singlerun_image_path = singlerun_csv_path.replace(".csv", ".png")
    elif type == "saddle_search":
        singlerun_name = f"{sim_name}={str(control_value)}, {secondary_variable}={str(secondary_value)}, process {str(proc_id)}"
        print(f"\n ------- Beginning simulation with {singlerun_name} ------- \n")
        singlerun_csv_path = csv_parent_dir + "/" + singlerun_name + ".csv"
        singlerun_image_path = singlerun_csv_path.replace(".csv", ".png")
    
    return singlerun_name, singlerun_csv_path, singlerun_image_path


def execute_single_run(
        params, control_variable, control_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv, color_dict,
        extra_parameters, dependent_var, opti_mode, statistic, type, init_csv, export_to_parent_csv=True,
        secondary_variable=None, secondary_value=None, no_runs_for_current_parameter_set=None, rerun_idx=None
    ):
    
    _params = params.copy()
    _params[control_variable] = control_value
    if secondary_variable:
        _params[secondary_variable] = secondary_value
    print("rerun idx:", rerun_idx)
    singlerun_name, singlerun_csv_path, singlerun_image_path = get_singlerun_name(
        type, sim_name, control_value, csv_parent_dir, secondary_variable=secondary_variable, secondary_value=secondary_value, proc_id=proc_id,
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


def execute_multiple_runs(params, primary_variable, primary_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv, color_dict, extra_parameters,
        dependent_var, opti_mode, statistic, type, init_csv, secondary_value=None, secondary_variable=None, dependent_result_range_mean=1, dependent_result_range_stdev=1
    ):
    
    params[secondary_variable] = secondary_value
    params["report_state"] = "False"
    _params = SimpleNamespace(**params)
    dependent_values = []
    for i in range(n_reruns):
        print(f"Beginning run {i+1} of {n_reruns}")
        dynamics, tree_cover_slope, largest_absolute_slope, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, _params = execute_single_run(
            params, primary_variable, primary_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv, color_dict, extra_parameters,
            dependent_var, opti_mode, statistic, type, init_csv, secondary_value=secondary_value, secondary_variable=secondary_variable, export_to_parent_csv=False
        )
        if dependent_var == "tree_cover_slope":
            dependent_values.append(tree_cover_slope)
        else:
            dependent_values.append(getattr(dynamics, dependent_var))
        print("Current {}: ".format(dependent_var), tree_cover_slope)
        if i < (n_reruns - 1):
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

    return dynamics, tree_cover_slope, largest_absolute_slope, result, initial_no_dispersals, stats.stdev(dependent_values), singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, _params


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
        params, primary_variable, primary_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv, color_dict,
        extra_parameters, dependent_var, opti_mode, statistic, secondary_variable, secondary_range, no_attempts
    ):
    no_attempts = int(no_attempts)
    print("\nExecuting saddle search with primary variable {} = {}".format(primary_variable, primary_value), "and secondary variable {} in {}".format(secondary_variable, secondary_range))
    init_csv = True
    
    # Perform the first two attempts
    argmin = secondary_range[0]
    argmax = secondary_range[1]
    dynamics, tree_cover_slope, largest_absolute_slope, argmin_result, initial_no_dispersals, argmin_results_stdev, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, _params = execute_multiple_runs(
        params, primary_variable, primary_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv,
        color_dict, extra_parameters, dependent_var, opti_mode, statistic, "saddle_search", init_csv, secondary_value=argmin, secondary_variable=secondary_variable
    )
    dynamics, tree_cover_slope, largest_absolute_slope, argmax_result, initial_no_dispersals, argmax_results_stdev, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, _params = execute_multiple_runs(
        params, primary_variable, primary_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv,
        color_dict, extra_parameters, dependent_var, opti_mode, statistic, "saddle_search", init_csv, secondary_value=argmax, secondary_variable=secondary_variable
    )
    secondary_value_trajectory = {argmin:argmin_result, argmax:argmax_result}
    assert argmin_result != argmax_result, f"Mean of dependent variable is equal (specifically: {argmin_result}) for both extremes of the given secondary parameter range ({argmin}, {argmax}). Check whether parameter range is realistic and max_timesteps (currently {params['max_timesteps']}) is sufficiently long."
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
        dynamics, tree_cover_slope, largest_absolute_slope, cur_secondary_value_result, cur_initial_no_dispersals, _, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, _params = execute_multiple_runs(
            params, primary_variable, primary_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv,
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
        
    
    return best_secondary_value, dynamics, tree_cover_slope, largest_absolute_slope, cur_secondary_value_result, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, _params


def get_secondary_value_if_applicable(secondary_variable, secondary_range):
    if secondary_variable is None:
        yield None
        return
    else:
        secondary_control_value = secondary_range[0]
        while secondary_control_value <= secondary_range[2]:
            yield secondary_control_value
            secondary_control_value += secondary_range[1] # New value is old value + step size


def yield_rerun_idx_if_not_saddle_search(type, n_reruns):
    if type == "saddle_search":
        yield None
        return
    else:
        if n_reruns == 0 or n_reruns == 1:
            yield 0
            return
        else:
            for rerun_idx in range(n_reruns):
                yield rerun_idx

            
def iterate_across_range(params, control_variable, control_range, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv, color_dict,
                        extra_parameters, type, dependent_var, opti_mode, statistic, secondary_variable, secondary_range, attempts=None):
    init_csv = True
    i = 0
    time_budget_per_run = 60 * 60
    no_runs_for_current_parameter_set = 0
    params["report_state"] = "False"
    control_value_minimum = control_range[0] + proc_id * (control_range[2] / (no_processes) )
    control_value = control_value_minimum
    
        
    for rerun_idx in yield_rerun_idx_if_not_saddle_search(type, n_reruns):
        control_value = control_value_minimum
        print("\n\n\n\n ------- Received rerun index: ", rerun_idx, "for batch type: ", type, "and n_reruns: ", n_reruns)
        while (control_value < control_range[1]) or (type == "constant"):    
            print("\n\n\n\n---new control val:", control_value)
            run_starttime = time.time()
            if (type == "range" or type == "constant"):
                for secondary_value in get_secondary_value_if_applicable(secondary_variable, secondary_range):
                    # Run the simulation and append its results to the total results csv
                    dynamics, tree_cover_slope, largest_absolute_slope, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, _params = execute_single_run(
                        params, control_variable, control_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv, color_dict, 
                        extra_parameters, dependent_var, opti_mode, statistic, type, init_csv, secondary_variable=secondary_variable, secondary_value=secondary_value, 
                        no_runs_for_current_parameter_set=1, rerun_idx=rerun_idx
                    )
                    params = vars(_params) # Update params in case they were modified during the run

                    _io.export_state(
                        dynamics, total_results_csv, init_csv, initial_no_dispersals=initial_no_dispersals, control_variable=control_variable, control_value=control_value,
                        tree_cover_slope=tree_cover_slope, secondary_variable=secondary_variable, secondary_value=secondary_value, extra_parameters=str(extra_parameters),
                        cfg=SimpleNamespace(**params)
                    )
            elif (type == "saddle_search"):
                secondary_value, dynamics, tree_cover_slope, largest_absolute_slope, guess_result, initial_no_dispersals, singlerun_name, singlerun_csv_path, singlerun_image_path, dependent_result_range_stdev, _params = execute_saddle_search(
                    params, control_variable, control_value, csv_parent_dir, proc_id, no_processes, n_reruns, sim_name, total_results_csv, color_dict, extra_parameters,
                    dependent_var, opti_mode, statistic, secondary_variable, secondary_range, attempts
                )
                params = vars(_params) # Update params in case they were modified during the run
                
                _io.export_state(
                    dynamics, total_results_csv, init_csv, control_variable=control_variable, control_value=control_value,
                    tree_cover_slope=tree_cover_slope, extra_parameters=str(extra_parameters), secondary_variable=secondary_variable, secondary_value=secondary_value,
                    dependent_var=dependent_var, dependent_val=guess_result, dependent_result_range_stdev=dependent_result_range_stdev, cfg=SimpleNamespace(**params)
                )
            else:
                raise RuntimeError("Invalid batch type '{}'. Exiting..".format(type))

            # Get a color image representation of the final state
            img = vis.get_image_from_grid(dynamics.state.grid, color_dict, collect_states=1)
            vis.save_image(img, singlerun_image_path, get_max(dynamics.state.grid.width, 1000))

            init_csv = False
    
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


def ensure_correct_arg_datatype(cfg):
    """Ensure that the arguments are of the correct datatype.

    cfg:
        cfg (dict): Dictionary containing the arguments and their values.

    Returns:
        dict: Dictionary with the arguments and their values, with the correct datatypes.
    """

    for key in cfg:
        cfg[key] = _ensure_correct_arg_datatype(cfg[key])
    
    return cfg
 

def read_config_from_file(_batch_cfg):
    with open(os.path.join(config.cfg.SOURCE_DIR, "batch_config", _batch_cfg.config), "r") as f:
        batch_cfg = json.load(f)     

    for k, v in vars(_batch_cfg).items():
        batch_cfg[k] = v

    return SimpleNamespace(**batch_cfg)


def create_color_dict(batch_cfg):
    no_colors = 100
    color_dict = batch_cfg.vis.get_color_dict(no_colors, begin=0.2, end=0.5)
    color_dict[0] = np.array((170, 255, 255), np.uint8)

    return color_dict


def set_random_seeds(batch_cfg):
    # Set random seed
    if (batch_cfg.random_seed == -999):
        batch_cfg.random_seed = random.randint(0, 100000000)
        batch_cfg.rng = np.random.default_rng(batch_cfg.random_seed)
        logger.info(f"Generated new global random seed ({batch_cfg.random_seed})")
    else:
        batch_cfg.rng = np.random.default_rng(batch_cfg.random_seed)
        logger.info(f"Using given global random seed ({batch_cfg.random_seed})")

    # Set random seed for fire frequency probability distribution. If -999 is given, a random seed will be generated. Otherwise, the given seed will be used.
    if batch_cfg.firefreq_random_seed == -999:
        batch_cfg.firefreq_random_seed = random.randint(0, 1000000)
        logger.info(f"Generated new random seed ({batch_cfg.firefreq_random_seed}) for fire frequency probability distribution.")
    else:
        batch_cfg.firefreq_random_seed = batch_cfg.firefreq_random_seed
        logger.info(f"Using given random seed ({batch_cfg.firefreq_random_seed}) for fire frequency probability distribution.")

    return batch_cfg


def determine_no_processes(batch_cfg):
    if batch_cfg.no_processes < 0:
        cpu_count = os.cpu_count()
        logger.info(f"No number of processes set. Detected {cpu_count} CPU cores, will therefore use {cpu_count-1} processes.")
        batch_cfg.no_processes = max(1, cpu_count - 1)

    return batch_cfg


def determine_total_results_csv(batch_cfg):
    batch_cfg.total_results_csv = batch_cfg.csv_parent_dir + "/{}_results.csv".format(batch_cfg.csv_parent_dir.split("state_data/")[1])
    return batch_cfg


def load_batch_config(batch_cfg):
    batch_cfg = read_config_from_file(batch_cfg)
    batch_cfg = parse_jobs(batch_cfg)
    batch_cfg = determine_csv_parent_dir(batch_cfg)
    batch_cfg = determine_total_results_csv(batch_cfg) 
    batch_cfg = set_random_seeds(batch_cfg)
    batch_cfg = determine_no_processes(batch_cfg)
    batch_cfg.vis = vis.Visualiser(batch_cfg)
    batch_cfg.color_dict = create_color_dict(batch_cfg)
    batch_cfg.headless = True
    batch_cfg.init_csv = True

    return batch_cfg


def run_batch(batch_cfg, proc_id, sim_counter):
    configure_logger(logname="batch.log", format="%(levelname)s %(processName)s: %(message)s", **vars(batch_cfg))

    while True:
        with sim_counter.get_lock():
            n_finished_simulations = sim_counter.value
            sim_counter.value += 1        

        job = batch_cfg.jobs.get(n_finished_simulations)
        if not job:
            break
        
        logger.debug(f"\n\nProcess {proc_id} starting simulation {n_finished_simulations+1}/{batch_cfg.jobs.count()} with parameters: {job}")
 
    return True


def parse_jobs(batch_cfg):
    batch_cfg.jobs = Jobs(**vars(batch_cfg))
    logger.info("Number of jobs to run (including reruns): {}".format(batch_cfg.jobs.count()))
    
    return batch_cfg


def configure_logger(logname=None, batch_verbosity=None, _format="%(levelname)s:%(name)s:%(message)s", **_):
    if batch_verbosity == "info":
        logging.basicConfig(filename=logname, level=logging.INFO, format=_format)
    elif batch_verbosity == "debug":
        logging.basicConfig(filename=logname, level=logging.DEBUG, format=_format)

    handler = get_stdout_logging_handler()
    logger.addHandler(handler)


def export_args(batch_cfg):
    args_json_path = batch_cfg.csv_parent_dir + "/parameters.json"
    with open(args_json_path, "w") as args_json:
        args = deepcopy(batch_cfg.jobs.get_specific_job(0, return_dict=True))
        for k, v in batch_cfg.arguments.items():
            args[k] = v
        json.dump(args, args_json)


def manage_processes(batch_cfg):
    procs = Processes()
    n_finished_simulations = Value("i", 0)
    logger.info("Starting processes...")
    for i in range(batch_cfg.no_processes):
        procs.add(Process(target=run_batch, args=(batch_cfg, i, n_finished_simulations)))
    logger.info("Finished starting processes.")

    while not procs.finished(print_progress=True):
        time.sleep(0.4)
    procs.join()

    logger.info("All processes have finished.")


def main(batch_cfg):
    os.remove("batch.log")
    configure_logger(logname="batch.log", **vars(batch_cfg))
    batch_cfg = load_batch_config(batch_cfg)
    export_args(batch_cfg)
    manage_processes(batch_cfg)

        
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-cfg', '--config', default="constant_batch", type=str)
    parser.add_argument('-np', '--no_processes', type=int, default=-1, help="Number of processes to run simultaneously.")
    parser.add_argument('-r', '--run', type=str, default=1, help="Unique run identifier in case of multiple runs per batch")
    parser.add_argument('-rs', '--random_seed', type=str, default=-999, help="Random seed for the batch. If seed != -999, each individual simulation may still get a unique random seed, but according to a reproducible sequence. If a random seed (!= -999) is provided in the batch config, this seed will be used for each simulation.")
    parser.add_argument('-rsf', '--firefreq_random_seed', type=str, default=-999, help="Random seed for fire frequencies for the batch. If != -999, each individual simulation may still get a unique random seed, but according to a reproducible sequence. If a random seed (!= -999) is provided in the batch config, this seed will be used for each simulation.")
    parser.add_argument('-vrb', '--batch_verbosity', type=str, default="info", help="Verbosity level for the batch runs. Set to 0 to suppress console output.")
    cfg = parser.parse_args()

    try:
        main(cfg)
    except Exception as e:
        logger.error("An error occurred during execution of the batch module.")
        determine_csv_parent_dir(cfg)
        csv_file = cfg.csv_parent_dir + "/Error_log.txt";
        with open(csv_file, 'a') as f:
            f.write(str(e))
            f.write(traceback.format_exc())
            
        traceback.print_exc()
        time.sleep(10)

