from ast import arg
import random
import traceback
import json
import sys
import logging
from types import SimpleNamespace
from copy import deepcopy
from argparse import ArgumentParser
from multiprocessing import Process, Value
import time
import os
import numpy as np

import config
import app
import file_handling as _io
import visualization as vis
import helpers as h

logger = logging.getLogger(__name__)


class Jobs:
    """Parses and stores the set of argument values of each simulation job as a namespace object."""

    def __init__(self, runs=None, arguments=None, rng=None, firefreq_rng=None, **args):
        # Set from given args
        self.n_runs = runs
        self.arg_changes = arguments
        self.rng = rng
        self.firefreq_rng = firefreq_rng
        
        # Set class defaults
        self.no_finished_rerun_cycles = 0
        self.sampling_mode = "regular" # alternative: "random"
        
        # Compute/derive
        self.defaults = self.get_defaults(args)
        self.update_default_values()
        self.jobs, self.job_indices = self.parse(arguments)

    def update_default_values(self):
        self.default_values = list(self.defaults.values())

    def generate_and_apply_random_seeds(self, job):
        # Apply new random seeds to a copy of the job, so that the original job remains unchanged (otherwise each rerun of a job will have the same seed)
        
        job_copy = deepcopy(job)
        if job.random_seed == -999:
            job_copy.random_seed = int(self.rng.integers(0, 100000000))
        if job.firefreq_random_seed == -999:
            job_copy.firefreq_random_seed = int(self.firefreq_rng.integers(0, 1000000))
        
        return job_copy

    def get_comprehensive_jobs_summary(self):
        summary = deepcopy(self.defaults)
        for key, arg_cfg in self.arg_changes.items():
            vec = self.get_vec(arg_cfg)
            summary[key] = {"value(s)": vec}
        
        return summary

    def get_control_variables(self, job):
        if isinstance(job, dict):
            job_ns = SimpleNamespace(**job)
        else: 
            job_ns = job
        
        ctrl_vars = {}
        for key, _val in self.arg_changes.items():
            # Only include control variables that vary across simulations (these are dictionaries).
            if type(_val) != dict:
                continue

            # Handle special arguments
            value = getattr(job_ns, key)
            if key == "treecover":
                key = "initial tree cover"
            elif _val.get("sub_arguments"): # Nested arguments
                sub_args = _val["sub_arguments"]
                for sub_arg_key in sub_args.keys():
                    sub_value = h.get_nested_dict_value(value, sub_arg_key)
                    ctrl_vars[key + ":" + sub_arg_key] = sub_value
                continue
            
            ctrl_vars[key] = value
        
        print("Control vars: ", ctrl_vars)
            
        # Add random seeds if they haven't been added yet
        for rs in ["random_seed", "firefreq_random_seed"]:
            if not rs in ctrl_vars:
                ctrl_vars[rs] = getattr(job_ns, rs)

        return ctrl_vars

    def get_defaults(self, args):
        defaults = config.get_all_defaults()
        defaults["headless"] = True
        defaults["verbosity"] = -1 # Suppress all non-critical print statements
        defaults["EXPORT_DIR"] = args["csv_parent_dir"]
        
        return defaults

    def get_key_idx(self, key):
        return list(self.defaults.keys()).index(key)

    def get_vec(self, arg_cfg):
        def round_to_significant_digits(vec, minim, maxim, stepsize):
            no_digits_after_comma = max(
                h.digits_after_decimal(minim),
                h.digits_after_decimal(maxim),
                h.digits_after_decimal(stepsize)
            )
            vec = [float(round(v, no_digits_after_comma)) for v in vec]
            return vec

        if type(arg_cfg) != dict:
            # If no dictionary is provided, we instead assume a constant value across all simulation jobs.
            value = arg_cfg
            vec = [value]
        elif arg_cfg["interpolation"] == "linear":
            minim, maxim, stepsize = arg_cfg["range"] + [arg_cfg["stepsize"]]
            vec = list(np.arange(minim, maxim + stepsize/2, stepsize))
            vec = round_to_significant_digits(vec, minim, maxim, stepsize)
        elif arg_cfg["interpolation"] == "powerten":
            minim, maxim, stepsize = arg_cfg["range"] + [arg_cfg["stepsize"]]
            exponents = list(np.arange(minim, maxim + stepsize/2, stepsize))
            vec = [10**p for p in exponents]
            vec = round_to_significant_digits(vec, minim, maxim, stepsize)
        elif arg_cfg["interpolation"] == "categorical":
            vec = arg_cfg["range"]
        elif arg_cfg["interpolation"] == "nested":
            # Create a dictionary containing the sub-arguments and their associated value ranges.
            vec = {}
            for sub_arg_key, sub_arg_cfg in arg_cfg["sub_arguments"].items():
                sub_vec = self.get_vec(sub_arg_cfg) # Parse the sub-arg config just like a normal arg config.
                vec[sub_arg_key] = sub_vec

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
        job_idx_generator = h.get_random_int_generator
        
        return jobs, job_idx_generator

    def do_post_processing(self, jobs, job_idx_generator, generator_cfg):
        job_indices = list(job_idx_generator(generator_cfg))
        jobs = [self.convert_value_set_to_namespace(job) for job in jobs]
        for job in jobs:
            h.check_cli_args(**vars(job))
        
        return jobs, job_indices

    def parse(self, arg_changes):
        job_count = self.derive_unique_job_count(arg_changes)
        logger.info("Expecting {} unique jobs.".format(job_count))
        generator_cfg = SimpleNamespace()
        if job_count < 1e6:
            jobs = self.parse_arg_values(arg_changes)
            generator_cfg.n = len(jobs)
            job_idx_generator = h.midpoint_gap_indices
        else:
            jobs, job_idx_generator = self.switch_to_random_sampling(job_count, generator_cfg)

        jobs, job_indices = self.do_post_processing(jobs, job_idx_generator, generator_cfg)
 
        return jobs, job_indices

    def parse_arg_values(self, arg_changes):
        value_sets = self.default_values.copy()
        for key, arg_cfg in arg_changes.items():
            vec = self.get_vec(arg_cfg)
            idx = self.get_key_idx(key)
            if type(vec) == dict:
                # Obtain base dictionary and update the defaults accordingly
                base_dict = arg_cfg["base"]
                is_single_value_set = type(value_sets[0]) != list
                self.apply_single_arg_change(value_sets, is_single_value_set, idx, base_dict)
                for sub_arg_key, sub_vec in vec.items():
                    value_sets = self.add_range(value_sets, idx, sub_vec, sub_keys=sub_arg_key)
            else:
                value_sets = self.add_range(value_sets, idx, vec)

        # If we're doing a multi-dimensional sensitivity analysis, we shuffle the order of the jobs
        if not type(value_sets[0]) == list:
            self.rng.shuffle(value_sets)

        return value_sets

    def get_specific_job(self, idx, return_dict=False):
        if return_dict:
            return vars(self.jobs[idx])
        else:
            return self.jobs[idx]

    def count(self):
        return self.n_runs * len(self.jobs)

    def get(self, n_started_simulations):
        if n_started_simulations >= self.count():
            return None

        idx = self.job_indices[n_started_simulations % len(self.jobs)]
        job = self.get_specific_job(idx)
        job_copy = self.generate_and_apply_random_seeds(job)

        return job_copy

    def apply_single_arg_change(self, value_sets, is_single_value_set, idx, value, sub_keys=None):
        if is_single_value_set:
            if sub_keys:
                main_value = deepcopy(value_sets[idx])
                sub_value = value
                h.set_nested_dict_value(main_value, sub_keys, sub_value)
            value_sets[idx] = value
        else:
            for value_set in value_sets:
                if sub_keys:
                    parent_dictionary = deepcopy(value_set[idx])
                    sub_value = value
                    h.set_nested_dict_value(parent_dictionary, sub_keys, sub_value)
                    value_set[idx] = parent_dictionary
                else:
                    value_set[idx] = value

        return value_sets

    def add_range(self, value_sets, idx, vec, sub_keys=None):
        """Add a range of values (i.e., a vector) to the argument sets at position idx.
        
        vec (list): Range of values to add.
        idx (int): Index into each argument list (='value_set') of the argument to modify.
        value_sets (list): list of argument lists (each sublist is a 'value_set'), or a single, flat value list if no arg changes have yet been made.
        sub_keys (str): Coded name of the sub-argument to modify in case of nested argument changes. Format: "<sub_key>:<sub_sub_key>:...:etc".
        """

        is_single_value_set = type(value_sets[0]) != list
        if is_single_value_set:
            expanded_value_sets = [value_sets.copy() for i in range(len(vec))]
        else:
            expanded_value_sets = []
            for i in range(len(vec)):
                expanded_value_sets += deepcopy(value_sets)

        stepsize = len(value_sets)
        for i, value in enumerate(vec):
            value_sets = self.apply_single_arg_change(deepcopy(value_sets), is_single_value_set, idx, value, sub_keys=sub_keys)
            begin = i * stepsize
            end = begin + stepsize
            if is_single_value_set:
                _value_sets = [value_sets]
            else:
                _value_sets = value_sets
                
            expanded_value_sets[begin:end] = _value_sets

        return expanded_value_sets


def get_new_batch_folder():
    batch_no = 1
    batch_folder = config.cfg.DATA_OUT_DIR + "/state_data/batch_000001"
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


def set_batch_random_seeds(batch_cfg):
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


def determine_n_processes(batch_cfg):
    if batch_cfg.n_processes < 0:
        cpu_count = os.cpu_count()
        logger.info(f"No number of processes set. Detected {cpu_count} CPU cores, will therefore use {cpu_count-1} processes.")
        batch_cfg.n_processes = max(1, cpu_count - 1)

    return batch_cfg


def determine_total_results_csv(batch_cfg):
    batch_cfg.total_results_csv = batch_cfg.csv_parent_dir + "/{}_results.csv".format(batch_cfg.csv_parent_dir.split("state_data/")[1])
    return batch_cfg


def load_batch_config(batch_cfg):
    batch_cfg = read_config_from_file(batch_cfg)
    batch_cfg = determine_csv_parent_dir(batch_cfg)
    batch_cfg = determine_total_results_csv(batch_cfg) 
    batch_cfg = set_batch_random_seeds(batch_cfg)
    batch_cfg = parse_jobs(batch_cfg)
    batch_cfg = determine_n_processes(batch_cfg)
    batch_cfg.vis = vis.Visualiser(batch_cfg)
    batch_cfg.color_dict = create_color_dict(batch_cfg)
    batch_cfg.headless = True

    return batch_cfg


def export_state(batch_cfg, dynamics, job, sim_cfg, init_csv):
    sim_cfg.control_vars = batch_cfg.jobs.get_control_variables(job)
    initialize_csv = False
    if init_csv.value == 1:
        with init_csv.get_lock(): # We use a lock to ensure the total results csv is only initialized once, by a single process.
            if init_csv.value == 1:
                init_csv.value = 0
                initialize_csv = True

    _io.export_state(
        dynamics, path=batch_cfg.total_results_csv, init_csv=initialize_csv, cfg=sim_cfg
    )


def run_sim(batch_cfg, job, init_csv):
    dynamics, sim_cfg = app.main(**vars(job))
    export_state(batch_cfg, dynamics, job, sim_cfg, init_csv)


def suppress_irrelevant_console_output():
    # Ignore a numpy warning that pops up during visualization
    np.seterr(invalid="ignore")
    
    # Redirect stdout to the null device, ensuring all console output of the child processes is suppressed.
    sys.stdout = open(os.devnull, "w")


def run_batch(batch_cfg, proc_id, sim_counter, finished_sim_counter, init_csv):
    # We don't want detailed information to pop up about every single simulation, so we suppress it.
    suppress_irrelevant_console_output()

    # Initialize the logger for this process.
    configure_logger(logname="batch.log", format="%(levelname)s %(processName)s: %(message)s", **vars(batch_cfg))

    while True:
        with sim_counter.get_lock():
            n_started_simulations = sim_counter.value
            sim_counter.value += 1        

        job = batch_cfg.jobs.get(n_started_simulations)
        if not job:
            break

        run_sim(batch_cfg, job, init_csv)
        with finished_sim_counter.get_lock():
            finished_sim_counter.value += 1
 
    return True


def parse_jobs(batch_cfg):
    batch_cfg.jobs = Jobs(**vars(batch_cfg))
    logger.info("Number of simulations to run (including re-runs): {}".format(batch_cfg.jobs.count()))
    
    return batch_cfg


def configure_logger(logname=None, batch_verbosity=None, _format="%(levelname)s:%(name)s:%(message)s", **_):
    if batch_verbosity == "info":
        logging.basicConfig(filename=logname, level=logging.INFO, format=_format)
    elif batch_verbosity == "debug":
        logging.basicConfig(filename=logname, level=logging.DEBUG, format=_format)

    handler = h.get_stdout_logging_handler()
    logger.addHandler(handler)


def remove_non_serializable_items(batch_cfg_copy):
    deletions = []
    for k, v in batch_cfg_copy.items():
        try:
            json.dumps(v)
        except TypeError:
            deletions.append(k)
    for k in deletions + ["help", "arguments_examples"]:
        del batch_cfg_copy[k]


def export_batch_cfg(batch_cfg):
    args_json_path = batch_cfg.csv_parent_dir + "/configuration.json"
    with open(args_json_path, "w") as args_json:
        # Get the arguments from one job as a representative set of arguments for the batch, and export these to json. 
        # This way, we have a record of the arguments used for the batch, without having to export the arguments 
        # of each individual simulation.
        batch_cfg_copy = vars(deepcopy(batch_cfg))
        args = batch_cfg.jobs.get_comprehensive_jobs_summary()
        batch_cfg_copy["arguments"] = args
        
        # Remove non-serializable items from config to avoid issues when exporting to json
        remove_non_serializable_items(batch_cfg_copy)
            
        # Export to json file
        json.dump(batch_cfg_copy, args_json, indent=4)


def report_sim_completions(batch_cfg, n_reported_completions, n_finished_simulations):
    while (n_finished_simulations.value > n_reported_completions) and (n_reported_completions < batch_cfg.jobs.count()):
        n_reported_completions += 1
        cfg_str = json.dumps(batch_cfg.jobs.get_control_variables(batch_cfg.jobs.get(n_reported_completions-1)), indent=4)
        print(f"\nFinished simulation {n_reported_completions}/{batch_cfg.jobs.count()} with arguments: \n{cfg_str}.")
    
    return n_reported_completions


def manage_processes(batch_cfg):
    # Initialize Processes-container and shared values
    procs = h.Processes()
    n_started_simulations = Value("i", 0)
    n_finished_simulations = Value("i", 0)
    init_csv = Value("i", 1)
    
    # Add processes to container and start them
    n_reported_completions = 0
    logger.info("Starting processes...")
    for i in range(batch_cfg.n_processes):
        procs.add(Process(target=run_batch, args=(batch_cfg, i, n_started_simulations, n_finished_simulations, init_csv)))
    logger.info("Finished starting processes.")

    # Wait for all processes to finish and report completed simulations in the meantime
    while not procs.finished(print_progress=True):
        time.sleep(0.5)
        n_reported_completions = report_sim_completions(batch_cfg, n_reported_completions, n_finished_simulations)
        
    procs.join()
    logger.info("Batch complete.")


def main(batch_cfg):
    if os.path.exists("batch.log"):
        os.remove("batch.log")
    configure_logger(logname="batch.log", **vars(batch_cfg))
    batch_cfg = load_batch_config(batch_cfg)
    export_batch_cfg(batch_cfg)
    manage_processes(batch_cfg)

        
if __name__ == "__main__":
    print("Starting batch process...")
    parser = ArgumentParser()
    parser.add_argument('-cfg', '--config', default="constant_batch", type=str)
    parser.add_argument('-np', '--n_processes', type=int, default=-1, help="Number of processes to run simultaneously.")
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

