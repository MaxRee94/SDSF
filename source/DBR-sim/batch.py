from ast import arg
from distutils.dep_util import newer_group
import random
from stat import FILE_ATTRIBUTE_ENCRYPTED
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

    def __init__(self, runs=None, arguments=None, rng=None, firefreq_rng=None, **batch_args):
        # Set from given args
        self.n_runs = runs
        self.arg_changes = arguments
        self.rng = rng
        self.firefreq_rng = firefreq_rng
        
        # Set class defaults
        self.no_finished_rerun_cycles = 0
        self.sampling_mode = "regular" # alternative: "random"
        self.first_keyframe_with_intersim_variation = {}
        
        # Compute/derive
        self.defaults = self.get_defaults(batch_args, arguments)
        self.batch_type = batch_args["type"]
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
            vec = self.get_vec(arg_cfg, key)
            summary[key] = {"value(s)": vec}
        
        return summary
    
    def get_first_keyframe_with_intersim_variation(self, arg_key, job):
        if self.first_keyframe_with_intersim_variation.get(arg_key) is None:
            return None
        else:
            keytime = self.first_keyframe_with_intersim_variation[arg_key]
            return job.keyframes[arg_key][keytime]

    def attach_control_variables(self, job):
        if isinstance(job, dict):
            job_ns = SimpleNamespace(**job)
        else: 
            job_ns = job
        
        ctrl_vars = {}
        control_config = deepcopy(self.arg_changes)
        control_config["keyframes"] = {} # Add keyframes field with '{}' as placeholder
        for key, _val in control_config.items():
            # Only include control variables that vary across simulations (these are dictionaries).
            if type(_val) != dict:
                continue
            
            # Ignore default settings of keyframed arguments. These will be overwritten by the corresponding keyframe list.
            if key in list(job_ns.keyframes.keys()):
                continue
            
            # Ignore other non-standard arguments (e.g., 'forest_suitability')
            if not self.is_standard_argument(key):
                continue

            # Get job-specific value of given control variable
            value = getattr(job_ns, key)
            
            # Handle special arguments
            if key == "treecover":
                key = "initial tree cover"
            if key == "keyframes":
                for arg_key, keyframe_vec in job_ns.keyframes.items():
                    # first_variable_keyframe = self.get_first_keyframe_with_intersim_variation(arg_key, job)
                    # if first_variable_keyframe is None:
                    #     first_keyframe_time = list(keyframe_vec.keys())[0]
                    #     keyframe_val = keyframe_vec[first_keyframe_time]
                    # else:
                    #     keyframe_val = first_variable_keyframe
                    ctrl_vars[arg_key] = keyframe_vec
                continue
            elif _val.get("sub_arguments"): # Nested arguments
                sub_args = _val["sub_arguments"]
                for sub_arg_key in sub_args.keys():
                    sub_value = h.get_nested_dict_value(value, sub_arg_key)
                    ctrl_vars[key + ":" + sub_arg_key] = sub_value
                continue
            
            ctrl_vars[key] = value
            
        # Add random seeds if they haven't been added yet
        for rs in ["random_seed", "firefreq_random_seed"]:
            if not rs in ctrl_vars:
                ctrl_vars[rs] = getattr(job_ns, rs)

        # Attach to job namespace
        job_ns.control_vars = ctrl_vars
 
        return job_ns

    def attach_sim_name(self, job):
        full_sim_name = _io.get_sim_name(job)
        job.sim_name = _io.get_sim_name(job, extra_short=True)

        # Ensure all data is written to sim directory
        job.EXPORT_DIR = os.path.join(job.EXPORT_DIR, full_sim_name)
        job.DATA_OUT_DIR = job.EXPORT_DIR
        
        return job

    def get_defaults(self, batch_args, arguments):
        defaults = config.get_all_defaults()
        
        # Overwrite some defaults that are required for a batch run.
        defaults["headless"] = True
        defaults["verbosity"] = -1 # Suppress all non-critical print statements
        defaults["EXPORT_DIR"] = batch_args["csv_parent_dir"]
        defaults["keyframes"] = {}
        for arg_key in arguments.keys():
            if arg_key not in defaults.keys():
                defaults[arg_key] = None # Add a default value of None to any missing argument keys
        
        return defaults

    def get_key_idx(self, key):
        return list(self.defaults.keys()).index(key)

    def get_keyframe_vec(self, arg_cfg, arg_key):
        # Expand keyframe vectors based on batch config settings.
        keyframe_vecs = {}
        longest_kfv = {}
        first_keyframe_with_intersim_variation = None
        for keyframe_cfg in arg_cfg["keyframes"]:
            time = keyframe_cfg["time"]
            keyframe_vec = self.get_vec(keyframe_cfg)
            keyframe_vecs[time] = keyframe_vec
            if len(keyframe_vec) > len(longest_kfv):
                longest_kfv = keyframe_vec
                self.first_keyframe_with_intersim_variation[arg_key] = time           
            
        # Duplicate constant keyframe values to construct a vector with a length equal to that of any interpolated keyframe vectors.
        # If no interpolated keyframe vectors exist, no duplication will occur and the constant keyframe vectors will simply be of length 1.
        keyframe_vecs_extended = {}
        for kf_time, kfv in keyframe_vecs.items():
            if len(kfv) < len(longest_kfv):
                kfv = kfv * len(longest_kfv)
            keyframe_vecs_extended[kf_time] = kfv
                
        # Check that each keyframe vector is of equal length now
        for _, kfv in keyframe_vecs_extended.items():
            assert(len(kfv) == len(longest_kfv)), "All non-constant keyframes should have an equal number of values. Check batch config JSON."
                
        # 'Transpose', to obtain keyframes for each individual simulation (rather than a list of values across all sims, for each keytime).
        vec = []
        for sim_idx in range(len(longest_kfv)):
            sim_keyframes = {}
            for keytime, kfv in keyframe_vecs_extended.items():
                sim_keyframes[keytime] = kfv[sim_idx] # The i-th value of each keyframe vector belongs to the corresponding sim.
            vec.append(sim_keyframes)

        return vec

    def get_vec(self, arg_cfg, arg_key=None):
        
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
        elif arg_cfg["interpolation"] == "constant":
            value = arg_cfg["value"]
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
        elif arg_cfg["interpolation"] == "keyframed":
            vec = self.get_keyframe_vec(arg_cfg, arg_key)

        return vec

    def convert_value_set_to_namespace(self, job_values):
        job = self.convert_job_to_dict(job_values)

        return SimpleNamespace(**job)

    def convert_job_to_dict(self, job_values):
        job = {}
        for i, (key, _) in enumerate(self.defaults.items()):
            job[key] = job_values[i]

        return job
    
    def arg_cfg_contains_keyframes(self, arg_cfg):
        if type(arg_cfg) == dict and arg_cfg.get("interpolation") == "keyframed":
            return True
        return False

    def vec_contains_keyframes(self, vec):
        if type(vec) == list and type(vec[0]) == dict:
            print("first")
            first_subkey_of_first_value = list(vec[0].keys())[0]
            if type(first_subkey_of_first_value) == int:
                return True
        if type(vec) == dict:
            print("second")
            for value in vec.values():
                if self.vec_contains_keyframes(value):
                    return True
        return False
    
    def contains_nested_args(self, vec):
        return type(vec) == dict

    def is_standard_argument(self, arg_key):
        if arg_key in ["forest_suitability"]:
            return False
        return True

    def derive_unique_job_count(self, arg_changes):
        job_count = 0
        if job_count == 0:
            job_count = 1
        for key, arg_cfg in arg_changes.items():
            vec = self.get_vec(arg_cfg, key)
            if self.arg_cfg_contains_keyframes(arg_cfg):
                job_count *= len(vec)
            elif self.contains_nested_args(vec):
                for sub_vec in vec.values():
                    job_count *= len(sub_vec)
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

    def do_post_processing(self, jobs, job_idx_generator, generator_cfg, arg_changes):
        job_indices = list(job_idx_generator(generator_cfg))
        jobs = [self.convert_value_set_to_namespace(job) for job in jobs]
        for job in jobs:
            h.check_cli_args(**vars(job))
            job = self.add_attachments(job, arg_changes)
        
        return jobs, job_indices
    
    def add_attachments(self, job, arg_changes):
        job = self.attach_control_variables(job)
        job = self.attach_sim_name(job)
        job.batch_type = self.batch_type
        
        return job
    
    def traverse_dict_and_apply_suitability_coefficients(self, thedict, coefficients, job):
        if thedict.get("relation_with_forest_suitability"):
            job_specific_relation = thedict["relation_with_forest_suitability"]
            for ck, cv in coefficients.items(): 
                job_specific_relation = job_specific_relation.replace(ck, str(cv))
                                
            thedict = "SUITABILITY-DERIVED:" + job_specific_relation
            return thedict
        else:
            for key, value in thedict.items():
                # Recurse
                if type(value) == dict:
                    thedict[key] = self.traverse_dict_and_apply_suitability_coefficients(value, coefficients, job)

        return thedict

    def apply_suitability_relation_coefficients_to_nested_arg(self, arg_key, arg_cfg, arg_changes, jobs):
        coefficient_keys = list(arg_cfg.get("sub_arguments").keys())
        for i, job in enumerate(jobs):
            coefficients = {k: getattr(job, arg_key).get(k) for k in coefficient_keys} # Get the coefficients for the current job.
            job_arg_dict = getattr(job, arg_key)
            job_arg_dict = self.traverse_dict_and_apply_suitability_coefficients(job_arg_dict, coefficients, job)
            setattr(job, arg_key, job_arg_dict)

    def apply_suitability_relation_coefficients(self, arg_changes, jobs):
        """In case of a bifurcation analysis with arguments whose value is determined by a relation with the 'forest_suitability' argument,
        ensure the coefficients in this relation are assigned the parsed values."""

        for arg_key, arg_cfg in arg_changes.items():
            if type(arg_cfg) == dict and arg_cfg.get("base", {}):
                is_flat_suitability_arg = arg_cfg["base"].get("relation_with_forest_suitability") # A flat suitability arg is a non-dict argument controlled by forest suitability.
                if is_flat_suitability_arg:
                    relation = arg_cfg["base"]["relation_with_forest_suitability"]
                    coefficient_keys = list(arg_cfg.get("sub_arguments").keys())
                    for job in jobs:
                        job_specific_relation = relation
                        for ck in coefficient_keys: 
                            coefficient_value = getattr(job, arg_key).get(ck)
                            job_specific_relation = job_specific_relation.replace(ck, str(coefficient_value))
                        setattr(job, arg_key, "SUITABILITY-DERIVED:" + job_specific_relation)
                else:
                    # Argument is a dictionary argument. We need to check for each sub-argument whether it is controlled by forest
                    # suitability, and if so, apply the corresponding coefficients.
                    self.apply_suitability_relation_coefficients_to_nested_arg(arg_key, arg_cfg, arg_changes, jobs)
        
        return jobs

    def expand_forest_suitability_arg(self, arg_changes, jobs):
        forest_suitability_vec = self.get_vec(arg_changes["forest_suitability"]["value"], "forest_suitability")
        
        for job in jobs:
            job.forest_suitability = forest_suitability_vec
        
        return jobs

    def parse(self, arg_changes):
        if self.batch_type == "bifurcation_analysis":
            jobs, job_indices = self.parse_sensitivity_analysis(arg_changes)
            jobs = self.apply_suitability_relation_coefficients(arg_changes, jobs)
            jobs = self.expand_forest_suitability_arg(arg_changes, jobs)
            return jobs, job_indices
        else:
            return self.parse_sensitivity_analysis(arg_changes)

    def parse_sensitivity_analysis(self, arg_changes):
        job_count = self.derive_unique_job_count(arg_changes)
        logger.info("Expecting {} unique jobs.".format(job_count))
        generator_cfg = SimpleNamespace()
        if job_count < 1e6:
            jobs = self.parse_arg_values(arg_changes)
            generator_cfg.n = len(jobs)
            job_idx_generator = h.midpoint_gap_indices
        else:
            jobs, job_idx_generator = self.switch_to_random_sampling(job_count, generator_cfg)

        jobs, job_indices = self.do_post_processing(jobs, job_idx_generator, generator_cfg, arg_changes)
 
        return jobs, job_indices

    def parse_arg_values(self, arg_changes):
        value_sets = self.default_values.copy()
        
        # Expand argument values based on batch config settings.
        expanded_arguments = {"keyframes": {"idx": self.get_key_idx("keyframes")}}
        for key, arg_cfg in arg_changes.items():
            vec = self.get_vec(arg_cfg, key)
            idx = self.get_key_idx(key)
            if self.arg_cfg_contains_keyframes(arg_cfg):
                if self.contains_nested_args(vec):
                    for subkey, sub_vec in deepcopy(vec).items():
                        if self.vec_contains_keyframes(sub_vec):
                            expanded_arguments["keyframes"][key + ":" + subkey] = sub_vec
                            vec[key + ":" + subkey] = ["KEYFRAMED"] # Placeholder will be overwritten by self.apply_first_keyframes()
                            del vec[subkey]
                else:
                    expanded_arguments["keyframes"][key] = vec
                    vec = ["KEYFRAMED"] # Placeholder will be overwritten by self.apply_first_keyframes()
            expanded_arguments[key] = {"vec": vec, "idx": idx, "arg_cfg": arg_cfg}
        
        # Add an additional dimension to the jobs tensor, and apply the expanded argument values across all existing jobs.
        for key, expanded_cfg in expanded_arguments.items():
            if key == "keyframes":
                # Create base dictionary and update the default dictionary accordingly
                base_dict = {subkey : {} for subkey in expanded_cfg.keys() if subkey != "idx"}    # Store an empty dictionary for each subkey.
                                                                                                                # A 'subkey' is the name of an argument, stored inside
                                                                                                                # the keyframes dict.
                is_single_value_set = type(value_sets[0]) != list
                self.apply_single_arg_change(value_sets, is_single_value_set, expanded_cfg["idx"], base_dict)
                
                # Add ranges of subkeys
                for subkey, keyframe_sets in expanded_cfg.items():
                    if subkey == "idx":
                        continue
                    value_sets = self.add_range(value_sets, expanded_cfg["idx"], keyframe_sets, subkey=subkey, has_nested_keyframes=len(subkey.split(":")) > 1)
            elif self.contains_nested_args(expanded_cfg["vec"]):
                # Obtain base dictionary and update the default dictionary accordingly
                base_dict = expanded_cfg["arg_cfg"]["base"]
                is_single_value_set = type(value_sets[0]) != list
                self.apply_single_arg_change(value_sets, is_single_value_set, expanded_cfg["idx"], base_dict)
                
                # Apply ranges of sub-arguments
                for sub_arg_key, sub_vec in expanded_cfg["vec"].items():
                    if "KEYFRAMED" in sub_vec:
                        continue # If there are nested args that are keyframed, skip these (will be handled by 'if key == "keyframes"' block).
                    value_sets = self.add_range(value_sets, expanded_cfg["idx"], sub_vec, subkey=sub_arg_key)
            else:
                value_sets = self.add_range(value_sets, expanded_cfg["idx"], expanded_cfg["vec"])

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
        self.apply_first_keyframes(job_copy)

        return job_copy

    def apply_first_keyframes(self, job):
        if job.keyframes:
            for arg_key, keyframes in job.keyframes.items():
                keytimes = sorted(keyframes) # Sort keyframes in ascending order of time
                first_keyframe_value = keyframes[keytimes[0]] # Get value of first keyframe (i.e., the keyframe with the earliest time)]
                if ":" in arg_key:
                    # Handle nested keyframes
                    arg_key_split = arg_key.split(":")
                    parent_key = arg_key_split[0]
                    subkey = ":".join(arg_key_split[1:])
                    arg_parent_dict = getattr(job, parent_key)
                    h.set_nested_dict_value(arg_parent_dict, subkey, first_keyframe_value)
                else:
                    setattr(job, arg_key, first_keyframe_value)

    def apply_single_arg_change(self, value_sets, is_single_value_set, idx, value, subkey=None, has_nested_keyframes=False):
        if is_single_value_set:
            if subkey:
                parent_dict = deepcopy(value_sets[idx]) # Obtain 'base' dict
                sub_value = value
                if has_nested_keyframes:
                    parent_dict[subkey] = sub_value
                else:
                    h.set_nested_dict_value(parent_dict, subkey, sub_value)
                value = parent_dict
            value_sets[idx] = value
        else:
            for value_set in value_sets:
                if subkey:
                    parent_dict = deepcopy(value_set[idx]) # Obtain 'base' dict
                    sub_value = value
                    if has_nested_keyframes:
                        parent_dict[subkey] = sub_value
                    else:
                        h.set_nested_dict_value(parent_dict, subkey, sub_value)
                    value_set[idx] = parent_dict
                else:
                    value_set[idx] = value

        return value_sets

    def add_range(self, value_sets, idx, vec, subkey=None, has_nested_keyframes=False):
        """Add a range of values (i.e., a vector) to the argument sets at position idx.
        
        vec (list): Range of values to add.
        idx (int): Index into each argument list (='value_set') of the argument to modify.
        value_sets (list): list of argument lists (each sublist is a 'value_set'), or a single, flat value list if no arg changes have yet been made.
        subkey (str): Coded name of the sub-argument to modify in case of nested argument changes. Format: "<subkey>:<sub_subkey>:...:etc".
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
            value_sets = self.apply_single_arg_change(
                deepcopy(value_sets), is_single_value_set, idx, value,
                subkey=subkey, has_nested_keyframes=has_nested_keyframes
            )
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
    initialize_csv = False
    if init_csv.value == 1:
        with init_csv.get_lock(): # We use a lock to ensure the total results csv is only initialized once, by a single process.
            if init_csv.value == 1:
                init_csv.value = 0
                initialize_csv = True

    _io.export_state(
        dynamics, path=batch_cfg.total_results_csv, init_csv=initialize_csv, cfg=sim_cfg, use_updated_ctrl_vars=False
    )


def init_sim_dir(job):
    if not os.path.isdir(job.EXPORT_DIR):
        os.makedirs(job.EXPORT_DIR)


def export_sim_cfg(sim_cfg):
    with open(sim_cfg.EXPORT_DIR + "/sim_configuration.json", "w") as f:
        cfg = remove_non_serializable_items(vars(sim_cfg))
        for key in ["patches", "color_dict"]:
            if key in cfg.keys():
                del cfg[key]
        json.dump(cfg, f)


def run_sim(batch_cfg, job, init_csv):
    init_sim_dir(job)
    dynamics, sim_cfg = app.main(**vars(job))
    export_state(batch_cfg, dynamics, job, sim_cfg, init_csv)
    export_sim_cfg(sim_cfg)


def run_batch(batch_cfg, proc_id, sim_counter, finished_sim_counter, init_csv):
    # We don't want detailed information to pop up about every single simulation, so we suppress it.
    h.suppress_irrelevant_console_output()

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


def remove_non_serializable_items(mydict):
    """Remove non serializable items from given dict."""

    new_dict = {}
    for k, v in mydict.items():
        try:
            json.dumps(v)
            new_dict[k] = v
        except TypeError:
            pass
    
    for k in ["help", "arguments_examples"]:
        if k in list(mydict.keys()):
            del new_dict[k]

    return new_dict


def export_batch_cfg(batch_cfg):
    args_json_path = batch_cfg.csv_parent_dir + "/batch_configuration.json"
    with open(args_json_path, "w") as args_json:
        # Get the arguments from one job as a representative set of arguments for the batch, and export these to json. 
        # This way, we have a record of the arguments used for the batch, without having to export the arguments 
        # of each individual simulation.
        batch_cfg_copy = vars(deepcopy(batch_cfg))
        args = batch_cfg.jobs.get_comprehensive_jobs_summary()
        batch_cfg_copy["arguments"] = args
        
        # Remove non-serializable items from config to avoid issues when exporting to json
        batch_cfg_copy = remove_non_serializable_items(batch_cfg_copy)
            
        # Export to json file
        json.dump(batch_cfg_copy, args_json, indent=4)


def report_sim_completions(batch_cfg, n_reported_completions, n_finished_simulations):
    while (n_finished_simulations.value > n_reported_completions) and (n_reported_completions < batch_cfg.jobs.count()):
        n_reported_completions += 1
        cfg_str = json.dumps(batch_cfg.jobs.get(n_reported_completions-1).control_vars, indent=4)
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

