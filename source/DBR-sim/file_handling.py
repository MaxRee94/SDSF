from cgitb import lookup
import csv
import os
from pathlib import Path
import datetime
from shutil import ExecError
import cv2
import time
from config import *
import helpers
import numpy as np

import sine_pattern_generator as spg
import simple_noise_generator as sng



def get_tree_sizes(dynamics, bins="initialize", cfg=None):
    if type(bins) == str and bins == "initialize":
        _tree_sizes = [i for i in range(int(cfg.max_dbh))]
        counts, bins = np.histogram(_tree_sizes, bins=5)
    _tree_sizes = dynamics.state.get_tree_sizes()
    counts, bins = np.histogram(_tree_sizes, bins=bins)
    return counts, bins


def create_output_dirs(main_output_dir):
    create_directory_if_not_exists(main_output_dir)
    create_directory_if_not_exists(os.path.join(main_output_dir, "state_data"))
    create_directory_if_not_exists(os.path.join(main_output_dir, "image_timeseries"))


def get_tree_ages(dynamics, bins="initialize"):
    if type(bins) == str and bins == "initialize":
        _tree_ages = [i for i in range(400)] # Harcoded max age of 400 years for now.
        counts, bins = np.histogram(_tree_ages, bins=5)
    _tree_ages = dynamics.state.get_tree_ages()
    counts, bins = np.histogram(_tree_ages, bins=bins)
    return counts, bins


def get_firefree_interval_stats(dynamics, _type):
    intervals = dynamics.get_firefree_intervals(_type)
    stdev = np.std(intervals)
    mean = np.mean(intervals)
    return mean, stdev


def load_json_config(file_path):
    """Load a JSON configuration file.
    
    Params:
        file_path (str): Path to the JSON file.
        
    Returns:
        config (dict): Loaded configuration dictionary.
    """
    with open(file_path, 'r') as json_file:
        config = json.load(json_file)
    return config


def set_heterogeneity_maps(dynamics, cfg):
    """Set input maps from cfg if provided.
    
    Params:
        dynamics (Dynamics): Dynamics object
        cfg(SimpleNamespace): User arguments
    """
    map_root_dir = f"{cfg.DATA_IN_DIR}/heterogeneity/"
    json_path = os.path.join(map_root_dir, cfg.heterogeneity)
    heterogeneity_map_cfg = load_json_config(json_path)
    map_types = {
        "grass_carrying_capacity": dynamics.state.grid.set_grass_carrying_capacity,
        "local_growth_multipliers": dynamics.state.grid.set_local_growth_multipliers,
    }
    for map_type, setter_func in map_types.items():
        m_cfg = heterogeneity_map_cfg[map_type]
        m_cfg = helpers.overwrite_from_global_arguments(m_cfg, vars(cfg))
        map_dir = os.path.join(map_root_dir, map_type)
        if m_cfg.get("filename"):
            impath = os.path.join(map_dir, m_cfg["filename"])
            image = cv2.imread(impath, cv2.IMREAD_GRAYSCALE)
            print(f"Read {map_type} map from image file:", impath)
        else:
            if m_cfg["type"] == "sine":
                image = cfg.spg.generate((dynamics.state.grid.width, dynamics.state.grid.width), **m_cfg)
            elif m_cfg["type"] == "noise":
                image = sng.generate(grid_width=dynamics.state.grid.width, cfg=cfg, **m_cfg)
            else:
                raise ValueError(f"Unknown {map_type} pattern type: {m_cfg['type']}")
            impath = os.path.join(map_dir, f"{map_type}.png")
            cv2.imwrite(impath, image)
            print(f"Generated {map_type} map saved to:", impath)

        image = cv2.resize(image, (dynamics.state.grid.width, dynamics.state.grid.width), interpolation=cv2.INTER_NEAREST)
        image = image / 255  # Normalize to 0-1
        setter_func(image)

    return dynamics, cfg


def create_directory_if_not_exists(_dir):
    if not os.path.exists(_dir):
        os.makedirs(_dir)
        if not os.path.exists(_dir):
            raise RuntimeError("Failed to create directory: " + _dir)
    else:
        raise RuntimeError("Directory already exists: " + _dir)


def export_state(
        dynamics, path="", init_csv=True, control_variable=None, control_value=None, tree_cover_slope=0,
        extra_parameters="", secondary_variable=None, secondary_value=None, dependent_var=None, dependent_val=None, initial_no_dispersals=None, dependent_result_range_stdev=None,
        cfg=None
    ):
    fieldnames = [
        "time", "treecover", "slope", "population_size", "#seeds_produced", "fires", "top_kills", "nonseedling_top_kills", "deaths",
        "extra_parameters", "mean tree size", "stdev tree size",
        "firefree_interval_mean", "firefree_interval_stdev", "firefree_interval_full_sim_mean", "firefree_interval_full_sim_stdev", "time_spent_moving",
        "shaded_out", "outcompeted_by_seedlings", "outcompeted_by_oldstems", "initial_no_dispersals",
        "recruits", "germination_attempts", "oldstem_competition_and_shading", "seedling_competition_and_shading", "perimeter_length", "perimeter-area_ratio",
        "basal_area"
    ]
    if dependent_var:
        fieldnames.insert(0, dependent_var)
    if secondary_variable:
        if secondary_variable == "treecover":
            secondary_variable = "initial_tree_cover"
        fieldnames.insert(0, secondary_variable)
    if initial_no_dispersals:
        fieldnames.insert(-3, "initial_number_of_dispersals")
    if control_variable:
        if control_variable == "treecover":
            control_variable = "initial tree cover"
        fieldnames.insert(0, control_variable)
    if dependent_result_range_stdev:
        fieldnames.insert(0, "dependent_result_range_stdev")
    if not cfg.rotate_randomly:
        fieldnames.insert(3, "global_rotation_offset")

    # Collect tree ages and sizes
    if init_csv:
        tree_sizes, cfg.treesize_bins = get_tree_sizes(dynamics, "initialize", cfg=cfg)
        tree_ages, cfg.tree_age_bins = get_tree_ages(dynamics, "initialize")
    else:
        tree_sizes, _ = get_tree_sizes(dynamics, cfg.treesize_bins)
        tree_ages, _ = get_tree_ages(dynamics, cfg.tree_age_bins)

    treesize_bins_strings = []
    tree_age_bins_strings = []
    for i in range(len(cfg.treesize_bins)-1):
        (size_lowb, size_highb) = (cfg.treesize_bins[i], cfg.treesize_bins[i+1])
        (age_lowb, age_highb) = (cfg.tree_age_bins[i], cfg.tree_age_bins[i+1])
        treesize_bins_strings.append("trees[dbh {0} - {1}]".format(round(size_lowb), round(size_highb)))
        tree_age_bins_strings.append("trees[age {0} - {1}]".format(round(age_lowb), round(age_highb)))
        fieldnames.insert(5, treesize_bins_strings[i])
        fieldnames.insert(6+i, tree_age_bins_strings[i])

    # Initialize CSV file with headers if it doesn't exist
    if init_csv and not os.path.exists(path):
        if path == "":
            print("\n\nExport location:", cfg.EXPORT_DIR)
            path = os.path.join(cfg.EXPORT_DIR, "Simulation_" + str(datetime.datetime.now()).replace(":", "-") + ".csv")
        with open(path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

    with open(path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        firefree_interval_mean, firefree_interval_stdev = get_firefree_interval_stats(dynamics, "current_iteration")
        firefree_interval_fullsim_mean, firefree_interval_fullsim_stdev = get_firefree_interval_stats(dynamics, "average")
        fires = dynamics.get_fires()
        fires = "|".join([str(fire) for fire in fires])
        result = {
            "time": str(dynamics.time),
            "treecover": str(dynamics.state.grid.get_tree_cover()), 
            "slope": str(tree_cover_slope),
            "population_size": str(dynamics.state.population.size()),
            "#seeds_produced": str(dynamics.seeds_produced),
            "fires": fires,
            treesize_bins_strings[0]: str(tree_sizes[0]),
            treesize_bins_strings[1]: str(tree_sizes[1]),
            treesize_bins_strings[2]: str(tree_sizes[2]),
            treesize_bins_strings[3]: str(tree_sizes[3]),
            treesize_bins_strings[4]: str(tree_sizes[4]),
            tree_age_bins_strings[0]: str(tree_ages[0]),
            tree_age_bins_strings[1]: str(tree_ages[1]),
            tree_age_bins_strings[2]: str(tree_ages[2]),
            tree_age_bins_strings[3]: str(tree_ages[3]),
            tree_age_bins_strings[4]: str(tree_ages[4]),
            "mean tree size": str(np.mean(dynamics.state.get_tree_sizes())),
            "stdev tree size": str(np.std(dynamics.state.get_tree_sizes())),
            "top_kills": str(dynamics.get_no_fire_induced_topkills()),
            "deaths": str(dynamics.get_no_fire_induced_deaths()),
            "extra_parameters": str(extra_parameters),
            "firefree_interval_mean": str(firefree_interval_mean),
            "firefree_interval_stdev": str(firefree_interval_stdev),
            "firefree_interval_full_sim_mean": str(firefree_interval_fullsim_mean),
            "firefree_interval_full_sim_stdev": str(firefree_interval_fullsim_stdev),
            "recruits": str(dynamics.get_no_recruits("all")),
            "initial_no_dispersals": str(dynamics.get_initial_no_dispersals()),
            "time_spent_moving": str(dynamics.get_fraction_time_spent_moving()),
            "shaded_out": str(dynamics.get_fraction_seedlings_dead_due_to_shade()),
            "nonseedling_top_kills": str(dynamics.get_no_fire_induced_nonseedling_topkills()),
            "outcompeted_by_seedlings": str(dynamics.get_fraction_seedlings_outcompeted()),
            "outcompeted_by_oldstems": str(dynamics.get_fraction_seedlings_outcompeted_by_older_trees()),
            "germination_attempts": str(dynamics.get_no_germination_attempts()),
            "oldstem_competition_and_shading": str(dynamics.get_fraction_cases_oldstem_competition_and_shading()),
            "seedling_competition_and_shading": str(dynamics.get_fraction_cases_seedling_competition_and_shading()),
            "perimeter-area_ratio": str(dynamics.get_perimeter_area_ratio()),
            "perimeter_length": str(dynamics.get_forest_perimeter_length()),
            "basal_area": str(dynamics.get_basal_area())
        }
        PAR_derived_area = float(result["perimeter_length"]) / float(result["perimeter-area_ratio"]) if float(result["perimeter-area_ratio"]) != 0 else 0
        treecover_derived_area = float(result["treecover"]) * cfg.grid_width**2 * cfg.cell_width**2
        #assert (PAR_derived_area == treecover_derived_area), f"Perimeter-length / Perimeter-area-ratio ({PAR_derived_area}) does not equal Tree-cover * Spatial domain area ({treecover_derived_area})."
        if not cfg.rotate_randomly:
            result["global_rotation_offset"] = str(cfg.global_rotation_offset)
        if control_variable:
            result[control_variable] = control_value
        if secondary_variable:
            result[secondary_variable] = secondary_value
        if dependent_var:
            result[dependent_var] = dependent_val
        if initial_no_dispersals:
            result["initial_number_of_dispersals"] = initial_no_dispersals
        if dependent_result_range_stdev:
            result["dependent_result_range_stdev"] = dependent_result_range_stdev
        writer.writerow(result)
    
    cfg.csv_path = path
    return cfg


def get_lookup_table(species, width, rcg_width):
    path = os.path.join(cfg.DATA_INTERNAL_DIR, f"lookup_table_{species}_width-{width}_rcg_width-{rcg_width}.npy")
    if os.path.exists(path):
        return np.load(path), path
    else:
        return None, path

def export_lookup_table(lookup_table, width, species):
    path = os.path.join(cfg.DATA_INTERNAL_DIR, f"lookup_table_{species}_width-{width}_rcg_width-{lookup_table.shape[0]}.npy")
    print("path: ", path)
    np.save(path, lookup_table)

def import_image(fpath):
    img = cv2.imread(fpath)
    return img

def save_numpy_array_to_file(array, path):
    np.save(path, array)

def save_state(dynamics):
    state_table = dynamics.state.get_state_table()
    time = str(dynamics.time).zfill(4)
    save_numpy_array_to_file(state_table, f"{cfg.DATA_OUT_DIR}/tree_positions/tree_positions_iter_{time}.npy")

def load_tree_dbh_values(_time):
    if os.path.exists(cfg.TREE_DBH_FILE):
        if _time == 1:
            os.remove(cfg.TREE_DBH_FILE)
            return {_time : {}}
    else:
        return {_time : {}}
    
     
    with open (cfg.TREE_DBH_FILE, "r") as dbh_json_file:
        tree_dbh_values = json.load(dbh_json_file)
            
    tree_dbh_values[_time] = {}
    
    return tree_dbh_values
    
def save_tree_dbh_values(tree_dbh_values):
    with open(cfg.TREE_DBH_FILE, "w") as dbh_json_file:
        json.dump(tree_dbh_values, dbh_json_file)

def update_state_report(dynamics):
    state_report_file = f"{cfg.DATA_OUT_DIR}/state_reports/state_report.npy"
    current_state = dynamics.state.get_state_table()
    germ_index = 3 
    death_index = germ_index + 1
    tree_dbh_values = load_tree_dbh_values(dynamics.time)

    if os.path.exists(state_report_file) and dynamics.time == 1:
        os.remove(state_report_file)
    elif os.path.exists(state_report_file):
        state_report = np.load(state_report_file)
        
        # Update the state report with newly germinated trees
        for i in range(current_state.shape[0]):
            _id = int(current_state[i, 0])
            if _id not in state_report[:, 0]:
                # Use tree format: id, x, y, time of germination, time of death
                new_tree = list(current_state[i])[:3] # Get the tree id, x and y coordinates
                new_tree = np.concatenate((new_tree, np.array([dynamics.time, 99999])), axis=0) # Add time of germination and placeholder for time of death
                state_report = np.vstack([state_report, new_tree])
            
        # Update dbh file
        for i in range(current_state.shape[0]):
            _id = int(current_state[i, 0])
            cur_dbh = float(current_state[i][3])
            tree_dbh_values[dynamics.time][_id] = cur_dbh
        
        # Update the state report with tree deaths
        for i in range(state_report.shape[0]):
            _id = int(state_report[i, 0])
            if (_id not in current_state[:, 0]) and (int(state_report[i][death_index]) == 99999):
                state_report[i][death_index] = dynamics.time - 1
    
    if dynamics.time == 1:
        state_report = np.zeros((current_state.shape[0], 5))
        state_report[:,:3] = current_state[:,:3]
        state_report[:,-1] = 99999

    save_numpy_array_to_file(state_report, state_report_file)
    save_tree_dbh_values(tree_dbh_values)




