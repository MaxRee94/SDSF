"""Create lookup table for batch simulation folders."""

import json
import os
from argparse import ArgumentParser
from config import *


def get_x_new_batch_folders(x):
    """
    Generate x new unique batch folder paths.

    Parameters
    ----------
    x : int
        Number of new batch folders to generate.

    Returns
    -------
    list
        List of x new batch folder paths.
    """
    batch_no = 1
    batch_folder = cfg.DATA_OUT_DIR + "/state_data/batch_000001"
    new_batch_folders = []
    while os.path.exists(batch_folder):
        batch_no += 1
        batch_folder = batch_folder.split("batch_")[0] + "batch_" + str(batch_no).zfill(6)
        
    for i in range(x):
        batch_folder = batch_folder.split("batch_")[0] + "batch_" + str(batch_no + i).zfill(6)
        new_batch_folders.append(batch_folder)

    return new_batch_folders


def create_batch_lookup_table(number_of_batchfolders=None):
    """
    Create a lookup table mapping batch indices to folder paths.

    Parameters
    ----------
    number_of_batchfolders : int, optional
        Number of batch folders to create. If None, defaults to cfg setting.

    Returns
    -------
    dict
        Dictionary mapping batch indices to folder paths.
    """
    new_batch_folders = get_x_new_batch_folders(number_of_batchfolders)
    batch_lookup_table = {}
    for i, folder in enumerate(new_batch_folders):
        batch_lookup_table[i] = folder

    return batch_lookup_table


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-no', '--number_of_batchfolders', type=int)
    args = parser.parse_args()
    batch_lookup_table = create_batch_lookup_table(**vars(args))
    output_path = os.path.normpath(os.path.join(cfg.DATA_OUT_DIR, "state_data/tmp/batchfolder_lookup_table.json"))
    with open(output_path, "w") as lookup_table_file:
        json.dump(batch_lookup_table, lookup_table_file)
    

