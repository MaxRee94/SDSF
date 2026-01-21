from dataclasses import dataclass
import re
import sys

import cv2
import numpy as np
import time
import os
from pathlib import Path
from matplotlib import pyplot as plt
import json
import scipy.stats as stats
import random
import math
from glob import glob
import shutil
from copy import deepcopy
import contextlib
import subprocess

import visualization as vis
import file_handling as io
from helpers import *
from config import *
import app



class OutputTests:
    def __init__(self, mode):
        self.mode = "generate_benchmark" if mode == "end2end_reset" else "evaluate"
        self.test_files = self.get_test_files()
        self.main_output_dir = self.create_main_output_dir()

    @staticmethod
    def get_test_files():
        test_files = []
        for f in glob(cfg.END2END_TESTCASE_DIR + "/*.json"):
            test_files.append(f)

        return test_files

    def create_main_output_dir(self):
        _dir = os.path.join(cfg.END2END_TEST_OUTPUT_DIR, time.strftime("%Y_%m_%d-%H,%M,%S"))
        io.create_directory_if_not_exists(_dir)
        
        return _dir

    def run_test(self, testfile, test_name, output_dir):
        stdout, stderr, returncode = self.run_test_process(testfile, output_dir, self.mode)
        
        self.write_cmdline_output_to_file(stdout, output_dir, "stdout.txt")
        if stderr:
            self.write_cmdline_output_to_file(stderr, output_dir, "stderr.txt")
    
        success = returncode == 0
        
        return success

    @staticmethod
    def run_test_process(testfile, output_dir, mode):
        p = subprocess.run(
            [sys.executable, "end2end_test.py", "--cfg_file", testfile, "--output_dir", output_dir, "--mode", mode],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=8
        )
        
        return p.stdout, p.stderr, p.returncode
    
    @staticmethod
    def write_cmdline_output_to_file(stdout, output_dir, filename):
        path = os.path.join(output_dir, filename)
        with open(path, "w") as f:
            f.write(stdout)

    def prepare_test(self, index, testfile):
        test_name = Path(testfile).stem
        output_dir = os.path.join(self.main_output_dir, test_name)
        io.create_output_dirs(output_dir)
        message = f"Running {index+1}/{len(self.test_files)}: {test_name}"
        print(message + "...", end="")        

        return test_name, output_dir, message

    def produce_test_report(self, message, success):
        sys.stdout.flush()
        success_basemessage = "    Passed test" if self.mode == "evaluate" else "    Generated benchmark for"
        failure_basemessage = "--> FAILED test" if self.mode == "evaluate" else "--> FAILED to generate benchmark for"
        if success:
            print("\r" + message.replace("Running", f"{success_basemessage} test") + "         ")
        else:
            print("\r" + message.replace("Running", f"{failure_basemessage} test") + "         ")

    def run_all(self):
        print("{} all end-to-end tests...".format("Running" if self.mode == "evaluate" else "Generating benchmarks for"))
        n_failed = 0
        for i, testfile in enumerate(self.test_files):
            test_name, output_dir, message = self.prepare_test(i, testfile)
            
            success = self.run_test(testfile, test_name, output_dir)
            n_failed += not success
            
            self.produce_test_report(message, success)

        return n_failed, len(self.test_files)


def move_subfolders_to_old(root_dir: str, old_name: str = "old") -> None:
    """
    Move all direct subfolders of `root_dir` into the existing subfolder `old_name`,
    except the `old_name` folder itself.

    Example layout before:
        end2end_test_output/
            run_001/
            run_002/
            old/

    After:
        end2end_test_output/
            old/
                run_001/
                run_002/
    """
    old_dir = os.path.join(root_dir, old_name)

    if not os.path.isdir(old_dir):
        raise FileNotFoundError(f"'old' directory does not exist: {old_dir}")

    for entry in os.listdir(root_dir):
        src = os.path.join(root_dir, entry)

        # Skip non-directories and the 'old' folder itself
        if entry == old_name or not os.path.isdir(src):
            continue

        dst = os.path.join(old_dir, entry)

        # shutil.move handles cross-filesystem moves as well
        shutil.move(src, dst)



def prepare(cfg):
    move_subfolders_to_old(cfg.END2END_TEST_OUTPUT_DIR)


def main(mode):
    print("Starting end-to-end tests...")

    prepare(cfg)
    tests = OutputTests(mode)
    n_failed, n_total = tests.run_all()
    
    if mode == "evaluate":
        print(f"\nEnd-to-end tests completed.", "All passed." if n_failed==0 else f"Tests FAILED in {n_failed} / {n_total} cases!")
    else:
        print(f"\nEnd-to-end benchmark creation complete.", "All succeeded." if n_failed==0 else f"Failure occurred in {n_failed} / {n_total} cases!")

