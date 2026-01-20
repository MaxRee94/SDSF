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
from copy import deepcopy
import contextlib
import subprocess

import visualization as vis
import file_handling as io
from helpers import *
from config import *
import app



class OutputTests:
    def __init__(self):
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
        stdout, stderr, returncode = self.run_test_process(testfile, output_dir)
        
        self.write_cmdline_output_to_file(stdout, output_dir, "stdout.txt")
        if stderr:
            self.write_cmdline_output_to_file(stderr, output_dir, "stderr.txt")
    
        if returncode == 0:
            success = self.parse_stdout(stdout)
        else:
            success = False
        
        return success

    @staticmethod
    def parse_stdout(stdout):
        success = False
        lines = stdout.splitlines()
        success = None
        for i in range (len(lines), 0, -1):
            line = lines[i-1]
            if "test_result" in line:
                result_json = json.loads(line.strip())
                if result_json["test_result"] == 1:
                    success = True
                else:
                    success = False
                break
        assert success is not None, "Could not parse test result from stdout."
        return success

    @staticmethod
    def run_test_process(testfile, output_dir):
        p = subprocess.run(
            [sys.executable, "end2end_test.py", "--cfg_file", testfile, "--output_dir", output_dir],
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
        if success:
            print("\r" + message.replace("Running", "    Passed test") + "         ")
        else:
            print("\r" + message.replace("Running", "--> FAILED test") + "         ")

    def run_all(self):
        print("Running all end-to-end tests...")
        n_failed = 0
        for i, testfile in enumerate(self.test_files):
            test_name, output_dir, message = self.prepare_test(i, testfile)
            
            success = self.run_test(testfile, test_name, output_dir)
            n_failed += not success
            
            self.produce_test_report(message, success)

        return n_failed, len(self.test_files)



def main():
    print("Starting end-to-end tests...")

    tests = OutputTests()
    n_failed, n_total = tests.run_all()
    
    print(f"\nEnd-to-end tests completed.", "All passed." if n_failed==0 else f"Tests FAILED in {n_failed} / {n_total} cases!")

