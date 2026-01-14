from dataclasses import dataclass
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

    @staticmethod
    def get_test_files():
        test_files = []
        for f in glob(cfg.END2END_TESTCASE_DIR + "/*.json"):
            test_files.append(f)

        return test_files

    def run_test(self, testfile, test_name):
        stdout, stderr = self.run_test_process(testfile)
        success = self.parse_stdout(stdout)
        
        self.write_stdout_to_file(stdout, test_name)
        
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
    def run_test_process(testfile):   
        p = subprocess.run(
            [sys.executable, "output_test.py", "--cfg_file", testfile],
            text=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=8
        )
        
        return p.stdout, p.stderr
    
    @staticmethod
    def write_stdout_to_file(stdout, test_name):
        timecode = time.strftime("%Y_%m_%d-%H,%M,%S")
        out_dir = os.path.join(cfg.END2END_TEST_OUTPUT_DIR, "commandline_output", timecode)
        os.makedirs(out_dir)
        with open(os.path.join(out_dir, test_name), "w") as f:
            f.write(stdout)

    def run_all(self):
        print("Running all end-to-end tests...")
        n_failed = 0
        for i, testfile in enumerate(self.test_files):
            test_name = Path(testfile).stem
            message = f"Running {i+1}/{len(self.test_files)}: {test_name}"
            print(message + "...", end="")
            
            success = self.run_test(testfile, test_name)
            n_failed += not success
            
            sys.stdout.flush()
            if success:
                print("\r" + message.replace("Running", ">>  Passed test") + "         ")
            else:
                print("\r" + message.replace("Running", ">>  FAILED test") + "         ")

        return n_failed, len(self.test_files)


def main():
    print("Starting end-to-end tests...")

    tests = OutputTests()
    n_failed, n_total = tests.run_all()
    
    print(f"\nEnd-to-end tests completed.", "All passed." if n_failed==0 else f"Tests FAILED in {n_failed} / {n_total} cases!")

