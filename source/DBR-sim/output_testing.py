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

import visualization as vis
import file_handling as io
from helpers import *
from config import *
import app



class Test:
    def __init__(self, cfg_file, default_args):
        self.args = self.load(cfg_file, default_args)
        self.name = Path(cfg_file).stem
        self.current_outputs = {}
        self.expected_outputs = {}
        
    def load(self, cfg_file, default_args):
        case_args = self.load_casefile(cfg_file)
        args = self.apply_case_args(case_args, default_args)
        args.headless = True

        return args

    @staticmethod
    def apply_case_args(case_args, default_args):
        args = default_args
        for key, value in vars(case_args).items():
            setattr(args, key, value)
        return args

    @staticmethod
    def load_casefile(cfg_file):
        with open(cfg_file, "r") as f:
            args = json.load(f)
        args = SimpleNamespace(**args)
        return args
     
    def run(self):
        app.main(**vars(self.args))
        time.sleep(1)

        return True


class OutputTests:
    def __init__(self):
        kwargs = get_all_defaults()
        kwargs = load_json_strings_if_any(kwargs)
        self.default_args = SimpleNamespace(**kwargs)
        self.tests = self.get_tests()

    def get_tests(self):
        test_files = self.get_test_files()
        tests = []
        for f in test_files:
            test = Test(f, self.default_args)
            tests.append(test)

        return tests

    @staticmethod
    def get_test_files():
        test_files = []
        print(cfg.OUTPUT_TESTCASE_DIR + "/*.json")
        for f in glob(cfg.OUTPUT_TESTCASE_DIR + "/*.json"):
            test_files.append(f)

        return test_files

    def run_all(self):
        print("Running all output tests...")
        for i, test in enumerate(self.tests):
            message = f"- Running test {i+1}/{len(self.tests)}: {test.name}"
            print(message + " ...", end="")
            success = test.run()
            sys.stdout.flush()
            if success:
                print("\r" + message.replace("Running", "Passed"))
            else:
                print("\r" + message.replace("Running", "FAILED"))
            
        print("All tests completed.")


def init_tests():
    print("Initializing tests...")

    tests = OutputTests()
    
    print("Tests initialized")
    
    return tests


def main():
    print("Starting output tests...")

    tests = init_tests()
    tests.run_all()
    
    print("\nOutput tests completed.")

