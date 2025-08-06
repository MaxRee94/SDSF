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

import visualization as vis
import file_handling as io
from helpers import *
from config import *
import controllable_pattern_generator as cpg

sys.path.append(BUILD_DIR)
#from x64.Debug import dbr_cpp as cpp
from x64.Release import dbr_cpp as cpp


class Test:
    def __init__(self):
        self.inputs = {}
        self.current_outputs = {}
        self.expected_outputs = {}
     
    def run(self, default_inputs):
        print("Running test...")

        print("Test completed successfully.")


class OutputTests:
    def __init__(self):
        self.tests = []
        self.default_inputs = get_all_defaults()

    def run_all(self):
        print("Running all output tests...")
        for i in enumerate(self.tests):
            print(f"Running test: {i}/{len(self.tests)}")
            self.tests[i].run(self.default_inputs)
            
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

