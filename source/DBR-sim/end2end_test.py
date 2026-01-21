
from calendar import c
from pathlib import Path
from copy import deepcopy
from argparse import ArgumentParser
import cv2
import shutil

import visualization as vis
import file_handling as io
from helpers import *
from config import *
import app

import time



class Test:
    def __init__(self, cfg):
        kwargs = get_all_defaults()
        kwargs = load_json_strings_if_any(kwargs)
        default_cfg = SimpleNamespace(**kwargs)
        self.output_dir = cfg.output_dir
        self.mode = cfg.mode
        self.cfg = self.load(cfg.cfg_file, default_cfg)
        self.name = Path(cfg.cfg_file).stem
        self.current_outputs = {}
        self.expected_outputs = {}
        
    def load(self, cfg_file, default_cfg):
        case_cfg = self.load_casefile(cfg_file)
        cfg = self.apply_case_cfg(case_cfg, default_cfg)
        cfg = self.overwrite_output_dir(cfg, self.output_dir)
        self.check_cfg(cfg)
        cfg.headless = True
        return cfg

    @staticmethod
    def overwrite_output_dir(cfg, output_dir):
        cfg.DATA_OUT_DIR = output_dir
        cfg = derive_output_dirs(cfg)
        return cfg

    @staticmethod
    def check_cfg(cfg):
        assert cfg.random_seed != -999, "End-to-end test cases must specify a random seed."

    @staticmethod
    def apply_case_cfg(case_cfg, default_cfg):
        cfg = deepcopy(default_cfg)
        for key, value in vars(case_cfg).items():
            setattr(cfg, key, value)
        return cfg

    @staticmethod
    def load_casefile(cfg_file):
        with open(cfg_file, "r") as f:
            cfg = json.load(f)
        cfg = SimpleNamespace(**cfg)
        return cfg

    def copy_cover_image(self):
        dst = os.path.join(self.output_dir, "cover.png")
        shutil.copyfile(self.cfg.cover_img_path, dst)

    def run(self):
        _, _, _, _, self.cfg = app.main(**vars(self.cfg))
        self.copy_cover_image()

        return True


class Evaluator:
    def __init__(self, testcase):
        self.testcase = testcase
        self.benchmark_dir = self.get_benchmark_dir()
    
    def get_benchmark_dir(self):
        benchmark_dir = os.path.join(cfg.END2END_TEST_BENCHMARK_DIR, self.testcase)
        io.create_directory_if_not_exists(benchmark_dir)
        return benchmark_dir

    def evaluate(self, test):
        
        return True
 

def evaluate(test, evaluator, args):
    result = evaluator.evaluate(test)
    print(json.dumps({"test_result": result}))


def replace_benchmark(test, evaluator, args):
    # Copy current outputs to expected outputs
    io.remove_dir_contents(evaluator.benchmark_dir)
    print("Copying file tree from {} to {}".format(test.output_dir, evaluator.benchmark_dir))
    io.copy_tree_if_needed(test.output_dir, evaluator.benchmark_dir)
    sys.stdout.flush()


def main():
    parser = ArgumentParser()
    parser.add_argument("--cfg_file", type=str, required=True, help="Path to the test configuration file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory.")
    parser.add_argument("--mode", type=str, required=True, help=(
        "Whether to generate a benchmark (use 'generate_benchmark') or evaluate based on existing benchmark (use 'evaluate').")
    )
    args = parser.parse_args()
    test = Test(args)
    
    evaluator = Evaluator(test.name)
    test.run()
    if args.mode == "evaluate":
        evaluate(test, evaluator, args)
    elif args.mode == "generate_benchmark":
        replace_benchmark(test, evaluator, args)
    else: 
        raise ValueError("Unknown mode: {}".format(args.mode))


if __name__ == "__main__":
    main()
    
