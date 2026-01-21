
from calendar import c
from pathlib import Path
from copy import deepcopy
from argparse import ArgumentParser
import cv2
import shutil
import time
import filecmp
from glob import glob

import visualization as vis
import file_handling as io
from helpers import *
from config import *
from end2end_evaluation import Evaluator
import app


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

    def rename_state_data_file(self):
        src = list(glob(os.path.join(self.output_dir, "state_data", "Simulation*.csv")))[0]
        dst = os.path.join(self.output_dir, "state_data", "Simulation.csv")
        os.rename(src, dst)

    def do_post_test_filemanagement(self):
        self.copy_cover_image()
        self.rename_state_data_file()

    def run(self):
        _, _, _, _, self.cfg = app.main(**vars(self.cfg))
        self.do_post_test_filemanagement()



def replace_benchmark(test, evaluator, args):
    # Replace benchmark with current output
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
    
    evaluator = Evaluator(test.name, cfg)
    test.run()
    if args.mode == "evaluate":
        evaluator.evaluate(test)
    elif args.mode == "generate_benchmark":
        replace_benchmark(test, evaluator, args)
    else: 
        raise ValueError("Unknown mode: {}".format(args.mode))


if __name__ == "__main__":
    main()
    
