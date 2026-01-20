
from calendar import c
from pathlib import Path
from copy import deepcopy
from argparse import ArgumentParser

import visualization as vis
import file_handling as io
from helpers import *
from config import *
import app

import time



class Test:
    def __init__(self, args):
        kwargs = get_all_defaults()
        kwargs = load_json_strings_if_any(kwargs)
        default_args = SimpleNamespace(**kwargs)
        self.output_dir = args.output_dir
        self.mode = args.mode
        self.args = self.load(args.cfg_file, default_args)
        self.name = Path(args.cfg_file).stem
        self.current_outputs = {}
        self.expected_outputs = {}
        
    def load(self, cfg_file, default_args):
        case_args = self.load_casefile(cfg_file)
        args = self.apply_case_args(case_args, default_args)
        args = self.overwrite_output_dir(args, self.output_dir)
        args.headless = True
        return args

    @staticmethod
    def overwrite_output_dir(args, output_dir):
        args.DATA_OUT_DIR = output_dir
        args = derive_output_dirs(args)
        return args

    @staticmethod
    def apply_case_args(case_args, default_args):
        args = deepcopy(default_args)
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
    
