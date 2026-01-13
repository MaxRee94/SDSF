
from pathlib import Path
from copy import deepcopy
from argparse import ArgumentParser

import visualization as vis
import file_handling as io
from helpers import *
from config import *
import app



class Test:
    def __init__(self, args):
        kwargs = get_all_defaults()
        kwargs = load_json_strings_if_any(kwargs)
        default_args = SimpleNamespace(**kwargs)
        self.args = self.load(args.cfg_file, default_args)
        self.name = Path(args.cfg_file).stem
        self.current_outputs = {}
        self.expected_outputs = {}
        
    def load(self, cfg_file, default_args):
        case_args = self.load_casefile(cfg_file)
        args = self.apply_case_args(case_args, default_args)
        args.headless = True

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
    

def main():
    parser = ArgumentParser()
    parser.add_argument("--cfg_file", type=str, required=True, help="Path to the test configuration file.")
    args = parser.parse_args()
    test = Test(args)
    result = 1 if test.run() else 0
    
    print(json.dumps({"test_result": result}))


if __name__ == "__main__":
    main()
    
