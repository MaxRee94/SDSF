import numpy as np
import sys
import json
import argparse
from config import ParameterConfig
import logging
import itertools
from collections import deque
import heapq
import warnings


def get_stdout_logging_handler():
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)s %(levelname)s:    %(message)s')
    handler.setFormatter(formatter)

    return handler


logger = logging.getLogger(__name__)
handler = get_stdout_logging_handler()
logger.addHandler(handler)


class Processes:
    def __init__(self):
        self.procs = []
        self.active_proc_count = 0

    def add(self, proc):
        self.procs.append(proc)
        self.procs[-1].start()
        self.active_proc_count += 1

    def __getitem__(self, idx):
        return self.procs[idx]

    def join(self):
        for proc in self.procs:
            proc.join()
        return True

    def finished(self, print_progress=False):
        _active_proc_count = sum(proc.is_alive() for proc in self.procs)
        if print_progress and _active_proc_count != self.active_proc_count:
            logger.info(f"    {len(self.procs) - _active_proc_count} out of {len(self.procs)} processes have finished.")
        self.active_proc_count = _active_proc_count
        return self.active_proc_count == 0


def suppress_warning(category=None, message=None):
    assert category is not None and message is not None, "Both category and warning must be provided to suppress a warning."
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=category, message=message)


def midpoint_gap_indices(cfg):
    """
    Yield indices in an order where each next index is the midpoint of the
    currently largest remaining gap between already-sampled indices.

    Starts with endpoints (0 and n-1), then repeatedly inserts midpoints of
    the largest unsplit interval.
    """
    if cfg.n <= 0:
        return
    if cfg.n == 1:
        yield 0
        return

    # yield endpoints first
    yield 0
    yield cfg.n - 1

    # max-heap of intervals: (-length, lo, hi)
    # We only split intervals that still contain unsampled indices (hi - lo >= 2).
    heap = []
    if cfg.n - 1 - 0 >= 2:
        heapq.heappush(heap, (-(cfg.n - 1 - 0), 0, cfg.n - 1))

    while heap:
        _, lo, hi = heapq.heappop(heap)
        if hi - lo < 2:
            continue

        mid = (lo + hi) // 2
        yield mid

        # push left and right sub-intervals if they still have room for new points
        if mid - lo >= 2:
            heapq.heappush(heap, (-(mid - lo), lo, mid))
        if hi - mid >= 2:
            heapq.heappush(heap, (-(hi - mid), mid, hi))


def get_random_int_generator(cfg):
    while True:
        yield cfg.rng.randint(cfg.low, cfg.high)


def digits_after_decimal(num):
    num = str(num)
    if '.' in num:
        return len(num.split('.', 1)[1])
    return 0


def add_kwargs(parser):
    """Dynamically add arguments using `ParameterConfig` defined in `config.py`.

    Arguments are dynamically loaded from a single config to reduce code duplication
    between cli, api and gui. Please view `config.py` > `ParameterConfig`
    for reference of all arguments and their configurations.

    Args:
        parser(argparse.ArgumentParser instance): The argument parser to add the arguments to.
    """
    for _, arg_config in ParameterConfig():
        keys = arg_config["keys"]["cli"]
        settings = arg_config["settings"]
        parser.add_argument(*keys, **settings)


def overwrite_from_global_arguments(args, global_args):
    """Overwrite arguments from global arguments.

    Args:
        args (dict): Dictionary of arguments to overwrite.
        global_args (dict): Dictionary of  global arguments.
    """
    for key, value in args.items():
        if global_args.get(key):
            args[key] = global_args[key]
    return args


def inject_unknown_args(args, unknown_args):
    # Process unknown arguments: assume form --key value
    it = iter(unknown_args)
    for token in it:
        if token.startswith("--"):
            key = token.lstrip("-")
            try:
                value = next(it)
                if value.startswith("--"):  # standalone flag
                    raise ValueError(f"Standalone flag detected ({value}, for key {key}) where value expected.")
                else:
                    if value.isnumeric():
                        value = int(value)
                    else:
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                    setattr(args, key, value)
            except StopIteration:
                setattr(args, key, True)
        else:
            raise ValueError(f"Unexpected token {token} in unknown arguments.")

    return args


def parse_args(parser=None):
    """Parse commandline arguments.

    Returns:
        kwargs(dict): Dictionary containing the cli args and their values.
    """
    if parser is None:
        # Create a new ArgumentParser instance if none is provided
        parser = argparse.ArgumentParser()
    
    add_kwargs(parser)
    args, unknown_args = parser.parse_known_args()

    # Handle unknown arguments (if any)
    if args.allow_unknown_args and unknown_args:
        kwargs = inject_unknown_args(args, unknown_args)
    else:
        # Re-start parser without unknown args allowed to raise error
        _ = parser.parse_args()
    
    kwargs = vars(args)
        
    return kwargs


def load_json_strings_if_any(kwargs):
    # Convert JSON strings to dicts
    for k, v in kwargs.items():
        if type(v) == str and "{" in v:
            kwargs[k] = json.loads(v) 
    return kwargs


class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)


def get_max(val1, val2):
    if val1 > val2:
        return val1
    else:
        return val2
    

def get_2d_dist(p1, p2):
    return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def apply_user_args_to_configuration(args, cfg):
    for key, value in vars(args).items():
        setattr(cfg, key, value)
    return cfg


def is_number(inputString):
    return all(char.isdigit() for char in inputString)
