"""Parse commandline arguments."""

import sys
import json
from argparse import ArgumentParser
from config import ParameterConfig
import app


def add_kwargs(parser):
    """Dynamically add arguments using `ParameterConfig` defined in `config.py`.

    Arguments are dynamically loaded from a single config to reduce code duplication
    between cli, api and gui. Please view `config.py` > `ParameterConfig`
    for reference of all arguments and their configurations.

    Args:
        parser(ArgumentParser instance): The argument parser to add the arguments to.
    """
    for _, arg_config in ParameterConfig():
        keys = arg_config["keys"]["cli"]
        settings = arg_config["settings"]
        parser.add_argument(*keys, **settings)


def parse_cli_args():
    """Parse commandline arguments.

    Returns:
        kwargs(dict): Dictionary containing the cli args and their values.
    """
    parser = ArgumentParser()
    add_kwargs(parser)
    kwargs = vars(parser.parse_args())
    
    # Convert JSON strings to dicts
    for k, v in kwargs.items():
        if type(v) == str and "{" in v:
            kwargs[k] = json.loads(v) 
        
    return kwargs


def check_cli_args(**user_args):
    assert (user_args["grid_width"] % user_args["resource_grid_width"] == 0), "Resource grid width must be a divisor of grid width."


def main():
    """Set up and launch DBR-sim  either directly (using cli args) or by first showing the GUI."""
    cli_args_provided = len(sys.argv) > 1
    if cli_args_provided:
        kwargs = parse_cli_args()
        check_cli_args(**kwargs)
        tests = kwargs.pop("tests")
        if tests == "unit_tests":
            import unit_testing as ut
            ut.main(**kwargs)
            return
        elif tests == "output_tests":
            import output_testing as ot
            ot.main()
            return
        else:
            app.main(**kwargs)
    else:
        raise NotImplementedError("No cli args were provided. Trying to fall back to GUI-mode, but GUI has not been implemented (yet).")
    

if __name__ == '__main__':
    main()
