"""Parse commandline arguments."""

import sys
import json
from argparse import ArgumentParser
from config import ParameterConfig
import app
from helpers import *


def check_cli_args(**user_args):
    assert (user_args["grid_width"] % user_args["resource_grid_width"] == 0), "Resource grid width must be a divisor of grid width."


def main():
    """Set up and launch DBR-sim  either directly (using cli args) or by first showing the GUI."""
    cli_args_provided = len(sys.argv) > 1
    if cli_args_provided:
        kwargs = parse_args()
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
