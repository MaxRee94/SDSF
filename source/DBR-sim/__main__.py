"""Parse commandline arguments."""

import sys
import json
from config import ParameterConfig
import app
from helpers import *



def main():
    """Set up and launch DBR-sim  either directly (using cli args) or by first showing the GUI."""
    cli_args_provided = len(sys.argv) > 1
    if cli_args_provided:
        kwargs = parse_args()
        kwargs = load_json_strings_if_any(kwargs)
        check_cli_args(**kwargs)
        tests = kwargs.pop("tests")
        if tests == "unit_tests":
            import unit_testing as ut
            ut.main(**kwargs)
            return
        elif tests == "end2end" or tests == "end2end_reset":
            import end2end_test_administrator as e2e
            e2e.main(tests)
            return
        else:
            app.main(**kwargs)
    else:
        raise NotImplementedError("No cli args were provided. Trying to fall back to GUI-mode, but GUI has not been implemented (yet).")
    

if __name__ == '__main__':
    main()
