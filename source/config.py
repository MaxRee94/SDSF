"""DBR-sim defaults and constants."""


_parameter_config = {
    "gridsize": {
        "keys": {
            "cli": ["--gridsize", "-gs"]
        },
        "settings": {
            "nargs": "*",
            "type": int,
            "help": (
                "The number of grid cells to be used along one axis of the spatial domain."
            ),
        },
        "default": 100,
    },
}


class ParameterConfig():
    """Wrapper class for _parameter_config."""

    def __init__(self):
        """Initialize and load data from _parameter_config."""
        self.data = {}
        self.load()

    def load(self):
        """Update `self.data` by loading parameter config data.

        Some parameters have a function stored as a default. We execute the stored
        function to get the concrete default value to be used. This is done for
        parameters with defaults that are project-specific. For example, the default
        for the parameter `template` is the function 'get_project_publish_template'.
        When executed, this function returns the publish template defined in the
        project entity (a string), which is then used to overwrite the default in the
        config data.
        """
        for name, cfg in _parameter_config.copy().items():
            _cfg = cfg.copy()
            default = _cfg.get("default")
            if default and callable(default):
                # Execute the stored function to get the concrete default specific to
                # the resolved project.
                _cfg["default"] = default()

            self.data[name] = _cfg

    def __iter__(self):
        """Yield all parameter names and -configurations.

        Yields:
            name(str): Name of the parameter.
            cfg(dict): Contains the configuration of the parameter.
        """
        for name, cfg in self.data.items():
            yield name, cfg

    def get(self, name):
        """Return the configuration associated with the given name of the parameter.

        Args:
            name(str): The name of the parameter to return the configuration for.
        Returns:
            (dict): The configuration associated with the name of the parameter.
        """
        return self.data[name]
