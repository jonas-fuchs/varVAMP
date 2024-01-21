"""
Imports configuration modules (default_config and custom modules).
Represents the imported configurations in a VarVAMPConfig object, does some
basic sanity checks along the way, and merges the default configuration with
custom configuration provided through the $VARVAMP_CONFIG environment variable.

The final VarVAMPConfig instance is what is used as "the" config throughout
varVAMP.
"""

import importlib.machinery
import importlib.util
import os

from varvamp.scripts import default_config


class VarVAMPConfig():
    """
    Parses, validates and stores imported config data.
    """

    __all__ = default_config.__all__

    def __init__(self, config, strict=True):
        defined_attrs = [attr for attr in dir(config) if attr[0] != '_']
        for attr in defined_attrs:
            if attr not in self.__all__:
                # either this is a corrupted default_config,
                # or a custom_config with a wrong param name,
                # or a wrong file was passed as a custom_config by the user.
                raise AttributeError(
                    f'{config.__name__} lists an unknown parameter: "{attr}".'
                    'Please check your VarVAMP config file!'
                )
        for config_option in self.__all__:
            if strict or hasattr(config, config_option):
                # for a custom_config strict will be set to False which will
                setattr(self, config_option, getattr(config, config_option))

    def update(self, other_config):
        for config_option in self.__all__:
            # let any param defined in a custom config overwrite the definition
            # from the default_config
            if hasattr(other_config, config_option):
                setattr(
                    self, config_option, getattr(other_config, config_option)
                )


def load_custom_config(custom_config_path):
    """
    Imports a custom config from a file path and returns its contents as a
    VarVAMPConfig object.
    """

    loader = importlib.machinery.SourceFileLoader(
        'custom_config', custom_config_path
    )
    spec = importlib.util.spec_from_loader('custom_config', loader)
    custom_config = importlib.util.module_from_spec(spec)
    loader.exec_module(custom_config)
    return VarVAMPConfig(custom_config, strict=False)


config = VarVAMPConfig(default_config)

custom_config_file = os.getenv('VARVAMP_CONFIG')
if custom_config_file:
    # try to:
    # 1. turn the value of $VARVAMP_CONFIG into an absolute path
    # 2. import from that path
    # 3. update the default configuration with the contents of the custom config
    config.update(
        load_custom_config(
            os.path.normpath(os.path.join(os.getcwd(), custom_config_file))
        )
    )
