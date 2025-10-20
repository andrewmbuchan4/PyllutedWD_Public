#!/usr/bin/env python

import argparse as ap
import configparser
import pickle
from pathlib import Path

import timescale_interpolator as ti


abbreviations = {
    # Necessary to keep file names below 100 characters!
    # (Fixed in v3.11 of multinest, but I use v3.10)
    "Earthlike": "EL",
    "NonEarthlike": "NEL",
    "MantleOnly": "MO",
    "Default": "D",
    "TestPrior": "TP",
    "HighPressure": "HP",
    "LowPressure": "LP",
    "RaisedPressure": "RP",
    "EarthMantle": "EM",
    "Meteorite": "M",
    "NELRevamp": "NELR",
}

hierarchy_abbreviations = {
    # Necessary to keep file names below 100 characters!
    # (Fixed in v3.11 of multinest, but I use v3.10)
    "Hierarchy_Basic": "HB",
    "Hierarchy_Default": "HD",
    "Hierarchy_Test": "HT",
    "Hierarchy_Test2": "HT2",
    "Hierarchy_Test3": "HT3",
    "Hierarchy_OM": "HOM",
    "Hierarchy_OM_reduced": "HOMr",
    "Hierarchy_Earth": "HE",
    "Hierarchy_Meteorite": "HM",
    "Hierarchy_Revamp": "HR",
}

true_values = ["true"]
# ^ A list of values that will be treated as truthy when parsing config, but all in
# lower case (since we will take the lower case of the input value)


# I assume that files will not be moved around!
def get_path_to_parent():
    # Path(__file__) is the path of this file
    # resolve() returns the absolute path
    # parents[1] is the path of the parent directory
    return f"{Path(__file__).resolve().parents[1]}/"


def get_path_to_data():
    return f"{get_path_to_parent()}data/"


def get_path_to_feni():
    return f"{get_path_to_data()}feni/"


def get_path_to_src():
    return f"{get_path_to_parent()}src/"


def get_path_to_default_graphs():
    return f"{get_path_to_src()}graphs/"


def get_path_to_utils():
    return f"{get_path_to_parent()}utils/"


class ConfigWrapper:

    def __init__(self, default_config="configuration.ini"):
        self.config = configparser.ConfigParser()
        self.default_config = default_config

    def get(self, section, variable):
        try:
            return self.config.get(section, variable)
        except configparser.NoSectionError:
            # This could be dangerous if someone tries to read a non existent section
            # and then ends up overwriting all their config settings.
            # Looks like this might be solvable by repeated calls to config.read?
            self.config.read(get_path_to_src() + self.default_config)
            self.config.sections()
            return self.config.get(section, variable)


def read_from_pickle(path_to_file):
    pkl_file = open(path_to_file, "rb")
    data = pickle.load(pkl_file)
    pkl_file.close()
    return data


config_wrapper = ConfigWrapper()


def get_path_to_output_base_dir():
    toret = config_wrapper.get("Paths", "output_dir")
    if not toret.endswith("/"):
        toret += "/"
    return toret


def get_path_to_da_pollution_tables_dir():
    toret = config_wrapper.get("Paths", "da_pollution_tables_dir")
    if not toret.endswith("/"):
        toret += "/"
    assert Path(toret).is_dir(), (
        f"Path to DA pollution tables directory does not exist: {toret}"
        + "\nMost likely you have not set this path in configuration.ini"
    )
    return toret


def get_path_to_pocomc_dir():
    toret = config_wrapper.get("Paths", "pocomc_dir")
    if not toret.endswith("/"):
        toret += "/"
    return toret


def get_timescale_types_to_use():
    toret = config_wrapper.get("Settings", "timescale_types")
    return sorted(
        list(set([ti.TimescaleType[tt_str.strip()] for tt_str in toret.split(",")]))
    )


def get_thermohaline_modes_to_use():
    toret = config_wrapper.get("Settings", "thermohaline_modes")
    return sorted(
        list(
            set([tm_str.strip().lower() in true_values for tm_str in toret.split(",")])
        )
    )


def get_suppress_graphical_output():
    toret = config_wrapper.get("Settings", "suppress_graphical_output")
    return toret.lower() in true_values


def get_live_points():
    toret = config_wrapper.get("Settings", "live_points").strip()
    return int(toret)


def get_default_logg():
    toret = config_wrapper.get("Settings", "default_logg").strip()
    return float(toret)


def get_default_mass():
    toret = config_wrapper.get("Settings", "default_mass").strip()
    return float(toret)


def get_default_ca():
    toret = config_wrapper.get("Settings", "default_ca").strip()
    return float(toret)


def get_verbose():
    toret = config_wrapper.get("Settings", "verbose").strip()
    return bool(toret)


def get_resume():
    toret = config_wrapper.get("Settings", "resume").strip()
    return bool(toret)


def get_seed():
    toret = config_wrapper.get("Settings", "seed").strip()
    if not toret:
        # Happens if string was empty - we catch this at manager level and replace with
        # a default value
        raise AttributeError
    return int(toret)


def get_differentiation_model():
    toret = config_wrapper.get("Settings", "differentiation_model").strip()
    return toret


def get_pollution_models_to_use():
    toret = config_wrapper.get("Settings", "pollution_models")
    return [pm_str.strip() for pm_str in toret.split(",")]


def get_wd_input_file():
    toret = config_wrapper.get("Files", "wd_input_file").strip()
    return toret


def get_stellar_compositions_file():
    toret = config_wrapper.get("Files", "stellar_compositions_file").strip()
    return toret


def get_path_to_pipeline_base_dir():
    return f"{get_path_to_output_base_dir()}pipeline/"


def get_path_to_pylluted_dir():
    # Keep this short to fit in 100 characters
    return f"{get_path_to_output_base_dir()}r/"


def get_path_to_historical_output_dir():
    # This will only be useful if you have the output from Harrison et al. 2021 in this
    # directory - I can send this, or it can be generated by running the PWDCode.py
    # script in an earlier version of the codebase
    return f"{get_path_to_output_base_dir()}output_for_harrison2021/"


def set_up_configuration():
    args = parse_command_line_arguments()
    config_wrapper.config.read(f"{get_path_to_src()}{args.config_file}")
    config_wrapper.config.sections()


def parse_command_line_arguments():
    parser = ap.ArgumentParser(description="Configuration Filename")
    parser.add_argument(
        dest="config_file",
        type=str,
        help="Name of configuration file (will look in "
        + get_path_to_src()
        + " for a file of this name)",
    )
    # parser.add_argument(
    #    dest='stellar_compositions_filename',
    #    type=str,
    #    help='Stellar compositions file name (will look in '
    #    + get_path_to_data()
    #    + ' for a file of this name)'
    # )
    # parser.add_argument(
    #    dest='n_live_points',
    #    type=int,
    #    help='Number of live points to run models with'
    # )
    # parser.add_argument(
    #    dest='enhancement_model',
    #    type=str,
    #    help='Enhancement model name'
    # )
    # parser.add_argument(
    #    dest='base_dir',
    #    type=str,
    #    help='Directory to store output'
    # )
    # parser.add_argument(
    #    '--seed',
    #    default=-1,
    #    dest='seed',
    #    type=int,
    #    help='Seed for random number generation. Should be set to -1 for purposes"
    #        +'other than testing'
    # )
    # parser.add_argument(
    #    dest='pollution_model_names',
    #    type=str,
    #    help='Pollution model names, separated by spaces',
    #    nargs='+'
    # )
    parsed_arguments = parser.parse_args()
    return parsed_arguments


def main():
    print(f"Source directory is: {get_path_to_parent()}")
    print(f"Looking for input data in: {get_path_to_data()}")
    print(f"Looking for partitioning data in: {get_path_to_feni()}")
    print(f"Looking for utility scripts in: {get_path_to_utils()}")
    print(f"Output directory: {get_path_to_output_base_dir()}")
    print(f"Will put synthetic_pipeline output in: {get_path_to_pipeline_base_dir()}")
    print(f"Old data in: {get_path_to_historical_output_dir()}")
    print(f"Thermohaline data in : {get_path_to_da_pollution_tables_dir()}")


if __name__ == "__main__":
    main()
