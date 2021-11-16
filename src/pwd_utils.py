#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse as ap

from pathlib import Path

abbreviations = {  # Necessary to keep file names below 100 characters! (Fixed in v3.11 of Pymultinest, but I use v3.10)
    'Earthlike': 'EL',
    'NonEarthlike': 'NEL',
    'MantleOnly': 'MO',
    'Default': 'D',
    'TestPrior': 'TP',
    'HighPressure': 'HP',
    'LowPressure': 'LP',
    'RaisedPressure': 'RP'
}

# I assume that files will not be moved around!
def get_path_to_parent():
    # Path(__file__) is the path of this file
    # resolve() returns the absolute path
    # parents[1] is the path of the parent directory
    return str(Path(__file__).resolve().parents[1])

def get_path_to_data():
    return get_path_to_parent() + '/data/'
    
def get_path_to_feni():
    return get_path_to_parent() + '/feni_src/feni/'
    
def get_path_to_utils():
    return get_path_to_parent() + '/utils/'
    
def get_path_to_original_src():
    return get_path_to_parent() + '/original_codebase/'

def parse_command_line_arguments():
    parser = ap.ArgumentParser(description='Manager Arguments')
    parser.add_argument(
        dest='wd_data_filename',
        type=str,
        help='White Dwarf abundances, errors and sinking timescales (will look in ' + get_path_to_data() + ' for a file of this name)'
    )
    parser.add_argument(
        dest='stellar_compositions_filename',
        type=str,
        help='Stellar compositions file name (will look in ' + get_path_to_data() + ' for a file of this name)'
    )
    parser.add_argument(
        dest='n_live_points',
        type=int,
        help='Number of live points to run models with'
    )
    parser.add_argument(
        dest='enhancement_model',
        type=str,
        help='Enhancement model name'
    )
    parser.add_argument(
        dest='base_dir',
        type=str,
        help='Directory to store output'
    )
    parser.add_argument(
        '--seed',
        default=-1,
        dest='seed',
        type=int,
        help='Seed for random number generation. Should be set to -1 for purposes other than testing'
    )
    parser.add_argument(
        dest='pollution_model_names',
        type=str,
        help='Pollution model names, separated by spaces',
        nargs='+'
    )
    parsed_arguments = parser.parse_args()
    return parsed_arguments
