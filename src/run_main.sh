#!/usr/bin/env bash
# TODO: Consider using something more sophisticated than $0. Could be overkill
main_dir=$(dirname "$0")
main_path=${main_dir}"/main.py"
#python -m cProfile -s cumtime $main_path BlouinConglomNewTimescales.csv StellarCompositionsSortFE.csv 2000 NonEarthlike <your_filepath_here> Hierarchy_Default
python $main_path BlouinConglomNewTimescales.csv StellarCompositionsSortFE.csv 2000 NonEarthlike <your_filepath_here> Hierarchy_Default
