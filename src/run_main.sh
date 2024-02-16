#!/usr/bin/env bash
#module restore PyllutedWD
main_dir=$(dirname "$0")
main_path=${main_dir}"/main.py"
#python -m cProfile -s cumtime $main_path WDInputData.csv StellarCompositionsSortFE.csv 2000 NonEarthlike Hierarchy_Default
#python $main_path WDInputData.csv SolarComposition.csv 2000 Meteorite Hierarchy_Meteorite
python $main_path WDInputData.csv StellarCompositionsSortFE.csv 2000 NonEarthlike Hierarchy_Default
