#!/usr/bin/env bash
#module restore PyllutedWD
main_dir=$(dirname "$0")
main_path=${main_dir}"/main.py"

python $main_path configuration.ini
