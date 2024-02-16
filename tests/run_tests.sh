#!/usr/bin/env bash
#module restore PyllutedWD
# TODO: Consider using something more sophisticated than $0. Could be overkill
# Also fetching the current directory name seems to not work on Windows!?
test_dir=$(dirname "$0")
test_path=${test_dir}"/unit_tests.py"
python $test_path
