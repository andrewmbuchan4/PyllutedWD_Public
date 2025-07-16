#!/usr/bin/env bash
test_dir=$(dirname "$0")
test_path=${test_dir}"/unit_tests.py"
python $test_path
