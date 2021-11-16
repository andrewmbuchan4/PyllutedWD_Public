# PyllutedWD
Python code for analysis of polluted white dwarfs

Before starting, there are a few files which have a placeholder filepath: <your_filepath_here>
These should be replaced appropriately

Affected files:
src/graph_factory.py
src/multi_hist.py
src/plot_multi_pressure_histogram.py
src/plot_retrieval_stats.py
src/run_main.sh
src/stats_extractor.py
tests/unit_tests.py

original_codebase/PWDCode.py is also affected, but can be ignored unless you specifically want to run it

To get started:

sh run_tests.sh (in the tests directory) runs unit tests. Run this first to check everything works

sh run_main.sh (in the src directory) runs the main python code.
The python entry point is main.py (which also specifies the system(s) to run the code on)
