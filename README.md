## PyllutedWD

Python code for analysis of polluted white dwarfs, available from [GitHub](https://github.com/andrewmbuchan4/PyllutedWD_Public)

As used in [Planets or asteroids? A geochemical method to constrain the masses of White Dwarf pollutants](https://arxiv.org/abs/2111.08779)

## To get started

Before starting, there are a few files which have a placeholder filepath: <your_filepath_here>
These should be replaced appropriately

<your_filepath_here> should be replaced with an output path for the main code in the following files:

- src/multi_hist.py
- src/plot_multi_pressure_histogram.py
- src/run_main.sh
- tests/unit_tests.py

The main code should create an output subdirectory within this path. <your_filepath_here> should be replaced with the path of the output subdirectory in the following files:

- src/plot_retrieval_stats.py
- src/stats_extractor.py

<your_filepath_here> should be replaced with a path for video output (if activated) in the following files:

- src/graph_factory.py

original_codebase/PWDCode.py is also affected, but can be ignored unless you specifically want to run it (in which case, replace it with an output file path)

To check the code works, navigate to the tests directory and run

```
sh run_tests.sh
```

To run the main python code, navigate to the src directory and run

```
sh run_main.sh
```

The python entry point is main.py, which also specifies the system(s) to run the code on according to their 0-indexed position in the white dwarf input file specified in run_main.sh (data/BlouinConglomNewTimescales.csv by default)

Running main.py on one or more systems creates raw output files and graphical output in the designated output path
Rerunning the code will skip the PyMultiNest calculation and instead load the output from the previous run
To re-run the PyMultiNest calculation (e.g. if code has been updated), delete the associated output files in the output/ subdirectory in the designated output path first

## Software requirements

- Python 3.x

Non-standard modules used:

- PyMultiNest (this code written using v3.10)
- Numba
- Matplotlib
- xlrd

## Contact

For help, please contact Andy Buchan at amb237@cam.ac.uk
