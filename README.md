## PyllutedWD

Python code for analysis of polluted white dwarfs, available from [GitHub](https://github.com/andrewmbuchan4/PyllutedWD_Public)

As used in:

[Planets or asteroids? A geochemical method to constrain the masses of White Dwarf pollutants](https://academic.oup.com/mnras/article/510/3/3512/6472245)<br />
[Asynchronous accretion can mimic diverse white dwarf pollutants I: core and mantle fragments](https://academic.oup.com/mnras/article/519/2/2646/6827923)<br />
[Asynchronous accretion can mimic diverse white dwarf pollutants II: water content](https://academic.oup.com/mnras/article/519/2/2663/6827922)<br />
[Rapid formation of exoplanetesimals revealed by white dwarfs](https://www.nature.com/articles/s41550-022-01815-8)<br />
[Seven white dwarfs with circumstellar gas discs II: tracing the composition of exoplanetary building blocks](https://academic.oup.com/mnras/article/532/4/3866/7697554)<br />
[White dwarf constraints on geological processes at the population level](https://academic.oup.com/mnras/article/532/2/2705/7701788)<br />
[Host star and exoplanet composition: Polluted white dwarf reveals depletion of moderately refractory elements in planetary material](https://www.aanda.org/articles/aa/full_html/2025/01/aa51621-24/aa51621-24.html)<br />
[Measurements of three exo-planetesimal compositions: a planetary core, a chondritic body, and an icy Kuiper belt analogue](https://academic.oup.com/mnras/article/541/2/1377/8174999)<br />
[White dwarfs as probes of extrasolar planet compositions and fundamental astrophysics](https://doi.org/10.48550/arXiv.2507.03029)

Earlier versions of this code were introduced in:

[Polluted white dwarfs: constraints on the origin and geology of exoplanetary material](https://academic.oup.com/mnras/article/479/3/3814/5046489)<br />
[Bayesian constraints on the origin and geology of exoplanetary material using a population of externally polluted white dwarfs](https://academic.oup.com/mnras/article/504/2/2853/6188376)

## To get started

My workflow to get this running in a new environment was roughly as follows:

Firstly, install additional python modules

```
pip install --user corner
pip install --user xlrd
pip install --user pymultinest
pip install --user xlsxwriter
pip install --user ternary-diagram
pip install --user rpy2
```

Next, install PyMultiNest, build the shared object (.so) file and add it to LD_LIBRARY_PATH

```
git clone https://github.com/JohannesBuchner/MultiNest
cd MultiNest/build/
cmake ..
make
export LD_LIBRARY_PATH=$HOME/path/to/MultiNest/lib:$LD_LIBRARY_PATH (replace path name as appropriate)
```

For any calculations involving thermohaline mixing, the code uses tables from Evan Bauer's DA_Pollution_Tables github repository. Download these to an appropriate location using:

```
git clone https://github.com/evbauer/DA_Pollution_Tables
```

Next, install and load relevant modules (may or may not be necessary, depending on whether you are working with the Environment Modules system).

```
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load SciPy-bundle/2022.05
module load numba/0.56.4-CUDA-11.7.0
module load Tkinter/3.10.4
module load texlive/2020
module load R/4.2.1
```

For convenience, you can save this set of modules using
```
module save PyllutedWD
```
so that in future they can all be loaded in one command:
```
module restore PyllutedWD
```
You may want to add the above command to run_main.sh and run_tests.sh so that you can load modules and run the code in one go.

The synthetic_pipeline script (specifically, the multivariate_tests module it imports) uses the R package 'cramer'.
To install this, assuming R is already installed, run the following command from the command line
```
R
```
And then, within the R interface, run the following command:
```
install.packages("cramer")
```

Finally, update src/configuration.ini with the path you would like to save output to, as well as the location of the DA pollution tables.

## Validation

To check the code works, navigate to the tests directory and run

```
sh run_tests.sh
```

All tests should pass, and the message 'OK' should be printed to the terminal.

## Running the code

To run the main python code, which models white dwarf pollution within a Bayesian framework, navigate to the src directory and run

```
sh run_main.sh
```

To run the synthetic pipeline code, which generates and models synthetic polluted white dwarfs, navigate to the src directory and run

```
python synthetic_pipeline.py
```

## The Bayesian code (main.py)

The entry point is main.py. A typical command line call to main.py can be found in run_main.sh, and looks like this:

```
python main.py configuration.ini
```

### Input

The only command line argument is the name of a configuration file (by default, it looks in configuration.ini), which contains the various parameters and settings. These are:
```
[Paths]
output_dir: Location to save output data
da_pollution_tables_dir: Location of the DA_Pollution_Tables directory
pocomc_dir: Install location of pocoMC (For experimentation - not necessary!)

[Files]
- wd_input_file: csv file containing data from polluted WDs to be modelled, relative to /data/ directory.
    Compatible with PEWDD.
    The abundance columns (e.g., 'log(O/{H,He})') should specify number abundances of the relevant element (O in this case) relative to the dominant atmospheric element (the element specified in 'atmosphere').
    A typical entry will just be a negative number (e.g., '-8.5') but the minus sign can be omitted for convenience.
    The corresponding error column (e.g., 'log(O/{H,He})e') should contain the 1 sigma error estimate on the abundance, or a -1 to indicate that the abundance is an upper bound (asymmetric errors are not yet supported).
    For abundances and errors, indicate no data by leaving the entry blank.
    (Timescale columns (e.g., 't_Al') are ignored in this version of the code).
- stellar_compositions_file: csv file containing stellar compositions, relative to /data/ directory.
    This file contains 11 columns, corresponding to the following (number) abundance ratios, on a linear scale (not log): Al/Mg, Ti/Mg, Ca/Mg, Ni/Mg, Fe/Mg, Cr/Mg, Si/Mg, Na/Mg, O/Mg, C/Mg, N/Mg.
    The rows can be sorted according to the application.
    The default file, StellarCompositionsSortFE.csv, is sorted by Fe abundance, which makes sense if iron core formation is considered to be the main (or one of the main) compositional variables.

[Settings]
- timescale_types: Set(s) of timescale grids to use, comma separated.  
    By default, configuration.ini lists all possible options.
- thermohaline_modes: Either 'True', 'False', or 'True, False', to indicate whether thermohaline mixing should be included, neglected, or both.
- suppress_graphical_output: If True, do not produce any graphical output. If False, graphs will be generated in post-processing.
- live_points: the number of live points that PyMultiNest will use.
    The higher the number, the slower the code will run but the better constrained the Bayesian evidence/posteriors will be - I typically use 2000, or 20 for testing purposes.
- differentiation_model: Prescription for how core-mantle differentiation is calculated.
    In normal usage, this should be either 'Earthlike' or 'NonEarthlike' (other options listed in src/enhancement_model.py).
    I typically use 'NonEarthlike' - this allows the model to use pressure/oxygen fugacity variables to explore a range of compositions.
- pollution_models: Parameter set to be explored.
    By default, 'Hierarchy_Default' is used, which means that multiple combinations of parameters will be explored using a hierarchy of lists of parameters, which is the typical use case.
    Multiple hierarchies exist, see src/model_parameters.py.
    It is also possible to specify multiple individual models (not hierarchies) by listing them one after another.
- seed: Seed used for random number generation in PyMultiNest.
    By default, this is -1 (which means it will use the system clock).
- default_logg: The value of log(g), in cgs units, to assume for white dwarfs if no value is specified.
    By default, this is 8.
- default_ca: The default abundance of Ca to assume for white dwarfs with no Ca detection, in log10 number abundance relative to H/He.
    This is only used when calculating sinking timescales with grids which use this as a parameter.
    The default value is -15.
- verbose: Whether to run PyMultiNest in verbose mode.
- resume: Whether to run PyMultiNest in resume mode.
    If True (default), this means that if a system has already been run (or partially run) for a given set of parameters, the code will continue where it left off, which is a useful timesaver.
    If False, the code will start from scratch, which is useful (and necessary!) if you are rerunning with edited abundance data or something like that.
```

By default, the code will run on all systems specified in the white dwarf data input file. To run on a subset of these systems, change the argument in the call to manager.run() in main.py. The argument should be a list containing integers specifying the row(s) in the white dwarf input file of the systems to run (the first row below the header is row 0).

### Output

The output data folders often have truncated names, to accommodate MultiNest's awkward 100-character path limit. The output directory is structured thus:

```
<output_dir>/
├── pipeline
└── r[esults]
    ├── <system1>
    :   ├── <timescale type abbrev.>
    :   :   └── {n[ot considering thermohaline mixing], t[hermohaline mixing considered]}
        :       ├── c
                |   └── <pymultinest output files>
                ├── <plot>.pdf
                ├── <plot>.pdf.txt
                ├── <plot>.png
                └── ..._stats.csv

```

Outputs are stored in `output_dir/r` based on the `output_dir` specified in the configuration file (see above).
Within `r/`, each system has its own subdirectory.
Within each system's directory, there is a folder corresponding to each of the timescale types specified in the config (see above); the abbreviations used are as follows (and also in `src/timescale_interpolator.py`):
- `3o` : `Bedard3DOvershoot`
- `3p` : `Bedard3DOvershootPatched`
- `bn` : `BedardNoOvershoot`
- `bo` : `BedardOvershoot`
- `kn` : `KoesterNoOvershoot`
- `ko` : `KoesterOvershoot`
- `m`  : `MWDD`
- `vo` : `BedardVariableOvershoot`

Within each of these is a folder for each of the thermohaline modes specified in the config (see above):
`n` if `thermohaline_modes` is `False`;
`t` if `True`;
both if `True, False`.

Within these folders is a folder `c/` (containing `pymultinest` outputs), any generated graphs (in pdf, png, and .pdf.txt forms), and a file called `..._stats.csv`. This csv file contains the important quantitative outputs for this particular system/(timescale type)/(thermohaline mode). Here I briefly summarise the key/potentially unclear outputs in the csv file:
- Near the top is a table listing all the models that were run (`Model` column), and the Bayesian evidence of each (`ln_Z_model` column). Higher (less negative) is better!
    - The third column in this table is called `Good fit?`. This column checks whether the best model is actually able to fit the data well.
- The next table shows sinking timescales for each element.
- Below is a line saying `Results from model:` followed by the name of a model. Below this point, the results refer to this model specifically (until you reach another such line).
    - The percentile values on each parameter are calculated from resampling randomly from the individual posteriors - forward modelling using these median values will not necessarily be the same as the median fit!
- The `Disc Composition` row contains the relative abundances of each element (specified in the `Elements` row) in the disc at the point of formation, i.e., the bulk composition of the pollutant's parent body
- The `Parent Core Number Fraction` entry specifies the predicted core fraction of the pollutant's parent body (based on pressure/oxygen fugacity), not the pollutant itself. If this is extremely low (i.e., << 0.01), it indicates that the differentiation model converged to an unphysical result.
- Similarly, the `Radius /km` and `Mass /M_Earth` entries refer to the parent body
- `delta time/Myrs` is the median value of t - t_event, where t is the time since accretion and t_event is accretion event lifetime. The 3rd and 4th columns are the upper and lower errors on this value (similarly elsewhere)
- The `Build Up %`, `Steady State %` and `Declining %` entries specify the posterior probability of accretion being in each of those phases
- The `Temperature /K` entry specifies the median temperature characterising the extent of volatile depletion (which can be interpreted as the temperature during formation)
- Among the various oxygen excess outputs, the key ones are at the bottom under `Excess Oxygen Semi-Sampling results`
    - `Sigma excess (default)` is the sigma significance of an oxygen excess or deficit, using the default oxidation scheme
    - `Median fractional excess (default)` is the (median value of) the fraction of oxygen which cannot be assigned to metal oxides, using the default oxidation scheme


## The synthetic WD code (`synthetic_pipeline.py`)

The entry point is synthetic_pipeline.py. There are no command line arguments.

Output will be stored in /path/to/output/pipeline where path/to/output is the path you specify as output_dir. Each pipeline setup will have its own subdirectory. Within this subdirectory are further subdirectories for each combination of population + observer + modeller, and graphical output for the pipeline as a whole.

Within the pipeline directory, another directory called popdumps is created. This stores a data dump as a csv for each population + observer + modeller combination. Each dump can be shared by multiple pipelines. Before running a pipeline, the code checks to see if the corresponding popdump file exists, and if so will load it instead of calculating from scratch. If the code has changed and you wish to rerun a combination, you will therefore need to delete any deprecated popdumps.

Several graphs show distributions of certain variables, with one of four labels. These labels are:
- Input: the parameter distribution as sampled directly from the input configuration
- Pollution: the distribution of a proxy element across the population. This is the 'true' pollution
- Observed: the distribution of a proxy element across the population after applying detection thresholds and random noise
- Modelled: the parameter distribution obtained when using the Observed data to retrieve the initial Input distribution

## Editing graphs

Across the whole codebase, any output graphs come with a corresponding .txt dump file describing the graph. You can replot most types of graph from this dump using the dict_plot_from_dump.py script in the utils directory (python dict_plot_from_dump.py /path/to/txt/file/text_file.txt). By editing the .txt file you can therefore edit a graph without rerunning everything

## Software requirements

- Python 3.x
- R (for some statistical tests - not essential)

Non-standard python modules used:

- PyMultiNest (this code written using v3.10)
- numba
- matplotlib
- xlrd
- corner
- xlsxwriter
- ternary-diagram
- rpy2

## Contact

For help, please contact Andy Buchan at andy.buchan@warwick.ac.uk
