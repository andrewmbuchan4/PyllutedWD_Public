## PyllutedWD

Python code for analysis of polluted white dwarfs, available from [GitHub](https://github.com/andrewmbuchan4/PyllutedWD_Public)

As used in [Planets or asteroids? A geochemical method to constrain the masses of White Dwarf pollutants](https://arxiv.org/abs/2111.08779)

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

Secondly, install PyMultiNest, build the shared object (.so) file and add it to LD_LIBRARY_PATH

```
git clone https://github.com/JohannesBuchner/MultiNest
cd MultiNest/build/
cmake ..
make
export LD_LIBRARY_PATH=$HOME/path/to/MultiNest/lib:$LD_LIBRARY_PATH (replace path name as appropriate)
```

Thirdly, install and load relevant modules (may or may not be necessary, depending on whether you are working with the Environment Modules system).

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

Finally, update src/configuration.ini with the path you would like to save output to.

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

## Input and output for the Bayesian code (main.py)

The entry point is main.py. A typical command line call to main.py can be found in run_main.sh, and looks like this:

```
python main.py WDInputData.csv StellarCompositionsSortFE.csv 2000 NonEarthlike Hierarchy_Default
```
The first command line argument, WDInputData.csv, specifies the file containing data from polluted white dwarfs to be modelled. The code will look for this file in the /data/ directory. The 'Name' column must not contain whitespace. The 'Type' column should be 'H' or 'He' for H- or He-dominated atmospheres, respectively. The abundance columns (e.g., 'log(Al/Hx)') should specify number abundances of the relevant element (Al in this case) relative to the dominant atmospheric element (the element specified in 'Type'). A typical entry will just be a negative number (e.g., '-8.5') but the minus sign can be omitted for convenience, and you can specify that the abundance is actually an upper bound by prefixing with < (e.g., '<-8.5'). The corresponding error column (e.g., 'error log(Al/Hx)') should contain the 1 sigma error estimate on the abundance. Asymmetric errors are not yet supported. This column is ignored if the abundance was an upper bound. For abundances and errors, indicate no data by leaving the entry blank. The final set of columns (e.g., 't_Al') specify sinking timescales for the corresponding element (Al in this case) in years. Setting one of these timescales to 0 will trigger the code to calculate the corresponding timescale itself. By default, I set all of them to 0.

The second command line argument, StellarCompositionsSortFE.csv, specifies the file containing stellar compositions. The code will look for this file in the /data/ directory. This file contains 11 columns, corresponding to the following (number) abundance ratios, on a linear scale (not log):
Al/Mg, Ti/Mg, Ca/Mg, Ni/Mg, Fe/Mg, Cr/Mg, Si/Mg, Na/Mg, O/Mg, C/Mg, N/Mg
The rows can be sorted according to the application. The file specified here, StellarCompositionsSortFE.csv, is sorted by Fe abundance, which makes sense if iron core formation is considered to be the main (or one of the main) compositional variables.

The third argument is the number of live points that PyMultiNest will use. The higher the number, the slower the code will run but the better constrained the Bayesian evidence/posteriors will be - I typically use 2000, or 20 for testing purposes.

The fourth argument is the 'enhancement model', which essentially sets the prescription for how core--mantle differentiation is calculated, either 'Earthlike' or 'NonEarthlike'. I typically use NonEarthlike - this allows the model to use pressure/oxygen fugacity variables to explore a range of compositions

The fifth argmuent is the name of the pollution model, which sets the parameters to be explored. In this example, 'Hierarchy_Default' is used, which means that multiple combinations of parameters will be explored using a hierarchy of lists of parameters, which is the typical use case. Multiple hierarchies exist, see model_parameters.py. It is also possible to specify multiple models (not hierarchies) in the command line by listing them one after another.

The seed for random number generation can optionally be specified using the --seed flag.

By default, the code will run on all systems specified in the white dwarf data input file. To run on a subset of these systems, change the argument in the call to manager.run() in main.py. The argument should be a list containing integers specifying the row(s) in the white dwarf input file of the systems to run (the first row below the header is row 0).

Output will be stored in /path/to/output/r where path/to/output is the path you specify in configuration.ini. Each system will have its own subdirectory. This subdirectory contains any generated graphs, a further subdirectory called c containing the PyMultiNest output, and a .csv file with a variety of output quantities. Here I briefly summarise the key/potentially unclear outputs in the csv file:
- Near the top is a table listing all the models that were run (column Model), and the Bayesian evidence of each (column ln_Z_model). Higher (less negative) is better!
- The final column in this table is called 'Good fit?'. Use this column to check whether the best model is actually able to fit the data well
- Below the sinking timescales should be a line saying 'Results from model:' followed by the name of a model. Below this point, until you reach another such line, the results refer to this model specifically
- The percentile values on each parameter are calculated from resampling randomly from the individual posteriors - forward modelling using these median values will not necessarily be the same as the median fit!
- The Disc Composition row contains the relative abundances of each element (specified in the Elements row) in the disc at the point of formation, i.e., the bulk composition of the pollutant's parent body
- The Parent Core Number Fraction entry specifies the predicted core fraction of the pollutant's parent body (based on pressure/oxygen fugacity), not the pollutant itself
- Similarly, the Radius and Mass entries refer to the parent body
- delta time is the median value of t - t_event, where t is the time since accretion and t_event is accretion event lifetime. The 3rd and 4th columns are the upper and lower errors on this value (similarly elsewhere)
- The Build Up, Steady State and Declining entries specify the posterior probability of accretion being in each of those phases
- The Temperature entry specifies the median temperature characterising the extent of volatile depletion (which can be interpreted as the temperature during formation)
- Among the various oxygen excess outputs, the key ones are at the bottom under Excess Oxygen Semi-Sampling results
- Sigma excess (default) is the sigma significance of an oxygen excess or deficit, using the default oxidation scheme
- Median fractional excess (default) is the (median value of) the fraction of oxygen which cannot be assigned to metal oxides, using the default oxidation scheme

Note that rerunning the code will skip the PyMultiNest calculation and instead load the output from the previous run (stored in the c subdirectory)
To re-run the PyMultiNest calculation (e.g. if code has been updated), delete the associated output files first

## Input and output for the synthetic white dwarf code (synthetic_pipeline.py)

The entry point is synthetic_pipeline.py. There are no command line arguments.

Output will be stored in /path/to/output/pipeline where path/to/output is the path you specify in configuration.ini. Each pipeline setup will have its own subdirectory. Within this subdirectory are further subdirectories for each combination of population + observer + modeller, and graphical output for the pipeline as a whole.

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
