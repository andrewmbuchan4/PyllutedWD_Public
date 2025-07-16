#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import csv
import itertools
import numpy as np
import os
import warnings

import chemistry_info as ci
import live_data as ld
import model_analyser as ma
import model_parameters as mp
import pollution_model as pm
import pwd_utils as pu
import stellar_composition as sc
import timescale_interpolator as ti
import white_dwarf as wd

class Manager:

    def __init__(self):
        self.wd_data_filename = pu.get_wd_input_file()
        self.stellar_compositions_filename = pu.get_stellar_compositions_file()
        self.enhancement_model = pu.get_differentiation_model()
        try:
            self.seed = pu.get_seed()
        except AttributeError:
            # Then seed was not supplied: assume -1
            self.seed = -1
        pollution_model_names = pu.get_pollution_models_to_use()
        if pollution_model_names is not None and len(pollution_model_names) == 1 and pollution_model_names[0].startswith('Hierarchy'):
            # Then we're running a hierarchy which will dynamically create models to run
            self.model_names = None
            self.use_hierarchy = True
            self.hierarchy_name = pollution_model_names[0]
            try:
                self.parameter_hierarchy = mp.hierarchy_definitions_dict[self.hierarchy_name]
            except KeyError as e:
                print('Error: Did not recognise hierarchy ' + pollution_model_names[0])
                print('Available hierarchies: ' + ', '.join(mp.hierarchy_definitions_dict.keys()))
                raise
        else:
            # Then we're running a specified batch of models
            self.model_names = pollution_model_names
            self.use_hierarchy = False
            self.parameter_hierarchy = None
            self.hierarchy_name = None

        self.timescale_interpolator = ti.TimescaleInterpolator()

        self.timescale_types_to_run = pu.get_timescale_types_to_use()
        self.thermohaline_regimes_to_run = pu.get_thermohaline_modes_to_use()
        self.suppress_graphical_output = pu.get_suppress_graphical_output()

        self.white_dwarfs = list()
        self.default_mass = pu.get_default_mass()
        self.default_logg = pu.get_default_logg()
        self.default_Ca_value = pu.get_default_ca()
        self.verbose = pu.get_verbose()
        self.resume = pu.get_resume()
        self.n_live_points = pu.get_live_points()
        self.stellar_compositions = None
        self.models = dict()

        # Think this next bit should be removed - the concept of a comparison prior actually makes no sense
        self.default_prior = 'Default'  # Use this prior except for running certain comparison priors at the end of a hierarchy execution. This placement isn't ideal since prior_functions.py also has this as its default prior
        self.comparison_priors = list()# ['HighPressure', 'LowPressure']

        self.load_global_data()
        if not self.use_hierarchy:
            # This is a bit ugly: we should maybe just do the logic of this function in self.run_system?
            self.load_models()
        # I assume there will only ever be one of these:
        self.analyser = ma.ModelAnalyser(None)
        self.executed_models = dict()

    def get_output_dir(self, wd_name, create=True):
        dirname = pu.get_path_to_pylluted_dir() + wd_name + '/'
        if create:
            self.create_dir(dirname)
        return dirname

    def get_timescale_output_dir(self, wd_name, timescale_type, consider_thermohaline=False, create=True):
        ct_string = 't' if consider_thermohaline else 'n'
        dirname = self.get_output_dir(wd_name, create) + timescale_type.short_str() + '/' + ct_string + '/'
        return dirname

    def get_chains_dir(self, wd_name, timescale_type, consider_thermohaline=False, create=True):
        dirname = self.get_timescale_output_dir(wd_name, timescale_type, consider_thermohaline, create) + 'c/' # Choosing a short name because of the 100 character limit
        if create:
            self.create_dir(dirname)
        return dirname

    def create_dir(self, dir_to_make):
        try:
            os.makedirs(dir_to_make)
        except FileExistsError:
            pass

    def get_white_dwarf_by_name(self, wd_name):
        for i, wd in enumerate(self.white_dwarfs):
            if wd.full_name() == wd_name:
                return i, wd
        print('Warning! Did not recognise white dwarf: ' + wd_name)

    def load_global_data(self):
        self.load_wd_data(self.wd_data_filename)
        self.load_compositions(self.stellar_compositions_filename)

    def load_wd_data(self, wd_data_filename):
        if wd_data_filename is None:
            return
        self.white_dwarfs = list()
        potential_elements_to_ignore = [ci.Element.C, ci.Element.N]
        subsolar_limits = {
            ci.Element.C: collections.OrderedDict({
                ci.Element.Mg: 0.95, # Mg and Si are very similar so just setting them the same
                ci.Element.Si: 0.95
            }),
            ci.Element.N: collections.OrderedDict({
                ci.Element.Mg: 0.35,
                ci.Element.Si: 0.35
            })
        }
        with open(pu.get_path_to_data() + wd_data_filename, encoding='utf-8') as wdcsv:
            for row in csv.DictReader(wdcsv):
                wd_name_raw = row['Star']
                wd_name = "".join(wd_name_raw.split())
                try:
                    wd_name_tag = row['Variant']
                except KeyError:
                    wd_name_tag = ''
                try:
                    wd_name_fullname = row['Identifier']
                except KeyError:
                    wd_name_fullname = ''
                if wd_name_fullname is not None and wd_name_fullname != '' and (wd_name_tag is None or wd_name_tag == ''):
                    try:
                        wd_name_identifier = wd_name_fullname.split(wd_name_raw)[1]
                    except IndexError:
                        print()
                        print(wd_name_raw)
                        print(wd_name_fullname)
                        wd_name_identifier = wd_name_fullname
                    wd_name_tag = "".join(wd_name_identifier.split())
                wd_name_tuple = (wd_name, wd_name_tag)
                wd_type_raw = row['atmosphere']
                try:
                    wd_type = ci.Element[wd_type_raw]
                except KeyError:
                    raise KeyError('Invalid atmospheric type: must be H or He, received ' + str(wd_type_raw))

                wd_abundance_data_raw = dict()

                for el in ci.metals: # i.e., not H and He
                    abundance_column = 'log(' + str(el) + '/H(e))'
                    error_column = 'log(' + str(el) + '/H(e))e'
                    try:
                        raw_error = row[error_column]
                        error = float(raw_error)
                    except KeyError:
                        error = None # occurs if this element is just absent
                    except ValueError:
                    # Occurs if error missing, e.g. if its an upper bound - so assume 0
                        error = 0
                    pewdd_ub = False
                    if error is not None and error < 0:
                        # In PEWDD, this is an indication that the data type is an upper bound
                        pewdd_ub = True
                        error = 0
                    try:
                        raw_value = row[abundance_column]
                        value = float(raw_value)
                        if pewdd_ub:
                            data_point_type = wd.WhiteDwarfDataPointType.upper_bound
                        else:
                            data_point_type = wd.WhiteDwarfDataPointType.measurement
                    except KeyError:
                        value = None # occurs if this element is just absent from the table
                    except ValueError:
                        # Occurs if, for example, this is an upper bound so the value starts with <
                        if row[abundance_column].startswith('<'):
                            # Then this is an upper bound
                            data_point_type = wd.WhiteDwarfDataPointType.upper_bound
                            value = float(row[abundance_column][1:])
                        elif row[abundance_column].startswith('>'):
                            # Then this is an lower bound
                            data_point_type = wd.WhiteDwarfDataPointType.lower_bound
                            value = float(row[abundance_column][1:])
                        elif row[abundance_column] == '':
                            value = None
                        else:
                            # Then we don't know what it is
                            raise ValueError('Unable to parse abundances for ' + wd_name)
                    included = True
                    if value is not None:
                        wd_abundance_data_raw[el] = wd.WhiteDwarfDataPoint(data_point_type, value, error, included)
                abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
                # Remove C and/or N if they are significantly sub-solar
                # There is a bit of an issue here: in principle the 'sub-solar' description refers to the composition
                # at the point of condensation. The numbers here are at the point of observation. Could be different.
                # No obvious fix without assuming the answer, other than perhaps making the ref_limits more generous
                for el in ci.usual_elements:
                    if el in potential_elements_to_ignore:
                        peti_abundance = abundance_data.get_abundance(el)
                        if peti_abundance is not None:
                            auto_ignore = True
                            allow = True
                            for ref_el, ref_limit in subsolar_limits[el].items():
                                ref_abundance = abundance_data.get_abundance(ref_el)
                                if ref_abundance is not None:
                                    auto_ignore = False
                                    ratio = peti_abundance.value - ref_abundance.value
                                    if ratio < ref_limit:
                                        allow = False
                            if auto_ignore or (not allow):
                                abundance_data.abundance_data_dict[el].included = False
                try:
                    wd_mass = float(row['mass'])
                except ValueError:
                    wd_mass = self.default_mass
                try:
                    wd_mass_error = float(row['mass_err'])
                except ValueError:
                    wd_mass_error = 0
                # If any timescales (or logq) are 0, then we'll try to calculate them ourselves:
                try:
                    Teff = float(row['T_eff'])
                except ValueError:
                    Teff = None # It's one thing to have a default logg but a default Teff doesn't make much sense to me
                try:
                    Teff_error = int(row['T_eff_err'])
                except ValueError:
                    Teff_error = 0
                try:
                    logg = float(row['logg'])
                except ValueError:
                    logg = self.default_logg
                try:
                    logg_error = float(row['logg_err'])
                except ValueError:
                    logg_error = 0
                wd_property_data_raw = {
                    mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, wd_type),
                    #mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, Teff, Teff_error),
                    mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, logg, logg_error),
                    mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, wd_mass, wd_mass_error)
                }
                if Teff is not None:
                    wd_property_data_raw[mp.WDParameter.temperature] = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, Teff, Teff_error)
                    if wd_type == ci.Element.He and abundance_data.get_abundance(ci.Element.Ca) is None:
                        # Then we need to intervene: technically the timescale_interpolator won't be able to calculate a timescale
                        # So we'll just use a very small amount of Ca
                        timescales = self.timescale_interpolator.get_wd_timescales(wd_type, logg, Teff, self.default_Ca_value)
                    else:
                        if wd_type == ci.Element.He:
                            timescales = self.timescale_interpolator.get_wd_timescales(wd_type, logg, Teff, abundance_data.get_abundance(ci.Element.Ca).value)
                        else:
                            timescales = self.timescale_interpolator.get_wd_timescales(wd_type, logg, Teff)
                property_data = wd.WhiteDwarfPropertyData(wd_property_data_raw)
                white_dwarf = wd.WhiteDwarf(wd_name_tuple, property_data, abundance_data, timescales)
                self.white_dwarfs.append(white_dwarf)
                try:
                    # If we have an ID column, we can reassure ourselves that the index in our list of white dwarfs will be what we expect it to be
                    assert self.white_dwarfs.index(white_dwarf) == int(row['ID'])
                except KeyError:
                    pass

    def load_compositions(self, input_filename=None):
        if (self.stellar_compositions is not None) and (input_filename is not None):
            warnings.warn('Warning! Attempting to reload compositions')
        if input_filename is not None:
            self.stellar_compositions_filename = input_filename
        if self.stellar_compositions_filename is None:
            return
        self.stellar_compositions = self.load_generic_float_data_csv(self.stellar_compositions_filename)

    def load_generic_float_data_csv(self, input_filename):
        with open(pu.get_path_to_data() + input_filename, encoding='utf-8') as generic_csv:
            generic_list = [row for row in csv.reader(generic_csv)]
            generic_array = np.asarray(generic_list)
        return generic_array.astype(np.float)

    def publish_live_data(self, N_wd, timescale_type):
        elements_to_model = ci.usual_elements
        all_wd_abundances = list()
        all_wd_errors = list()
        all_wd_timescales = list()
        elements_present = list()

        white_dwarf = self.white_dwarfs[N_wd]

        wd_name = white_dwarf.full_name()
        all_wd_timescales = white_dwarf.get_timescales_as_array(timescale_type, elements_to_model)
        if all_wd_timescales is None:
            print('No timescales of type ' + str(timescale_type) + ' for ' + wd_name + ', will skip')
            return

        ld._live_white_dwarf = white_dwarf
        abundances, errors, upper_bounds, lower_bounds = white_dwarf.get_abundance_arrays(elements_to_model)
        ld._live_all_wd_abundances = abundances
        ld._live_all_wd_errors = errors
        ld._live_all_wd_timescales = all_wd_timescales
        ld._live_upper_bounds = upper_bounds
        ld._live_lower_bounds = lower_bounds
        ld._live_timescale_type = timescale_type

        ld._live_t_mg = all_wd_timescales[elements_to_model.index(ci.Element.Mg)]
        ld._live_stellar_compositions = self.stellar_compositions
        ld._live_q = white_dwarf.get_logq(timescale_type).value
        ld._live_mass = white_dwarf.get_mass().value
        ld._live_type = white_dwarf.get_atmospheric_type().value
        ld._live_elements_present = white_dwarf.get_elements_present()
        ld._live_Hx = white_dwarf.get_atmospheric_type().value
        ld._live_M_cvz = white_dwarf.get_logq_in_solar_masses(ld._live_timescale_type)
        ld._live_teff = white_dwarf.get_teff().value
        ld._live_logg = white_dwarf.get_logg().value

    def publish_live_model(self, model_name, consider_thermohaline, prior_name='Default'):
        ld._live_model = model_name
        ld._live_prior = prior_name
        ld._live_enhancement_model = self.enhancement_model
        ld._live_consider_thermohaline = consider_thermohaline

    def load_models(self, model_names=None):
        # If we're in this function, then we're going to run a predetermined set of models for each timescale type,
        # so we're going to set self.models to be a dict where the keys are the timescale types/therm flags and the values are the predetermined set
        if model_names is not None:
            self.model_names = model_names
        for timescale_type in self.timescale_types_to_run:
            self.models[timescale_type] = collections.OrderedDict()
            for consider_thermohaline in self.thermohaline_regimes_to_run:
                self.models[timescale_type][consider_thermohaline] = collections.OrderedDict()
                if self.model_names is not None:
                    for model_name in self.model_names:
                        self.models[timescale_type][consider_thermohaline][model_name] = pm.PollutionModel(
                            model_name,
                            timescale_type,
                            self.enhancement_model,
                            self.default_prior,
                            consider_thermohaline,
                            self.n_live_points,
                            self.seed,
                            self.verbose,
                            self.resume
                        )

    def execute_models(self, N_wd, timescale_type, consider_thermohaline):
        white_dwarf = self.white_dwarfs[N_wd]
        if self.models.get(timescale_type, dict()).get(consider_thermohaline) is not None:
            self.publish_live_data(N_wd, timescale_type)
            if ld._live_all_wd_timescales is None:
                print('No timescales of type ' + str(timescale_type) + ' for the following system, will skip')
                print(white_dwarf)
            max_filename_overflow = self.get_max_filename_overflow(white_dwarf, timescale_type, consider_thermohaline)
            for model_name, model in self.models[timescale_type][consider_thermohaline].items():
                if model_name not in self.executed_models.get((N_wd, timescale_type, consider_thermohaline), list()):
                    # Then we haven't already run this one, so must do so. Otherwise, we can ignore.
                    self.publish_live_model(model_name, consider_thermohaline, model.prior_name)
                    print('About to run model ' + str(model_name) + ' with therm: ' + str(consider_thermohaline) + ' and with prior ' + str(model.prior_name))
                    practical_output_dir = self.get_chains_dir(white_dwarf.abbreviated_name(max_filename_overflow), timescale_type, consider_thermohaline)
                    model.execute(practical_output_dir)
                    if (N_wd, timescale_type, consider_thermohaline) not in self.executed_models.keys():
                        self.executed_models[(N_wd, timescale_type, consider_thermohaline)] = list()
                    self.executed_models[(N_wd, timescale_type, consider_thermohaline)].append(model_name)

    def get_max_filename_overflow(self, white_dwarf, timescale_type, consider_thermohaline):
        max_filename_overflow = 0
        # NB: The following logic is a bit fragile: It only works if assuming that all other timescale items/thermohaline flags that
        # we're going to run have the same length as this current one (otherwise we'll find that one of those other combinations
        # may have a different length, so will get abbreviated differently, and end up in a different output directory)
        for model_name, model in self.models[timescale_type][consider_thermohaline].items():
            ideal_output_dir = self.get_chains_dir(white_dwarf.full_name(), timescale_type, consider_thermohaline, False)
            filename_overflow = model.get_file_name_overflow(ideal_output_dir)
            if filename_overflow > max_filename_overflow:
                max_filename_overflow = filename_overflow
        return max_filename_overflow

    def compare(self, N_wd, timescale_type, consider_thermohaline):
        wd_name = self.white_dwarfs[N_wd].full_name()
        print('Comparing models for ' + wd_name + ' (' + str(timescale_type) + ', Thermohaline = ' + str(consider_thermohaline) + ')')
        self.publish_live_data(N_wd, timescale_type)
        number_of_data_points = len(self.white_dwarfs[N_wd].get_elements_present())
        self.analyser.compare(self.models[timescale_type][consider_thermohaline], N_wd, number_of_data_points)
        max_ln_Z = None
        max_ln_Z_name = None
        base_model_name = None

        for model_name, model in self.models[timescale_type][consider_thermohaline].items():  # should be an OrderedDict with the base model in first position
            if base_model_name is None:
                base_model_name = model_name
                max_ln_Z = model.comparison[base_model_name]['ln_Z_model']
                max_ln_Z_name = model_name
            else:
                if model.comparison[max_ln_Z_name]['ln_Z_model'] > max_ln_Z:
                    max_ln_Z = model.comparison[max_ln_Z_name]['ln_Z_model']
                    max_ln_Z_name = model_name
        print('Best model was ' + max_ln_Z_name)
        for model_name, model in self.models[timescale_type][consider_thermohaline].items():
            if model_name == max_ln_Z_name:
                model.best_model = True
            else:
                model.best_model = False
        return max_ln_Z_name

    def run(self, systems_to_run=None):
        if systems_to_run is None:  # Run the whole input file by default
            systems_to_run = range(len(self.white_dwarfs))
        for system in systems_to_run:
            self.run_system(system)

    def build_model(self, hierarchy_level_list):
        # Make a pollution model consisting of the parameters in the relevent hierarchy levels of the parameter hierarchy
        model_name = pu.hierarchy_abbreviations[self.hierarchy_name]
        parameters_to_use = list()
        for hl in hierarchy_level_list:
            parameters_to_use += self.parameter_hierarchy[hl]
            model_name += str(hl)
        mp.model_definitions_dict[model_name] = dict()
        for potential_param in mp.ModelParameter:
            mp.model_definitions_dict[model_name][potential_param] = potential_param in parameters_to_use
        return model_name

    def register_models_with_altered_priors(self, best_model_name, wd_name, timescale_type, consider_thermohaline):
        toret = list()
        for prior in self.comparison_priors:
            try:
                new_model_name = best_model_name + '_' + pu.abbreviations[prior]
            except KeyError as e:
                message = 'Warning! Could not run using unrecognised prior ' + prior + '. Ensure it is added to abbreviations in pwd_utils.py'
                print(message)
                raise KeyError(message) from e
            required_params = mp.model_definitions_dict[best_model_name]  # Making a copy here because of potential for changing contents of original
            mp.model_definitions_dict[new_model_name] = required_params
            self.models[timescale_type][consider_thermohaline][new_model_name] = pm.PollutionModel(
                new_model_name,
                timescale_type,
                self.enhancement_model,
                prior,
                consider_thermohaline,
                self.n_live_points,
                self.seed,
                self.verbose,
                self.resume
            )
            toret.append(new_model_name)
        return toret

    def run_system(self, system_to_run):
        white_dwarf = self.white_dwarfs[system_to_run]
        for timescale_type in self.timescale_types_to_run:
            if self.use_hierarchy:
                self.models[timescale_type] = collections.OrderedDict()
            for consider_thermohaline in self.thermohaline_regimes_to_run:
                if self.use_hierarchy:
                    self.models[timescale_type][consider_thermohaline] = collections.OrderedDict()
                successfully_ran_timescale_type = False
                if self.use_hierarchy:
                    best_model = None
                    if len(list(self.parameter_hierarchy.keys())) > 4:  # Too many combinations to try every one
                        accepted_hierarchy_levels = [0]
                        max_hierarchy_level = 0
                        while self.parameter_hierarchy.get(max_hierarchy_level+1) is not None:
                            base_model_name = self.build_model(accepted_hierarchy_levels)
                            comparison_model_name = self.build_model(accepted_hierarchy_levels + [max_hierarchy_level+1])
                            print('Attempting to register ' + base_model_name + ' and ' + comparison_model_name)
                            self.models[timescale_type][consider_thermohaline][base_model_name] = pm.PollutionModel(
                                base_model_name,
                                timescale_type,
                                self.enhancement_model,
                                self.default_prior,
                                consider_thermohaline,
                                self.n_live_points,
                                self.seed,
                                self.verbose,
                                self.resume
                            )
                            self.models[timescale_type][consider_thermohaline][comparison_model_name] = pm.PollutionModel(
                                comparison_model_name,
                                timescale_type,
                                self.enhancement_model,
                                self.default_prior,
                                consider_thermohaline,
                                self.n_live_points,
                                self.seed,
                                self.verbose,
                                self.resume
                            )
                            successfully_ran_timescale_type = self.execute_models(system_to_run, timescale_type, consider_thermohaline) # We don't do anything with this output though?
                            best_model = self.compare(system_to_run, timescale_type, consider_thermohaline)
                            max_hierarchy_level += 1
                            if comparison_model_name == best_model:
                                accepted_hierarchy_levels += [max_hierarchy_level]
                    else:
                        non_zero_levels = list(self.parameter_hierarchy.keys())
                        non_zero_levels.remove(0)
                        all_combos = list()
                        for length in range(0, len(non_zero_levels)+1):
                            for subset in itertools.combinations(non_zero_levels, length):
                                all_combos.append([0] + list(subset))
                        for combo in all_combos:
                            model_name = self.build_model(combo)
                            print('Attempting to register ' + model_name)
                            self.models[timescale_type][consider_thermohaline][model_name] = pm.PollutionModel(
                                model_name,
                                timescale_type,
                                self.enhancement_model,
                                self.default_prior,
                                consider_thermohaline,
                                self.n_live_points,
                                self.seed,
                                self.verbose,
                                self.resume
                            )
                        self.execute_models(system_to_run, timescale_type, consider_thermohaline)
                        best_model = self.compare(system_to_run, timescale_type, consider_thermohaline)
                        #self.final_hierarchy = [int(c) for c in best_model if c.isdigit()]
                    ap_models = self.register_models_with_altered_priors(best_model, white_dwarf.full_name(), timescale_type, consider_thermohaline)
                    self.execute_models(system_to_run, timescale_type, consider_thermohaline)
                    for ap_model in ap_models:
                        self.analyser.compare_two_models(self.models[timescale_type][consider_thermohaline][ap_model], self.models[timescale_type][consider_thermohaline][best_model], output_dir, system_to_run, len(white_dwarf.get_elements_present()))
                        for other_ap_model in ap_models:
                            if ap_model is not other_ap_model:  # Using 'is not' rather than '!=' checks that they are (not) actually the same object instance, rather than (not) objects of equal value, which is want we want here
                                self.analyser.compare_two_models(self.models[timescale_type][consider_thermohaline][ap_model], self.models[timescale_type][consider_thermohaline][other_ap_model], output_dir, system_to_run, len(white_dwarf.get_elements_present()))
                    self.analyse_models(system_to_run, timescale_type, consider_thermohaline)
                else:
                    self.execute_models(system_to_run, timescale_type, consider_thermohaline)
                    self.compare(system_to_run, timescale_type, consider_thermohaline)
                    self.analyse_models(system_to_run, timescale_type, consider_thermohaline)

    def analyse_models(self, N_wd, timescale_type, consider_thermohaline):
        if self.suppress_graphical_output:
            print('Graphical output is being suppressed')
        white_dwarf = self.white_dwarfs[N_wd]
        #wd_name = white_dwarf.full_name()
        self.create_dir(self.get_timescale_output_dir(white_dwarf.abbreviated_name_as_used, timescale_type, consider_thermohaline))
        self.analyser.update_graph_dir(self.get_timescale_output_dir(white_dwarf.abbreviated_name_as_used, timescale_type, consider_thermohaline))
        self.publish_live_data(N_wd, timescale_type)
        hierarchy_name = pu.hierarchy_abbreviations[self.hierarchy_name] if self.model_names is None else '_'.join(self.model_names)

        prior_name = None

        for model_name, model in self.models[timescale_type][consider_thermohaline].items():
            if prior_name is not None and model.prior_name != prior_name:
                raise ValueError('Did not foresee a situation where models are being run with different priors - need to make this more flexible!')
            prior_name = model.prior_name

        # This needs to mimic the logic from pollution_model - should probably actually use that logic rather than copy it
        stats_file = self.get_timescale_output_dir(white_dwarf.abbreviated_name_as_used, timescale_type, consider_thermohaline) + white_dwarf.full_name() + '_' + timescale_type.short_str()
        if consider_thermohaline:
            stats_file += '_t'
        stats_file += '_p' + str(self.n_live_points) + '_' + hierarchy_name + '_' + pu.abbreviations[self.enhancement_model] + '_' + pu.abbreviations[prior_name] + '_stats.csv'

        print('Writing to file: ' + stats_file)
        with open(stats_file, 'w', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow(['System Name:', white_dwarf.system_name])
            to_write.writerow(['Variant Name:', white_dwarf.variant])
            to_write.writerow(['Full Name:', white_dwarf.full_name()])
            to_write.writerow(['T_eff:', str(white_dwarf.get_teff())])
            to_write.writerow(['log(g):', str(white_dwarf.get_logg())])
            to_write.writerow(['Timescale Type:', str(timescale_type)])
            to_write.writerow(['Thermohaline:', str(consider_thermohaline)])
            to_write.writerow(['Input file used:', self.wd_data_filename])
            to_write.writerow(['System ID in input file:', str(N_wd)])
            to_write.writerow([
                'Model',
                'Best model?',
                'Good fit?',
                'Base model',
                'ln_Z_model',
                'ln_Z_base',
                'Bayes_factor_model_base',
                'n_sigma_model_base',
                'chi_model',
                'chi_base',
                'chi_model_per_data_point',
                'chi_base_per_data_point',
                'Model params'
            ])

        for model_name, model in self.models[timescale_type][consider_thermohaline].items():
            self.publish_live_model(model_name, consider_thermohaline)
            self.analyser.dump_model_stats(
                model_name,
                model,
                stats_file
            )

        #self.analyser.find_diff_sigma(self.models, N_wd, len(ld._live_non_zero_wd_abundances), stats_file)


        # TODO: This stuff still gets printed in the csv file even when the model didn't run (i.e., the 3DOvershoot case for He dominated systems) with no indication that the output is not valid!
        self.analyser.find_parameter_sigma('Differentiation', mp.ModelParameter.fragment_core_frac, self.models[timescale_type][consider_thermohaline], N_wd, len(white_dwarf.get_elements_present()), stats_file)
        self.analyser.find_parameter_sigma('Heating', mp.ModelParameter.formation_distance, self.models[timescale_type][consider_thermohaline], N_wd, len(white_dwarf.get_elements_present()), stats_file)
        self.analyser.find_best_heated_model(self.models[timescale_type][consider_thermohaline], N_wd, len(white_dwarf.get_elements_present()), stats_file)
        plot_best_nondiff = False
        plot_best_heated = False
        bonus_plots = list()
        bonus_plot_names = dict()#{
        #    'HM02' : 'Best model without heating'
        #}
        bonus_fits = dict()
        bonus_error_lows = dict()
        bonus_error_highs = dict()

        with open(stats_file, 'a', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow([])
            to_write.writerow(['Interpolated Quantity', 'Timescale or Value'])
            try:
                to_write.writerow(['log(q)', str(white_dwarf.get_logq(timescale_type))])
                for element, value in white_dwarf.get_timescale_values_dict(timescale_type).items():
                    to_write.writerow([str(element), str(value)])
            except KeyError:
                # Occurs if there were no timescales for this timescale_type, so there is not going to be any other output from here on!
                to_write.writerow(['None available, None available'])
                return

        models_to_dump = list()
        for model_name, model in self.models[timescale_type][consider_thermohaline].items():
            if model.best_model:
                models_to_dump = [model_name] + models_to_dump
            elif (model.best_model_without_parameter.get(mp.ModelParameter.fragment_core_frac, False) and plot_best_nondiff) or (model.best_heated_model and plot_best_heated) or model_name in bonus_plots:
                models_to_dump.append(model_name)
            else:
                pass
        for model_name in models_to_dump:
            model = self.models[timescale_type][consider_thermohaline][model_name]
            self.publish_live_model(model_name, consider_thermohaline)
            self.analyser.make_plots_and_dump_fit(
                white_dwarf,
                N_wd,
                timescale_type,
                model.actual_output_dir,
                stats_file,
                model_name,
                model,
                self.enhancement_model,
                consider_thermohaline,
                bonus_fits,
                bonus_error_lows,
                bonus_error_highs,
                self.suppress_graphical_output
            )

