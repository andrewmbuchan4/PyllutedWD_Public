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

    def __init__(self, args):
        self.wd_data_filename = args.wd_data_filename
        self.stellar_compositions_filename = args.stellar_compositions_filename
        self.enhancement_model = args.enhancement_model
        try:
            self.seed = args.seed
        except AttributeError:
            # Then seed was not supplied: assume -1
            self.seed = -1
        if args.pollution_model_names is not None and len(args.pollution_model_names) == 1 and args.pollution_model_names[0].startswith('Hierarchy'):
            # Then we're running a hierarchy which will dynamically create models to run
            self.model_names = None
            self.use_hierarchy = True
            self.hierarchy_name = args.pollution_model_names[0]
            try:
                self.parameter_hierarchy = mp.hierarchy_definitions_dict[self.hierarchy_name]
            except KeyError as e:
                print('Error: Did not recognise hierarchy ' + args.pollution_model_names[0])
                print('Available hierarchies: ' + ', '.join(mp.hierarchy_definitions_dict.keys()))
                raise
        else:
            # Then we're running a specified batch of models
            self.model_names = args.pollution_model_names
            self.use_hierarchy = False
            self.parameter_hierarchy = None
            self.hierarchy_name = None
        self.timescale_interpolator = ti.TimescaleInterpolator()
        self.final_hierarchy = None
        self.white_dwarfs = list()
        self.default_logg = 8
        self.default_Ca_value = -15
        self.n_live_points = args.n_live_points
        self.stellar_compositions = None
        self.models = collections.OrderedDict()
        self.default_prior = 'Default'  # Use this prior except for running certain comparison priors at the end of a hierarchy execution. This placement isn't ideal since prior_functions.py also has this as its default prior
        self.comparison_priors = list()# ['HighPressure', 'LowPressure']
        self.load_global_data()
        if self.model_names is not None:
            self.load_models()
        # I assume there will only ever be one of these:
        self.analyser = ma.ModelAnalyser(None)
        self.executed_models = dict()

    def get_output_dir(self, wd_name):
        filename = pu.get_path_to_pylluted_dir() + wd_name + '/'
        self.create_dir(filename)
        return filename

    def get_chains_dir(self, wd_name):
        filename = self.get_output_dir(wd_name) + 'c/' # Choosing a short name because of the 100 character limit
        self.create_dir(filename)
        return filename

    def create_dir(self, dir_to_make):
        try:
            os.makedirs(dir_to_make)
        except FileExistsError:
            pass

    def load_global_data(self):
        self.load_wd_data(self.wd_data_filename)
        self.load_compositions(self.stellar_compositions_filename)

    #def load_wd_data(self, wd_data_filename):
    #    # Name	Type	Mass	Teff     log(g)       log(q)	log(Al/X)	error log(Al/X)	log(Ti/X)	  error log(Ti/X)	log(Ca/X)	error log(Ca/X)	log(Ni/X)   error log(Ni/X)	log(Fe/X)	error log(Fe/X)	log(Cr/X)	error log(Cr/X)	log(Mg/X)	error log(Mg/X)	log(Si/X)	error log(Si/X)	log(Na/X)	error log(Na/X)	log(O/X)	error log(O/X)	log(C/X)	error log(C/X)	log(N/X)	error log(N/X)	t_Al	t_Ti	t_Ca	t_Ni	t_Fe	t_Cr	t_Mg	t_Si	t_Na	t_O	t_C	t_N
    #    # 0      1        2       3         4           5               6            7            8              9            10             11            12              13         14              15           16             17          18               19        20               21          22              23        24                25        26             27           28            29             30      31      32      33      34      35      36    37       38      39   40  41
    #    if (self.wd_names is not None) and (wd_data_filename is not None):
    #        warnings.warn('Warning! Attempting to reload WD data')
    #    if wd_data_filename is None:
    #        return
    #    self.wd_names = list()
    #    self.wd_types = dict()
    #    self.wd_masses = dict()
    #    self.wd_qs = dict()
    #    self.wd_timescales = dict()
    #    self.wd_abundances = dict()
    #    self.wd_abundance_upper_bounds = dict()
    #    self.wd_abundance_lower_bounds = dict()
    #    self.wd_errors = dict()
    #    self.wd_excluded_abundances = dict()
    #    self.wd_excluded_abundance_upper_bounds = dict()
    #    self.wd_excluded_abundance_lower_bounds = dict()
    #    self.wd_excluded_errors = dict()
    #    input_abundance_indices = {
    #        ci.Element.Al: 6,
    #        ci.Element.Ti: 8,
    #        ci.Element.Ca: 10,
    #        ci.Element.Ni: 12,
    #        ci.Element.Fe: 14,
    #        ci.Element.Cr: 16,
    #        ci.Element.Mg: 18,
    #        ci.Element.Si: 20,
    #        ci.Element.Na: 22,
    #        ci.Element.O: 24,
    #        ci.Element.C: 26,
    #        ci.Element.N: 28
    #    }
    #    input_timescale_indices = {
    #        ci.Element.Al: 30,
    #        ci.Element.Ti: 31,
    #        ci.Element.Ca: 32,
    #        ci.Element.Ni: 33,
    #        ci.Element.Fe: 34,
    #        ci.Element.Cr: 35,
    #        ci.Element.Mg: 36,
    #        ci.Element.Si: 37,
    #        ci.Element.Na: 38,
    #        ci.Element.O: 39,
    #        ci.Element.C: 40,
    #        ci.Element.N: 41
    #    }
    #    #potential_elements_to_ignore = list()
    #    potential_elements_to_ignore = [ci.Element.C, ci.Element.N]
    #    subsolar_limits = {
    #        ci.Element.C: collections.OrderedDict({
    #            ci.Element.Mg: 0.95, # Mg and Si are very similar so just setting them the same
    #            ci.Element.Si: 0.95
    #        }),
    #        ci.Element.N: collections.OrderedDict({
    #            ci.Element.Mg: 0.35,
    #            ci.Element.Si: 0.35
    #        })
    #    }
    #    with open(pu.get_path_to_data() + wd_data_filename, encoding='utf-8') as wdcsv:
    #        row_count = 0
    #        for row in csv.reader(wdcsv):
    #            if row_count == 0: # Heading row
    #                pass  # TODO: Find column indices dynamically using the heading row
    #            else:
    #                wd_name = row[0]
    #                self.wd_names.append(wd_name)
    #                wd_type_raw = row[1]
    #                wd_type = ci.Element[wd_type_raw]  # TODO: Maybe put this in a try/except
    #                self.wd_types[wd_name] = wd_type
    #                self.wd_masses[wd_name] = float(row[2])
    #                self.wd_abundances[wd_name] = dict()
    #                self.wd_abundance_upper_bounds[wd_name] = dict()
    #                self.wd_abundance_lower_bounds[wd_name] = dict()
    #                self.wd_errors[wd_name] = dict()
    #                self.wd_excluded_abundances[wd_name] = dict()
    #                self.wd_excluded_abundance_upper_bounds[wd_name] = dict()
    #                self.wd_excluded_abundance_lower_bounds[wd_name] = dict()
    #                self.wd_excluded_errors[wd_name] = dict()
    #                self.wd_timescales[wd_name] = dict()
    #                for el in ci.usual_elements:
    #                    input_abundance_index = input_abundance_indices[el]
    #                    input_error_index = input_abundance_index + 1
    #                    try:
    #                        self.wd_abundances[wd_name][el] = float(row[input_abundance_index])
    #                        self.wd_abundance_upper_bounds[wd_name][el] = None
    #                        self.wd_abundance_lower_bounds[wd_name][el] = None
    #                    except ValueError:
    #                        # Occurs if, for example, this is an upper bound so the value starts with <
    #                        if row[input_abundance_index].startswith('<'):
    #                            # Then this is an upper bound
    #                            self.wd_abundance_upper_bounds[wd_name][el] = float(row[input_abundance_index][1:])
    #                            self.wd_abundances[wd_name][el] = 0.0
    #                            self.wd_abundance_lower_bounds[wd_name][el] = None
    #                        elif row[input_abundance_index].startswith('>'):
    #                            # Then this is an lower bound
    #                            self.wd_abundance_lower_bounds[wd_name][el] = float(row[input_abundance_index][1:])
    #                            self.wd_abundances[wd_name][el] = 0.0
    #                            self.wd_abundance_upper_bounds[wd_name][el] = None
    #                        else:
    #                            # Then we don't know what it is
    #                            raise ValueError('Unable to parse abundances for ' + wd_name)
    #                    self.wd_errors[wd_name][el] = float(row[input_error_index])
    #                    input_timescale_index = input_timescale_indices[el]
    #                    self.wd_timescales[wd_name][el] = float(row[input_timescale_index])
    #                # Remove C and/or N if they are significantly sub-solar
    #                # There is a bit of an issue here: in principle the 'sub-solar' description refers to the composition
    #                # at the point of condensation. The numbers here are at the point of observation. Could be different.
    #                # No obvious fix without assuming the answer, other than perhaps making the ref_limits more generous
    #                for peti in potential_elements_to_ignore:
    #                    peti_value = self.wd_abundances[wd_name][peti]
    #                    if peti_value == 0.0:
    #                        peti_value = self.wd_abundance_upper_bounds[wd_name][peti]
    #                    if peti_value not in [0.0, None]:
    #                        auto_ignore = True
    #                        allow = True
    #                        for ref_el, ref_limit in subsolar_limits[peti].items():
    #                            ref_abundance = self.wd_abundances[wd_name][ref_el]
    #                            if ref_abundance != 0.0:
    #                                auto_ignore = False
    #                                ratio = peti_value - ref_abundance
    #                                if ratio < ref_limit:
    #                                    allow = False
    #                                break
    #                        if auto_ignore or (not allow):
    #                            if self.wd_abundances[wd_name][peti] != 0.0:
    #                                self.wd_excluded_abundances[wd_name][peti] = self.wd_abundances[wd_name][peti]
    #                                self.wd_excluded_errors[wd_name][peti] = self.wd_errors[wd_name][peti]
    #                            if self.wd_abundance_upper_bounds[wd_name].get(peti) is not None:
    #                                self.wd_excluded_abundance_upper_bounds[wd_name][peti] = self.wd_abundance_upper_bounds[wd_name][peti]
    #                            if self.wd_abundance_lower_bounds[wd_name].get(peti) is not None:
    #                                self.wd_excluded_abundance_lower_bounds[wd_name][peti] = self.wd_abundance_lower_bounds[wd_name][peti]
    #                            self.wd_abundances[wd_name][peti] = 0.0
    #                            self.wd_errors[wd_name][peti] = 0.0
    #                            self.wd_abundance_upper_bounds[wd_name][peti] = None
    #                            self.wd_abundance_lower_bounds[wd_name][peti] = None
    #                # If any timescales (or logq) are 0, then we'll try to calculate them ourselves:
    #                try:
    #                    logg = float(row[4])
    #                except ValueError:
    #                    logg = self.default_logg
    #                Teff = int(row[3])
    #                self.wd_loggs[wd_name] = logg
    #                self.wd_teffs[wd_name] = Teff
    #                if wd_type == ci.Element.He and self.wd_abundances[wd_name][ci.Element.Ca] == 0:
    #                    # Then we need to intervene: technically the timescale_interpolator won't be able to calculate a timescale
    #                    # So we'll just use a very small amount of Ca
    #                    backup_timescales = self.timescale_interpolator.get_wd_timescales(str(wd_type), logg, Teff, self.default_Ca_value)
    #                else:
    #                    backup_timescales = self.timescale_interpolator.get_wd_timescales(str(wd_type), logg, Teff, self.wd_abundances[wd_name][ci.Element.Ca])
    #                for el in ci.usual_elements:
    #                    if (self.wd_timescales[wd_name][el] is None or self.wd_timescales[wd_name][el] == 0) and backup_timescales is not None:
    #                        self.wd_timescales[wd_name][el] = backup_timescales[el]
    #                logq = float(row[5])
    #                if (logq is None or logq == 0) and backup_timescales is not None:
    #                    logq = backup_timescales['logq']
    #                self.wd_qs[wd_name] = logq
    #            row_count += 1

    def load_wd_data(self, wd_data_filename):
        # Name	Type	Mass	Teff     log(g)       log(q)	log(Al/Hx)	error log(Al/Hx)	log(Ti/Hx)	  error log(Ti/Hx)	log(Ca/Hx)	error log(Ca/Hx)	log(Ni/Hx)   error log(Ni/Hx)	log(Fe/Hx)	error log(Fe/Hx)	log(Cr/Hx)	error log(Cr/Hx)	log(Mg/Hx)	error log(Mg/Hx)	log(Si/Hx)	error log(Si/Hx)	log(Na/Hx)	error log(Na/Hx)	log(O/Hx)	error log(O/Hx)	log(C/Hx)	error log(C/Hx)	log(N/Hx)	error log(N/Hx)	t_Al	t_Ti	t_Ca	t_Ni	t_Fe	t_Cr	t_Mg	t_Si	t_Na	t_O	t_C	t_N
        # 0      1        2       3         4           5               6            7            8              9            10             11            12              13         14              15           16             17          18               19        20               21          22              23        24                25        26             27           28            29             30      31      32      33      34      35      36    37       38      39   40  41
        if wd_data_filename is None:
            return
        input_abundance_indices = {
            ci.Element.Al: 9,
            ci.Element.Ti: 11,
            ci.Element.Ca: 13,
            ci.Element.Ni: 15,
            ci.Element.Fe: 17,
            ci.Element.Cr: 19,
            ci.Element.Mg: 21,
            ci.Element.Si: 23,
            ci.Element.Na: 25,
            ci.Element.O: 27,
            ci.Element.C: 29,
            ci.Element.N: 31
        }
        input_timescale_indices = {
            ci.Element.Al: 33,
            ci.Element.Ti: 34,
            ci.Element.Ca: 35,
            ci.Element.Ni: 36,
            ci.Element.Fe: 37,
            ci.Element.Cr: 38,
            ci.Element.Mg: 39,
            ci.Element.Si: 40,
            ci.Element.Na: 41,
            ci.Element.O: 42,
            ci.Element.C: 43,
            ci.Element.N: 44
        }
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
            row_count = 0
            for row in csv.reader(wdcsv):
                if row_count == 0: # Heading row
                    pass  # TODO: Find column indices dynamically using the heading row
                else:
                    wd_name = row[0]
                    wd_type_raw = row[1]
                    try:
                        wd_type = ci.Element[wd_type_raw]  # TODO: Maybe put this in a try/except
                    except KeyError:
                        raise KeyError('Invalid atmospheric type: must be H or He, received ' + str(wd_type_raw))

                    wd_abundance_data_raw = dict()
                    wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc')
                    timescale_dict = dict()

                    for el in ci.usual_elements:
                        input_abundance_index = input_abundance_indices[el]
                        input_error_index = input_abundance_index + 1
                        try:
                            data_point_type = wd.WhiteDwarfDataPointType.measurement
                            value = float(row[input_abundance_index])
                        except ValueError:
                            # Occurs if, for example, this is an upper bound so the value starts with <
                            if row[input_abundance_index].startswith('<'):
                                # Then this is an upper bound
                                data_point_type = wd.WhiteDwarfDataPointType.upper_bound
                                value = float(row[input_abundance_index][1:])
                            elif row[input_abundance_index].startswith('>'):
                                # Then this is an lower bound
                                data_point_type = wd.WhiteDwarfDataPointType.lower_bound
                                value = float(row[input_abundance_index][1:])
                            elif row[input_abundance_index] == '':
                                value = None
                            else:
                                # Then we don't know what it is
                                raise ValueError('Unable to parse abundances for ' + wd_name)
                        input_timescale_index = input_timescale_indices[el]
                        timescale_dict[el] = float(row[input_timescale_index])
                        included = True
                        if value is not None:
                            try:
                                error = float(row[input_error_index])
                            except ValueError:
                            # Occurs if error missing, e.g. if its an upper bound - so assume 0
                                error = 0
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

                    wd_mass = float(row[2])
                    try:
                        wd_mass_error = float(row[3])
                    except ValueError:
                        wd_mass_error = 0
                    # If any timescales (or logq) are 0, then we'll try to calculate them ourselves:
                    Teff = int(row[4])
                    try:
                        Teff_error = int(row[5])
                    except ValueError:
                        Teff_error = 0
                    try:
                        logg = float(row[6])
                    except ValueError:
                        logg = self.default_logg
                    try:
                        logg_error = float(row[7])
                    except ValueError:
                        logg_error = 0
                    wd_property_data_raw = {
                        mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, wd_type),
                        mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, Teff, Teff_error),
                        mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, logg, logg_error),
                        mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, wd_mass, wd_mass_error)
                    }
                    if wd_type == ci.Element.He and abundance_data.get_abundance(ci.Element.Ca) is None:
                        # Then we need to intervene: technically the timescale_interpolator won't be able to calculate a timescale
                        # So we'll just use a very small amount of Ca
                        backup_timescales = self.timescale_interpolator.get_wd_timescales(str(wd_type), logg, Teff, self.default_Ca_value)
                    else:
                        if wd_type == ci.Element.He:
                            backup_timescales = self.timescale_interpolator.get_wd_timescales(str(wd_type), logg, Teff, abundance_data.get_abundance(ci.Element.Ca).value)
                        else:
                            backup_timescales = self.timescale_interpolator.get_wd_timescales(str(wd_type), logg, Teff)
                    for el in ci.usual_elements:
                        if timescale_dict[el] in [None, 0, 0.0, np.nan] and backup_timescales is not None:
                            timescale_dict[el] = backup_timescales[el]
                    logq = float(row[8])
                    if (logq is None or logq == 0) and backup_timescales is not None:
                        logq = backup_timescales['logq']
                    wd_property_data_raw = {
                        mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, wd_type),
                        mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, Teff, Teff_error),
                        mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, logg, logg_error),
                        mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, wd_mass, wd_mass_error),
                        mp.WDParameter.logq: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, logq, 0)
                    }
                    property_data = wd.WhiteDwarfPropertyData(wd_property_data_raw)
                    white_dwarf = wd.WhiteDwarf(wd_name, property_data, abundance_data)
                    white_dwarf.timescale_dict = timescale_dict
                    self.white_dwarfs.append(white_dwarf)
                row_count += 1

    def load_compositions(self, input_filename=None):
        if (self.stellar_compositions is not None) and (input_filename is not None):
            warnings.warn('Warning! Attempting to reload compositions')
        if input_filename is not None:
            self.stellar_compositions_filename = input_filename
        if self.stellar_compositions_filename is None:
            return
        self.stellar_compositions = self.load_generic_float_data_csv(self.stellar_compositions_filename)

    def load_generic_float_data_csv(self, input_filename):
        # TODO: Change this to a with clause to prevent ResourceWarnings
        with open(pu.get_path_to_data() + input_filename, encoding='utf-8') as generic_csv:
            generic_list = [row for row in csv.reader(generic_csv)]
            generic_array = np.asarray(generic_list)
        return generic_array.astype(np.float)

    # TODO Maybe calculate live data for all runs and index in later?
    def publish_live_data(self, N_wd):
        # TODO: at some point I would like to simplify this to just ld._live_white_dwarf = self.white_dwarfs[N_wd]

        #non_zero_wd_abundances = list()
        #non_zero_wd_errors = list()
        #non_zero_wd_timescales = list()
        all_wd_abundances = list()
        all_wd_errors = list()
        all_wd_timescales = list()
        elements_present = list()

        white_dwarf = self.white_dwarfs[N_wd]
        ld._live_white_dwarf = white_dwarf
        wd_name = white_dwarf.name

        all_wd_timescales = white_dwarf.get_timescales_as_array(ci.usual_elements)
        abundances, errors, upper_bounds, lower_bounds = white_dwarf.get_abundance_arrays(ci.usual_elements)
        ld._live_all_wd_abundances = abundances
        ld._live_all_wd_errors = errors
        ld._live_all_wd_timescales = all_wd_timescales
        ld._live_upper_bounds = upper_bounds
        ld._live_lower_bounds = lower_bounds

        #for el in ci.usual_elements:
        #    wd_abundance = self.wd_abundances[wd_name][el]
        #    all_wd_abundances.append(wd_abundance)
        #    all_wd_errors.append(self.wd_errors[wd_name][el])
        #    all_wd_timescales.append(self.wd_timescales[wd_name][el])
        #    if wd_abundance != 0:
        #        non_zero_wd_abundances.append(wd_abundance)
        #        non_zero_wd_errors.append(self.wd_errors[wd_name][el])
        #        non_zero_wd_timescales.append(self.wd_timescales[wd_name][el])
        #        elements_present.append(el)

        #ld._live_non_zero_wd_abundances = np.transpose(non_zero_wd_abundances)
        #ld._live_non_zero_wd_errors = np.transpose(non_zero_wd_errors)
        #ld._live_non_zero_wd_timescales = np.transpose(non_zero_wd_timescales)

        #ld._live_all_wd_abundances = np.transpose(all_wd_abundances)
        #ld._live_all_wd_errors = np.transpose(all_wd_errors)
        #ld._live_all_wd_timescales = np.transpose(all_wd_timescales)
        #ld._live_upper_bounds = self.wd_abundance_upper_bounds[wd_name]
        #ld._live_lower_bounds = self.wd_abundance_lower_bounds[wd_name]

        ld._live_t_mg = white_dwarf.timescale_dict[ci.Element.Mg]
        ld._live_stellar_compositions = self.stellar_compositions  # TODO: This doesn't need to be updated every time
        ld._live_q = white_dwarf.get_logq().value
        ld._live_mass = white_dwarf.get_mass().value
        ld._live_type = white_dwarf.get_atmospheric_type().value
        ld._live_elements_present = white_dwarf.get_elements_present()

    def publish_live_model(self, model_name, prior_name='Default'):
        ld._live_model = model_name
        ld._enhancement_model = self.enhancement_model
        ld._live_prior = prior_name

    def load_models(self, model_names=None):
        models_already_loaded = self.models is not None and len(self.models.keys()) > 0
        if models_already_loaded and model_names is not None:
            warnings.warn('Warning! Attempting to reload models')
        if model_names is not None:
            self.model_names = model_names
        if self.model_names is None:
            return
        for model_name in self.model_names:
            self.models[model_name] = pm.PollutionModel(model_name, self.enhancement_model, self.default_prior, self.n_live_points, self.seed)

    def execute_models(self, N_wd):
        white_dwarf = self.white_dwarfs[N_wd]
        if self.models is not None:
            self.publish_live_data(N_wd)
            for model_name, model in self.models.items():
                if model_name not in self.executed_models.get(N_wd, list()):
                    # Then we haven't already run this one, so must do so. Otherwise, we can ignore.
                    self.publish_live_model(model_name, model.prior_name)
                    print('About to run model ' + model_name + ' with prior ' + model.prior_name)
                    model.execute(N_wd, self.get_chains_dir(white_dwarf.name))
                    if N_wd not in self.executed_models.keys():
                        self.executed_models[N_wd] = list()
                    self.executed_models[N_wd].append(model_name)

    def compare(self, N_wd):
        wd_name = self.white_dwarfs[N_wd].name
        print('Comparing models for ' + wd_name)
        self.publish_live_data(N_wd)
        number_of_data_points = len(self.white_dwarfs[N_wd].get_elements_present())
        self.analyser.compare(self.models, self.get_chains_dir(wd_name), N_wd, number_of_data_points)
        max_ln_Z = None
        max_ln_Z_name = None
        base_model_name = None

        for model_name, model in self.models.items():  # models should be an OrderedDict with the base model in first position
            if base_model_name is None:
                base_model_name = model_name
                max_ln_Z = model.comparison[base_model_name]['ln_Z_model']
                max_ln_Z_name = model_name
            else:
                if model.comparison[max_ln_Z_name]['ln_Z_model'] > max_ln_Z:
                    max_ln_Z = model.comparison[max_ln_Z_name]['ln_Z_model']
                    max_ln_Z_name = model_name
        print('Best model was ' + max_ln_Z_name)
        for model_name, model in self.models.items():
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

    def register_models_with_altered_priors(self, best_model_name, wd_name):
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
            self.models[new_model_name] = pm.PollutionModel(new_model_name, self.enhancement_model, prior, self.n_live_points, self.seed)
            toret.append(new_model_name)
        return toret

    def run_system(self, system_to_run):
        white_dwarf = self.white_dwarfs[system_to_run]
        if self.use_hierarchy:
            best_model = None
            if len(list(self.parameter_hierarchy.keys())) > 4:  # Too many combinations to try every one
                accepted_hierarchy_levels = [0]
                max_hierarchy_level = 0
                self.models = collections.OrderedDict()
                while self.parameter_hierarchy.get(max_hierarchy_level+1) is not None:
                    base_model_name = self.build_model(accepted_hierarchy_levels)
                    comparison_model_name = self.build_model(accepted_hierarchy_levels + [max_hierarchy_level+1])
                    print('Attempting to register ' + base_model_name + ' and ' + comparison_model_name)
                    self.models[base_model_name] = pm.PollutionModel(base_model_name, self.enhancement_model, self.default_prior, self.n_live_points, self.seed)
                    self.models[comparison_model_name] = pm.PollutionModel(comparison_model_name, self.enhancement_model, self.default_prior, self.n_live_points, self.seed)
                    self.execute_models(system_to_run)
                    best_model = self.compare(system_to_run)
                    max_hierarchy_level += 1
                    if comparison_model_name == best_model:
                        accepted_hierarchy_levels += [max_hierarchy_level]
                self.final_hierarchy = accepted_hierarchy_levels
            else:
                self.models = collections.OrderedDict()
                non_zero_levels = list(self.parameter_hierarchy.keys())
                non_zero_levels.remove(0)
                all_combos = list()
                for length in range(0, len(non_zero_levels)+1):
                    for subset in itertools.combinations(non_zero_levels, length):
                        all_combos.append([0] + list(subset))
                for combo in all_combos:
                    model_name = self.build_model(combo)
                    print('Attempting to register ' + model_name)
                    self.models[model_name] = pm.PollutionModel(model_name, self.enhancement_model, self.default_prior, self.n_live_points, self.seed)
                self.execute_models(system_to_run)
                best_model = self.compare(system_to_run)
                self.final_hierarchy = [int(c) for c in best_model if c.isdigit()]
            ap_models = self.register_models_with_altered_priors(best_model, white_dwarf.name)
            self.execute_models(system_to_run)
            for ap_model in ap_models:
                self.analyser.compare_two_models(self.models[ap_model], self.models[best_model], output_dir, system_to_run, len(white_dwarf.get_elements_present()))
                for other_ap_model in ap_models:
                    if ap_model is not other_ap_model:  # Using 'is not' rather than '!=' checks that they are (not) actually the same object instance, rather than (not) objects of equal value, which is want we want here
                        self.analyser.compare_two_models(self.models[ap_model], self.models[other_ap_model], output_dir, system_to_run, len(white_dwarf.get_elements_present()))
            self.analyse_models(system_to_run)

        else:
            self.execute_models(system_to_run)
            self.compare(system_to_run)
            self.analyse_models(system_to_run)

    def analyse_models(self, N_wd):
        suppress_graphical_output = False # This flag is a workaround (It exists because I found that if running a lot of systems in one go, the graphical output can cause weird crashes)
        if suppress_graphical_output:
            print('Warning! Graphical output is being suppressed. Switch it on in manager.py -> Manager class -> analyse_models')
        white_dwarf = self.white_dwarfs[N_wd]
        wd_name = white_dwarf.name
        self.create_dir(self.get_output_dir(wd_name))
        self.analyser.update_graph_dir(self.get_output_dir(wd_name))
        self.publish_live_data(N_wd)
        hierarchy_name = pu.hierarchy_abbreviations[self.hierarchy_name] if self.model_names is None else '_'.join(self.model_names)
        stats_file = self.get_output_dir(wd_name) + wd_name + '_p' + str(self.n_live_points) + '_' + hierarchy_name + '_' + pu.abbreviations[self.enhancement_model] + '_stats.csv'
        print('Writing to file: ' + stats_file)
        with open(stats_file, 'w', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow(['System Name:', str(wd_name)])
            to_write.writerow(['Input file used:', self.wd_data_filename])
            to_write.writerow(['System ID in input file:', str(N_wd)])
            to_write.writerow([
                'Model',
                'Base model',
                'ln_Z_model',
                'ln_Z_base',
                'Bayes_factor_model_base',
                'n_sigma_model_base',
                'chi_model',
                'chi_base',
                'chi_model_per_data_point',
                'chi_base_per_data_point',
                'Model params',
                'Best model?',
                'Good fit?'
            ])

        for model_name, model in self.models.items():
            self.publish_live_model(model_name)
            self.analyser.dump_model_stats(
                model_name,
                model,
                stats_file
            )

        #self.analyser.find_diff_sigma(self.models, N_wd, len(ld._live_non_zero_wd_abundances), stats_file)

        self.analyser.find_parameter_sigma('Differentiation', mp.ModelParameter.fragment_core_frac, self.models, N_wd, len(white_dwarf.get_elements_present()), stats_file, self.get_chains_dir(wd_name))
        self.analyser.find_parameter_sigma('Heating', mp.ModelParameter.formation_distance, self.models, N_wd, len(white_dwarf.get_elements_present()), stats_file, self.get_chains_dir(wd_name))
        self.analyser.find_best_heated_model(self.models, N_wd, len(white_dwarf.get_elements_present()), stats_file, self.get_chains_dir(wd_name))
        plot_best_nondiff = False
        plot_best_heated = False
        bonus_plots = list()#['Hierarchy_Default_levels_01']
        bonus_plot_names = {
        #    'Hierarchy_Default_levels_0' : 'dummy1',
        #    'Hierarchy_Default_levels_01' : 'dummy2',
        #    'Hierarchy_Default_levels_012' : 'Max Z',
        #    'Hierarchy_Default_levels_0123' : 'Best fit'
        }
        bonus_fits = dict()
        bonus_error_lows = dict()
        bonus_error_highs = dict()

        with open(stats_file, 'a', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow([])
            to_write.writerow(['Interpolated Quantity', 'Timescale or Value'])
            to_write.writerow(['log(q)', str(white_dwarf.get_logq())])
            for element, value in white_dwarf.timescale_dict.items():
                to_write.writerow([str(element), str(value)])

        models_to_dump = list()
        for model_name, model in self.models.items():
            if model.best_model:
                models_to_dump = [model_name] + models_to_dump
            elif (model.best_model_without_parameter.get(mp.ModelParameter.fragment_core_frac, False) and plot_best_nondiff) or (model.best_heated_model and plot_best_heated) or model_name in bonus_plots:
                models_to_dump.append(model_name)
            else:
                pass
        for model_name in models_to_dump:
            model = self.models[model_name]
            self.publish_live_model(model_name)
            self.analyser.make_plots_and_dump_fit(
                white_dwarf,
                N_wd,
                self.get_chains_dir(wd_name),
                stats_file,
                model_name,
                model,
                self.enhancement_model,
                bonus_fits,
                bonus_error_lows,
                bonus_error_highs,
                suppress_graphical_output
            )

            #temp_stats, mass_stats, comp_fits_and_errors, eo_samples, semisampled_eo_dict = self.analyser.make_all_plots(
            #    white_dwarf,
            #    self.get_chains_dir(wd_name),
            #    model_name,
            #    model,
            #    N_wd,
            #    #self.wd_timescales[wd_name],
            #    #wd_name,
            #    #self.wd_abundances[wd_name],
            #    #self.wd_abundance_upper_bounds[wd_name],
            #    #self.wd_abundance_lower_bounds[wd_name],
            #    #self.wd_errors[wd_name],
            #    #self.wd_excluded_abundances[wd_name],
            #    #self.wd_excluded_abundance_upper_bounds[wd_name],
            #    #self.wd_excluded_abundance_lower_bounds[wd_name],
            #    #self.wd_excluded_errors[wd_name],
            #    self.enhancement_model,
            #    #self.wd_types[wd_name],
            #    bonus_fits,
            #    bonus_error_lows,
            #    bonus_error_highs,
            #    suppress_graphical_output
            #)
            ##if bonus_plots != []:
            ##    for key, fit in mk2_fits_and_errors[0].items():
            ##        bonus_fits[bonus_plot_names[model_name] + ' ' + key] = fit
            ##    for key, error_lows in mk2_fits_and_errors[1].items():
            ##        bonus_error_lows[bonus_plot_names[model_name] + ' ' + key] = error_lows
            ##    for key, error_highs in mk2_fits_and_errors[2].items():
            ##        bonus_error_highs[bonus_plot_names[model_name] + ' ' + key] = error_highs
            #self.analyser.dump_model_fit(
            #    white_dwarf,
            #    self.get_chains_dir(wd_name),
            #    model,
            #    N_wd,
            #    stats_file,
            #    comp_fits_and_errors,
            #    temp_stats,
            #    mass_stats,
            #    #wd_name,
            #    model_name,
            #    #self.wd_abundances[wd_name],
            #    #self.wd_abundance_upper_bounds[wd_name],
            #    #self.wd_abundance_lower_bounds[wd_name],
            #    #self.wd_errors[wd_name],
            #    #self.wd_excluded_abundances[wd_name],
            #    #self.wd_excluded_abundance_upper_bounds[wd_name],
            #    #self.wd_excluded_abundance_lower_bounds[wd_name],
            #    #self.wd_excluded_errors[wd_name],
            #    eo_samples,
            #    #self.wd_timescales[wd_name],
            #    semisampled_eo_dict,
            #    suppress_graphical_output
            #)
#            if (model.best_model_without_parameter.get(mp.ModelParameter.fragment_core_frac, False) and plot_best_nondiff) or (model.best_heated_model and plot_best_heated) or model_name in bonus_plots:
#                temp_stats, mass_stats, mk2_fits_and_errors, semisampled_eo_dict = self.analyser.make_all_plots(
#                    self.get_chains_dir(wd_name),
#                    model_name,
#                    model,
#                    N_wd,
#                    self.wd_timescales[wd_name],
#                    wd_name,
#                    self.wd_abundances[wd_name],
#                    self.wd_abundance_upper_bounds[wd_name],
#                    self.wd_abundance_lower_bounds[wd_name],
#                    self.wd_errors[wd_name],
#                    self.enhancement_model,
#                    self.wd_types[wd_name],
#                    bonus_fits,
#                    bonus_error_lows,
#                    bonus_error_highs
#                )
#                if bonus_plots != []:
#                    for key, fit in mk2_fits_and_errors[0].items():
#                        bonus_fits[bonus_plot_names[model_name] + ' ' + key] = fit
#                    for key, error_lows in mk2_fits_and_errors[1].items():
#                        bonus_error_lows[bonus_plot_names[model_name] + ' ' + key] = error_lows
#                    for key, error_highs in mk2_fits_and_errors[2].items():
#                        bonus_error_highs[bonus_plot_names[model_name] + ' ' + key] = error_highs
#                self.analyser.dump_model_fit(
#                    self.get_chains_dir(wd_name),
#                    model,
#                    N_wd,
#                    stats_file,
#                    temp_stats,
#                    mass_stats,
#                    wd_name,
#                    model_name,
#                    self.wd_abundances[wd_name],
#                    self.wd_abundance_upper_bounds[wd_name],
#                    self.wd_abundance_lower_bounds[wd_name],
#                    self.wd_errors[wd_name],
#                    mk2_fits_and_errors[3]
#                )

