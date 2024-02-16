#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum
from scipy.stats import binom
import csv
import numpy as np
import os
import random
import scipy.integrate as scint
import scipy.interpolate as si


import chemistry_info as ci
import complete_model as cm
import hollands_abundances as ha
import live_data as ld
import model_parameters as mp
import pwd_utils as pu
import synthetic_bandpass as sb
import synthetic_configurations as sc
import timescale_interpolator as ti

class PopulationParameter(Enum):
    size = 0
    wd_config = 1
    pollution_config = 2

    def __str__(self):
        return self.name

class SyntheticSystem:

    def __init__(self, wd_properties, pollution_properties, pollution_abundances, observed=None, observed_abundances=None, bandpass_magnitudes=None, modelled_properties=None, id_no=None):
        self.wd_properties = wd_properties
        self.pollution_properties = pollution_properties
        self.pollution_abundances = pollution_abundances
        self.observed_abundances = observed_abundances
        self.bandpass_magnitudes = bandpass_magnitudes
        self.modelled_properties = modelled_properties
        self.id = id_no
        self.observed = observed

    def __repr__(self):
        if self.id is None:
            toret = '\nWD [No ID]:'
        else:
            toret = '\nWD ' + str(self.id) + ':'
        try:
            toret += '\n' + ', '.join([str(parameter) + ': ' + str(value) for parameter, value in self.wd_properties.items()])
        except AttributeError:  # If wd_properties is None
            toret += '\n[No WD properties]'
        try:
            toret += '\nPollution properties:\n' + ', '.join([str(parameter) + ': ' + str(value) for parameter, value in self.pollution_properties.items()])
        except AttributeError:
            toret += '\nPollution properties:\n[No Pollution properties]'
        try:
            toret += '\nPollution:\n' + ', '.join([str(el) + ': ' + str(abundance) for el, abundance in self.pollution_abundances.items()])
        except AttributeError:
            toret += '\nPollution:\n[No Pollution values]'
        toret += '\nObserved:\n' + str(self.observed)
        try:
            toret += '\nObserved Abundances:\n' + ', '.join([str(el) + ': ' + str(abundance) for el, abundance in self.observed_abundances.items()])
        except AttributeError:
            toret += '\nObserved Abundances:\n[No observed abundance values]'
        try:
            toret += '\nBandpass Magnitudes:\n' + ', '.join([str(bandpass) + ': ' + str(magnitude) for bandpass, magnitude in self.bandpass_magnitudes.items()])
        except AttributeError:
            toret += '\nBandpass Magnitudes:\n[No bandpass magnitude values]'
        try:
            toret += '\nModelled:\n' + ', '.join([str(parameter) + ': ' + str(value) for parameter, value in self.modelled_properties.items()])
        except AttributeError:
            toret += '\nModelled:\n[No Modelled properties]'
        toret += '\n'
        return toret

    def compare_dicts_ignoring_None(self, dict1, dict2):
        if dict1 is None:
            return dict2 is None
        if dict2 is None:
            return dict1 is None
        for key, val1 in dict1.items():
            val2 = dict2.get(key, None)
            if val1 is None and val2 is None:
                pass
            elif val1 == val2:
                pass
            else:
                return False
        for key, val2 in dict2.items():
            if val2 is not None and key not in dict1:
                return False
        return True


    def __eq__(self, other):
        if isinstance(other, self.__class__):
            wd_properties_match = self.compare_dicts_ignoring_None(self.wd_properties, other.wd_properties)
            pollution_properties_match = self.compare_dicts_ignoring_None(self.pollution_properties, other.pollution_properties)
            pollution_abundances_match = self.compare_dicts_ignoring_None(self.pollution_abundances, other.pollution_abundances)
            observed_abundances_match = self.compare_dicts_ignoring_None(self.observed_abundances, other.observed_abundances)
            bandpass_magnitudes_match = self.compare_dicts_ignoring_None(self.bandpass_magnitudes, other.bandpass_magnitudes)
            modelled_properties_match = self.compare_dicts_ignoring_None(self.modelled_properties, other.modelled_properties)
            observed_matches = self.observed == other.observed
            return wd_properties_match and pollution_properties_match and pollution_abundances_match and observed_abundances_match and bandpass_magnitudes_match and modelled_properties_match and observed_matches
        else:
            return False

    def to_csv_row(self, wd_property_keys, pollution_property_keys, pollution_abundance_keys, print_observed, observed_abundance_keys, bandpass_magnitude_keys, modelled_property_keys):
        toret = list()
        toret.append(str(self.id))
        for wpk in wd_property_keys:
            try:
                toret.append(str(self.wd_properties.get(wpk, 'None')))
            except AttributeError:
                toret.append('None')
        for ppk in pollution_property_keys:
            try:
                toret.append(str(self.pollution_properties.get(ppk, 'None')))
            except AttributeError:
                toret.append('None')
        for pak in pollution_abundance_keys:
            try:
                toret.append(str(self.pollution_abundances.get(pak, 'None')))
            except AttributeError:
                toret.append('None')
        if print_observed:
            toret.append(str(self.observed))
        for oak in observed_abundance_keys:
            try:
                toret.append(str(self.observed_abundances.get(oak, 'None')))
            except AttributeError:
                toret.append('None')
        for bmk in bandpass_magnitude_keys:
            try:
                toret.append(str(self.bandpass_magnitudes.get(bmk, 'None')))
            except AttributeError:
                toret.append('None')
        for mpk in modelled_property_keys:
            try:
                toret.append(str(self.modelled_properties.get(mpk, 'None')))
            except AttributeError:
                toret.append('None')
        return toret

    def set_bandpass_magnitudes(self, bandpass_magnitudes):
        self.bandpass_magnitudes = bandpass_magnitudes

    def set_modelled_properties(self, modelled_properties):
        self.modelled_properties = modelled_properties

    def set_observed_abundances(self, observed_abundances):
        self.observed_abundances = observed_abundances

    def pollution_abundance_log_ratio(self, element1, element2):
        el1 = self.pollution_abundances.get(element1)
        el2 = self.pollution_abundances.get(element2)
        try:
            ratio = el1 - el2
        except TypeError:
            ratio = None
        return ratio

    def observed_abundance_log_ratio(self, element1, element2):
        el1 = self.observed_abundances.get(element1)
        el2 = self.observed_abundances.get(element2)
        try:
            ratio = el1 - el2
        except TypeError:
            ratio = None
        return ratio

    def check_elements_detected(self, list_of_elements):
        # Returns True if all elements in list_of_elements were detected
        if not self.observed:
            return False
        for element in list_of_elements:
            if self.observed_abundances.get(element, None) is None:
                return False
        return True

    def check_properties_modelled(self, list_of_properties):
        # Returns True if all properties in list_of_properties were modelled (with no Nones or nans anywhere to be seen!)
        if self.modelled_properties is None:
            return False
        forbidden_values = [None, np.nan]
        for m_property in list_of_properties:
            if self.modelled_properties.get(m_property, None) in forbidden_values:
                return False
            if isinstance(self.modelled_properties.get(m_property, None), list):
                for m_value in self.modelled_properties.get(m_property, None):
                    if m_value in forbidden_values:
                        return False
        return True

    def distance_to(self, other, distance_elements=[ci.Element.Ca, ci.Element.Mg, ci.Element.Fe]):
        # Calculate "distance" to another SyntheticSystem in element space
        # For present purposes, only care about Ca, Mg and Fe. I'm using a Euclidian in log space - (linear space is bad because Ca is naturally underabundant relative to Mg and Fe so would be ignored)
        # But I normalise to the total quantities of these elements since I actually only care about relative abundances
        self_total_camgfe = 0
        other_total_camgfe = 0
        for element in distance_elements:
            self_total_camgfe += 10**self.observed_abundances[element]
            other_total_camgfe += 10**other.observed_abundances[element]
        self_total_camgfe = np.log10(self_total_camgfe)
        other_total_camgfe = np.log10(other_total_camgfe)
        normalisation_offset = self_total_camgfe - other_total_camgfe
        distance_list = list()
        for element in distance_elements:
            distance = self.observed_abundances[element] - other.observed_abundances[element]
            distance_mod = distance - normalisation_offset
            distance_list.append(distance_mod)
        total_distance = 0
        for distance in distance_list:
            total_distance += distance**2
        total_distance = np.sqrt(total_distance)
        return total_distance

class SyntheticPopulation:

    def __init__(self, population_size, wd_config_to_use, pollution_config_to_use, output_filename=None, pre_made_population=None, timescale_interpolator=None):
        self.wd_config_to_use = wd_config_to_use
        self.stellar_compositions = None
        self.stellar_compositions_filename = 'StellarCompositionsSortFE.csv' # This should really be an argument (similarly in the modeller)
        self.pollution_config_to_use = pollution_config_to_use
        if timescale_interpolator is None:
            self.timescale_interpolator = ti.TimescaleInterpolator()
        else:
            self.timescale_interpolator = timescale_interpolator
        self.inverse_cdf_dict = dict()
        self.io_ks_test_results = dict()
        self.pop_ks_test_results = dict()
        self.spectral_type_dict = {
            'DA': ci.Element.H,
            'DB': ci.Element.He
        }
        self.wd_configurations = sc.wd_configurations
        self.pollution_configurations = sc.pollution_configurations
        self.id_string = 'ID'
        self.output_filename = output_filename
        if self.wd_config_to_use is not None and self.wd_config_to_use not in self.wd_configurations:
            raise ValueError('Unrecognised WD configuration ' + str(wd_config_to_use))
        if self.pollution_config_to_use is not None and self.pollution_config_to_use not in self.pollution_configurations:
            raise ValueError('Unrecognised pollution configuration ' + str(pollution_config_to_use))
        if pre_made_population is None:
            if self.output_filename is not None and os.path.isfile(self.output_filename):
                print(self.output_filename + ' already exists, loading population from file')
                # Then assume we are making a population that has already been made (and outputted) - so just load it up
                self.load_from_csv(self.output_filename, population_size)
            else:
                if output_filename is None:
                    print('No output filename supplied. Creating new population')
                else:
                    print(self.output_filename + ' does not already exist, creating new population')
                self.create_population(population_size)
        else:
            if isinstance(pre_made_population, str):
                print('Pre-made population set to ' + pre_made_population + ', will attempt to load from that file')
                # Then try to load it from a file called pre_made_population
                self.load_from_csv(pre_made_population, population_size)
            else:
                print('Pre-made population of ' + str(len(pre_made_population)) + ' systems was supplied. Will re-use it.')
                # Then assume this IS the pre made population, and they already have ids
                self.population = pre_made_population
        self.unsampled_indices = list(range(len(self)))

    def __repr__(self):
        toret = 'Synthetic Population of ' + str(len(self.population)) + ' systems:'
        system_count = 0
        for system in self.population:
            toret += '\nSystem ' + str(system_count)
            toret += str(system)
            system_count += 1
        return toret

    def convert_property_to_readable_str(self, system_property):
        if isinstance(system_property, ci.Element):
            return 'log(' + str(system_property) + '/Hx)'
        elif isinstance(system_property, mp.WDParameter):
            if system_property == mp.WDParameter.logg:
                return 'WD Log(g)'
            else:
                return 'WD ' + str(system_property)
        else:
            return str(system_property)

    def convert_readable_str_to_property(self, header_string):
        if header_string.startswith('WD'):
            new_str = header_string[3:].replace(' ', '_').replace('(', '').replace(')', '').lower()
            return mp.WDParameter[new_str]
        elif header_string.startswith('Input'):
            new_str = header_string[6:].replace(' ', '_').lower()
            if new_str == 'time_since_accretion':
                return mp.ModelParameter.t_sinceaccretion
            if new_str.endswith('fraction'):
                new_str = new_str[:-4]
                new_str = new_str.replace('_number_', '_')
            return mp.ModelParameter[new_str]
        elif header_string.startswith('Output'):
            new_str = header_string[7:].replace(' ', '_').lower()
            if new_str == 'time_since_accretion':
                return mp.ModelParameter.t_sinceaccretion
            if new_str.endswith('fraction'):
                new_str = new_str[:-4]
                new_str = new_str.replace('_number_', '_')
            return mp.ModelParameter[new_str]
        elif header_string.startswith('True'):
            new_str = header_string.split('(')[1].split('/')[0]
            return ci.Element[new_str]
        elif header_string == 'Observed?':
            return 'Observed'
        elif header_string.startswith('Observed'):
            new_str = header_string.split('(')[1].split('/')[0]
            return ci.Element[new_str]
        elif header_string == self.id_string:
            return self.id_string
        elif len(header_string) == 1:
            return sb.Bandpass[header_string]
        else:
            raise ValueError('Unrecognised header string: ' + header_string)

    def cast_string_to_float_or_int(self, input_string):
        if input_string == '' or input_string.isspace():
            return None
        if input_string == '-inf':
            return float('-inf')
        if input_string == 'None':
            return None
        if input_string in ['nan', 'np.nan']:
            return np.nan
        if input_string.startswith('[') and input_string.endswith(']'):
            # Then it's a list
            list_of_str = input_string[1:-1].split(',')
            toret = [self.cast_string_to_float_or_int(foi.strip()) for foi in list_of_str]
            return toret
        if '.' in input_string:
            return float(input_string)
        else:
            try:
                return int(input_string)
            except ValueError:
                return input_string

    def dump_to_csv(self, output_file=None, only_raw=False):
        pollution_properties = list() if only_raw else [pp for pp in mp.ModelParameter]
        pollution_abundance_keys = ci.usual_elements
        observed_abundance_keys = list() if only_raw else ci.usual_elements
        bandpasses = list() if only_raw else [b for b in sb.Bandpass]
        modelled_properties = list() if only_raw else [pp for pp in mp.ModelParameter]
        header = [self.id_string]
        header += [self.convert_property_to_readable_str(item) for item in mp.wd_parameters_for_synthesis]
        header += ['Input ' + self.convert_property_to_readable_str(item) for item in pollution_properties]
        header += ['True ' + self.convert_property_to_readable_str(item) for item in pollution_abundance_keys]
        if not only_raw:
            header += ['Observed?']
        header += ['Observed ' + self.convert_property_to_readable_str(item) for item in observed_abundance_keys]
        header += [self.convert_property_to_readable_str(item) for item in bandpasses]
        header += ['Output ' + self.convert_property_to_readable_str(item) for item in modelled_properties]
        if output_file is None:
            output_file = self.output_filename
        if output_file is None:
            print('Warning! No output file supplied, could not dump')
            return
        if os.path.isfile(output_file):
            preexisting_pop_size = sum(1 for line in open(output_file)) - 1  # Remembering to subtract header line
            if preexisting_pop_size >= len(self.population):
                # Then don't dump: this is a preexisiting population, which we presumably loaded up a subset from, and we're about to delete a bunch of data!
                print('File ' + output_file + ' already exists, aborting dump')
                return
            else:
                print('Warning: ' + output_file + ' already exists, but is smaller than generated population, proceeding with dump')
        else:
            print('Dumping to file ' + output_file)
        if '/' in output_file:
            storage_dir = output_file[:output_file.rfind('/')]
            os.makedirs(storage_dir, exist_ok=True)
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow(header)
            for system in self.population:
                to_write.writerow(system.to_csv_row(mp.wd_parameters_for_synthesis, pollution_properties, pollution_abundance_keys, not only_raw, observed_abundance_keys, bandpasses, modelled_properties))

    def load_from_csv(self, input_file, original_requested_pop_size=np.inf):
        if original_requested_pop_size is None:
            requested_pop_size = np.inf
        else:
            requested_pop_size = original_requested_pop_size
        header_dict = dict()
        header_row = None
        self.population = list()
        with open(input_file, encoding='utf-8') as input_csv:
            row_number = 0
            for row in csv.reader(input_csv):
                if row_number == 0:
                    # We need to sort out what the headings mean
                    header_row = row
                    for i, header in enumerate(row):
                        header_dict[i] = self.convert_readable_str_to_property(header)
                else:
                    wd_properties = dict()
                    pollution_properties = dict()
                    pollution_abundances = dict()
                    bandpass_magnitudes = dict()
                    modelled_properties = dict()
                    id_no = None
                    observed = None
                    observed_abundances = dict()
                    for i, value in enumerate(row):
                        required_property = header_dict[i]
                        if isinstance(required_property, mp.WDParameter):
                            wd_properties[required_property] = self.cast_string_to_float_or_int(value)
                        elif isinstance(required_property, mp.ModelParameter):
                            if header_row[i].startswith('Input'):
                                pollution_properties[required_property] = self.cast_string_to_float_or_int(value)
                            elif header_row[i].startswith('Output'):
                                if value != 'None':
                                    modelled_properties[required_property] = self.cast_string_to_float_or_int(value)
                            else:
                                raise ValueError('Unrecognised header:' + header_row[i])
                        elif isinstance(required_property, ci.Element):
                            if header_row[i].startswith('True'):
                                pollution_abundances[required_property] = self.cast_string_to_float_or_int(value)
                            elif header_row[i].startswith('Observed'):
                                if value != 'None':
                                    observed_abundances[required_property] = self.cast_string_to_float_or_int(value)
                            else:
                                raise ValueError('Unrecognised header:' + header_row[i])
                        elif isinstance(required_property, sb.Bandpass):
                            if value != 'None':
                                bandpass_magnitudes[required_property] = self.cast_string_to_float_or_int(value)
                        elif required_property == self.id_string:
                            id_no = int(value)
                        elif required_property == 'Observed':
                            if value == 'True':
                                observed = True
                            elif value == 'False':
                                observed = False
                            else:
                                pass
                        else:
                            raise ValueError('Unrecognised header:' + header_row[i])
                    if bandpass_magnitudes == dict():
                        bandpass_magnitudes = None
                    if observed_abundances == dict():
                        observed_abundances = None
                    if modelled_properties == dict():
                        modelled_properties = None
                    new_system = SyntheticSystem(wd_properties, pollution_properties, pollution_abundances, observed, observed_abundances, bandpass_magnitudes, modelled_properties, id_no)
                    self.population.append(new_system)
                    if len(self.population) >= requested_pop_size:
                        # If we only want to load up 5 systems, but the file contains 50, we should return early
                        return
                row_number += 1

        if requested_pop_size > len(self.population):
            print('Warning! Requested population size ' + str(requested_pop_size) + ' was greater than systems available to load (' + str(len(self.population)) + ')')
            if np.isfinite(requested_pop_size):
                print('Will create the remaining systems')
                self.create_population(requested_pop_size, False)

    def __len__(self):
        return len(self.population)

    def __iter__(self):
        return iter(self.population)

    def __next__(self):
        return next(self.population)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.population == other.population  # This could potentially be edited to ignore order?
        else:
            return False

    def __getitem__(self, index):
        return self.population[index]

    def load_compositions(self):
        if self.stellar_compositions is None: # So we only do this once
            stellar_compositions = self.load_generic_float_data_csv(self.stellar_compositions_filename)
            ld._live_stellar_compositions = stellar_compositions
            self.stellar_compositions = stellar_compositions

    def load_generic_float_data_csv(self, input_filename):
        with open(pu.get_path_to_data() + input_filename, encoding='utf-8') as generic_csv:
            generic_list = [row for row in csv.reader(generic_csv)]
            generic_array = np.asarray(generic_list)
        return generic_array.astype(np.float)

    def get_next_id(self):
        id_list = list()
        for system in self.population:
            id_list.append(system.id)
        if len(id_list) < 1:
            return 0
        highest_id = max(id_list)
        return highest_id + 1

    def extract_raw_population(self):
        toret = list()
        for system in self.population:
            toret.append(SyntheticSystem(system.wd_properties, system.pollution_properties, system.pollution_abundances, None, None, None, None, system.id))
        return toret

    def check_sinking_timescales(self, reference_timescale=1E5):
        total_count = 0
        long_timescale_count = 0
        for system in self.population:
            wd_timescales = self.timescale_interpolator.get_wd_timescales(
                str(self.spectral_type_dict[system.wd_properties[mp.WDParameter.spectral_type]]),
                system.wd_properties[mp.WDParameter.logg],
                system.wd_properties[mp.WDParameter.temperature],
                system.pollution_abundances[ci.Element.Ca]
            )
            if wd_timescales[ci.Element.Mg] > reference_timescale:
                long_timescale_count += 1
            total_count += 1
        print(str(long_timescale_count) + '/' + str(total_count) + ' systems exceeded t_Mg = ' + str(reference_timescale))

    def create_population(self, population_size, reset_population=True):
        print('Generating synthetic population of ' + str(population_size) + ' systems')
        if reset_population:
            self.population = list()
        pop_count = len(self.population)
        excluded_count = 0
        while pop_count < population_size:
            new_system = self.create_system()
            print(new_system)
            if new_system is None:
                excluded_count += 1
            else:
                all_nans = True
                for element, abundance in new_system.pollution_abundances.items():
                    if not np.isnan(abundance):
                        all_nans = False
                if not all_nans:
                    self.population.append(new_system)
                    pop_count += 1
                else:
                    raise
                    excluded_count += 1
            if pop_count % 100 == 0:
                print('Made ' + str(pop_count) + ' systems so far, excluded ' + str(excluded_count))

    def create_system(self, use_iterative_approach=False):
        # use_iterative_approach: True is now deprecated but keeping it as an option just in case I want to reproduce older results
        self.load_compositions()
        enhancement_model = 'Earthlike'
        t_formation = 1.5
        normalise_abundances = True
        snapshot_wd_atm = True
        wd_config_dict = self.wd_configurations[self.wd_config_to_use]
        wd_dict = dict()
        for parameter in mp.wd_parameters_for_synthesis:
            wd_dict[parameter] = self.draw_variable_from_distribution(wd_config_dict[parameter])
        pollution_config_dict = self.pollution_configurations[self.pollution_config_to_use]
        input_dict = dict()
        for parameter in mp.ModelParameter:
            input_dict[parameter] = self.draw_variable_from_distribution(pollution_config_dict[parameter])
        if isinstance(input_dict[mp.ModelParameter.formation_distance], (np.integer, int)):
            input_dict[mp.ModelParameter.formation_distance] = float(input_dict[mp.ModelParameter.formation_distance])
        CaHe = -10 # initial guess, only relevant for DBs
        #modified_pollution_frac = None
        modified_pollution_frac = input_dict[mp.ModelParameter.pollution_frac]
        HorHe = self.spectral_type_dict[wd_dict[mp.WDParameter.spectral_type]]
        HorHe_str = str(HorHe)
        if HorHe == ci.Element.H:
            wd_timescales = self.timescale_interpolator.get_wd_timescales(HorHe_str, wd_dict[mp.WDParameter.logg], wd_dict[mp.WDParameter.temperature], CaHe)
            modified_pollution_frac = self.modify_pollution_frac(
                input_dict[mp.ModelParameter.pollution_frac],
                input_dict[mp.ModelParameter.t_sinceaccretion],
                input_dict[mp.ModelParameter.accretion_timescale],
                wd_timescales
            )
            ld._live_all_wd_timescales = [wd_timescales[el] for el in ci.usual_elements] # TODO: We also shouldn't need to do this twice
            ld._live_non_zero_wd_timescales = [wd_timescales[el] for el in ci.usual_elements]
            output_dict, ignore = cm.complete_model_calculation(
                input_dict[mp.ModelParameter.metallicity],
                input_dict[mp.ModelParameter.t_sinceaccretion],
                input_dict[mp.ModelParameter.formation_distance],
                input_dict[mp.ModelParameter.feeding_zone_size],
                input_dict[mp.ModelParameter.parent_core_frac],
                input_dict[mp.ModelParameter.parent_crust_frac],
                input_dict[mp.ModelParameter.fragment_core_frac],
                input_dict[mp.ModelParameter.fragment_crust_frac],
                modified_pollution_frac,
                input_dict[mp.ModelParameter.accretion_timescale],
                input_dict[mp.ModelParameter.pressure],
                input_dict[mp.ModelParameter.oxygen_fugacity],
                enhancement_model,
                t_formation,
                normalise_abundances,
                snapshot_wd_atm
            )
        elif HorHe == ci.Element.He:
            if use_iterative_approach:
                scaling_element = ci.Element.Mg
                tolerance = 0.00001
                rel_diff = tolerance + 1
                iterations = 0
                max_iterations = 20
                # max_iterations = -1 # For testing purposes (skipping the self-consistency check)
                initial_cahe = None
                initial_mghe = None
                initial_fehe = None
                while rel_diff > tolerance:
                    wd_timescales = self.timescale_interpolator.get_wd_timescales(HorHe_str, wd_dict[mp.WDParameter.logg], wd_dict[mp.WDParameter.temperature], CaHe)
                    modified_pollution_frac = self.modify_pollution_frac(
                        input_dict[mp.ModelParameter.pollution_frac],
                        input_dict[mp.ModelParameter.t_sinceaccretion],
                        input_dict[mp.ModelParameter.accretion_timescale],
                        wd_timescales,
                        scaling_element
                    )
                    ld._live_all_wd_timescales = [wd_timescales[el] for el in ci.usual_elements] # TODO: We also shouldn't need to do this twice
                    ld._live_non_zero_wd_timescales = [wd_timescales[el] for el in ci.usual_elements]
                    output_dict, ignore = cm.complete_model_calculation(
                        input_dict[mp.ModelParameter.metallicity],
                        input_dict[mp.ModelParameter.t_sinceaccretion],
                        input_dict[mp.ModelParameter.formation_distance],
                        input_dict[mp.ModelParameter.feeding_zone_size],
                        input_dict[mp.ModelParameter.parent_core_frac],
                        input_dict[mp.ModelParameter.parent_crust_frac],
                        input_dict[mp.ModelParameter.fragment_core_frac],
                        input_dict[mp.ModelParameter.fragment_crust_frac],
                        modified_pollution_frac,
                        input_dict[mp.ModelParameter.accretion_timescale],
                        input_dict[mp.ModelParameter.pressure],
                        input_dict[mp.ModelParameter.oxygen_fugacity],
                        enhancement_model,
                        t_formation,
                        normalise_abundances,
                        snapshot_wd_atm
                    )
                    if output_dict is None:
                        print('Warning! WD creation had None result, returning None')
                        return None
                    if iterations == 0:
                        initial_cahe = output_dict[ci.Element.Ca]
                        initial_mghe = output_dict[ci.Element.Mg]
                        initial_fehe = output_dict[ci.Element.Fe]
                    diff = output_dict[ci.Element.Ca] - CaHe
                    rel_diff = abs(diff/CaHe)
                    if iterations % 20 == 0 and iterations > 0:
                        print()
                        print('Making WD (iteration ' + str(iterations) + ') ')
                        print('Request Pfrac: ' + str(input_dict[mp.ModelParameter.pollution_frac]))
                        print('Modified Pfrac: ' + str(modified_pollution_frac))
                        print('Input CaHe: ' + str(CaHe))
                        print('Initial CaHe: ' + str(initial_cahe))
                        print('Initial MgHe: ' + str(initial_mghe))
                        print('Initial FeHe: ' + str(initial_fehe))
                        print('Initial CaFe: ' + str(initial_cahe - initial_fehe))
                        print('Initial MgFe: ' + str(initial_mghe - initial_fehe))
                        print('Output CaHe: ' + str(output_dict[ci.Element.Ca]))
                        print('Output MgHe: ' + str(output_dict[ci.Element.Mg]))
                        print('Output FeHe: ' + str(output_dict[ci.Element.Fe]))
                        print('Output CaFe: ' + str(output_dict[ci.Element.Ca] - output_dict[ci.Element.Fe]))
                        print('Output MgFe: ' + str(output_dict[ci.Element.Mg] - output_dict[ci.Element.Fe]))
                        print('rel_diff: ' + str(rel_diff))
                    if iterations > max_iterations:
                        print('Warning! WD creation exceeded max iterations (' + str(max_iterations) + '), using most recent attempt')
                        rel_diff = 0
                    CaHe = output_dict[ci.Element.Ca]
                    iterations += 1
            else:
                CaHe = -9 # Going slightly higher than -10. Based on this being kind of in the middle of the range of timescale variation, and a good compromise between high = less common but more observable
                scaling_element = ci.Element.Ca # Switched this - Ca is more representative than Mg
                wd_timescales = self.timescale_interpolator.get_wd_timescales(HorHe_str, wd_dict[mp.WDParameter.logg], wd_dict[mp.WDParameter.temperature], CaHe)
                modified_pollution_frac = self.modify_pollution_frac(
                    input_dict[mp.ModelParameter.pollution_frac],
                    input_dict[mp.ModelParameter.t_sinceaccretion],
                    input_dict[mp.ModelParameter.accretion_timescale],
                    wd_timescales,
                    scaling_element
                )
                ld._live_all_wd_timescales = np.array([wd_timescales[el] for el in ci.usual_elements])
                output_dict, ignore = cm.complete_model_calculation(
                    input_dict[mp.ModelParameter.metallicity],
                    input_dict[mp.ModelParameter.t_sinceaccretion],
                    input_dict[mp.ModelParameter.formation_distance],
                    input_dict[mp.ModelParameter.feeding_zone_size],
                    input_dict[mp.ModelParameter.parent_core_frac],
                    input_dict[mp.ModelParameter.parent_crust_frac],
                    input_dict[mp.ModelParameter.fragment_core_frac],
                    input_dict[mp.ModelParameter.fragment_crust_frac],
                    modified_pollution_frac,
                    input_dict[mp.ModelParameter.accretion_timescale],
                    input_dict[mp.ModelParameter.pressure],
                    input_dict[mp.ModelParameter.oxygen_fugacity],
                    enhancement_model,
                    t_formation,
                    normalise_abundances,
                    snapshot_wd_atm
                )
                if output_dict is None:
                    print('Warning! WD creation had None result, returning None')
                    return None
        else:
            raise ValueError('Unrecognised WD atmosphere type: ' + HorHe_str + ', must be H or He')
        input_dict[mp.ModelParameter.pollution_frac] = modified_pollution_frac
        new_system = SyntheticSystem(wd_dict, input_dict, output_dict, None, None, None, None, self.get_next_id())
        return new_system

    def modify_pollution_frac(self, pollution_frac, t_sinceaccretion, accretion_timescale, wd_timescales, scaling_element=ci.Element.Ca):
        #t_sinceaccretion is in Myr
        #accretion_timescale is in yr
        t_sinceaccretion_yrs = 1000000*t_sinceaccretion
        t_scaled_declining = (t_sinceaccretion_yrs - accretion_timescale)/wd_timescales[scaling_element]  # Mg is a somewhat arbitrary choice - Ca is actually more representative, no?
        if t_scaled_declining > 0:
            # Then we are in declining phase, and should reduce pollution fraction to prevent certain silly situations where pollution fraction is like -6 and it's all oxygen
            correction = t_scaled_declining/np.log(10)
            # This is essentially saying that the pollution fraction we inputted was actually the steady state pollution fraction
            return pollution_frac - correction
        else:
            return pollution_frac

    def draw_variable_from_distribution(self, distribution_description):
        distribution_type = distribution_description[0]
        distribution_args = distribution_description[1]
        if distribution_type == sc.Distribution.Uniform:
            # distribution_args = [min, max]
            assert len(distribution_args) == 2
            return np.random.uniform(distribution_args[0], distribution_args[1])
        elif distribution_type == sc.Distribution.Normal:
            # distribution_args = [mu, sigma]
            assert len(distribution_args) == 2
            return np.random.normal(distribution_args[0], distribution_args[1])
        elif distribution_type == sc.Distribution.Delta:
            # distribution_args = [list of possible values to choose from]
            return np.random.choice(distribution_args)
        elif distribution_type == sc.Distribution.Triangle:
            # distribution_args = [left, mode, right]
            assert len(distribution_args) == 3
            return np.random.triangular(distribution_args[0], distribution_args[1], distribution_args[2])
        elif distribution_type == sc.Distribution.Slope:
            # distribution_args = [left_edge_x, left_edge_y relative to right_edge_y, right_edge_x]
            assert len(distribution_args) == 3
            left_edge_x = distribution_args[0]
            left_edge_y = distribution_args[1]
            right_edge_x = distribution_args[2]
            right_edge_y = 1
            inv_gradient = (right_edge_x - left_edge_x)/(right_edge_y - left_edge_y)
            x_intercept = left_edge_x - (inv_gradient*left_edge_y)
            if left_edge_y < right_edge_y:
                sample = left_edge_x - 1
                while sample < left_edge_x:
                    sample = np.random.triangular(x_intercept, right_edge_x, right_edge_x)
                return sample
            elif left_edge_y > right_edge_y:
                sample = right_edge_x + 1
                while sample > right_edge_x:
                    sample = np.random.triangular(left_edge_x, left_edge_x, x_intercept)
                return sample
            else:
                return np.random.uniform(left_edge_x, right_edge_x)
        elif distribution_type == sc.Distribution.CustomFunction:
            # distribution_args = [function name, function, min, max]
            # First time we use this, we should cache some information about the function, i.e. the cdf
            if distribution_args[0] not in self.inverse_cdf_dict:
                self.cache_inverse_cdf(distribution_args[0], distribution_args[1], distribution_args[2], distribution_args[3])
            return float(self.inverse_cdf_dict[distribution_args[0]](np.random.uniform(0, 1)))
        elif distribution_type == sc.Distribution.CustomDistribution:
            # distribution_args = [function name, bins, counts]
            # First time we use this, we should cache some information about the function, i.e. the cdf
            if distribution_args[0] not in self.inverse_cdf_dict:
                self.cache_inverse_cdf_from_table(distribution_args[0], distribution_args[1], distribution_args[2])
            return float(self.inverse_cdf_dict[distribution_args[0]](np.random.uniform(0, 1)))
        else:
            raise ValueError('Unrecognised distribution type: ' + repr(distribution_type))

    def cache_inverse_cdf(self, function_name, function, min_val, max_val):
        sample_count = 100000
        areas = np.zeros(sample_count)
        points_to_sample = np.linspace(min_val, max_val, sample_count)
        for i, pts in enumerate(points_to_sample):
            area, error = scint.quad(function, min_val, pts)
            if area > 0 and error/area > 0.1:
                print('Warning! Integration error was: ' + str(error))
            areas[i] = area
        total_area = areas[-1]
        normalised_areas = np.zeros(sample_count)
        for i, area in enumerate(areas):
            normalised_areas[i] = area/total_area
        self.inverse_cdf_dict[function_name] = si.interp1d(normalised_areas, points_to_sample)

    def cache_inverse_cdf_from_table(self, function_name, bins, counts):
        #Assume bins is a list of (left edge, right edge) tuples
        bin_boundaries = np.zeros(len(counts) + 1)
        for i, bin_tuple in enumerate(bins):
            bin_boundaries[i] = bin_tuple[0]
        bin_boundaries[-1] = bins[-1][1]
        total_counts = 0
        for c in counts:
            total_counts += c
        normalised_counts = np.zeros(len(counts) + 1)
        normalised_counts[0] = 0
        count_so_far = 0
        for i, c in enumerate(counts):
            count_so_far += c
            normalised_counts[i + 1] = count_so_far/total_counts
        # WARNING: This caching seems to have a strange bug where it's possible to sample above the maximum input value
        # if you supply it with bins that go above that value (even if there's nothing in those bins). Workaround:
        # Make sure that your bins extend as far as the values and no further!
        self.inverse_cdf_dict[function_name] = si.interp1d(normalised_counts, bin_boundaries)

    def input_values(self, parameter):
        if isinstance(parameter, mp.ModelParameter):
            return self.pollution_input_values(parameter)
        else: # Then it must be a WDParameter
            return self.wd_values(parameter)

    def wd_values(self, wd_parameter):
        toret = list()
        for system in self.population:
            if wd_parameter == mp.WDParameter.spectral_type:
                toret.append(0 if system.wd_properties[wd_parameter] == 'DA' else 1)
            else:
                toret.append(system.wd_properties.get(wd_parameter))
        return toret

    def pollution_input_values(self, input_parameter):
        toret = list()
        for system in self.population:
            if system.pollution_properties is None:
                toret.append(None)
            else:
                toret.append(system.pollution_properties.get(input_parameter))
        return toret

    def pollution_abundance_values(self, element):
        toret = list()
        for system in self.population:
            toret.append(system.pollution_abundances.get(element))
        return toret

    def get_common_el_values(self, list_of_elements):
        toret = list()
        for element in list_of_elements:
            toret.append(list())
        for system in self.population:
            if system.check_elements_detected(list_of_elements):
                for j, element in enumerate(list_of_elements):
                    toret[j].append(system.observed_abundances[element])
        return toret

    def observed_abundances(self, element):
        toret = list()
        for system in self.population:
            if (not system.observed) or (system.observed_abundances is None):
                toret.append(None)
            else:
                toret.append(system.observed_abundances.get(element))
        return toret

    def pollution_abundance_log_ratios(self, element1, element2):
        toret = list()
        for system in self.population:
            ratio = system.pollution_abundance_log_ratio(element1, element2)
            el1 = system.pollution_abundances.get(element1)
            el2 = system.pollution_abundances.get(element2)
            try:
                ratio = el1 - el2
            except TypeError:
                ratio = None
            toret.append(ratio)
        return toret

    def observed_abundance_log_ratios(self, element1, element2):
        toret = list()
        for system in self.population:
            if (not system.observed) or (system.observed_abundances is None):
                toret.append(None)
            else:
                ratio = system.observed_abundance_log_ratio(element1, element2)
                el1 = system.observed_abundances.get(element1)
                el2 = system.observed_abundances.get(element2)
                try:
                    ratio = el1 - el2
                except TypeError:
                    ratio = None
                toret.append(ratio)
        return toret

    def observed_abundance_bools(self, element1, element2):
        # Returns a list of (bool, bool) tuples
        # The bools indicate whether element1/element 2 were observed respectively
        toret = list()
        for system in self.population:
            if (not system.observed) or (system.observed_abundances is None):
                toret.append((False, False))
            else:
                el1 = system.observed_abundances.get(element1)
                el2 = system.observed_abundances.get(element2)
                if el1 is None and el2 is None:
                    toret.append((False, False))
                elif el1 is None and el2 is not None:
                    toret.append((False, True))
                elif el1 is not None and el2 is None:
                    toret.append((True, False))
                else:
                    toret.append((True, True))
        return toret

    def modelled_values(self, parameter, collapse=False):
        toret = list()
        for system in self.population:
            if system.modelled_properties is None:
                toret.append(None)
            else:
                if system.modelled_properties.get(parameter) is None:
                    toret.append(None)
                else:
                    if collapse:
                        if isinstance(system.modelled_properties.get(parameter), list):
                            for sampled_val in system.modelled_properties.get(parameter):
                                toret.append(sampled_val)
                        else:
                            toret.append(system.modelled_properties.get(parameter))
                    else:
                        toret.append(system.modelled_properties.get(parameter))
        return toret

    def modelled_io_pairs(self, parameter):
        i_vals = list()
        o_vals = list()
        unmodelled_i_vals = list()
        for system in self.population:
            input_val = system.pollution_properties.get(parameter)
            output_val = None
            if system.modelled_properties is not None:
                output_val = system.modelled_properties.get(parameter)
            if output_val in [np.nan, None]:
                unmodelled_i_vals.append(input_val)
            elif isinstance(output_val, list) and all(ov in [np.nan, None] for ov in output_val):
                unmodelled_i_vals.append(input_val)
            else:
                i_vals.append(input_val)
                o_vals.append(output_val)
        assert len(i_vals) == len(o_vals)
        return i_vals, o_vals, unmodelled_i_vals

    def find_systems_with_abundances(self, element_dictionary):
        toret = list()
        for system in self.population:
            match = True
            for element, value in element_dictionary.items():
                if system.pollution_abundances.get(element, value) != value:
                    match = False
            if match:
                toret.append(system)
        return toret

    def get_all_ids(self):
        toret = list()
        for system in self.population:
            toret.append(system.id)
        return toret

    def find_systems_with_id(self, id_no):
        toret = list()
        for system in self.population:
            if system.id == id_no:
                toret.append(system)
        if len(toret) == 0:
            return None
        if len(toret) > 1:
            print('Warning! Found ' + str(len(toret)) + ' systems with id ' + str(id_no) + ', returning the first one')
        return toret[0]

    def assess_core_mantle_trustworthiness(self):
        list_of_cr_tuples = list()
        mantle_drop_outs = 0
        core_drop_outs = 0
        for system in self.population:
            is_core_rich = system.pollution_properties[mp.ModelParameter.fragment_core_frac] > 0.17
            appears_core_rich = None
            if system.modelled_properties is not None:
                min_fcnf = min(system.modelled_properties[mp.ModelParameter.fragment_core_frac])
                max_fcnf = max(system.modelled_properties[mp.ModelParameter.fragment_core_frac])
                if min_fcnf is None or max_fcnf is None:
                    pass
                else:
                    if min_fcnf > 0.17:
                        appears_core_rich = True
                    elif max_fcnf < 0.17:
                        appears_core_rich = False
                    else:
                        pass
            if appears_core_rich is None:
                # Then the system droppped out
                if is_core_rich:
                    core_drop_outs += 1
                else:
                    mantle_drop_outs += 1
            cr_tuple = (is_core_rich, appears_core_rich)
            list_of_cr_tuples.append(cr_tuple)
        true_positives = 0
        false_positives = 0
        true_negatives = 0
        false_negatives = 0
        for cr_tuple in list_of_cr_tuples:
            # Avoiding 'else' in this section because there might be Nones
            if cr_tuple[1] is True:
                # It looked core-rich...
                if cr_tuple[0] is True:
                    # ... and it was
                    true_positives += 1
                if cr_tuple[0] is False:
                    # ... but it wasn't
                    false_positives += 1
            if cr_tuple[1] is False:
                # It looked mantle-rich...
                if cr_tuple[0] is True:
                    # ... but it wasn't
                    false_negatives += 1
                if cr_tuple[0] is False:
                    # ... and it was
                    true_negatives += 1
        print()
        print('Testing for Core-Richness')
        print('True Positives: ' + str(true_positives))
        print('False Positives: ' + str(false_positives))
        print('Total Positives: ' + str(false_positives+true_positives))
        print('True Negatives: ' + str(true_negatives))
        print('False Negatives: ' + str(false_negatives))
        print('Total Negatives: ' + str(true_negatives+false_negatives))
        try:
            print('Core-Rich Trustworthiness: ' + str(true_positives/(true_positives+false_positives)))
            print('Mantle-Rich Trustworthiness: ' + str(true_negatives/(true_negatives+false_negatives)))
            print('Apparent Core-rich proportion: ' + str((false_positives+true_positives)/(false_positives+true_positives+true_negatives+false_negatives)))
        except ZeroDivisionError:
            print('Could not calculate Core-Rich/Mantle-Rich Trustworthiness')
        print('Core-rich dropouts: ' + str(core_drop_outs))
        print('Mantle-rich dropouts: ' + str(mantle_drop_outs))

    def get_observed_subset(self):
        subset = list()
        for system in self.population:
            if system.observed:
                subset.append(system)
        return SyntheticPopulation(None, self.wd_config_to_use, self.pollution_config_to_use, None, subset, self.timescale_interpolator)

    def get_subset_with_detected_elements(self, list_of_elements):
        # Pick out only systems with detections of all elements listed in list_of_elements
        subset = list()
        hack_to_exclude_systems_in_dec = False
        for system in self.population:
            #if system.check_elements_detected(list_of_elements):
                #subset.append(system)
            if system.check_elements_detected(list_of_elements):
                if hack_to_exclude_systems_in_dec:
                    if system.pollution_properties[mp.ModelParameter.t_sinceaccretion]*1000000 < (system.pollution_properties[mp.ModelParameter.accretion_timescale] + 3000000):
                        subset.append(system)
                else:
                    subset.append(system)
        if hack_to_exclude_systems_in_dec:
            print('Warning! Current sample selection excludes systems > 3 Myr into declining phase!')
        return SyntheticPopulation(None, self.wd_config_to_use, self.pollution_config_to_use, None, subset, self.timescale_interpolator)

    def get_subset_with_modelled_properties(self, list_of_properties):
        # Pick out only systems with detections of all elements listed in list_of_elements
        subset = list()
        for system in self.population:
            if system.check_properties_modelled(list_of_properties):
                subset.append(system)
        return SyntheticPopulation(None, self.wd_config_to_use, self.pollution_config_to_use, None, subset, self.timescale_interpolator)

    def get_ca_noise_threshold(self, pollution_frac):
        # Leftover from an earlier, more complicated function...
        return 0

    def get_most_polluted_subset(self, subset_size=None, element=None, impose_noise_pol_frac_correlation=False):
        # By default, we'll measure pollution by the (true) pollution fraction. (Ideally, should probably calculate would you would infer the pollution fraction to be based on the observed elements)
        # Supplying an element as the element argument means we'll instead look for the highest observed values of that element
        if subset_size is None:
            subset_size = np.inf
        subset_size = min(subset_size, len(self))
        if element is None:
            sorted_pop = sorted(self.population, key=lambda x: x.pollution_properties[mp.ModelParameter.pollution_frac], reverse=True)
        else:
            sorted_pop = sorted(self.population, key=lambda x: x.observed_abundances.get(element, -np.inf), reverse=True)
        if impose_noise_pol_frac_correlation:
            # Experimental! In order to try and reproduce the Hollands Ca/Fe dist, filter out systems where the Ca was overestimated
            # This is a bit of a hacky way to do it post hoc - a better way would of course be to impose this at the point of synthesis
            # The name impose_noise_pol_frac_correlation is a misnomer now
            sorted_pop = list(filter(lambda x: x.observed_abundances.get(ci.Element.Ca, np.inf) - x.pollution_abundances[ci.Element.Ca] < self.get_ca_noise_threshold(x.pollution_properties[mp.ModelParameter.pollution_frac]), sorted_pop))
        subset = sorted_pop[0:subset_size]
        return SyntheticPopulation(None, self.wd_config_to_use, self.pollution_config_to_use, None, subset, self.timescale_interpolator)

    def get_random_subset(self, subset_size=None, prevent_reuse=False):
        if subset_size is None:
            subset_size = np.inf
        subset_size = min(subset_size, len(self))
        #We're going round the houses (i.e. sampling the indices, rather than the population itself) so that we can maintain order
        if prevent_reuse:
            indices = self.unsampled_indices
        else:
            indices = range(len(self))
        if len(indices) < subset_size:
            print('Could not extract subset, requested size (' + str(subset_size) + ') > available systems (' + str(len(indices)) + ')')
            return None
        random_sample = random.sample(indices, k=subset_size)
        subset = list()
        i = 0
        for system in self.population:
            if i in random_sample:
                subset.append(system)
            i += 1
        self.unsampled_indices = [i for i in self.unsampled_indices if i not in random_sample]  # For now, let's always keep track of this regardless of prevent_reuse
        return SyntheticPopulation(None, self.wd_config_to_use, self.pollution_config_to_use, None, subset, self.timescale_interpolator)

    def reset_unsampled_indices(self):
        self.unsampled_indices = list(range(len(self)))

    def knn_comparison_against_other_pop(self, other, nearest_neighbours=3):
        correct_classification = 0
        incorrect_classification = 0
        assert len(self.population) == len(other.population)
        for pop_index, comparison_pop in enumerate([self.population, other.population]):
            for system in comparison_pop:
                if system.id is None:
                    raise ValueError('System ID must not be None')
                pop_distance_tuples = list() # [(bool, numerical) ... ] where bool is whether it was from 'self', numerical is the distance
                if pop_index == 0:
                    id_to_ignore = random.choice(other.get_all_ids()) # Pick an id to ignore randomly to even up the sample sizes (since we cannot compare a data point to itself)
                else:
                    id_to_ignore = random.choice(self.get_all_ids())
                for other_system in self.population:
                    if pop_index == 0 and system.id == other_system.id:
                        pass
                    elif pop_index == 1 and other_system.id == id_to_ignore:
                        pass
                    else:
                        distance = system.distance_to(other_system)
                        pop_distance_tuples.append((True, distance))
                for other_system in other.population:
                    if pop_index == 1 and system.id == other_system.id:
                        pass
                    elif pop_index == 0 and other_system.id == id_to_ignore:
                        pass
                    else:
                        distance = system.distance_to(other_system)
                        pop_distance_tuples.append((False, distance))
                random.shuffle(pop_distance_tuples) # Just in case of any ties when sorting, the outcome of the tie breaker is now random rather than going in the order added
                nearest_id_distance_tuples = sorted(pop_distance_tuples, key=lambda x: x[1])[0:nearest_neighbours] # Sort by distance, take closest ones
                votes_for_self = 0
                votes_for_other = 0
                for nidt in nearest_id_distance_tuples:
                    if nidt[0]:
                        votes_for_self += 1
                    else:
                        votes_for_other += 1
                if votes_for_self > votes_for_other: # nearest_neighbours should be odd, so no need to worry about draws
                    if pop_index == 0:
                        correct_classification += 1
                    else:
                        incorrect_classification += 1
                else:
                    if pop_index == 0:
                        incorrect_classification += 1
                    else:
                        correct_classification += 1
        return correct_classification, incorrect_classification, SyntheticPopulation.get_knn_p_value(correct_classification, incorrect_classification)

    @staticmethod
    def get_knn_p_value(correct_classification, incorrect_classification):
        N = correct_classification + incorrect_classification
        if correct_classification < 0 or incorrect_classification < 0 or N <= 0:
            return None
        # If this p_value is less than 0.5, it suggests that the two populations were different to some extent (since we were able to tell which group a randomly selected point belonged to (better than guessing))
        # If p_value > 0.5, the classification is mostly incorrect which is a bit odd - it suggests that the populations are non-randomly interleaved to some extent
        p_value = 1 - binom.cdf(correct_classification - 1, N, 0.5)
        return p_value

def generate_hollands_cool_dz_sample():
    list_of_systems = list()
    wd_properties_list, abundances_list = ha.get_hollands_all_el_values()
    assert len(wd_properties_list) == len(abundances_list)
    for wd_index, wd_properties in enumerate(wd_properties_list):
        abundances = abundances_list[wd_index]
        list_of_systems.append(
            SyntheticSystem(
                wd_properties,
                None,
                None,
                True,
                abundances,
                None,
                None,
                wd_index
            )
        )
    synth_pop = SyntheticPopulation(
        None,
        None,
        None,
        None,
        list_of_systems
    )
    return synth_pop

def demonstrate_synthesis():
    synth_pop = SyntheticPopulation(
        3,
        'HollandsDBs',
        'DBDeltaFcf'
    )
    print(synth_pop)

def main():
    demonstrate_synthesis()

if __name__ == '__main__':
    main()
