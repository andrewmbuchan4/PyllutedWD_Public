#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np
import os

import graph_factory as gf
import pwd_utils as pu

def find_stats_files(path, system_prefixes):
    stats_files = list()
    for system_prefix in system_prefixes:
        for directory_or_file in os.listdir(path):
            if os.path.isdir(path + directory_or_file) and directory_or_file.startswith(system_prefix):
                for filename in os.listdir(path + directory_or_file):
                    if filename.endswith('stats.csv'):
                        stats_files.append(path + directory_or_file + '/' + filename)
    stats_files.sort()
    return stats_files

def extract_retrieval_data(path, system_prefixes): # This function could do with some revision
    stats_files = find_stats_files(path, system_prefixes)
    data_dict = dict()
    for stat_file in stats_files:
        system = None
        pressure_5 = None  # = 5th percentile, i.e. 2 sigma lower limit
        pressure_16 = None  # = 16th percentile, i.e. 1 sigma lower limit
        pressure_50 = None  # = 50th percentile, i.e. median
        pressure_84 = None  # = 84th percentile, i.e. 1 sigma upper limit
        pressure_95 = None  # = 95th percentile, i.e. 2 sigma upper limit
        pressure_index = None
        fcf_5 = None  # = 5th percentile, i.e. 2 sigma lower limit
        fcf_16 = None  # = 16th percentile, i.e. 1 sigma lower limit
        fcf_50 = None  # = 50th percentile, i.e. median
        fcf_84 = None  # = 84th percentile, i.e. 1 sigma upper limit
        fcf_95 = None  # = 95th percentile, i.e. 2 sigma upper limit
        fcf_index = None
        fO2_5 = None  # = 5th percentile, i.e. 2 sigma lower limit
        fO2_16 = None  # = 16th percentile, i.e. 1 sigma lower limit
        fO2_50 = None  # = 50th percentile, i.e. median
        fO2_84 = None  # = 84th percentile, i.e. 1 sigma upper limit
        fO2_95 = None  # = 95th percentile, i.e. 2 sigma upper limit
        fO2_index = None
        distance_5 = None  # = 5th percentile, i.e. 2 sigma lower limit
        distance_16 = None  # = 16th percentile, i.e. 1 sigma lower limit
        distance_50 = None  # = 50th percentile, i.e. median
        distance_84 = None  # = 84th percentile, i.e. 1 sigma upper limit
        distance_95 = None  # = 95th percentile, i.e. 2 sigma upper limit
        distance_index = None
        pcf_50 = None
        with open(stat_file, encoding='utf-8') as output_csv:
            print()
            print('Reading ' + stat_file)
            reading_from_correct_section = False
            for row in csv.reader(output_csv):
                if len(row) > 0:
                    if row[0] == 'System Name:':
                        system = row[1]
                    if row[0] == 'Best model name:':
                        best_model_name = row[1]
                        print('Best model was ' + best_model_name)
            output_csv.seek(0)
            for row in csv.reader(output_csv):
                if len(row) > 0:
                    if row[0] == best_model_name:
                        good_fit = row[12] == 'True'
                        if not good_fit:
                            print('WARNING! Fit was not good')
                        else:
                            print('(Good fit!)')
                    if row[0] == 'Results from model:':
                        if row[1] == best_model_name:
                            reading_from_correct_section = True
                        else:
                            reading_from_correct_section = False
                    if reading_from_correct_section:
                        if row[0] == 'Parameter:':
                            try:
                                fcf_index = row.index('Fragment Core Fraction')
                            except ValueError:
                                fcf_index = None
                            try:
                                fO2_index = row.index('Oxygen Fugacity /ΔIW')
                            except ValueError:
                                fO2_index = None
                            try:
                                pressure_index = row.index('Pressure /GPa')
                            except ValueError:
                                pressure_index = None
                            try:
                                distance_index = row.index('log(Formation Distance/AU)')
                            except ValueError:
                                distance_index = None
                        elif row[0] == '5th percentile:':
                            if fcf_index is not None:
                                fcf_5 = float(row[fcf_index])
                            if fO2_index is not None:
                                fO2_5 = float(row[fO2_index])
                            if pressure_index is not None:
                                pressure_5 = float(row[pressure_index])
                            if distance_index is not None:
                                distance_5 = float(row[distance_index])
                        elif row[0] == '16th percentile:':
                            if fcf_index is not None:
                                fcf_16 = float(row[fcf_index])
                            if fO2_index is not None:
                                fO2_16 = float(row[fO2_index])
                            if pressure_index is not None:
                                pressure_16 = float(row[pressure_index])
                            if distance_index is not None:
                                distance_16 = float(row[distance_index])
                        elif row[0] == 'Median:':
                            if fcf_index is not None:
                                fcf_50 = float(row[fcf_index])
                            if fO2_index is not None:
                                fO2_50 = float(row[fO2_index])
                            if pressure_index is not None:
                                pressure_50 = float(row[pressure_index])
                            if distance_index is not None:
                                distance_50 = float(row[distance_index])
                        elif row[0] == '84th percentile:':
                            if fcf_index is not None:
                                fcf_84 = float(row[fcf_index])
                            if fO2_index is not None:
                                fO2_84 = float(row[fO2_index])
                            if pressure_index is not None:
                                pressure_84 = float(row[pressure_index])
                            if distance_index is not None:
                                distance_84 = float(row[distance_index])
                        elif row[0] == '95th percentile:':
                            if fcf_index is not None:
                                fcf_95 = float(row[fcf_index])
                            if fO2_index is not None:
                                fO2_95 = float(row[fO2_index])
                            if pressure_index is not None:
                                pressure_95 = float(row[pressure_index])
                            if distance_index is not None:
                                distance_95 = float(row[distance_index])
                        elif row[0] == 'Parent Core Number Fraction:':
                            try:
                                pcf_50 = float(row[1])
                            except IndexError:
                                pass  # There was no differentiation
        if system is not None:
            data_dict[system] = {
                'Pressure': {'5': pressure_5, '16': pressure_16, '50': pressure_50, '84': pressure_84, '95': pressure_95},
                'Fragment Core Fraction': {'5': fcf_5, '16': fcf_16, '50': fcf_50, '84': fcf_84, '95': fcf_95},
                'Oxygen Fugacity': {'5': fO2_5, '16': fO2_16, '50': fO2_50, '84': fO2_84, '95': fO2_95},
                'Formation Distance': {'5': distance_5, '16': distance_16, '50': distance_50, '84': distance_84, '95': distance_95},
                'Parent Core Fraction': {'50': pcf_50}
            }
    return data_dict

def preprocess_data(data, systems=['Earthfcf', 'Marsfcf']):
    # Extract the actual data points to plot
    to_plot = dict()
    for system in systems:
        to_plot[system] = dict()
        to_plot[system]['fcfs'] = list()
        variables = ['Pressure', 'Fragment Core Fraction', 'Oxygen Fugacity', 'Formation Distance', 'Parent Core Fraction']
        for variable in variables:
            to_plot[system][variable] = {'1s_upper': list(), '1s_lower': list(), '2s_upper': list(), '2s_lower': list(), 'medians': list()}
    for system, system_data in data.items():
        for test_system in systems:
            print()
            print(139)
            print(test_system)
            if system.startswith(test_system):
                # Then this is one we care about
                fragment_core_fraction_str = system.split(test_system + 'Fcf')[1] # 0p4 or something like that
                if 'p' in fragment_core_fraction_str:
                    fragment_core_fraction = float(fragment_core_fraction_str.split('p')[0] + '.' + fragment_core_fraction_str.split('p')[1])
                else:
                    fragment_core_fraction = float(fragment_core_fraction_str)
                print(fragment_core_fraction)
                to_plot[test_system]['fcfs'].append(fragment_core_fraction)
                for variable in variables:
                    print(variable)
                    value_5 = system_data[variable].get('5')
                    value_16 = system_data[variable].get('16')
                    value_50 = system_data[variable].get('50')
                    value_84 = system_data[variable].get('84')
                    value_95 = system_data[variable].get('95')
                    if value_50 is None:
                        to_plot[test_system][variable]['medians'].append(np.NaN)
                    else:
                        to_plot[test_system][variable]['medians'].append(value_50)
                    print(to_plot[test_system][variable]['medians'])
                    print()
                    print(to_plot[test_system]['fcfs'])
                    print()
                    if value_84 is not None and value_50 is not None:
                        to_plot[test_system][variable]['1s_upper'].append(value_84 - value_50)
                    else:
                        to_plot[test_system][variable]['1s_upper'].append(np.NaN)
                    if value_16 is not None and value_50 is not None:
                        to_plot[test_system][variable]['1s_lower'].append(value_50 - value_16)
                    else:
                        to_plot[test_system][variable]['1s_lower'].append(np.NaN)
                    if value_95 is not None and value_50 is not None:
                        to_plot[test_system][variable]['2s_upper'].append(value_95 - value_50)
                    else:
                        to_plot[test_system][variable]['2s_upper'].append(np.NaN)
                    if value_5 is not None and value_50 is not None:
                        to_plot[test_system][variable]['2s_lower'].append(value_50 - value_5)
                    else:
                        to_plot[test_system][variable]['2s_lower'].append(np.NaN)
    return to_plot

def plot_retrieved_stats():
    #systems = ['RealEarthMix', 'RealMarsMix']
    systems = ['SyntheticEarthMix', 'SyntheticMarsMix']
    data = extract_retrieval_data(pu.get_path_to_pylluted_dir(), systems)
    preprocessed_data = preprocess_data(data, systems)
    graph_fac = gf.GraphFactory()
    for variable in ['Pressure', 'Fragment Core Fraction', 'Oxygen Fugacity', 'Formation Distance', 'Parent Core Fraction']:
        earth_plot = graph_fac.plot_retrieved_variable_against_fcf(
            'Earth',
            variable,
            preprocessed_data[systems[0]]['fcfs'],
            preprocessed_data[systems[0]][variable]['medians'],
            preprocessed_data[systems[0]][variable]['1s_upper'],
            preprocessed_data[systems[0]][variable]['1s_lower'],
            preprocessed_data[systems[0]][variable]['2s_upper'],
            preprocessed_data[systems[0]][variable]['2s_lower'],
            'Synthetic'
        )
        earth_plot['retrieved_variable_plot']['subplots']['subplot1']['y_hide_ticks'] = [0]
        
        mars_plot = graph_fac.plot_retrieved_variable_against_fcf(
            'Mars',
            variable,
            preprocessed_data[systems[1]]['fcfs'],
            preprocessed_data[systems[1]][variable]['medians'],
            preprocessed_data[systems[1]][variable]['1s_upper'],
            preprocessed_data[systems[1]][variable]['1s_lower'],
            preprocessed_data[systems[1]][variable]['2s_upper'],
            preprocessed_data[systems[1]][variable]['2s_lower'],
            'Synthetic'
        )
        mars_plot['retrieved_variable_plot']['subplots']['subplot1']['legend'] = False
        
        graph_fac.multipanelise([earth_plot, mars_plot], 2, 1, 'retrieved_' + variable + '_Synthetic_multipanel.pdf', 15, 10, 0.07, 0)

def main():
    plot_retrieved_stats()

if __name__ == '__main__':
    main()
