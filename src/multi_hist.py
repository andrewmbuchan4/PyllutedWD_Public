#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np
import os
import xlrd

import graph_factory as gf

#map external model numbers to internal model numbers
model_map_dict = {
    1: 2,
    4: 4,
    6: 27,
    2: 25,
    5: 5,
    3: 26,
    9: 21,
    8: 28
}

def get_wd_names():
    wdsdcsv = open('../original_codebase/wd_data_1112.csv')
    namelist =  [row[0] for row in csv.reader(wdsdcsv)]
    return namelist

def get_best_fits(path):
    wdsdcsv = open(path + 'best_fits_hb20.csv')
    toret = dict()
    i = 0
    for row in csv.reader(wdsdcsv):
        if i > 0:
            toret[row[0]] = row[1]
        i += 1
    return toret

def find_xlsx_files(path=None):
    xlsx_files = list()
    for filename in os.listdir(path):
        if filename.endswith('.xlsx'):
            xlsx_files.append(filename)
    xlsx_files.sort()
    return xlsx_files

def find_ewp_files(path, xlsx_files, wd_names, best_fits):
    # Translate each xlsx file into the relevant .dat file
    # Template: obs number + 'model' + model number + 'post_equal_weights.dat'
    ewp_files = dict()
    for filename in xlsx_files:
        wd_name = filename.split('PWD')[0]
        try:
            obs_number = wd_names.index(wd_name)
        except ValueError:
            obs_number = None
        if obs_number is None:
            continue # Not a Hollands WD
        if (obs_number > 200) and (obs_number not in [249, 250]):
            continue # Not a Hollands WD
        if obs_number == 85:
            continue # This is J1055, the spectroscopic binary
        print(filename)
        best_fit = best_fits[wd_name]
        external_model_number = int(best_fit.split(' =')[0].split('M')[1])
        internal_model_number = model_map_dict[external_model_number]
        ewp_files[wd_name] = path + str(obs_number) + 'model' + str(internal_model_number) + 'post_equal_weights.dat'
    print('Analysing ' + str(len(ewp_files)) + ' systems')
    assert len(ewp_files) == 202  # There are 202 WDs in our sample
    return ewp_files

def get_binned_temperature_stats(xlsx_files, wd_names, path=None):
    temp_stats = dict()
    for filename in xlsx_files:
        wd_name = filename.split('PWD')[0]
        try:
            obs_number = wd_names.index(wd_name)
        except ValueError:
            obs_number = None
        if obs_number is None:
            continue # Not a Hollands WD
        if obs_number > 200:
            continue # Not a Hollands WD
        if obs_number == 85:
            continue # This is J1055, the spectroscopic binary
        print(filename)
        workbook = xlrd.open_workbook(path + filename)
        sheet = workbook.sheet_by_index(0)
        first_bin_value = sheet.cell_value(12, 3) # D13
        if first_bin_value == '':
            print('No temp vals detected for ' + wd_name)
            continue  # No temp vals for this WD
        bin_centres_raw = sheet.col_values(3)
        temp_vals_raw = sheet.col_values(4)
        
        bin_heading_index = bin_centres_raw.index('Bin Value Temperature/K')
        bin_centres = [be for be in bin_centres_raw[bin_heading_index+1:] if be != '']

        temp_heading_index = temp_vals_raw.index('Bin Count Temperature/K')
        temp_vals = [tv for tv in temp_vals_raw[temp_heading_index+1:] if tv != '']
        temp_stats[wd_name] = (bin_centres, temp_vals)
    return temp_stats

def get_all_accretion_lifetime_stats(files):
    all_values = dict()
    for wd_name, file_name in files.items():
        print(file_name)
        post_equal_weights = np.loadtxt(file_name, ndmin=2)
        system_vals = post_equal_weights[:, -2] #  I want the second to last column
        all_values[wd_name] = system_vals#(10**system_vals)/1000000
    #for wd, system_vals in all_values.items():
    #    for test_val in system_vals:
    #        if test_val > 8 or test_val < 0:
    #            raise # sanity check
    return all_values

def make_all_lifetime_hist_plot(all_values, combine=True):

    half_bin_size = 0.05#0.5
    bins, x_bar_centres = generate_bins_and_bar_centres(0, 8, half_bin_size)
    #bins, x_bar_centres = generate_bins_and_bar_centres(0, 100, half_bin_size)
    graph_fac = gf.GraphFactory()
    all_heights = list()
    wd_names = list()
    for wd_name, system_vals in all_values.items():
        heights, bins2 = np.histogram(
            system_vals,
            bins,
            density=True
        )
        wd_names.append(wd_name)
        all_heights.append(heights)
    if combine:
        #superlist = np.array([])
        #for wd_name, wd_values in all_values.items():
        #    print(wd_name)
        #    for value in wd_values:
        #        superlist = np.append(superlist, value)
        #    print(len(superlist))
        #print(superlist)
        #superlist = np.concatenate([10**(key_value_pair[1]) for key_value_pair in all_values.items()])
        #print(superlist)
        #print(len(superlist))
        #exp_mean = np.mean(superlist)
        #print(exp_mean)
        #print(np.log10(exp_mean))
        #exp_median = np.percentile(superlist, 50)
        #exp_sigma_upper = np.percentile(superlist, 84)
        #exp_sigma_lower = np.percentile(superlist, 16)
        #median = np.log10(exp_median)
        #print(median)
        #errorplus = np.log10(exp_sigma_upper - exp_median)
        #errorminus = np.log10(exp_median - exp_sigma_lower)
        all_exp_means = np.array([np.mean(10**key_value_pair[1]) for key_value_pair in all_values.items()])
        assert len(all_exp_means) == 202
        exp_median = np.percentile(all_exp_means, 50)
        exp_sigma_upper = np.percentile(all_exp_means, 84)
        exp_sigma_lower = np.percentile(all_exp_means, 16)
        print(exp_median)
        print(exp_sigma_upper)
        print(exp_sigma_lower)
        median = np.log10(exp_median)
        errorplus = np.log10(exp_sigma_upper) - np.log10(exp_median)
        errorminus = np.log10(exp_median) - np.log10(exp_sigma_lower)
        print(median)
        print(errorplus)
        print(errorminus)
        rounded_median_str = '%.2f' % median
        rounded_errorplus_str = '%.2f' % errorplus
        rounded_errorminus_str = '%.2f' % errorminus
        text_dict = {
            'median_text': {
                'x_pos': 7.9,
                'text_string': 'log(Accretion Event Lifetime /Yrs) = $' + rounded_median_str + ' ^{+' + rounded_errorplus_str + '}_{-' + rounded_errorminus_str + '}$',
                'horizontalalignment': 'right'
            }
        }
        averaged_heights = np.mean(all_heights, axis=0)  # This logic hopefully weights all the WDs equally
        graph_fac.make_histogram(x_bar_centres, [averaged_heights], ['Hollands et al. 2017 data'], 'bestmodel', half_bin_size*2, 1.1, 'log(Accretion Event Lifetime /Yrs)', '_acc_lifetime_agg_dist', text_dict, None)
        #graph_fac.make_histogram(x_bar_centres, [averaged_heights], ['Hollands et al. 2017 data'], 'bestmodel', half_bin_size*2, 1.1, 'Accretion Event Lifetime /Myr', '_acc_lifetime_agg_dist_lin', None, None)
    else:
        graph_fac.make_histogram(x_bar_centres, all_heights[0:4], wd_names[0:4], 'bestmodel', half_bin_size*2, 1.1, 'log(Accretion Event Lifetime /Yrs)', '_acc_lifetime_sep_dist', None, None)
        #graph_fac.make_histogram(x_bar_centres, [all_heights[0:4]], [wd_names[0:4]], 'bestmodel', half_bin_size*2, 1.1, 'Accretion Event Lifetime /Myr', '_acc_lifetime_sep_dist_lin', None, None)


def make_all_temperature_hist_plot(temp_stats, combine=True):
    # temp_stats structure: {WD name: (bins, values)}
    all_heights = list()
    bin_centres = None
    for wd_name, wd_stats in temp_stats.items():
        all_heights.append(wd_stats[1])
        bin_centres = wd_stats[0]  # All systems should have been binned the same way
    bin_size = bin_centres[1]-bin_centres[0]
    graph_fac = gf.GraphFactory()
    text_dict = {
        'icy_text': {
            'x_pos': 0,
            'text_string': 'Volatile\nrich',
            'horizontalalignment': 'center'
        },
        'dry_text': {
            'x_pos': 700,
            'text_string': 'Depleted in Volatiles',
            'horizontalalignment': 'center'
        },
        'eh_text': {
            'x_pos': 2175,
            'text_string': 'Depleted in Moderate Volatiles',
            'horizontalalignment': 'center'
        }
    }
    line_dict = {
        'vline1': {
            'x_start': 220
        },
        'vline2': {
            'x_start': 1250
        }
    }
    if combine:
        averaged_heights = np.mean(all_heights, axis=0)
        graph_fac.make_histogram(bin_centres, [averaged_heights], ['Hollands et al. 2017 data'], 'bestmodel', bin_size, 1.1, 'Temperature /K', '_temp_agg_dist', text_dict, line_dict)
    else:
        #wd_names = list(temp_stats.keys())
        selected_heights = list()
        wd_names = [
            #'SDSSJ1405+1549', # DV
            #'SDSSJ1340+2702', # DV
            #'SDSSJ1336+3547', # DV
            'SDSSJ0116+2050', #Not in Amy's list  # DV #B
            #'SDSSJ0047+1628', # DV
            #'SDSSJ0807+4930', #Not in Amy's list  H   Icy + DV
            'SDSSJ1234+5208', # H DV
            #'SDSSJ1024+1014', # DV into DMV
            'SDSSJ0916+2540', # H  Big DMV peak
            'SDSSJ1040+2407', # H  DMV (slightly less)
            #'SDSSJ0736+4118',  # DV into DMV
            #'SDSSJ1038-0036', # DMV
            #'SDSSJ1149+0519', #Not in Amy's list   Icy + DMV
            #'SDSSJ1405+1549', # DV
            #'SDSSJ1411+3410'  # DV into DMV
        ]
        for wd_name in wd_names:
            selected_heights.append(temp_stats[wd_name][1])
        graph_fac.make_histogram(bin_centres, selected_heights, wd_names, 'bestmodel', bin_size, 1.1, 'Temperature /K', '_temp_sep_dist', text_dict, line_dict)

def generate_bins_and_bar_centres(min_bin_edge, max_bin_edge, half_bin_size):
    bin_size = 2 * half_bin_size
    bins = np.arange(min_bin_edge, max_bin_edge + bin_size, bin_size)
    bar_centres = [x + half_bin_size for x in bins]
    bar_centres.pop()
    return bins, bar_centres

def main():
    path = pu.get_path_to_historical_output_dir()
    wd_names = get_wd_names()
    best_fits = get_best_fits(path)
    xlsx_files = find_xlsx_files(path)
    make_temp_plot = False
    make_tacc_plot = True
    if make_temp_plot:
        temp_stats = get_binned_temperature_stats(xlsx_files, wd_names, path)
        make_all_temperature_hist_plot(temp_stats)
        make_all_temperature_hist_plot(temp_stats, False)
    if make_tacc_plot:
        ewp_files = find_ewp_files(path, xlsx_files, wd_names, best_fits)
        stats = get_all_accretion_lifetime_stats(ewp_files)
        make_all_lifetime_hist_plot(stats)
        make_all_lifetime_hist_plot(stats, False)

if __name__ == '__main__':
    main()
