#!/usr/bin/env python
# -*- coding: utf-8 -*-

from argparse import Namespace

import csv
import json
import os
import numpy as np
import pymultinest as pn

import chemistry_info as ci
import manager as mn
import model_parameters as mp
import pwd_utils as pu

def load_manager():
    manager = mn.Manager(  #This is just to load the stellar compositions
        Namespace(
            wd_data_filename='WDInputData.csv',
            stellar_compositions_filename='StellarCompositionsSortFE.csv',
            n_live_points = 0,
            pollution_model_names=['Model_24'],
            enhancement_model='Earthlike',
            base_dir=pu.get_path_to_data()
        )
    )
    return manager

def get_hollands_ids():
    hollands_ids = list(range(0, 201))
    hollands_ids.remove(85)
    hollands_ids.extend([249, 250])
    return hollands_ids

def get_hollands_names(manager=None):
    if manager is None:
        manager = load_manager()
    hollands_names = list()
    i = 0
    hollands_ids = get_hollands_ids()
    for wd in manager.white_dwarfs:
        if i in hollands_ids:
            hollands_names.append(wd.name)
        i += 1
    return hollands_names

def get_hollands_abundance_values(element, manager=None):
    if manager is None:
        manager = load_manager()
    hollands_abundance_values = list()
    i = 0
    hollands_ids = get_hollands_ids()
    for wd in manager.white_dwarfs:
        if i in hollands_ids:
            abundance_info = wd.get_abundance(element)
            if abundance_info is not None and abundance_info.included:
                hollands_abundance_values.append(abundance_info.value)
        i += 1
    return hollands_abundance_values

def get_hollands_abundances_for_dt_plot(element, manager=None):
    if manager is None:
        manager = load_manager()
    hollands_abundance_tuples = list()
    i = 0
    hollands_ids = get_hollands_ids()
    for wd in manager.white_dwarfs:
        if i in hollands_ids:
            abundance_info = wd.get_abundance(element)
            if abundance_info is not None and abundance_info.included:
                hollands_abundance_tuples.append((wd.get_teff().value, abundance_info.value))
        i += 1
    return hollands_abundance_tuples, None

def find_json_params_files(wd_index):
    json_files = list()
    for filename in os.listdir(pu.get_path_to_historical_output_dir()):
        if filename.endswith('.json') and filename.startswith(str(wd_index) + 'model'):
            json_files.append(filename)
    json_files.sort()
    return json_files

def extract_posteriors(wd_name):
    # Step 1: find out what the index of that wd_name is (in the wd_data_1112.csv file)
    # Step 2: Read the best_fits file to find out what the parameters for the best model was
    # Step 3: Look up what that model would have been called in PWDCode.py (nb: not necessarily the same as the model name listed in the best fits file) by finding the matching params file
    # Step 4: Read from file <wd index>model<model number>post_equal_weights.dat, with the n_params set to the model dimensions
    wd_index = None
    with open(pu.get_path_to_original_src() + 'wd_data_1112.csv', encoding='utf-8') as wdcsv:
        row_count = 0
        for row in csv.reader(wdcsv):
            if row[0] == wd_name:
                wd_index = row_count
                break
            row_count += 1
    model_parameters = list()
    with open(pu.get_path_to_historical_output_dir() + 'best_fits_hb20.csv', encoding='utf-8') as bfcsv:
        for row in csv.reader(bfcsv):
            if row[0] == wd_name:
                # The last column has the parameters
                model_parameters_as_text = row[-1]
                model_parameters_raw = model_parameters_as_text[1:-1].split(',')
                model_parameters = [model_parameter.replace('\'','').lstrip().rstrip() for model_parameter in model_parameters_raw]
                break
    json_files = find_json_params_files(wd_index)
    matching_json_files = list()
    matching_models = list()
    true_model_parameters = None
    for json_file in json_files:
        with open(pu.get_path_to_historical_output_dir() + json_file) as f:
            params = json.load(f)
            if sorted(params) == sorted(model_parameters):
                matching_json_files.append(json_file)
                true_model_parameters = params
    for matching_json_file in matching_json_files:
        model_no = int(matching_json_file.split('model')[1].split('params')[0].replace('bu', '').replace('ss', '').replace('h', '').replace('i', ''))
        if model_no not in matching_models:
            matching_models.append(model_no)
    assert len(matching_models) == 1  # Not sure what to do if this isn't true
    model_basename = str(wd_index) + 'model' + str(matching_models[0])
    n_dims = len(true_model_parameters)
    print()
    print(wd_name)
    print(model_basename)
    print(n_dims)
    a = pn.Analyzer(n_params = n_dims, outputfiles_basename = pu.get_path_to_historical_output_dir() + model_basename)
    stats = a.get_stats()
    weightpost = a.get_equal_weighted_posterior()[:, 0:n_dims]  # This is basically just excluding the final column of the ...post_equal_weights.dat file
    best_fit = a.get_best_fit()
    return stats, weightpost, best_fit['log_likelihood'], best_fit['parameters'], true_model_parameters, matching_models[0]

def get_median_t_t_event_myr(wd_name):
    stats, weightpost, log_likelihood, bes_fit_parameter_values, model_parameters, model_number = extract_posteriors(wd_name)
    t_sinceaccretion_index = model_parameters.index('Time since Accretion/Myrs')
    accretion_timescale_index = model_parameters.index('log(Accretion Event Timescale/Yrs)')
    t_sinceaccretion = np.percentile(weightpost[:,t_sinceaccretion_index], 50)
    accretion_timescale = (10**np.percentile(weightpost[:,accretion_timescale_index], 50))/1000000
    return t_sinceaccretion, accretion_timescale

def get_hollands_abundances_and_timescales_dict(manager=None):
    if manager is None:
        manager = load_manager()
    hollands_abundance_dict = dict()
    i = 0
    hollands_ids = get_hollands_ids()
    for wd in manager.white_dwarfs:
        if i in hollands_ids:
            t_sinceaccretion, accretion_timescale = get_median_t_t_event_myr(wd.name)
            hollands_abundance_dict[wd.name] = {
                'Median time since accretion /Myr': t_sinceaccretion,
                'Median accretion timescale /Myr': accretion_timescale
            }
            for el, abundance in wd.get_abundance_values_dict().items():
                hollands_abundance_dict[wd.name][el] = abundance
        i += 1
    return hollands_abundance_dict

def get_hollands_pol_fracs(manager=None):
    if manager is None:
        manager = load_manager()
    hollands_pol_frac_dict = dict()
    i = 0
    hollands_ids = get_hollands_ids()
    for wd in manager.white_dwarfs:
        if i in hollands_ids:
            hollands_pol_frac_dict[wd.name] = wd.estimate_minimum_pollution_fraction()
        i += 1
    return hollands_pol_frac_dict

def get_hollands_el_el_values(element1, element2, manager=None):
    if manager is None:
        manager = load_manager()
    hollands_el_el_values = list()
    i = 0
    hollands_ids = get_hollands_ids()
    for wd in manager.white_dwarfs:
        if i in hollands_ids:
            hollands_el_el_values.append(wd.get_log_ratio(element1, element2))
        i += 1
    return hollands_el_el_values

def get_hollands_common_el_values(list_of_elements, manager=None):
    # Returns lists of element abundances for all wds that have ALL the specified elements, such that the length of each list is the same
    if manager is None:
        manager = load_manager()
    toret = list()
    for element in list_of_elements:
        toret.append(list())
    i = 0
    hollands_ids = get_hollands_ids()
    for wd in manager.white_dwarfs:
        if i in hollands_ids:
            elements_present = wd.get_elements_present()
            include = True
            for element in list_of_elements:
                if element not in elements_present:
                    include = False
            if include:
                for j, element in enumerate(list_of_elements):
                    toret[j].append(wd.get_abundance(element).value)
        i += 1
    return toret

def get_hollands_all_el_values(manager=None):
    # Returns dict of element abundances for all wds
    if manager is None:
        manager = load_manager()
    wd_properties = list()
    abundances = list()
    i = 0
    hollands_ids = get_hollands_ids()
    for wd in manager.white_dwarfs:
        if i in hollands_ids:
            abundances.append(wd.get_abundance_values_dict())
            wd_properties.append({
                mp.WDParameter.spectral_type: wd.get_spectral_type().value,
                mp.WDParameter.temperature: wd.get_teff().value,
                mp.WDParameter.logg: wd.get_logg().value
            })
        i += 1
    return wd_properties, abundances

def main():
    manager = load_manager()
    names = get_hollands_names(manager)
    print(names)
    print(len(names))
    fe_values = get_hollands_abundance_values(ci.Element.Fe, manager)
    print(fe_values)
    print(max(fe_values))
    ca_values = get_hollands_abundance_values(ci.Element.Ca, manager)
    print(ca_values)
    print(max(ca_values))
    mg_values = get_hollands_abundance_values(ci.Element.Mg, manager)
    print(mg_values)
    print(max(mg_values))
    max_pol_frac = -100
    pol_frac_dict = get_hollands_pol_fracs()
    for k, v in pol_frac_dict.items():
        if v > max_pol_frac:
            max_pol_frac = v
    print(max_pol_frac)
    toprint = get_hollands_el_el_values(ci.Element.Ca, ci.Element.Fe, manager)
    print(toprint)
    print(len(toprint))
    toprint2 = get_hollands_el_el_values(ci.Element.Mg, ci.Element.Fe, manager)
    print(toprint2)
    print(len(toprint2))
    cr_values = get_hollands_abundance_values(ci.Element.Cr, manager)
    print(cr_values)
    print(len(cr_values))
    ni_values = get_hollands_abundance_values(ci.Element.Ni, manager)
    print(ni_values)
    print(len(ni_values))
    print(get_hollands_common_el_values([ci.Element.Ca, ci.Element.Fe, ci.Element.Mg, ci.Element.Cr]))

if __name__ == '__main__':
    main()
