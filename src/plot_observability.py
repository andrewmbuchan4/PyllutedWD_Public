#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np

import chemistry_info as ci
import complete_model as cm
import graph_factory as gf
import manager as mn

from argparse import Namespace

pressure_args = {
    ci.Element.Cr:{
        #(0, 0): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 0, -2],
        #(0, 30): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 30, -2],
        #(0, 60): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 60, -2],
        #(0.1, 0): [458, 1.5, 0, 0.05, None, None, 0.1, None, -4, 6.5, 0, -2],
        #(0.1, 30): [458, 1.5, 0, 0.05, None, None, 0.1, None, -4, 6.5, 30, -2],
        #(0.1, 60): [458, 1.5, 0, 0.05, None, None, 0.1, None, -4, 6.5, 60, -2],
        #(0.99, 0): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 0, -2],
        #(0.99, 30): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 30, -2],
        #(0.99, 60): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 60, -2]
        
        (0.99, 0): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 0, -2],
        #(0.99, 15): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 15, -2],
        (0.99, 60): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 60, -2]
    },
    ci.Element.Ni:{
        (0, 0): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 0, -2],
        (0, 30): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 30, -2],
        (0, 60): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 60, -2]
    },
    ci.Element.Si:{
        (0.99, 0): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 0, -2],
        (0.99, 30): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 30, -2],
        (0.99, 60): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 60, -2]
    }
}
    
fcf_args = {
    ci.Element.Cr:{
    #(0, 0): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 0, -2],
    #(0, 30): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 30, -2],
    #(0, 60): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 60, -2],
    #(0.1, 0): [458, 1.5, 0, 0.05, None, None, 0.1, None, -4, 6.5, 0, -2],
    #(0.1, 30): [458, 1.5, 0, 0.05, None, None, 0.1, None, -4, 6.5, 30, -2],
    #(0.1, 60): [458, 1.5, 0, 0.05, None, None, 0.1, None, -4, 6.5, 60, -2],
    #(0.99, 0): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 0, -2],
    #(0.99, 30): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 30, -2],
    #(0.99, 60): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 60, -2]
        
        (0, 30): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 30, -2],
        (0.75, 30): [458, 1.5, 0, 0.05, None, None, 0.75, None, -4, 6.5, 30, -2],
        (0.99, 30): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 30, -2]
    },
    ci.Element.Ni:{
        (0, 30): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 30, -2],
        (0.1, 30): [458, 1.5, 0, 0.05, None, None, 0.1, None, -4, 6.5, 30, -2],
        (0.99, 30): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 30, -2]
    },
    ci.Element.Si:{
        (0, 30): [458, 1.5, 0, 0.05, None, None, 0, None, -4, 6.5, 30, -2],
        #(0.3, 0): [458, 1.5, 0, 0.05, None, None, 0.3, None, -4, 6.5, 0, -2],
        (0.75, 30): [458, 1.5, 0, 0.05, None, None, 0.75, None, -4, 6.5, 30, -2],
        (0.99, 30): [458, 1.5, 0, 0.05, None, None, 0.99, None, -4, 6.5, 30, -2]
    }
}

contour_args = {
    ci.Element.Cr: [458, 20, 0, 0.05, None, None, 'N/A', None, -4, 8, 'N/A', -2],
    ci.Element.Ni: [458, 20, 0, 0.05, None, None, 'N/A', None, -4, 8, 'N/A', -2],
    ci.Element.Si: [458, 20, 0, 0.05, None, None, 'N/A', None, -4, 8, 'N/A', -2]
}

def generate_pollution_fraction_values():
    return [1, 0, -4, -7]

def generate_pressure_values():
    #return [0, 30, 60]
    #return [0, 20, 40, 60]
    return [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60]

def generate_fcf_values():
    #return [0, 0.5, 1]
    #return [0, 0.25, 0.5, 0.75, 0.99]
    return [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99]

def get_element_key(element):
    return str(element) + '/Hx'

def import_MWDD_data(file_name):
    toret = dict()
    row_i = 0
    with open(file_name, encoding='utf-8') as mwddcsv:
        for row in csv.reader(mwddcsv):
            assert len(row) == 8
            if row_i == 0:
                pass # header
            else:
                # wdid,teff,logsihe,logsih,logcrhe,logcrh,lognihe,lognih
                toret[row[0]] = {'Teff': row[1], 'Si/He': row[2], 'Si/H': row[3], 'Cr/He': row[4], 'Cr/H': row[5], 'Ni/He': row[6], 'Ni/H': row[7]}
            row_i += 1
    print(str(len(toret)) + ' WDs imported')
    return toret

def generate_data(test_args_all):
    manager = mn.Manager(
        Namespace(
            wd_data_filename='WDInputData.csv',
            stellar_compositions_filename='StellarCompositionsSortFE.csv',
            n_live_points = 0, # This argument shouldn't matter, in fact we only use the manager to access observational data so nothing else matters
            enhancement_model = 'NonEarthlike',
            base_dir = '.',
            pollution_model_names=['Model_24']
        )
    )
    manager.publish_live_data(1)  # Using 1 as a dummy observation
    
    data_to_plot = dict()
    for element, test_args in test_args_all.items():
        data_to_plot[element] = dict()
        for arg_name, arg_set in test_args.items():
            data_to_plot[element][arg_name] = dict()
            for el in ci.usual_elements:
                data_to_plot[element][arg_name][get_element_key(el)] = list()
            for pol_frac in generate_pollution_fraction_values():
                test_result = cm.complete_model_calculation(
                    arg_set[0], # metallicity
                    arg_set[1], # t_since
                    arg_set[2], # distance
                    arg_set[3], # feeding zone
                    arg_set[4], # parent core frac
                    arg_set[5], # parent crust frac
                    arg_set[6], # fragment core frac
                    arg_set[7], # fragment crust frac
                    pol_frac, # pollution frac
                    10**(arg_set[9]), # accretion timescale
                    arg_set[10], # pressure
                    arg_set[11], # fO2
                    'NonEarthlike'
                )
                for el in ci.usual_elements:
                    data_to_plot[element][arg_name][get_element_key(el)].append(test_result[0][el])
    return data_to_plot

def generate_contour_data(test_args_all):
    manager = mn.Manager(
        Namespace(
            wd_data_filename='WDInputData.csv',
            stellar_compositions_filename='StellarCompositionsSortFE.csv',
            n_live_points = 0, # This argument shouldn't matter, in fact we only use the manager to access observational data so nothing else matters
            enhancement_model = 'NonEarthlike',
            base_dir = '.',
            pollution_model_names=['Model_24']
        )
    )
    manager.publish_live_data(1)  # Using 1 as a dummy observation
    
    data_to_plot = dict()
    pol_frac_sample = list()
    pol_frac_wds = dict()
    for element, arg_set in test_args_all.items():
        data_to_plot[element] = dict()
        for p in generate_pressure_values():
            for fcf in generate_fcf_values():
                test_result = cm.complete_model_calculation(
                    arg_set[0], # metallicity
                    arg_set[1], # t_since
                    arg_set[2], # distance
                    arg_set[3], # feeding zone
                    arg_set[4], # parent core frac
                    arg_set[5], # parent crust frac
                    fcf, # fragment core frac
                    arg_set[7], # fragment crust frac
                    arg_set[8], # pollution frac
                    10**(arg_set[9]), # accretion timescale
                    p, # pressure
                    arg_set[11], # fO2
                    'NonEarthlike'
                )
                data_to_plot[element][(p, fcf)] = test_result[0]
                total_p = 0
                for el, val in test_result[0].items():
                    total_p += 10**val
                total_p = np.log10(total_p)
                pol_frac_sample.append(total_p)
    for wd, abundances in manager.wd_abundances.items():
        total_p = 0
        for el, val in abundances.items():
            if val != 0:
                total_p += 10**val
        total_p = np.log10(total_p)
        pol_frac_wds[wd] = total_p
    #print('Pollution fractions in sample:')
    #print(pol_frac_sample)
    #print('Pollution fractions in wds:')
    #print(pol_frac_wds)
    return data_to_plot
    
def tabulate_contour_data(contour_data_to_plot, path=None):
    # structure: contour_data_to_plot[(P, fcf)] = {element: el/Hx}
    prev_args = None
    last_element = None
    for element, args in contour_args.items():
        print(element)
        print(prev_args)
        print(args)
        if prev_args is not None:
            if prev_args != args:
                print('Could not tabulate data because args were inconsistent between elements')
                return
        prev_args = args
        last_element = element
    path_to_use = 'observability_table.tex'
    Ps_to_tabulate = [0, 10, 20, 30, 40, 50, 60]
    fcfs_to_tabulate = [0, 0.2, 0.4, 0.6, 0.8, 0.99]
    pol_frac_hack = -1 # This is a hack. This number needs to be the same as the one in plot_contour_observability_data labelled THIS ONE
    if path is not None:
        path_to_use = path + path_to_use
    with open(path_to_use, 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            '%Fragment Core Fraction',
            'Pressure /GPa',
            'log(Mg/Fe)',
            'log(Ca/Fe)',
            'log(Cr/Hx)',
            'log(Ni/Hx)',
            'log(Si/Hx)'
        ])
        list_of_all_rows = list()
        for P_fcf_key, el_abundances in contour_data_to_plot[last_element].items(): # arbitrarily use the last element. We checked earlier that they're all the same
            list_of_vals = [
                str(P_fcf_key[1]),
                str(P_fcf_key[0]),
                str(round(el_abundances[ci.Element.Mg] - el_abundances[ci.Element.Fe], 2)),
                str(round(el_abundances[ci.Element.Ca] - el_abundances[ci.Element.Fe], 2)),
                str(round(el_abundances[ci.Element.Cr] + pol_frac_hack, 2)),
                str(round(el_abundances[ci.Element.Ni] + pol_frac_hack, 2)),
                str(round(el_abundances[ci.Element.Si] + pol_frac_hack, 2))
            ]
            if P_fcf_key[0] in Ps_to_tabulate and P_fcf_key[1] in fcfs_to_tabulate:
                list_of_all_rows.append([' & '.join(list_of_vals) + ' \\\\'])
        for row in sorted(list_of_all_rows):
            to_write.writerow(row)

def plot_observability_data(fcf_data_to_plot, pressure_data_to_plot):
    graph_fac = gf.GraphFactory()
    
    cr_fcf_plot = graph_fac.plot_observability(fcf_data_to_plot[ci.Element.Cr], ci.Element.Ca, ci.Element.Cr, ci.Element.Ca, ci.Element.Fe)
    cr_pressure_plot = graph_fac.plot_observability(pressure_data_to_plot[ci.Element.Cr], ci.Element.Ca, ci.Element.Cr, ci.Element.Ca, ci.Element.Fe)
    
    ni_fcf_plot = graph_fac.plot_observability(fcf_data_to_plot[ci.Element.Ni], ci.Element.Ca, ci.Element.Ni, ci.Element.Ca, ci.Element.Fe)
    ni_pressure_plot = graph_fac.plot_observability(pressure_data_to_plot[ci.Element.Ni], ci.Element.Ca, ci.Element.Ni, ci.Element.Ca, ci.Element.Fe)
    
    si_fcf_plot = graph_fac.plot_observability(fcf_data_to_plot[ci.Element.Si], ci.Element.Ca, ci.Element.Si, ci.Element.Ca, ci.Element.Fe)
    si_pressure_plot = graph_fac.plot_observability(pressure_data_to_plot[ci.Element.Si], ci.Element.Ca, ci.Element.Si, ci.Element.Ca, ci.Element.Fe)
    si_pressure_plot['observability_plot']['subplots']['subplot1']['x_hide_ticks'] = [0]
    
    graph_fac.multipanelise([
        cr_fcf_plot, cr_pressure_plot,
        ], 1, 2, ['observability_multipanel_cr_only.png'], 5, 10, 0, 0)
        
    graph_fac.multipanelise([
        ni_fcf_plot, ni_pressure_plot,
        ], 1, 2, ['observability_multipanel_ni_only.png'], 5, 10, 0, 0)
        
    graph_fac.multipanelise([
        si_fcf_plot, si_pressure_plot,
        ], 1, 2, ['observability_multipanel_si_only.png'], 5, 10, 0, 0)
    
    graph_fac.multipanelise([
        cr_fcf_plot, cr_pressure_plot,
        ni_fcf_plot, ni_pressure_plot,
        si_fcf_plot, si_pressure_plot
        ], 3, 2, ['observability_multipanel.pdf'], 15, 10, 0, 0)
        
def extract_el_el(contour_data_to_plot_el, numerator_el, denom_el, pol_frac=-6):
    fcf_index = 0
    toret_array = np.zeros((len(generate_fcf_values()), len(generate_pressure_values())))
    toret_hx_array = np.zeros((len(generate_fcf_values()), len(generate_pressure_values())))
    while fcf_index < len(generate_fcf_values()):
        p_index = 0
        while p_index < len(generate_pressure_values()):
            p = generate_pressure_values()[p_index]
            fcf = generate_fcf_values()[fcf_index]
            toret_array[fcf_index][p_index] = contour_data_to_plot_el[(p, fcf)][numerator_el] - contour_data_to_plot_el[(p, fcf)][denom_el]
            toret_hx_array[fcf_index][p_index] = contour_data_to_plot_el[(p, fcf)][numerator_el] + pol_frac
            p_index += 1
        fcf_index += 1
    return toret_array, toret_hx_array

def build_guideline(contour_data_to_plot_el, guideline_num_el, guideline_denom_el, path_to_follow):
    toret = dict()
    toret['num_el'] = guideline_num_el
    toret['denom_el'] = guideline_denom_el
    toret['x_data'] = list()
    toret['y_data'] = list()
    toret['z_data'] = list()
    for point in path_to_follow:
        # Assume these are (p, fcf) tuples (i.e. x, y)
        toret['x_data'].append(point[0])
        toret['y_data'].append(point[1])
        toret['z_data'].append("{:.2f}".format(contour_data_to_plot_el[point][guideline_num_el] - contour_data_to_plot_el[point][guideline_denom_el]))
    return toret
    
def plot_contour_observability_data(contour_data_to_plot): # Structure: contour_data_to_plot_el[el_of_that_run][(p, fcf)][element_of_interest]
    graph_fac = gf.GraphFactory()
    
    comparison_el = ci.Element.Ca
    guideline_num_el = ci.Element.Mg
    guideline_denom_el = ci.Element.Fe
    pol_frac = -1    # THIS ONE
    path_to_follow = [(10, 0), (10, 0.2), (10, 0.4), (10, 0.6), (10, 0.8), (10, 0.99), (30, 0.99), (50, 0.99)]

    processed_data = dict()
    processed_hx_data = dict()
    guideline_hx = dict()
    #data_to_tabulate = dict() #  [(P, fcf)] = {element: el/Hx}
    for element in contour_data_to_plot.keys():
        processed_data_tuple = extract_el_el(contour_data_to_plot[element], element, comparison_el, pol_frac)
        processed_data[element] = processed_data_tuple[0]
        processed_hx_data[element] = processed_data_tuple[1]
        guideline_hx[element] = build_guideline(contour_data_to_plot[element], guideline_num_el, guideline_denom_el, path_to_follow)
    all_plots = list()
    all_hx_plots = list()
    for el in contour_args.keys():
        all_plots.append(graph_fac.plot_3d_log_elel_ratio(
            processed_data[el],
            generate_pressure_values(),
            generate_fcf_values(),
            'Observability',
            el,
            comparison_el,
            None,
            None,
            None,
            'contour',
            'Pressure',
            'fcf',
            'Dummy3',
            True
        ))
        all_hx_plots.append(graph_fac.plot_3d_log_elel_ratio(
            processed_hx_data[el],
            generate_pressure_values(),
            generate_fcf_values(),
            'Observability',
            el,
            ci.Element.Placeholder,
            None,
            None,
            None,
            'contour',
            'Pressure',
            'fcf',
            'Dummy3',
            True,
            [guideline_hx[el]]
        ))
    graph_fac.multipanelise(all_plots, 3, 1, ['observability_contour_multipanel_ss.pdf'], 15, 10, 0, None)
    graph_fac.multipanelise(all_hx_plots, 3, 1, ['observability_contour_hx_multipanel_ss.pdf'], 15, 10, 0, None)

def plot_el_v_teff_data(el_v_teff_data):
    for element in [ci.Element.Ni, ci.Element.Cr, ci.Element.Si]:
        teff_vals = {'H': list(), 'He': list()}
        teff_upper_bound_vals = {'H': list(), 'He': list()}
        el_vals = {'H': list(), 'He': list()}
        el_upper_bound_vals = {'H': list(), 'He': list()}
        He_key = str(element) + '/He'
        H_key = str(element) + '/H'
        for WD, WD_data in el_v_teff_data.items():
            if WD_data[He_key] != '':
                upper_bound = False
                try:
                    val = float(WD_data[He_key])
                except ValueError:
                    val = float(WD_data[He_key][1:])
                    upper_bound = True
                try:
                    teff = float(WD_data['Teff'])
                except ValueError:
                    teff = None
                if teff is not None:
                    if upper_bound:
                        teff_upper_bound_vals['He'].append(teff)
                        el_upper_bound_vals['He'].append(val)
                    else:
                        teff_vals['He'].append(teff)
                        el_vals['He'].append(val)
            elif WD_data[H_key] != '':
                upper_bound = False
                try:
                    val = float(WD_data[H_key])
                except ValueError:
                    val = float(WD_data[H_key][1:])
                    upper_bound = True
                try:
                    teff = float(WD_data['Teff'])
                except ValueError:
                    teff = None
                if teff is not None:
                    if upper_bound:
                        teff_upper_bound_vals['H'].append(teff)
                        el_upper_bound_vals['H'].append(val)
                    else:
                        teff_vals['H'].append(teff)
                        el_vals['H'].append(val)
            else:
                pass
        graph_fac = gf.GraphFactory()
        print(teff_vals)
        print(teff_upper_bound_vals)
        print(el_vals)
        print(el_upper_bound_vals)
        graph_fac.plot_el_v_Teff(element, el_vals, teff_vals, el_upper_bound_vals, teff_upper_bound_vals)

def plot_original_observability(pressure_args, fcf_args):
    pressure_data_to_plot = generate_data(pressure_args)
    fcf_data_to_plot = generate_data(fcf_args)
    plot_observability_data(fcf_data_to_plot, pressure_data_to_plot)
    
def plot_contour_observability(contour_args):
    contour_data_to_plot = generate_contour_data(contour_args)
    tabulate_contour_data(contour_data_to_plot)
    plot_contour_observability_data(contour_data_to_plot)

def plot_el_v_teff(file_name):
    el_v_teff_data = import_MWDD_data(file_name)
    plot_el_v_teff_data(el_v_teff_data)
    
def main():
    #plot_original_observability(pressure_args, fcf_args)
    plot_contour_observability(contour_args)
    #plot_el_v_teff('../data/MWDD-export_CrNiSi_v_Teff.csv')

if __name__ == '__main__':
    main()
