#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chemistry_info as ci
import geology_info as gi
import graph_factory as gf

import numpy as np

def convert_number_composition_to_wd_fit(composition_dict):
    geo_model_for_solar_abundances = gi.GeologyModel()
    toret = dict()
    for comp_name, comp in composition_dict.items():
        if comp[ci.Element.Mg] != 0:
            toret[comp_name] = list()
            for element, el_abundance in comp.items():
                if element != ci.Element.Mg:
                    XMg = np.log10(el_abundance/comp[ci.Element.Mg])
                    XMg_normalised = XMg - geo_model_for_solar_abundances.solar_abundances[element]
                    toret[comp_name].append(XMg_normalised)
    return toret

def main():
    graph_fac = gf.GraphFactory('graphs')
    system = 'Earth'
    file_prefix = 'layer'
    test_composition = {
        'Bulk': {
            ci.Element.Al: 0.0153,
            ci.Element.Ti: 0.00044,
            ci.Element.Ca: 0.011099,
            ci.Element.Ni: 0.008066,
            ci.Element.Fe: 0.149057,
            ci.Element.Cr: 0.002351,
            ci.Element.Mg: 0.16482,
            ci.Element.Si: 0.149118,
            ci.Element.Na: 0.002037,
            ci.Element.O: 0.482879,
            ci.Element.C: 0.001581,
            ci.Element.N: 0.000046429
        },
        'Core': {
            ci.Element.Al: 0, # meant to be 0
            ci.Element.Ti: 0, # meant to be 0
            ci.Element.Ca: 0, # meant to be 0
            ci.Element.Ni: 0.0444,
            ci.Element.Fe: 0.7676,
            ci.Element.Cr: 0.00868,
            ci.Element.Mg: 0, # meant to be 0
            ci.Element.Si: 0.1071,
            ci.Element.Na: 0.00001, # meant to be 0
            ci.Element.O: 0, # meant to be 0
            ci.Element.C: 0.0083482,
            ci.Element.N: 0.0002684
        },
        'Mantle': {
            ci.Element.Al: 0.018,
            ci.Element.Ti: (1200/1000000)/2, #division by 2 is an approximate correction for weight -> number
            ci.Element.Ca: 0.013,
            ci.Element.Ni: 0.001,
            ci.Element.Fe: 0.024,
            ci.Element.Cr: 0.001,
            ci.Element.Mg: 0.198,
            ci.Element.Si: 0.158,
            ci.Element.Na: 0.002,
            ci.Element.O: 0.581,
            ci.Element.C: 0,
            ci.Element.N: 0.001
        },
        'Crust': {
            ci.Element.Al: 0.0641,
            ci.Element.Ti: 0.0042,
            ci.Element.Ca: 0.04452,
            ci.Element.Ni: 0.0000371,
            ci.Element.Fe: 0.0314,
            ci.Element.Cr: 0.000139,
            ci.Element.Mg: 0.04167,
            ci.Element.Si: 0.181509,
            ci.Element.Na: 0.01771,
            ci.Element.O: 0.6011,
            ci.Element.C: 0,
            ci.Element.N: 0
        }
    }
    impure_core_frac = 0.9
    test_composition['Impure Core'] = dict()
    for el in ci.usual_elements:
        test_composition['Impure Core'][el] = (impure_core_frac*test_composition['Core'][el]) + ((1-impure_core_frac)*test_composition['Mantle'][el])
    all_wd_abundances = {
        ci.Element.Al: -6,
        ci.Element.Ti: -6,
        ci.Element.Ca: -6,
        ci.Element.Ni: -6,
        ci.Element.Fe: -6,
        ci.Element.Cr: -6,
        ci.Element.Mg: -6,
        ci.Element.Si: -6,
        ci.Element.Na: -6,
        ci.Element.O: -6,
        ci.Element.C: -6,
        ci.Element.N: -6
    }
    all_wd_abundance_errors = {
        ci.Element.Al: 0.1,
        ci.Element.Ti: 0.1,
        ci.Element.Ca: 0.1,
        ci.Element.Ni: 0.1,
        ci.Element.Fe: 0.1,
        ci.Element.Cr: 0.1,
        ci.Element.Mg: 0.1,
        ci.Element.Si: 0.1,
        ci.Element.Na: 0.1,
        ci.Element.O: 0.1,
        ci.Element.C: 0.1,
        ci.Element.N: 0.1
    }
    fit_dict = convert_number_composition_to_wd_fit(test_composition)
    error_low_dict = None
    error_high_dict = None
    upper_limits = None
    lower_limits = None
    video = False
    hack_legend_to_only_show_pressure = False
    excluded_wd_abundances = None
    excluded_wd_abundance_errors = None
    excluded_wd_upper_bounds = None
    excluded_wd_lower_bounds = None
    graph_fac.make_composition_plot_mk2(
        system,
        file_prefix,
        all_wd_abundances,
        all_wd_abundance_errors,
        fit_dict,
        error_low_dict,
        error_high_dict,
        upper_limits,
        lower_limits,
        video,
        hack_legend_to_only_show_pressure,
        excluded_wd_abundances,
        excluded_wd_abundance_errors,
        excluded_wd_upper_bounds,
        excluded_wd_lower_bounds
    )

if __name__ == '__main__':
    main()
