#!/usr/bin/env python
# -*- coding: utf-8 -*-

from argparse import Namespace

import chemistry_info as ci
import detection_thresholds as dt
import hollands_abundances as ha
import mwdd_abundances as mwdd
import graph_factory as gf
import pwd_utils as pu

def build_wd_data(element, DA=True):
    if DA:
        wd_data = {'DA': list(), 'DB': None}
        wd_upper_bounds = {'DA': list(), 'DB': None}
        wd_data['DA'], wd_upper_bounds['DA'] = mwdd.get_mwdd_abundances(element, ci.Element.H, True, True)
    else:
        wd_data = {'DB': list(), 'DA': None}
        wd_upper_bounds = {'DB': list(), 'DA': None}
        wd_data['DB'], wd_upper_bounds['DB'] = ha.get_hollands_abundances_for_dt_plot(element)
    return wd_data, wd_upper_bounds

def plot_detection_thresholds(DA=True):
    # Construct dict by element containing gradients/intercepts
    to_plot = dict()
    if DA:
        detection_threshold_whitelist = ['Default']
    else:
        detection_threshold_whitelist = ['Hollands']
    #manager = mn.Manager(  #This is just to load the stellar compositions
    #    Namespace(
    #        wd_data_filename='WDMSextract.csv',
    #        stellar_compositions_filename='StellarCompositionsSortFE.csv',
    #        n_live_points = 0,
    #        pollution_model_names=['Model_24'],
    #        enhancement_model='Earthlike',
    #        base_dir=pu.get_path_to_data()
    #    )
    #)
    for threshold_set_name, threshold_set in dt.threshold_bank.items():
        if threshold_set_name in detection_threshold_whitelist:
            for spectral_type, threshold_dict in threshold_set.items():
                for element, definition in threshold_dict.items():
                    if element not in to_plot:
                        to_plot[element] = dict()
                    if spectral_type not in to_plot[element]:
                        to_plot[element][spectral_type] = dict()
                    to_plot[element][spectral_type][threshold_set_name] = definition
    graph_fac = gf.GraphFactory()
    list_of_plots = list()
    if DA:
        min_teff = 4000
        max_teff = 30000
        elements_to_plot = [ci.Element.Ca, ci.Element.Fe, ci.Element.Mg, ci.Element.Ni, ci.Element.Ti, ci.Element.Cr, ci.Element.Al, ci.Element.Na, ci.Element.O]
    else:
        min_teff = 4000
        max_teff = 10000
        elements_to_plot = [ci.Element.Ca, ci.Element.Fe, ci.Element.Mg, ci.Element.Cr]
    for i, element in enumerate(elements_to_plot):
        wd_data, wd_upper_bounds = build_wd_data(element, DA)
        new_plot = graph_fac.make_detection_thresholds_plot(to_plot, element, wd_data, wd_upper_bounds, min_teff, max_teff, i == len(elements_to_plot) - 1)
        list_of_plots.append(new_plot)
    if DA:
        graph_fac.multipanelise(
            list_of_plots,
            5,
            2,
            ['detection_thresholds.pdf', 'detection_thresholds.png'],
            13,
            10,
            0.3,
            0,
            False,
            True
        )
    else:
        graph_fac.multipanelise(
            list_of_plots,
            2,
            2,
            ['detection_thresholds_dz.pdf', 'detection_thresholds_dz.png'],
            7,
            10,
            0.3,
            0,
            False,
            True
        )

def main():
    #plot_detection_thresholds(True)
    plot_detection_thresholds(False)

if __name__ == '__main__':
    main()
