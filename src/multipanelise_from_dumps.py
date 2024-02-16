#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import csv
#import numpy as np
#import os

import graph_factory as gf
import pwd_utils as pu

dump_config1 = {
    'y_dimension': 3,
    'x_dimension': 3,
    'filenames': ['hollands_ksp_multiplots.pdf', 'hollands_ksp_multiplots.png'],
    'fig_height': 10,
    'fig_width': 15,
    'gridspec_wspace': 0,
    'gridspec_hspace': 0,
    'sharey_axes': True,
    'sharex_axes': True,
    'list_of_dumps': [
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsTidal_Hollands_Ca_Fe_fragment_core_frac.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsTidal_Hollands_Mg_Fe_fragment_core_frac.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsTidal_Hollands_Ca_Mg_fragment_core_frac.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsCollisional_Hollands_Ca_Fe_fragment_core_frac.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsCollisional_Hollands_Mg_Fe_fragment_core_frac.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsCollisional_Hollands_Ca_Mg_fragment_core_frac.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsDeltaPop_Hollands_Ca_Fe_fragment_core_frac.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsDeltaPop_Hollands_Mg_Fe_fragment_core_frac.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/pipeline_ksp_comparison_SyntheticHollandsDeltaPop_Hollands_Ca_Mg_fragment_core_frac.pdf.txt'
    ]
}

dump_config2 = {
    'y_dimension': 1,
    'x_dimension': 3,
    'filenames': ['HollandsDelta0p1_CaFe_MgFe_CaMg.pdf', 'HollandsDelta0p1_CaFe_MgFe_CaMg.png'],
    'fig_height': 10,
    'fig_width': 15,
    'gridspec_wspace': 0,
    'gridspec_hspace': 0,
    'sharey_axes': True,
    'sharex_axes': False,
    'list_of_dumps': [
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/SyntheticHollandsDeltaPop_RealisticObservererr0p1_StandardModeller/fcf_Observed_pdf_Ca_Fe.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/SyntheticHollandsDeltaPop_RealisticObservererr0p1_StandardModeller/fcf_Observed_pdf_Mg_Fe.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'deltavarerr/SyntheticHollandsDeltaPop_RealisticObservererr0p1_StandardModeller/fcf_Observed_pdf_Ca_Mg.pdf.txt'
    ]
}

dump_config3 = {
    'y_dimension': 3,
    'x_dimension': 3,
    'filenames': ['ternaryplots.pdf', 'ternaryplots.png'],
    'fig_height': 10,
    'fig_width': 15,
    'gridspec_wspace': 0.2,
    'gridspec_hspace': 0.2,
    'sharey_axes': False,
    'sharex_axes': False,
    'list_of_dumps': [
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsDeltaPop_HollandsObservererr0_NullModeller/ternary_plot_sampled_SyntheticHollandsDeltaPop_HollandsObservererr0.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsDeltaPop_HollandsObservererr0p2_NullModeller/ternary_plot_sampled_SyntheticHollandsDeltaPop_HollandsObservererr0p2.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsDeltaPop_HollandsObservererr0p4_NullModeller/ternary_plot_sampled_SyntheticHollandsDeltaPop_HollandsObservererr0p4.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsCollisional_HollandsObservererr0_NullModeller/ternary_plot_sampled_SyntheticHollandsCollisional_HollandsObservererr0.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsCollisional_HollandsObservererr0p2_NullModeller/ternary_plot_sampled_SyntheticHollandsCollisional_HollandsObservererr0p2.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsCollisional_HollandsObservererr0p4_NullModeller/ternary_plot_sampled_SyntheticHollandsCollisional_HollandsObservererr0p4.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsTidal_HollandsObservererr0_NullModeller/ternary_plot_sampled_SyntheticHollandsTidal_HollandsObservererr0.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsTidal_HollandsObservererr0p2_NullModeller/ternary_plot_sampled_SyntheticHollandsTidal_HollandsObservererr0p2.pdf.txt',
        pu.get_path_to_pipeline_base_dir() + 'dzcomparison/SyntheticHollandsTidal_HollandsObservererr0p4_NullModeller/ternary_plot_sampled_SyntheticHollandsTidal_HollandsObservererr0p4.pdf.txt'
    ]
}

def multipanelise_from_dumps(list_of_dump_configs):
    graph_fac = gf.GraphFactory()
    for dump_config in list_of_dump_configs:
        graph_fac.multipanelise_from_dumps(
            dump_config['list_of_dumps'],
            dump_config['y_dimension'],
            dump_config['x_dimension'],
            dump_config['filenames'],
            dump_config['fig_height'],
            dump_config['fig_width'],
            dump_config['gridspec_wspace'],
            dump_config['gridspec_hspace'],
            dump_config['sharey_axes'],
            dump_config['sharex_axes']
        )

def main():
    list_of_dump_configs = [dump_config3]
    multipanelise_from_dumps(list_of_dump_configs)

if __name__ == '__main__':
    main()
