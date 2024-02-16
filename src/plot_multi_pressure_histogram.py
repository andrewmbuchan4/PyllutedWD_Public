#!/usr/bin/env python
# -*- coding: utf-8 -*-

import model_analyser as ma
import graph_factory as gf
import pwd_utils as pu

def main():
    analyser = ma.ModelAnalyser(pu.get_path_to_output_graphs_dir())
    graph_fac = gf.GraphFactory()
    
    text_dict_dict = {
        'HPM0': {
            'hpm_fcf0_text': {
                'x_pos': 30, # Change to 0 for plotting separate dists...,
                'text_string': 'High pressure, 0\% core',
                'horizontalalignment': 'center',
                'fontsize': 24
            },
            'title_text': 'Pure material'
        },
        'HPM2': {
            'hpm_fcf2_text': {
                'x_pos': 30, # Change to 0 for plotting separate dists...,
                'text_string': 'High pressure, 2\% core',
                'horizontalalignment': 'center',
                'fontsize': 24
            },
            'title_text': 'Effect of impurity'
        },
        'LPM0': {
            'lpm_fcf0_text': {
                'x_pos': 30, # Change to 0 for plotting separate dists...,
                'text_string': 'Low pressure, 0\% core',
                'horizontalalignment': 'center',
                'fontsize': 24
            }
        },
        'LPM2': {
            'lpm_fcf2_text': {
                'x_pos': 30, # Change to 0 for plotting separate dists...,
                'text_string': 'Low pressure, 2\% core',
                'horizontalalignment': 'center',
                'fontsize': 24
            }
        },
        'HPC75': {
            'hpc_fcf75_text': {
                'x_pos': 30, # Change to 0 for plotting separate dists...,
                'text_string': 'High pressure, 75\% core',
                'horizontalalignment': 'center',
                'fontsize': 24
            }
        },
        'HPC99': {
            'hpc_fcf99_text': {
                'x_pos': 30, # Change to 0 for plotting separate dists...,
                'text_string': 'High pressure, 99\% core',
                'horizontalalignment': 'center',
                'fontsize': 24
            }
        },
        'LPC75': {
            'lpc_fcf75_text': {
                'x_pos': 30, # Change to 0 for plotting separate dists...,
                'text_string': 'Low pressure, 75\% core',
                'horizontalalignment': 'center',
                'fontsize': 24
            }
        },
        'LPC99': {
            'lpc_fcf99_text': {
                'x_pos': 30, # Change to 0 for plotting separate dists...,
                'text_string': 'Low pressure, 99\% core',
                'horizontalalignment': 'center',
                'fontsize': 24
            }
        }
    }
    
    #HPM
    hpm_fcf0_plot_dict = analyser.make_multisystem_pressure_plot([381, 383, 385], ['0.05 dex error', '0.1 dex error', '0.2 dex error'], text_dict_dict['HPM0'])  # fcf = 0
    hpm_fcf2_plot_dict = analyser.make_multisystem_pressure_plot([382, 384], ['0.05 dex error', '0.1 dex error'], text_dict_dict['HPM2'])# fcf = 0.02
    graph_fac.multipanelise([hpm_fcf0_plot_dict, hpm_fcf2_plot_dict], 1, 2, ['IdealHPM_multipressure_dist.pdf', 'IdealHPM_multipressure_dist.png'], 7, 20, None, 0, False, True)
    #LPM
    lpm_fcf0_plot_dict = analyser.make_multisystem_pressure_plot([389, 391, 393, 395], ['0.05 dex error', '0.1 dex error', '0.2 dex error', '0.4 dex error'], text_dict_dict['LPM0'])# fcf = 0
    lpm_fcf2_plot_dict = analyser.make_multisystem_pressure_plot([390, 392, 394], ['0.05 dex error', '0.1 dex error', '0.2 dex error'], text_dict_dict['LPM2'])# fcf = 0.02
    graph_fac.multipanelise([lpm_fcf0_plot_dict, lpm_fcf2_plot_dict], 1, 2, ['IdealLPM_multipressure_dist.pdf', 'IdealLPM_multipressure_dist.png'], 7, 20, None, 0, False, True)
    #HPC
    hpc_fcf75_plot_dict = analyser.make_multisystem_pressure_plot([397, 399, 401, 403], ['0.05 dex error', '0.1 dex error', '0.2 dex error', '0.4 dex error'], text_dict_dict['HPC75'])  # fcf = 0.75
    hpc_fcf99_plot_dict = analyser.make_multisystem_pressure_plot([398, 400, 402, 404], ['0.05 dex error', '0.1 dex error', '0.2 dex error', '0.4 dex error'], text_dict_dict['HPC99'])# fcf = 0.99
    graph_fac.multipanelise([hpc_fcf75_plot_dict, hpc_fcf99_plot_dict], 1, 2, ['IdealHPC_multipressure_dist.pdf', 'IdealHPC_multipressure_dist.png'], 7, 20, None, 0, False, True)
    #LPC
    lpc_fcf75_plot_dict = analyser.make_multisystem_pressure_plot([405, 407, 409, 411], ['0.05 dex error', '0.1 dex error', '0.2 dex error', '0.4 dex error'], text_dict_dict['LPC75'])  # fcf = 0.75
    lpc_fcf99_plot_dict = analyser.make_multisystem_pressure_plot([406, 408, 410, 412], ['0.05 dex error', '0.1 dex error', '0.2 dex error', '0.4 dex error'], text_dict_dict['LPC99'])# fcf = 0.99
    graph_fac.multipanelise([lpc_fcf75_plot_dict, lpc_fcf99_plot_dict], 1, 2, ['IdealLPC_multipressure_dist.pdf', 'IdealLPC_multipressure_dist.png'], 7, 20, None, 0, False, True)
    
    all_plots = [
        hpm_fcf0_plot_dict,
        hpm_fcf2_plot_dict,
        
        lpm_fcf0_plot_dict,
        lpm_fcf2_plot_dict,
        
        hpc_fcf99_plot_dict,
        hpc_fcf75_plot_dict,

        lpc_fcf99_plot_dict,
        lpc_fcf75_plot_dict,
        ]
    
    plots_to_remove_legend = [
        hpm_fcf0_plot_dict,
        hpm_fcf2_plot_dict,
        
        #lpm_fcf0_plot_dict,
        lpm_fcf2_plot_dict,
        
        hpc_fcf99_plot_dict,
        hpc_fcf75_plot_dict,

        lpc_fcf99_plot_dict,
        lpc_fcf75_plot_dict,
    ]
    
    for plottrl in all_plots:
        if plottrl in plots_to_remove_legend:
            plottrl['hist_plot']['subplots']['subplot1']['legend'] = False
        else:
            plottrl['hist_plot']['subplots']['subplot1']['legend_loc'] = 'lower right'
    
    graph_fac.multipanelise(all_plots, 4, 2, ['Ideal_allsystem_multipressure_dist.pdf'], 18, 15, None, 0, False, True)

if __name__ == '__main__':
    main()
