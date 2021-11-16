import matplotlib.pyplot as plt
import numpy as np
import shutil, os
import sys

import model_parameters as mp
import pwd_utils as pu

sys.path.append(pu.get_path_to_utils())

import dict_plotter as dp
    
def load_points_files(file_prefix, output_dir='output'):
    phys_live_file = output_dir + '/' + file_prefix + '_phys_live.points'  # Active points. N + 2 columns: <params>, log(Z), node number
    ev_file = output_dir + '/' + file_prefix + '_ev.dat'  # Rejected points. N + 3 columns: <params>, log(Z), log(prior mass), node number
    active_points = np.loadtxt(phys_live_file, ndmin=2)
    rejected_points = np.loadtxt(ev_file, ndmin=2)
    model_name = file_prefix.split('_lp')[0]  # Assume _lp doesn't appear in the model name...
    live_points = int(file_prefix.split('_lp')[1].split('_')[0])
    param_list = list()
    if model_name.startswith('Hierarchy'):
        hierarchy_name = model_name.split('_levels')[0] # Assume _levels or levels_ doesn't appear in the hierarchy name...
        levels = model_name.split('levels_')[1]
        hierarchy_definition = mp.hierarchy_definitions_dict[hierarchy_name]
        for level in levels:
            param_list += hierarchy_definition[int(level)]
    else:
        model_definition = mp.model_definitions_dict[model_name]
        for param, present in model_definition.items():
            if present:
                param_list.append(param)
    sorted_param_list = list()
    for param in mp.ModelParameter:
        if param in param_list:
            sorted_param_list.append(param)
    assert len(active_points[0]) == len(param_list) + 2
    assert len(rejected_points[0]) == len(param_list) + 3
    return active_points, rejected_points, sorted_param_list, live_points

def plot_ellipses(file_prefix, active_points, rejected_points, sorted_param_list, live_points):
    param_index_1 = 0
    param_index_2 = 0
    while param_index_2 <= param_index_1:
        if param_index_1 == param_index_2:
            param_index_1 += 1
            param_index_2 = 0
        if param_index_1 == len(sorted_param_list):
            return
        print('Plotting ellipses for indices ' + str(param_index_1) + ' and ' + str(param_index_2))
        plot_dict = {
            'ellipses': {
                'show': False,
                'filenames': [file_prefix + '_ellipses' + str(param_index_1) + str(param_index_2) + '.pdf'],
                'fig_height': 15,
                'fig_width': 10,
                'subplots': {
                    'Active': {
                        'subplot_region': 311,
                        'legend': False,
                        'legend_loc': 'best',
                        'legend_text_size': 10,
                        'title_text': r'Active ellipse ' + file_prefix,
                        'title_fontsize': 10,
                        'xlabel_text': str(sorted_param_list[param_index_1]),
                        'xlabel_fontsize': 8,
                        'x_tick_fontsize': 8,
                        'ylabel_text': str(sorted_param_list[param_index_2]),
                        'ylabel_fontsize': 8,
                        'y_tick_fontsize': 8,
                        'font': 'STIXGeneral',
                        'series': {
                            'active': {
                                'type': dp.PlotType.scatter_3d,
                                'x_data': active_points[:,param_index_1],
                                'y_data': active_points[:,param_index_2],
                                'z_data': active_points[:,len(sorted_param_list)],
                                'fill': False,
                                'legend': False,
                                'cbar_label': 'log(Z)'
                            }
                        }
                    },
                    'LastRejected': {
                        'subplot_region': 312,
                        'legend': False,
                        'legend_loc': 'best',
                        'legend_text_size': 10,
                        'title_text': r'Last rejected points ' + file_prefix,
                        'title_fontsize': 10,
                        'xlabel_text': str(sorted_param_list[param_index_1]),
                        'xlabel_fontsize': 8,
                        'x_tick_fontsize': 8,
                        'ylabel_text': str(sorted_param_list[param_index_2]),
                        'ylabel_fontsize': 8,
                        'y_tick_fontsize': 8,
                        'font': 'STIXGeneral',
                        'series': {
                            'last_rejected': {
                                'type': dp.PlotType.scatter_3d,
                                'x_data': rejected_points[-live_points:,param_index_1],
                                'y_data': rejected_points[-live_points:,param_index_2],
                                'z_data': rejected_points[-live_points:,len(sorted_param_list)],
                                'fill': False,
                                'legend': False,
                                'cbar_label': 'log(Z)'
                            }
                        }
                    },
                    'FirstRejected': {
                        'subplot_region': 313,
                        'legend': False,
                        'legend_loc': 'best',
                        'legend_text_size': 10,
                        'title_text': r'First rejected points ' + file_prefix,
                        'title_fontsize': 10,
                        'xlabel_text': str(sorted_param_list[param_index_1]),
                        'xlabel_fontsize': 8,
                        'x_tick_fontsize': 8,
                        'ylabel_text': str(sorted_param_list[param_index_2]),
                        'ylabel_fontsize': 8,
                        'y_tick_fontsize': 8,
                        'font': 'STIXGeneral',
                        'series': {
                            'first_rejected': {
                                'type': dp.PlotType.scatter_3d,
                                'x_data': rejected_points[:live_points,param_index_1],
                                'y_data': rejected_points[:live_points,param_index_2],
                                'z_data': rejected_points[:live_points,len(sorted_param_list)],
                                'fill': False,
                                'legend': False,
                                'cbar_label': 'log(Z)'
                            }
                        }
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output()
        param_index_2 += 1

def plot_final_dists(file_prefix, active_points, sorted_param_list):
    plot_dict = {
        'final_dists': {
            'show': False,
            'filenames': [file_prefix + '_final_dists.pdf'],
            'fig_height': 15,
            'fig_width': 10,
            'subplots': dict()
        }
    }
    
    for p, param in enumerate(sorted_param_list):
        plot_dict['final_dists']['subplots'][str(param)] = {
            'subplot_region': (len(sorted_param_list), 1, p + 1),
            'legend': False,
            'legend_loc': 'best',
            'legend_text_size': 10,
            'title_text': r'Final distribution of active points for run ' + file_prefix,
            'title_fontsize': 10,
            'xlabel_text': str(param),
            'xlabel_fontsize': 8,
            'x_tick_fontsize': 8,
            'ylabel_text': 'log(Z)',
            'ylabel_fontsize': 8,
            'y_tick_fontsize': 8,
            'font': 'STIXGeneral',
            'series': {
                'active': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': active_points[:,p],
                    'y_data': active_points[:,len(sorted_param_list)],
                    'line_color': 'b',
                    'line_marker': 'x',
                    'line_markersize': 1,
                    'line_style': 'None',
                    'zorder': 1,
                    'legend': True
                }
            }
        }
        
    plotter = dp.DictPlotter(plot_dict)
    plotter.draw()
    plotter.yield_output()
        
def plot_system(file_prefix):
    active_points, rejected_points, sorted_param_list, live_points = load_points_files(file_prefix)
    plot_final_dists(file_prefix, active_points, sorted_param_list)
    plot_ellipses(file_prefix, active_points, rejected_points, sorted_param_list, live_points)

def main():
    #progress_plotter_nonlive('Hierarchy_Default_levels_0_lp1000_obs95_phys_live.points', 4)
    plot_system('Hierarchy_Default_levels_0_lp1000_obs95')
    plot_system('Model_4_equiv_lp20_obs0')
    #load_points_files('Hierarchy_Default_levels_03_lp1000_obs160')
    
if __name__ == '__main__':
    main()
