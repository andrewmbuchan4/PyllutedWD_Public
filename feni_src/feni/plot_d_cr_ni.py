#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import kd as kd
import constants as constants
import numpy as np

# Assume that files will not be moved around!
def get_path_to_parent():
    # Path(__file__) is the path of this file
    # resolve() returns the absolute path
    # parents[1] is the path of the parent directory
    return str(Path(__file__).resolve().parents[2])

def get_path_to_utils():
    return get_path_to_parent() + '/utils/'

sys.path.append(get_path_to_utils())

import dict_plotter as dp

def D_to_sf(D, D_earth):
    return (D * (1 + D_earth)) / (D_earth * (1 + D))

def generate_data():
    D = kd.D('data/part_param_fischer_update.dat', 'data/int_param.dat', 'data/composition.dat')
    
    Pmin = 0.1
    Pmax = 100
    Pinc = 0.1
    
    Tmin = 1400
    Tmax = 4400
    
    dIW = -2 # As in Fischer et al 2015
    
    ele_set = ['hf','u','ta','pb','nb','si','mn','zn','ga','v','cr','cu','ti','fe','w','p','co','ni','o']
    
    T_dict = {
        'run1': lambda P: 1400,
        'run2': lambda P: 2900,
        'run3': lambda P: 4400,
        'run4': lambda P: constants.Tpdliq(P),
    }
    
    p_data = dict()
    cr_data = dict()
    cr_sd_up_data = dict()
    cr_sd_down_data = dict()
    ni_data = dict()
    ni_sd_up_data = dict()
    ni_sd_down_data = dict()
    
    sf_cr_data = dict()
    sf_cr_sd_up_data = dict()
    sf_cr_sd_down_data = dict()
    sf_ni_data = dict()
    sf_ni_sd_up_data = dict()
    sf_ni_sd_down_data = dict()
    
    p_range = np.linspace(Pmin, Pmax, 30)
    t_range = np.linspace(Tmin, Tmax, 30)
    p_mesh, t_mesh = np.meshgrid(p_range, t_range)
    d_cr_array = np.zeros((len(t_range), len(p_range)))
    d_ni_array = np.zeros((len(t_range), len(p_range)))
    
    fo2_range = np.linspace(-7, 2, 30)
    d_cr_fo2_array = list()
    d_ni_fo2_array = list()
    d_cr_fo2_sd_up_array = list()
    d_cr_fo2_sd_down_array = list()
    d_ni_fo2_sd_up_array = list()
    d_ni_fo2_sd_down_array = list()
    
    earth_pressure = 54 # As in Fischer et al 2015
    ds_earth, ds_sd_earth = D.mkd(earth_pressure, constants.Tpdliq(earth_pressure), dIW, constants.nbot, constants.gammaFe_sil, ele_set)
    d_ni_earth = ds_earth['ni']
    d_cr_earth = ds_earth['cr']
    
    for run_name, T_func in T_dict.items():
        p_data[run_name] = list()
        cr_data[run_name] = list()
        cr_sd_up_data[run_name] = list()
        cr_sd_down_data[run_name] = list()
        ni_data[run_name] = list()
        ni_sd_up_data[run_name] = list()
        ni_sd_down_data[run_name] = list()
        sf_cr_data[run_name] = list()
        sf_cr_sd_up_data[run_name] = list()
        sf_cr_sd_down_data[run_name] = list()
        sf_ni_data[run_name] = list()
        sf_ni_sd_up_data[run_name] = list()
        sf_ni_sd_down_data[run_name] = list()
        
        for P in p_range:
            T = T_func(P)
            Ds, Ds_sd = D.mkd(P, T, dIW, constants.nbot, constants.gammaFe_sil, ele_set)
            p_data[run_name].append(P)
            cr_data[run_name].append(Ds['cr'])
            cr_sd_up_data[run_name].append(Ds_sd['cr'][0])
            cr_sd_down_data[run_name].append(Ds_sd['cr'][1])
            ni_data[run_name].append(Ds['ni'])
            ni_sd_up_data[run_name].append(Ds_sd['ni'][0])
            ni_sd_down_data[run_name].append(Ds_sd['ni'][1])
            
            sf_cr_data[run_name].append(D_to_sf(Ds['cr'], d_cr_earth))
            sf_cr_sd_up_data[run_name].append(D_to_sf(Ds_sd['cr'][0], d_cr_earth))
            sf_cr_sd_down_data[run_name].append(D_to_sf(Ds_sd['cr'][1], d_cr_earth))
            sf_ni_data[run_name].append(D_to_sf(Ds['ni'], d_ni_earth))
            sf_ni_sd_up_data[run_name].append(D_to_sf(Ds_sd['ni'][0], d_ni_earth))
            sf_ni_sd_down_data[run_name].append(D_to_sf(Ds_sd['ni'][1], d_ni_earth))
    
    p_index = 0
    for P in p_range:
        t_index = 0
        for T in t_range:
            Ds, Ds_sd = D.mkd(P, T, dIW, constants.nbot, constants.gammaFe_sil, ele_set)
            d_cr_array[t_index][p_index] = Ds['cr']
            d_ni_array[t_index][p_index] = Ds['ni']
            t_index += 1
        p_index += 1
    
    P_fo2 = 60
    T_fo2 = constants.Tpdliq(P_fo2)
    for fo2 in fo2_range:
        Ds, Ds_sd = D.mkd(P_fo2, T_fo2, fo2, constants.nbot, constants.gammaFe_sil, ele_set)
        d_cr_fo2_array.append(Ds['cr'])
        d_cr_fo2_sd_up_array.append(Ds_sd['cr'][0])
        d_cr_fo2_sd_down_array.append(Ds_sd['cr'][1])
        d_ni_fo2_array.append(Ds['ni'])
        d_ni_fo2_sd_up_array.append(Ds_sd['ni'][0])
        d_ni_fo2_sd_down_array.append(Ds_sd['ni'][1])
    
    print(d_cr_fo2_array)
    print(d_cr_fo2_sd_up_array)
    print(d_cr_fo2_sd_down_array)
    print(d_ni_fo2_array)
    print(d_ni_fo2_sd_up_array)
    print(d_ni_fo2_sd_down_array)
    return (
        p_range,
        t_range,
        fo2_range,
        p_data,
        cr_data,
        cr_sd_up_data,
        cr_sd_down_data,
        ni_data,
        ni_sd_up_data,
        ni_sd_down_data,
        d_cr_array,
        d_ni_array,
        d_cr_fo2_array,
        d_cr_fo2_sd_up_array,
        d_cr_fo2_sd_down_array,
        d_ni_fo2_array,
        d_ni_fo2_sd_up_array,
        d_ni_fo2_sd_down_array,
        sf_cr_data,
        sf_cr_sd_up_data,
        sf_cr_sd_down_data,
        sf_ni_data,
        sf_ni_sd_up_data,
        sf_ni_sd_down_data,
    )

def main():
    (
        p_range,
        t_range,
        fo2_range,
        p_data,
        cr_data,
        cr_sd_up_data,
        cr_sd_down_data,
        ni_data,
        ni_sd_up_data,
        ni_sd_down_data,
        d_cr_array,
        d_ni_array,
        d_cr_fo2_array,
        d_cr_fo2_sd_up_array,
        d_cr_fo2_sd_down_array,
        d_ni_fo2_array,
        d_ni_fo2_sd_up_array,
        d_ni_fo2_sd_down_array,
        sf_cr_data,
        sf_cr_sd_up_data,
        sf_cr_sd_down_data,
        sf_ni_data,
        sf_ni_sd_up_data,
        sf_ni_sd_down_data
    ) = generate_data()
    
    plot_dict = {
        'plot1': {
            'show': False,
            'filenames': ['d_cr.pdf'],
            'subplots': {
                'subplot1': {
                    'subplot_region': 111,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'title_text': r'Variation of $D_{Cr}$ with Pressure',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$D_{Cr}$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 100,
                    'y_min': 0,
                    'y_max': 40,
                    'font': 'STIXGeneral',
                    'series': {
                        #'T = 1400K': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run1'],
                        #    'y_data': cr_data['run1'],
                        #    'line_type': 'r-',
                        #    'line_linewidth': 1,
                        #    'legend': True,
                        #},
                        #'T = 1400K (+1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run1'],
                        #    'y_data': cr_sd_up_data['run1'],
                        #    'line_type': 'r--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = 1400K (+1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run1'],
                        #    'y_data': cr_data['run1'],
                        #    'y_shade_data': cr_sd_up_data['run1'],
                        #    'shade_colour': 'r',
                        #    'shade_alpha': 0.3,
                        #},
                        #'T = 1400K (-1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run1'],
                        #    'y_data': cr_sd_down_data['run1'],
                        #    'line_type': 'r--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = 1400K (-1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run1'],
                        #    'y_data': cr_data['run1'],
                        #    'y_shade_data': cr_sd_down_data['run1'],
                        #    'shade_colour': 'r',
                        #    'shade_alpha': 0.3,
                        #},
                        #'T = 2900K': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run2'],
                        #    'y_data': cr_data['run2'],
                        #    'line_type': 'g-',
                        #    'line_linewidth': 1,
                        #    'legend': True,
                        #},
                        #'T = 2900K (+1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run2'],
                        #    'y_data': cr_sd_up_data['run2'],
                        #    'line_type': 'g--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = 2900K (+1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run2'],
                        #    'y_data': cr_data['run2'],
                        #    'y_shade_data': cr_sd_up_data['run2'],
                        #    'shade_colour': 'g',
                        #    'shade_alpha': 0.3,
                        #},
                        #'T = 2900K (-1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run2'],
                        #    'y_data': cr_sd_down_data['run2'],
                        #    'line_type': 'g--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = 2900K (-1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run2'],
                        #    'y_data': cr_data['run2'],
                        #    'y_shade_data': cr_sd_down_data['run2'],
                        #    'shade_colour': 'g',
                        #    'shade_alpha': 0.3,
                        #},
                        #'T = 4400K': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run3'],
                        #    'y_data': cr_data['run3'],
                        #    'line_type': 'b-',
                        #    'line_linewidth': 1,
                        #    'legend': True,
                        #},
                        #'T = 4400K (+1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run3'],
                        #    'y_data': cr_sd_up_data['run3'],
                        #    'line_type': 'b--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = 4400K (+1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run3'],
                        #    'y_data': cr_data['run3'],
                        #    'y_shade_data': cr_sd_up_data['run3'],
                        #    'shade_colour': 'b',
                        #    'shade_alpha': 0.3,
                        #},
                        #'T = 4400K (-1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run3'],
                        #    'y_data': cr_sd_down_data['run3'],
                        #    'line_type': 'b--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = 4400K (-1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run3'],
                        #    'y_data': cr_data['run3'],
                        #    'y_shade_data': cr_sd_down_data['run3'],
                        #    'shade_colour': 'b',
                        #    'shade_alpha': 0.3,
                        #},
                        'T = T(Peridotite liquidus)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': cr_data['run4'],
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'T = T(Peridotite liquidus) (+1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': cr_sd_up_data['run4'],
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'T = T(Peridotite liquidus) (+1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': p_data['run4'],
                            'y_data': cr_data['run4'],
                            'y_shade_data': cr_sd_up_data['run4'],
                            'shade_colour': 'r',
                            'shade_alpha': 0.3,
                        },
                        'T = T(Peridotite liquidus) (-1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': cr_sd_down_data['run4'],
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'T = T(Peridotite liquidus) (-1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': p_data['run4'],
                            'y_data': cr_data['run4'],
                            'y_shade_data': cr_sd_down_data['run4'],
                            'shade_colour': 'r',
                            'shade_alpha': 0.3,
                        },
                        'P_c marker': {
                            'type': dp.PlotType.vline,
                            'x_start': (1825. - 1420.)/104.42,
                            'y_min': 0,
                            'y_max': 110,
                            'line_linewidth': 0.5,
                        },
                        'P_c text': {
                            'type': dp.PlotType.text,
                            'x_pos': 3.4,
                            'y_pos': -5.2,
                            'text_string': '$P_C$',
                        }
                    }
                }
            }
        },
        'plot2': {
            'show': True,
            'filenames': ['d_cr_3d.pdf'],
            'subplots': {
                'subplot1': {
                    'subplot_region': 111,
                    'title_text': r'Variation of $D_{Cr}$ with Pressure and Temperature',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 100,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'font': 'STIXGeneral',
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            'z_data': d_cr_array,
                            'cbar_label': r'$D_{Cr}$'
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': constants.Tpdliq,
                            'x_start': 0,
                            'x_end': 100,
                            'x_points': 1000,
                            'line_type': 'k--',
                            'line_markersize': 1,
                            'line_linewidth': 2,
                            'legend': True
                        }
                    }
                }
            }
        },
        'plot3': {
            'show': False,
            'filenames': ['d_ni.pdf'],
            'subplots': {
                'subplot1': {
                    'subplot_region': 111,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'title_text': r'Variation of $D_{Ni}$ with Pressure',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$D_{Ni}$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 100,
                    'y_min': 10,
                    'y_max': 100000,
                    'y_scale': 'log',
                    'font': 'STIXGeneral',
                    'series': {
                    #    'T = 1400K': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run1'],
                    #        'y_data': ni_data['run1'],
                    #        'line_type': 'r-',
                    #        'line_linewidth': 1,
                    #        'legend': True,
                    #    },
                    #    'T = 1400K (+1 sigma)': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run1'],
                    #        'y_data': ni_sd_up_data['run1'],
                    #        'line_type': 'r--',
                    #        'line_linewidth': 0.5
                    #    },
                    #    'T = 1400K (+1 sigma) shade': {
                    #        'type': dp.PlotType.shade,
                    #        'x_data': p_data['run1'],
                    #        'y_data': ni_data['run1'],
                    #        'y_shade_data': ni_sd_up_data['run1'],
                    #        'shade_colour': 'r',
                    #        'shade_alpha': 0.3,
                    #    },
                    #    'T = 1400K (-1 sigma)': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run1'],
                    #        'y_data': ni_sd_down_data['run1'],
                    #        'line_type': 'r--',
                    #        'line_linewidth': 0.5
                    #    },
                    #    'T = 1400K (-1 sigma) shade': {
                    #        'type': dp.PlotType.shade,
                    #        'x_data': p_data['run1'],
                    #        'y_data': ni_data['run1'],
                    #        'y_shade_data': ni_sd_down_data['run1'],
                    #        'shade_colour': 'r',
                    #        'shade_alpha': 0.3,
                    #    },
                    #    'T = 2900K': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run2'],
                    #        'y_data': ni_data['run2'],
                    #        'line_type': 'g-',
                    #        'line_linewidth': 1,
                    #        'legend': True,
                    #    },
                    #    'T = 2900K (+1 sigma)': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run2'],
                    #        'y_data': ni_sd_up_data['run2'],
                    #        'line_type': 'g--',
                    #        'line_linewidth': 0.5
                    #    },
                    #    'T = 2900K (+1 sigma) shade': {
                    #        'type': dp.PlotType.shade,
                    #        'x_data': p_data['run2'],
                    #        'y_data': ni_data['run2'],
                    #        'y_shade_data': ni_sd_up_data['run2'],
                    #        'shade_colour': 'g',
                    #        'shade_alpha': 0.3,
                    #    },
                    #    'T = 2900K (-1 sigma)': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run2'],
                    #        'y_data': ni_sd_down_data['run2'],
                    #        'line_type': 'g--',
                    #        'line_linewidth': 0.5
                    #    },
                    #    'T = 2900K (-1 sigma) shade': {
                    #        'type': dp.PlotType.shade,
                    #        'x_data': p_data['run2'],
                    #        'y_data': ni_data['run2'],
                    #        'y_shade_data': ni_sd_down_data['run2'],
                    #        'shade_colour': 'g',
                    #        'shade_alpha': 0.3,
                    #    },
                    #    'T = 4400K': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run3'],
                    #        'y_data': ni_data['run3'],
                    #        'line_type': 'b-',
                    #        'line_linewidth': 1,
                    #        'legend': True,
                    #    },
                    #    'T = 4400K (+1 sigma)': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run3'],
                    #        'y_data': ni_sd_up_data['run3'],
                    #        'line_type': 'b--',
                    #        'line_linewidth': 0.5
                    #    },
                    #    'T = 4400K (+1 sigma) shade': {
                    #        'type': dp.PlotType.shade,
                    #        'x_data': p_data['run3'],
                    #        'y_data': ni_data['run3'],
                    #        'y_shade_data': ni_sd_up_data['run3'],
                    #        'shade_colour': 'b',
                    #        'shade_alpha': 0.3,
                    #    },
                    #    'T = 4400K (-1 sigma)': {
                    #        'type': dp.PlotType.scatter_2d,
                    #        'x_data': p_data['run3'],
                    #        'y_data': ni_sd_down_data['run3'],
                    #        'line_type': 'b--',
                    #        'line_linewidth': 0.5
                    #    },
                    #    'T = 4400K (-1 sigma) shade': {
                    #        'type': dp.PlotType.shade,
                    #        'x_data': p_data['run3'],
                    #        'y_data': ni_data['run3'],
                    #        'y_shade_data': ni_sd_down_data['run3'],
                    #        'shade_colour': 'b',
                    #        'shade_alpha': 0.3,
                    #    },
                        'T = T(Peridotite liquidus)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': ni_data['run4'],
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'T = T(Peridotite liquidus) (+1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': ni_sd_up_data['run4'],
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'T = T(Peridotite liquidus) (+1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': p_data['run4'],
                            'y_data': ni_data['run4'],
                            'y_shade_data': ni_sd_up_data['run4'],
                            'shade_colour': 'r',
                            'shade_alpha': 0.3,
                        },
                        'T = T(Peridotite liquidus) (-1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': ni_sd_down_data['run4'],
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'T = T(Peridotite liquidus) (-1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': p_data['run4'],
                            'y_data': ni_data['run4'],
                            'y_shade_data': ni_sd_down_data['run4'],
                            'shade_colour': 'r',
                            'shade_alpha': 0.3,
                        },
                        'P_c marker': {
                            'type': dp.PlotType.vline,
                            'x_start': (1825. - 1420.)/104.42,
                            'y_min': 0,
                            'y_max': 100000,
                            'line_linewidth': 0.5,
                        },
                        'P_c text': {
                            'type': dp.PlotType.text,
                            'x_pos': 3.2,
                            'y_pos': 6,
                            'text_string': '$P_C$',
                        }
                    }
                }
            }
        },
        'plot4': {
            'show': False,
            'filenames': ['d_ni_3d.pdf'],
            'subplots': {
                'subplot1': {
                    'subplot_region': 111,
                    'title_text': r'Variation of $D_{Ni}$ with Pressure and Temperature',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 100,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'font': 'STIXGeneral',
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            'z_data': d_ni_array,
                            'cbar_label': r'$D_{Ni}$',
                            'z_scale': 'log'
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': constants.Tpdliq,
                            'x_start': 0,
                            'x_end': 100,
                            'x_points': 1000,
                            'line_type': 'k--',
                            'line_markersize': 1,
                            'line_linewidth': 2,
                            'legend': True
                        }
                    }
                }
            }
        },
        'plot5': {
            'show': False,
            'filenames': ['d_cr_ni_fo2.pdf'],
            'subplots': {
                'subplot1': {
                    'subplot_region': 111,
                    'legend': True,
                    'title_text': r'Variation of $D_{Cr}$ and $D_{Ni}$ with Oxygen Fugacity $f(O_2)$',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': '$f(O_2)$ relative to IW buffer',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$D$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'font': 'STIXGeneral',
                    'x_min': -7,
                    'x_max': 2,
                    'y_scale': 'log',
                    'series': {
                        'Cr': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_cr_fo2_array,
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'Cr v fO2 (+1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_cr_fo2_sd_up_array,
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'Cr v fO2 (+1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': fo2_range,
                            'y_data': d_cr_fo2_array,
                            'y_shade_data': d_cr_fo2_sd_up_array,
                            'shade_colour': 'r',
                            'shade_alpha': 0.3
                        },
                        'Cr v fO2 (-1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_cr_fo2_sd_down_array,
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'Cr v fO2 (-1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': fo2_range,
                            'y_data': d_cr_fo2_array,
                            'y_shade_data': d_cr_fo2_sd_down_array,
                            'shade_colour': 'r',
                            'shade_alpha': 0.3
                        },
                        'Ni': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_ni_fo2_array,
                            'line_type': 'k-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'Ni v fO2 (+1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_ni_fo2_sd_up_array,
                            'line_type': 'k--',
                            'line_linewidth': 0.5
                        },
                        'Ni v fO2 (+1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': fo2_range,
                            'y_data': d_ni_fo2_array,
                            'y_shade_data': d_ni_fo2_sd_up_array,
                            'shade_colour': 'k',
                            'shade_alpha': 0.3
                        },
                        'Ni v fO2 (-1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_ni_fo2_sd_down_array,
                            'line_type': 'k--',
                            'line_linewidth': 0.5
                        },
                        'Ni v fO2 (-1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': fo2_range,
                            'y_data': d_ni_fo2_array,
                            'y_shade_data': d_ni_fo2_sd_down_array,
                            'shade_colour': 'k',
                            'shade_alpha': 0.3
                        }
                    }
                }
            }
        },
        'plot6': {
            'show': False,
            'filenames': ['sf_cr.pdf'],
            'subplots': {
                'subplot1': {
                    'subplot_region': 111,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'title_text': r'Variation of $S_{Cr}$ with Pressure',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$S_{Cr}$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 100,
                    #'y_min': 0,
                    #'y_max': 100,
                    'font': 'STIXGeneral',
                    'series': {
                        'T = T(Peridotite liquidus)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': sf_cr_data['run4'],
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'T = T(Peridotite liquidus) (+1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': sf_cr_sd_up_data['run4'],
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'T = T(Peridotite liquidus) (+1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': p_data['run4'],
                            'y_data': sf_cr_data['run4'],
                            'y_shade_data': sf_cr_sd_up_data['run4'],
                            'shade_colour': 'r',
                            'shade_alpha': 0.3,
                        },
                        'T = T(Peridotite liquidus) (-1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': sf_cr_sd_down_data['run4'],
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'T = T(Peridotite liquidus) (-1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': p_data['run4'],
                            'y_data': sf_cr_data['run4'],
                            'y_shade_data': sf_cr_sd_down_data['run4'],
                            'shade_colour': 'r',
                            'shade_alpha': 0.3,
                        },
                        #'P_c marker': {
                        #    'type': dp.PlotType.vline,
                        #    'x_start': (1825. - 1420.)/104.42,
                        #    'y_min': 0,
                        #    'y_max': 110,
                        #    'line_linewidth': 0.5,
                        #},
                        #'P_c text': {
                        #    'type': dp.PlotType.text,
                        #    'x_pos': 3.4,
                        #    'y_pos': -5.2,
                        #    'text_string': '$P_C$',
                        #}
                    }
                }
            }
        },
        'plot7': {
            'show': False,
            'filenames': ['sf_ni.pdf'],
            'subplots': {
                'subplot1': {
                    'subplot_region': 111,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'title_text': r'Variation of $S_{Ni}$ with Pressure',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$S_{Ni}$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 100,
                    #'y_min': 0,
                    #'y_max': 100,
                    'font': 'STIXGeneral',
                    'series': {
                        'T = T(Peridotite liquidus)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': sf_ni_data['run4'],
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'T = T(Peridotite liquidus) (+1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': sf_ni_sd_up_data['run4'],
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'T = T(Peridotite liquidus) (+1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': p_data['run4'],
                            'y_data': sf_ni_data['run4'],
                            'y_shade_data': sf_ni_sd_up_data['run4'],
                            'shade_colour': 'r',
                            'shade_alpha': 0.3,
                        },
                        'T = T(Peridotite liquidus) (-1 sigma)': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': sf_ni_sd_down_data['run4'],
                            'line_type': 'r--',
                            'line_linewidth': 0.5
                        },
                        'T = T(Peridotite liquidus) (-1 sigma) shade': {
                            'type': dp.PlotType.shade,
                            'x_data': p_data['run4'],
                            'y_data': sf_ni_data['run4'],
                            'y_shade_data': sf_ni_sd_down_data['run4'],
                            'shade_colour': 'r',
                            'shade_alpha': 0.3,
                        },
                        #'P_c marker': {
                        #    'type': dp.PlotType.vline,
                        #    'x_start': (1825. - 1420.)/104.42,
                        #    'y_min': 0,
                        #    'y_max': 110,
                        #    'line_linewidth': 0.5,
                        #},
                        #'P_c text': {
                        #    'type': dp.PlotType.text,
                        #    'x_pos': 3.4,
                        #    'y_pos': -5.2,
                        #    'text_string': '$P_C$',
                        #}
                    }
                }
            }
        }
    }
    
    plotter = dp.DictPlotter(plot_dict)
    plotter.draw()
    plotter.yield_output()

if __name__ == '__main__':
    main()
