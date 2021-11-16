#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys

import chemistry_info as ci
import geology_info as gi
import graph_factory as gf
import pwd_utils as pu

sys.path.append(pu.get_path_to_utils())

import dict_plotter as dp

def generate_data():
    geo_model = gi.GeologyModel()
    
    Pmin = 0.1
    Pmax = 60
    
    Tmin = 2000
    Tmax = 4000
    
    T_dict = {
        'run1': 1400,
        'run2': 2900,
        'run3': 4400,
        'run4': None # ie use the liquidus in PAMELA
    }
    
    p_data = dict()
    cr_data = dict()
    cr_sd_up_data = dict()
    cr_sd_down_data = dict()
    ni_data = dict()
    ni_sd_up_data = dict()
    ni_sd_down_data = dict()
    
    #TODO: Just have a dict with multiple keys instead of multiple dicts
    fe_data = dict()
    c_data = dict()
    o_data = dict()
    si_data = dict()
    
    sf_cr_data = dict()
    sf_cr_sd_up_data = dict()
    sf_cr_sd_down_data = dict()
    sf_ni_data = dict()
    sf_ni_sd_up_data = dict()
    sf_ni_sd_down_data = dict()
    
    p_range = np.linspace(Pmin, Pmax, 200)
    t_range = np.linspace(Tmin, Tmax, 200)
    p_mesh, t_mesh = np.meshgrid(p_range, t_range)
    d_cr_array = np.zeros((len(t_range), len(p_range)))
    d_ni_array = np.zeros((len(t_range), len(p_range)))
    d_si_array = np.zeros((len(t_range), len(p_range)))
    kd_deviation_ni_array = np.zeros((len(t_range), len(p_range)))
    
    fo2_range = np.linspace(-3, -1, 50)
    d_cr_fo2_array = list()
    d_ni_fo2_array = list()
    d_cr_fo2_sd_up_array = list()
    d_cr_fo2_sd_down_array = list()
    d_ni_fo2_sd_up_array = list()
    d_ni_fo2_sd_down_array = list()
    
    #TODO: Just have a dict with multiple keys instead of multiple dicts
    d_fe_fo2_array = list()
    d_c_fo2_array = list()
    d_o_fo2_array = list()
    d_si_fo2_array = list()
    
    earth_pressure = geo_model.get_earth_differentiation_pressure()
    earth_fO2 = geo_model.get_earth_oxygen_fugacity()
    abundances_earth, cnf_earth, ds_earth, dummy = geo_model.form_a_planet_iteratively(earth_pressure, earth_fO2)
    print('Earth numbers')
    print(abundances_earth)
    print(ds_earth)
    print(cnf_earth)
    
    for run_name, T in T_dict.items():
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
        fe_data[run_name] = list()
        c_data[run_name] = list()
        o_data[run_name] = list()
        si_data[run_name] = list()
        
        for P in p_range:
            abundances, cnf, Ds, dummy = geo_model.form_a_planet_iteratively(P, earth_fO2, T)
            p_data[run_name].append(P)
            cr_data[run_name].append(Ds[ci.Element.Cr])
            #cr_sd_up_data[run_name].append(Ds_sd[ci.Element.Cr][0])
            #cr_sd_down_data[run_name].append(Ds_sd[ci.Element.Cr][1])
            ni_data[run_name].append(Ds[ci.Element.Ni])
            #ni_sd_up_data[run_name].append(Ds_sd[ci.Element.Ni][0])
            #ni_sd_down_data[run_name].append(Ds_sd[ci.Element.Ni][1])
            fe_data[run_name].append(Ds[ci.Element.Fe])
            c_data[run_name].append(Ds[ci.Element.C])
            o_data[run_name].append(Ds[ci.Element.O])
            si_data[run_name].append(Ds[ci.Element.Si])
    p_index = 0
    ni_valence = 2  # Could get this from PAMELA?
    log_kd_ni_earth = 0.3  # From Fischer+ 2015
    for P in p_range:
        t_index = 0
        for T in t_range:
            abundances, cnf, Ds, dummy = geo_model.form_a_planet_iteratively(P, earth_fO2, T)
            try:
                d_cr_array[t_index][p_index] = Ds[ci.Element.Cr]
                d_ni_array[t_index][p_index] = Ds[ci.Element.Ni]
                d_si_array[t_index][p_index] = Ds[ci.Element.Si]
                kd_ni = Ds[ci.Element.Ni]/(Ds[ci.Element.Fe]**(ni_valence/2))  # This should be independent of fO2!
                kd_deviation_ni_array[t_index][p_index] = min(abs(np.log10(kd_ni) - log_kd_ni_earth)/log_kd_ni_earth, 1)
            except TypeError:  # Then Ds is None
                d_cr_array[t_index][p_index] = None
                d_ni_array[t_index][p_index] = None
                d_si_array[t_index][p_index] = None
                kd_deviation_ni_array[t_index][p_index] = None
            t_index += 1
        p_index += 1
    print(kd_deviation_ni_array)
    
    for fo2 in fo2_range:
        abundances, cnf, Ds, dummy = geo_model.form_a_planet_iteratively(earth_pressure, fo2)
        try:
            d_cr_fo2_array.append(Ds[ci.Element.Cr])
            #d_cr_fo2_sd_up_array.append(Ds_sd[ci.Element.Cr][0])
            #d_cr_fo2_sd_down_array.append(Ds_sd[ci.Element.Cr][1])
        except TypeError:
            d_cr_fo2_array.append(None)
            #d_cr_fo2_sd_up_array.append(None)
            #d_cr_fo2_sd_down_array.append(None)
        try:
            d_ni_fo2_array.append(Ds[ci.Element.Ni])
            #d_ni_fo2_sd_up_array.append(Ds_sd[ci.Element.Ni][0])
            #d_ni_fo2_sd_down_array.append(Ds_sd[ci.Element.Ni][1])
        except TypeError:
            d_ni_fo2_array.append(None)
            #d_ni_fo2_sd_up_array.append(None)
            #d_ni_fo2_sd_down_array.append(None)
        try:
            d_fe_fo2_array.append(Ds[ci.Element.Fe])
        except TypeError:
            d_fe_fo2_array.append(None)
        try:
            d_c_fo2_array.append(Ds[ci.Element.C])
        except TypeError:
            d_c_fo2_array.append(None)
        try:
            d_o_fo2_array.append(Ds[ci.Element.O])
        except TypeError:
            d_o_fo2_array.append(None)
        try:
            d_si_fo2_array.append(Ds[ci.Element.Si])
        except TypeError:
            d_si_fo2_array.append(None)
    
    print(d_cr_fo2_array)
    print(d_cr_fo2_sd_up_array)
    print(d_cr_fo2_sd_down_array)
    print(d_ni_fo2_array)
    print(d_ni_fo2_sd_up_array)
    print(d_ni_fo2_sd_down_array)
    return (
        geo_model.config_name,
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
        d_si_array,
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
        kd_deviation_ni_array,
        fe_data,
        c_data,
        o_data,
        si_data,
        d_fe_fo2_array,
        d_c_fo2_array,
        d_o_fo2_array,
        d_si_fo2_array,
        ds_earth
    )

def main():
    (
        run_tag,
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
        d_si_array,
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
        kd_deviation_ni_array,
        fe_data,
        c_data,
        o_data,
        si_data,
        d_fe_fo2_array,
        d_c_fo2_array,
        d_o_fo2_array,
        d_si_fo2_array,
        ds_earth
    ) = generate_data()
    geo_model = gi.GeologyModel()
    plot_dict = {
        'plot1': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_cr.pdf'],
            'subplots': {
                'd_cr': {
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
                            #'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        #'T = T(Peridotite liquidus) (+1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run4'],
                        #    'y_data': cr_sd_up_data['run4'],
                        #    'line_type': 'r--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = T(Peridotite liquidus) (+1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run4'],
                        #    'y_data': cr_data['run4'],
                        #    'y_shade_data': cr_sd_up_data['run4'],
                        #    'shade_colour': 'r',
                        #    'shade_alpha': 0.3,
                        #},
                        #'T = T(Peridotite liquidus) (-1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run4'],
                        #    'y_data': cr_sd_down_data['run4'],
                        #    'line_type': 'r--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = T(Peridotite liquidus) (-1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run4'],
                        #    'y_data': cr_data['run4'],
                        #    'y_shade_data': cr_sd_down_data['run4'],
                        #    'shade_colour': 'r',
                        #    'shade_alpha': 0.3,
                        #},
                        'P_c marker': {
                            'type': dp.PlotType.vline,
                            'x_start': (1825. - 1420.)/104.42,
                            'y_min': 0,
                            'y_max': 110,
                            'line_linewidth': 0.5,
                            'line_styles': 'dashed'
                        },
                        'P_c text': {
                            'type': dp.PlotType.text,
                            'x_pos': 3.4,
                            'y_pos': -5.2,
                            'text_string': '$P_C$',
                            'horizontalalignment': 'center'
                        }
                    }
                }
            }
        },
        'plot2': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_cr_3d.png'],
            'dpi': 150,
            'subplots': {
                'd_cr_3d': {
                    'subplot_region': 111,
                    'title_text': r'Variation of $D_\textrm{Cr}$ with Pressure and Temperature',
                    'title_fontsize': 16,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 16,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 16,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 60,
                    'y_max': 4000,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 14,
                    'font': 'STIXGeneral',
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            'z_data': d_cr_array,
                            'fill': True,
                            'cbar_label': r'$D_{Cr}$',
                            'levels': np.linspace(d_cr_array.min(), d_cr_array.max(), 20),
                            'cbar_ticks': np.linspace(1.5, 6.5, 6),
                            'cbar_labelfontsize': 12,
                            'cbar_shrink': 0.9,
                            'cbar_labelpad': 25
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': geo_model.pamela.peridotite_liquidus,
                            'x_start': 0,
                            'x_end': 100,
                            'x_points': 1000,
                            'line_color': 'r',
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
            'filenames': ['graphs/' + run_tag + '_d_ni.pdf'],
            'subplots': {
                'd_ni': {
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
                        #'T = T(Peridotite liquidus) (+1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run4'],
                        #    'y_data': ni_sd_up_data['run4'],
                        #    'line_type': 'r--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = T(Peridotite liquidus) (+1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run4'],
                        #    'y_data': ni_data['run4'],
                        #    'y_shade_data': ni_sd_up_data['run4'],
                        #    'shade_colour': 'r',
                        #    'shade_alpha': 0.3,
                        #},
                        #'T = T(Peridotite liquidus) (-1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': p_data['run4'],
                        #    'y_data': ni_sd_down_data['run4'],
                        #    'line_type': 'r--',
                        #    'line_linewidth': 0.5
                        #},
                        #'T = T(Peridotite liquidus) (-1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': p_data['run4'],
                        #    'y_data': ni_data['run4'],
                        #    'y_shade_data': ni_sd_down_data['run4'],
                        #    'shade_colour': 'r',
                        #    'shade_alpha': 0.3,
                        #},
                        'P_c marker': {
                            'type': dp.PlotType.vline,
                            'x_start': (1825. - 1420.)/104.42,
                            'y_min': 0,
                            'y_max': 100000,
                            'line_linewidth': 0.5
                            #'line_styles': '--'
                        },
                        'P_c text': {
                            'type': dp.PlotType.text,
                            'x_pos': 3.2,
                            'y_pos': 6,
                            'text_string': '$P_C$',
                            'horizontalalignment': 'center'
                        }
                    }
                }
            }
        },
        'plot4': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_ni_3d.png'],
            'dpi': 150,
            'subplots': {
                'd_ni_3d': {
                    'subplot_region': 111,
                    'title_text': r'Variation of $D_\textrm{Ni}$ with Pressure and Temperature',
                    'title_fontsize': 16,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 16,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 16,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 60,
                    'y_max': 4000,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 14,
                    'font': 'STIXGeneral',
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            'z_data': np.log10(d_ni_array),
                            'fill': True,
                            'cbar_label': r'$log(D_{Ni})$',
                            #'z_scale': 'log',
                            'levels': np.linspace(np.log10(d_ni_array.min()), np.log10(d_ni_array.max()), 20),
                            'cbar_ticks': np.linspace(1.5, 3.5, 5),
                            'cbar_labelfontsize': 12,
                            'cbar_shrink': 0.9,
                            'cbar_labelpad': 25
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': geo_model.pamela.peridotite_liquidus,
                            'x_start': 0,
                            'x_end': 100,
                            'x_points': 1000,
                            'line_color': 'r',
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
            'filenames': ['graphs/' + run_tag + '_d_cr_ni_fo2.pdf'],
            'subplots': {
                'd_cr_ni_fo2': {
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
                        #'Cr v fO2 (+1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': fo2_range,
                        #    'y_data': d_cr_fo2_sd_up_array,
                        #    'line_type': 'r--',
                        #    'line_linewidth': 0.5
                        #},
                        #'Cr v fO2 (+1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': fo2_range,
                        #    'y_data': d_cr_fo2_array,
                        #    'y_shade_data': d_cr_fo2_sd_up_array,
                        #    'shade_colour': 'r',
                        #    'shade_alpha': 0.3
                        #},
                        #'Cr v fO2 (-1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': fo2_range,
                        #    'y_data': d_cr_fo2_sd_down_array,
                        #    'line_type': 'r--',
                        #    'line_linewidth': 0.5
                        #},
                        #'Cr v fO2 (-1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': fo2_range,
                        #    'y_data': d_cr_fo2_array,
                        #    'y_shade_data': d_cr_fo2_sd_down_array,
                        #    'shade_colour': 'r',
                        #    'shade_alpha': 0.3
                        #},
                        'Ni': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_ni_fo2_array,
                            'line_type': 'k-',
                            'line_linewidth': 1,
                            'legend': True
                        }
                        #'Ni v fO2 (+1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': fo2_range,
                        #    'y_data': d_ni_fo2_sd_up_array,
                        #    'line_type': 'k--',
                        #    'line_linewidth': 0.5
                        #},
                        #'Ni v fO2 (+1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': fo2_range,
                        #    'y_data': d_ni_fo2_array,
                        #    'y_shade_data': d_ni_fo2_sd_up_array,
                        #    'shade_colour': 'k',
                        #    'shade_alpha': 0.3
                        #},
                        #'Ni v fO2 (-1 sigma)': {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': fo2_range,
                        #    'y_data': d_ni_fo2_sd_down_array,
                        #    'line_type': 'k--',
                        #    'line_linewidth': 0.5
                        #},
                        #'Ni v fO2 (-1 sigma) shade': {
                        #    'type': dp.PlotType.shade,
                        #    'x_data': fo2_range,
                        #    'y_data': d_ni_fo2_array,
                        #    'y_shade_data': d_ni_fo2_sd_down_array,
                        #    'shade_colour': 'k',
                        #    'shade_alpha': 0.3
                        #}
                    }
                }
            }
        },
        'plot6': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_kd_deviation_ni_3d.pdf'],
            'subplots': {
                'd_ni_3d': {
                    'subplot_region': 111,
                    'title_text': r'Relative error in $log(KD_{Ni})$, capped at 1',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'y_min': 2000,
                    'y_max': 4000,
                    'x_min': 35,
                    'x_max': 75,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'font': 'STIXGeneral',
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            'z_data': kd_deviation_ni_array,
                            'fill': True,
                            'cbar_label': r'Error',
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': geo_model.pamela.peridotite_liquidus,
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
        'plot7': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_v_p_mult.pdf'],
            'subplots': {
                'd_mult': {
                    'subplot_region': 111,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'title_text': r'Variation of $D_{element}$ with Pressure at fO2 = IW - 2',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$D$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 100,
                    'y_min': 0.00001,
                    'y_max': 1000000,
                    'y_scale': 'log',
                    'font': 'STIXGeneral',
                    'series': {
                        'Ni': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': ni_data['run4'],
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'Cr': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': cr_data['run4'],
                            'line_type': 'b-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'Fe': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': fe_data['run4'],
                            'line_type': 'k-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'C': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': c_data['run4'],
                            'line_type': 'c-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'O': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': o_data['run4'],
                            'line_type': 'y-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'Si': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': si_data['run4'],
                            'line_type': 'g-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'P_c marker': {
                            'type': dp.PlotType.vline,
                            'x_start': (1825. - 1420.)/104.42,
                            'y_min': 0,
                            'y_max': 1000000,
                            'line_linewidth': 0.5,
                        },
                        'P_c text': {
                            'type': dp.PlotType.text,
                            'x_pos': 5.5,
                            'y_pos': 0.00004,
                            'text_string': '$P_C$',
                            'horizontalalignment': 'center'
                        }
                    }
                }
            }
        },
        'plot8': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_v_fO2_mult.pdf'],
            'subplots': {
                'd_mult': {
                    'subplot_region': 111,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'title_text': r'Variation of $D_{element}$ with Oxygen Fugacity',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Oxygen Fugacity ( ' + r'$\Delta$' + 'IW)',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$D$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    #'x_min': 0,
                    #'x_max': 100,
                    #'y_min': 10,
                    #'y_max': 100000,
                    'y_scale': 'log',
                    'font': 'STIXGeneral',
                    'series': {
                        'Ni': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_ni_fo2_array,
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'Cr': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_cr_fo2_array,
                            'line_type': 'b-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'Fe': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_fe_fo2_array,
                            'line_type': 'k-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'C': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_c_fo2_array,
                            'line_type': 'c-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'O': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_o_fo2_array,
                            'line_type': 'y-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'Si': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': d_si_fo2_array,
                            'line_type': 'g-',
                            'line_linewidth': 1,
                            'legend': True
                        }
                    }
                }
            }
        },
        'plot9': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_v_p_rel_mult.pdf'],
            'subplots': {
                'd_mult': {
                    'subplot_region': 111,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'title_text': r'Relative variation of $D_{element}$ with Pressure at fO2 = IW - 2',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$D$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 65,
                    'y_min': 0.00001,
                    'y_max': 100,
                    'y_scale': 'log',
                    'font': 'STIXGeneral',
                    'series': {
                        'Ni': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': [p/ds_earth[ci.Element.Ni] for p in ni_data['run4']],
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'Cr': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': [p/ds_earth[ci.Element.Cr] for p in cr_data['run4']],
                            'line_type': 'b-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'Fe': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': [p/ds_earth[ci.Element.Fe] for p in fe_data['run4']],
                            'line_type': 'k-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'C': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': [p/ds_earth[ci.Element.C] for p in c_data['run4']],
                            'line_type': 'c-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'O': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': [p/ds_earth[ci.Element.O] for p in o_data['run4']],
                            'line_type': 'y-',
                            'line_linewidth': 1,
                            'legend': True,
                        },
                        'Si': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': p_data['run4'],
                            'y_data': [p/ds_earth[ci.Element.Si] for p in si_data['run4']],
                            'line_type': 'g-',
                            'line_linewidth': 1,
                            'legend': True,
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
                            'horizontalalignment': 'center'
                        }
                    }
                }
            }
        },
        'plot10': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_v_fO2_rel_mult.pdf'],
            'subplots': {
                'd_mult': {
                    'subplot_region': 111,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'title_text': r'Variation of $D_{element}$ with Oxygen Fugacity',
                    'title_fontsize': 12,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Oxygen Fugacity (' + r'$\Delta$' + 'IW)',
                    'xlabel_fontsize': 10,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': r'$D$',
                    'ylabel_fontsize': 10,
                    'ylabel_fontweight': 'bold',
                    'x_min': -3,
                    'x_max': -1.5,
                    'y_min': 0.3,
                    'y_max': 4,
                    'y_scale': 'log',
                    'font': 'STIXGeneral',
                    'series': {
                        'Ni': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': [p/ds_earth[ci.Element.Ni] for p in d_ni_fo2_array],
                            'line_type': 'r-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'Cr': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': [p/ds_earth[ci.Element.Cr] for p in d_cr_fo2_array],
                            'line_type': 'b-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'Fe': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': [p/ds_earth[ci.Element.Fe] for p in d_fe_fo2_array],
                            'line_type': 'k-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'C': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': [p/ds_earth[ci.Element.C] for p in d_c_fo2_array],
                            'line_type': 'c-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'O': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': [p/ds_earth[ci.Element.O] for p in d_o_fo2_array],
                            'line_type': 'y-',
                            'line_linewidth': 1,
                            'legend': True
                        },
                        'Si': {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': fo2_range,
                            'y_data': [p/ds_earth[ci.Element.Si] for p in d_si_fo2_array],
                            'line_type': 'g-',
                            'line_linewidth': 1,
                            'legend': True
                        }
                    }
                }
            }
        },
        'plot11': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_si_3d.png'],
            'dpi': 150,
            'subplots': {
                'd_si_3d': {
                    'subplot_region': 111,
                    'title_text': r'Variation of $D_\textrm{Si}$ with Pressure and Temperature',
                    'title_fontsize': 16,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 16,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 16,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 60,
                    'y_max': 4000,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 14,
                    'font': 'STIXGeneral',
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            #'z_data': d_si_array,
                            'z_data': np.log10(d_si_array),
                            'fill': True,
                            'cbar_label': r'$log(D_{Si})$',
                            #'z_scale': 'log',
                            'levels': np.linspace(np.log10(d_si_array.min()), np.log10(d_si_array.max()), 20),
                            'cbar_ticks': np.linspace(-2.8, 0, 5),
                            'cbar_labelfontsize': 12,
                            'cbar_shrink': 0.9,
                            'cbar_labelpad': 25
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': geo_model.pamela.peridotite_liquidus,
                            'x_start': 0,
                            'x_end': 100,
                            'x_points': 1000,
                            'line_color': 'r',
                            'line_markersize': 1,
                            'line_linewidth': 2,
                            'legend': True
                        }
                    }
                }
            }
        },
        'plot12': {
            'show': False,
            'filenames': ['graphs/' + run_tag + '_d_cr_ni_si_3d.pdf'],
            'dpi': 150,
            'fig_width': 7,
            'fig_height': 12,
            'gridspec_y_x': (3, 1),
            'gridspec_hspace': 0,
            'subplots': {
                'd_cr_3d': {
                    'subplot_region': 111,
                    #'title_text': r'Variation of $D_\textrm{Cr}$ with Pressure and Temperature',
                    'title_fontsize': 16,
                    'title_fontweight': 'bold',
                    #'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 16,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 16,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 60,
                    'y_max': 4000,
                    'x_tick_fontsize': 14,
                    'y_tick_fontsize': 13,
                    'legend': True,
                    'legend_loc': 'best',
                    'legend_text_size': 14,
                    'font': 'STIXGeneral',
                    'gridspec_index': 0,
                    'sharex_subplot': 'd_ni_3d',
                    'y_hide_ticks': [0],
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            'z_data': d_cr_array,
                            'fill': True,
                            'cbar_label': r'$D_\textrm{Cr}$',
                            'levels': np.linspace(d_cr_array.min(), d_cr_array.max(), 20),
                            'cbar_ticks': np.linspace(1.5, 6.5, 6),
                            'cbar_labelfontsize': 16,
                            'cbar_tickfontsize': 14,
                            'cbar_shrink': 0.9,
                            'cbar_labelpad': 25,
                            'colour_map_colours': ['#fde725','#5ec962','#21918c', '#3b528b', '#440154'],
                            'colour_map_min': d_cr_array.min(),
                            'colour_map_max': d_cr_array.max()
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': geo_model.pamela.peridotite_liquidus,
                            'x_start': 0,
                            'x_end': 100,
                            'x_points': 1000,
                            'line_color': 'r',
                            'line_markersize': 1,
                            'line_linewidth': 2,
                            'legend': True
                        }
                    }
                },
                'd_ni_3d': {
                    'subplot_region': 111,
                    #'title_text': r'Variation of $D_\textrm{Ni}$ with Pressure and Temperature',
                    'title_fontsize': 16,
                    'title_fontweight': 'bold',
                    #'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 16,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 16,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 60,
                    'y_max': 4000,
                    'x_tick_fontsize': 14,
                    'y_tick_fontsize': 13,
                    'legend': False,
                    'legend_loc': 'best',
                    'legend_text_size': 14,
                    'font': 'STIXGeneral',
                    'gridspec_index': 1,
                    'sharex_subplot': 'd_si_3d',
                    'y_hide_ticks': [0],
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            'z_data': np.log10(d_ni_array),
                            'fill': True,
                            'cbar_label': r'$\textrm{log}_{10}(D_\textrm{Ni})$',
                            #'z_scale': 'log',
                            #'levels': np.linspace(np.log10(d_ni_array.min()), np.log10(d_ni_array.max()), 20),
                            'levels': np.linspace(np.log10(d_ni_array.min()), np.log10(d_ni_array.max()), 20),
                            'cbar_ticks': np.linspace(1.5, 3.5, 5),
                            'cbar_labelfontsize': 16,
                            'cbar_tickfontsize': 14,
                            'cbar_shrink': 0.9,
                            'cbar_labelpad': 25,
                            'cbar_labelrotation': 90,
                            'colour_map_colours': ['#fde725','#5ec962','#21918c', '#3b528b', '#440154'],
                            'colour_map_min': np.log10(d_ni_array.min()),
                            'colour_map_max': np.log10(d_ni_array.min()) + (np.log10(d_si_array.max()) - np.log10(d_si_array.min()))
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': geo_model.pamela.peridotite_liquidus,
                            'x_start': 0,
                            'x_end': 100,
                            'x_points': 1000,
                            'line_color': 'r',
                            'line_markersize': 1,
                            'line_linewidth': 2,
                            'legend': True
                        }
                    }
                },
                'd_si_3d': {
                    'subplot_region': 111,
                    #'title_text': r'Variation of $D_\textrm{Si}$ with Pressure and Temperature',
                    'title_fontsize': 16,
                    'title_fontweight': 'bold',
                    'xlabel_text': 'Pressure / GPa',
                    'xlabel_fontsize': 16,
                    'xlabel_fontweight': 'bold',
                    'ylabel_text': 'Temperature / K',
                    'ylabel_fontsize': 16,
                    'ylabel_fontweight': 'bold',
                    'x_min': 0,
                    'x_max': 60,
                    'y_max': 4000,
                    'x_tick_fontsize': 14,
                    'y_tick_fontsize': 13,
                    'legend': False,
                    'legend_loc': 'best',
                    'legend_text_size': 14,
                    'font': 'STIXGeneral',
                    'gridspec_index': 2,
                    'series': {
                        '3dscatter': {
                            'type': dp.PlotType.contour_scatter,
                            'x_data': p_range,
                            'y_data': t_range,
                            #'z_data': d_si_array,
                            'z_data': np.log10(d_si_array),
                            'fill': True,
                            'cbar_label': r'$\textrm{log}_{10}(D_\textrm{Si})$',
                            #'z_scale': 'log',
                            'levels': np.linspace(np.log10(d_si_array.min()), np.log10(d_si_array.max()), 20),
                            'cbar_ticks': np.linspace(-2.8, 0, 5),
                            'cbar_labelfontsize': 16,
                            'cbar_tickfontsize': 14,
                            'cbar_shrink': 0.9,
                            'cbar_labelpad': 25,
                            'cbar_labelrotation': 90,
                            'colour_map_colours': ['#fde725','#5ec962','#21918c', '#3b528b', '#440154'],
                            'colour_map_min': np.log10(d_si_array.min()),
                            'colour_map_max': np.log10(d_si_array.max())
                        },
                        'Peridotite liquidus': {
                            'type': dp.PlotType.function_2d,
                            'function': geo_model.pamela.peridotite_liquidus,
                            'x_start': 0,
                            'x_end': 100,
                            'x_points': 1000,
                            'line_color': 'r',
                            'line_markersize': 1,
                            'line_linewidth': 2,
                            'legend': True
                        }
                    }
                }
            }
        }
    }
    
    plotter = dp.DictPlotter(plot_dict)
    plotter.draw()
    plotter.yield_output()
    print(ds_earth[ci.Element.Ni])
    print(ni_data['run4'])
    print([p/ds_earth[ci.Element.Ni] for p in ni_data['run4']])
    
if __name__ == '__main__':
    main()
