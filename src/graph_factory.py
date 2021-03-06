#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections as cn
import matplotlib.pyplot as plt
import math
import numpy as np
import sys

import chemistry_info as ci
import geology_info as gi
import live_data as ld
import pwd_utils as pu

sys.path.append(pu.get_path_to_utils())

import dict_plotter as dp

class GraphFactory:
    
    def __init__(self, output_dir='graphs'):
        self.output_dir = output_dir
        
        # This should contain info such as the system we're looking at
        self.unit_dict = {
            'Pressure': 'GPa',
            'fO2': 'ΔIW',
            'fcf': '',
            'Sulfur Content': ''
        }
        self.symbol_dict = {
            'Pressure': 'P',
            'fO2': 'fO2',
            'fcf': 'fcf',
            'Sulfur Content': 'X_S'
        }
        
        self.colour_list = [
            '#000000',
            #'#FF0000',
            '#1355D2',
            '#DF4443',
            '#CD853F',
            '#24C024',
            '#FFD700',
            '#87CEFA',
            '#32CD32',
            '#BC8F8F',
            '#808000',
            '#FF7F50',
            '#bbbaaa',
            '#7FFFD4',
            '#5F9EA0',
            '#191970',
            '#FFA500',
            '#FFA07A',
            '#C0C0C0',
            '#800000',
            '#BA55D3',
            '#66CDAA',
            '#FF0FD5',
            '#B0E0E6',
            '#8A2BE2',
            '#98FB98',
            '#F0E68C',
            '#006400',
            '#8B0000',
            '#DAA520',
            '#D8BFD8'
        ]
        
        self.rainbow = [
            '#DF0000',
            '#DFA020',
            #'#BFBF00',
            '#00A444',
            '#0000FF',
            '#4B0082',
            '#EE82EE'
        ]
        
        self.colour_dict = {
            ci.Element.Ca: '#191970',
            ci.Element.Al: '#66CDAA',
            ci.Element.Fe: '#C0C0C0',
            ci.Element.Si: '#DF2E20',
            ci.Element.Mg: '#87CEFA',
            ci.Element.O: '#FFD700',
            ci.Element.C: '#7FFFD4',
            ci.Element.Cr: '#808000',
            ci.Element.Ni: '#22CD22',
            ci.Element.Placeholder: '#000000'
        }
        
    def colour_mixer(self, c1, c2, ratio):
        c1 = c1.lstrip('#')
        c2 = c2.lstrip('#')
        c1_rgb = [int(c1[i:i+2], 16) for i in [0, 2, 4]]
        c2_rgb = [int(c2[i:i+2], 16) for i in [0, 2, 4]]
        toret_rgb = [int(((1-ratio)*c1_rgb[i]) + (ratio*c2_rgb[i])) for i in [0, 1, 2]]
        print(toret_rgb)
        return '#%02x%02x%02x' % (toret_rgb[0], toret_rgb[1], toret_rgb[2])
    
    def discrete_cmap(self, N, base_cmap=None):
        base = plt.cm.get_cmap(base_cmap)
        color_list = base(np.linspace(0, 1, N))
        cmap_name = base.name + str(N)
        return base.from_list(cmap_name, color_list, N)
    
    def convert_to_relative_solar_abundance(self, original_abundances, original_errors=None):
        toret = list()
        toret_errs = list()
        elements_present = list()
        geo_model = gi.GeologyModel()
        for element, value in original_abundances.items():
            if value is np.nan or value is None:
                pass
            else:
                toret.append(original_abundances[element] - original_abundances[ci.Element.Mg] - ld._geo_model.solar_abundances[element])  # Technically this should also bear stellar Mg in mind. Can afford to ignore here because it's always 0 (since we report relative to Mg)
                elements_present.append(element)
                if original_errors is not None:
                    # Absolute errors should add in quadriture for addition. TODO: Add error on solar abundance
                    if element == ci.Element.Mg:
                        toret_errs.append(original_errors[ci.Element.Mg])
                    else:
                        toret_errs.append(np.sqrt((original_errors[element])**2 + (original_errors[ci.Element.Mg])**2))
        return toret, toret_errs, elements_present
    
    # Assume all abundances are logarithmic relative to some reference X. So if we are interested in an element E
    # log[(E/Mg)/(E/Mg)_solar] = log(E/X) - log(Mg/X) - log(E/Mg)_solar = observed_abundances[E] - observed_abundances[Mg] - some constant
    # Assume inputs are dicts with ci.Elements as keys. Values can be np.nans
    # The fitted_abundances_list and fitted_abundance_errors are a dict of those dicts (with the key as an arbitrary string label), because we may want to plot multiple fits
    def make_composition_plot(self, file_prefix, observed_abundances, observed_abundance_errors, fitted_abundances_dict, fitted_abundance_errors_dict=None):
        observations, errors, x_tick_labels = self.convert_to_relative_solar_abundance(observed_abundances, observed_abundance_errors)
        fits = dict()
        i = 0
        for key, fa in fitted_abundances_dict.items():
            fits[key] = self.convert_to_relative_solar_abundance(fa)[0]  # discard errors and element labels. TODO: Verify here that the fitted abundances have the same elements present!
            i += 1
        series_dict = {
            'observations': {
                'type': dp.PlotType.scatter_2d_error,
                'x_data': range(0, len(observations)),
                'y_data': observations,
                'y_error_data': errors,
                'line_type': 'ko',
                'line_linewidth': 1,
                'legend': True,
            }
        }
        colours = ['r', 'g', 'b', 'c', 'm', 'y'] # TODO: Expand on this and loop around if necessary
        fit_index = 0
        for series_name, series in fits.items():
            series_dict[series_name] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': range(0, len(series)),
                'y_data': series,
                'line_type': colours[fit_index % len(colours)] + '-',
                'line_linewidth': 1,
                'legend': True,
            }
            fit_index += 1
        plot_dict = {
            'testplot': {
                'show': False,
                'filenames': ['composition_' + file_prefix + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': r'Composition of ' + file_prefix + ' pollutant',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Element',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': '$\mathrm{log((X/Mg)/(X/Mg)}_{\mathrm{mean \, stellar}})$',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'x_tick_labels': x_tick_labels,
                        #'x_min': 0,
                        #'x_max': 100,
                        #'y_min': 0,
                        #'y_max': 100,
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        
    def make_composition_plot_mk2(self, system, file_prefix, all_wd_abundances, all_wd_abundance_errors, fit_dict, error_low_dict=None, error_high_dict=None, upper_limits=None, lower_limits=None, video=False, hack_legend_to_only_show_pressure=False):
        x_axis = ['Al','Ti','Ca','Ni','Fe','Cr','Si','Na','O','C','N']
        whitedwarfabundancesplot = list()
        whitedwarferrorsplot = list()
        geo_model_for_solar_abundances = gi.GeologyModel()  # We're only using this to access the solar abundances so no arguments needed
        upper_limits_for_plot = list()  # a list of bools to track which of the abundances are actually upper/lower bounds
        lower_limits_for_plot = list()
        upper_abundances_for_plot = list()
        rel_solar_abundances_for_plot = list()
        lower_abundances_for_plot = list()
        for el, abundance in all_wd_abundances.items():
            if el == ci.Element.Mg:
                # Skip magnesium
                pass
            elif np.isnan(abundance) or abundance == 0:
                upper_abundances_for_plot.append(geo_model_for_solar_abundances.upper_ratiod_to_solar[el])
                rel_solar_abundances_for_plot.append(0)
                lower_abundances_for_plot.append(geo_model_for_solar_abundances.lower_ratiod_to_solar[el])
                if upper_limits is not None and upper_limits.get(el) is not None:
                    whitedwarfabundancesplot.append(upper_limits[el] - all_wd_abundances[ci.Element.Mg] - geo_model_for_solar_abundances.solar_abundances[el])
                    whitedwarferrorsplot.append(0.1) # This sets arrow size
                    upper_limits_for_plot.append(True)
                    lower_limits_for_plot.append(False)
                elif lower_limits is not None and lower_limits.get(el) is not None:
                    whitedwarfabundancesplot.append(lower_limits[el] - all_wd_abundances[ci.Element.Mg] - geo_model_for_solar_abundances.solar_abundances[el])
                    whitedwarferrorsplot.append(0.1) # This sets arrow size
                    upper_limits_for_plot.append(False)
                    lower_limits_for_plot.append(True)
                else:
                    whitedwarfabundancesplot.append(np.nan)
                    whitedwarferrorsplot.append(0)
                    upper_limits_for_plot.append(False)
                    lower_limits_for_plot.append(False)
            else:
                upper_abundances_for_plot.append(geo_model_for_solar_abundances.upper_ratiod_to_solar[el])
                rel_solar_abundances_for_plot.append(0)
                lower_abundances_for_plot.append(geo_model_for_solar_abundances.lower_ratiod_to_solar[el])
                whitedwarfabundancesplot.append(all_wd_abundances[el] - all_wd_abundances[ci.Element.Mg] - geo_model_for_solar_abundances.solar_abundances[el])
                whitedwarferrorsplot.append(np.sqrt((all_wd_abundance_errors[el]**2) + (all_wd_abundance_errors[ci.Element.Mg]**2)))
                upper_limits_for_plot.append(False)
                lower_limits_for_plot.append(False)
        print(whitedwarfabundancesplot)
        print(whitedwarferrorsplot)
        

        obs_data_series_name = self.strip_system_suffix(system)
        obs_data_series_name += ' Observational Data'
        
        series_dict = {
            obs_data_series_name: {
                'type': dp.PlotType.scatter_2d_error,
                #'x_data': range(0, len(observations)),
                'x_data': x_axis,
                'y_data': whitedwarfabundancesplot,
                'y_error_data': whitedwarferrorsplot,
                'line_type': 'ko',
                'line_linewidth': 1,
                'capsize': 2,
                'capthick': 1,
                'legend': True,
                'upper_limits': upper_limits_for_plot,
                'lower_limits': lower_limits_for_plot
            }
        }
        colours = self.colour_list
        colour_index = 0
        min_fit_value = None
        for fit_key, fit_data in fit_dict.items():
            if min_fit_value is None:
                min_fit_value = min(fit_data)
            else:
                min_fit_value = min(min_fit_value, min(fit_data))  # Potential improvement: only check values for which there is an observation
            if video or hack_legend_to_only_show_pressure:
                if fit_key == 'LP run':
                    fit_key_to_use = 'Median model (? GPa)'  # This is a hack
                else:
                    print(fit_key)
                    print(fit_key[0])
                    print(fit_key[1])
                    print(fit_key[2])
                    len_pressure = len(str(fit_key[1]))
                    padding_space = ' '*(5 - len_pressure)
                    print(fit_key[1])
                    print(len_pressure)
                    print(padding_space)
                    print(str(fit_key[1]) + padding_space + ' GPa')
                    fit_key_to_use = 'Fit at ' + str(fit_key[1]) + padding_space + ' GPa'
                print(fit_key_to_use)
            else:
                fit_key_to_use = fit_key
            series_dict[fit_key_to_use] = {
                'type': dp.PlotType.scatter_2d,
                #'x_data': range(0, len(observations)),
                'x_data': x_axis,
                'y_data': fit_data,
                'line_color': colours[colour_index],
                'line_style': '-',
                'line_linewidth': 2,
                'zorder': 2,
                'legend': True
            }
            if error_low_dict is not None and error_low_dict.get(fit_key) is not None and error_high_dict is not None and error_high_dict.get(fit_key) is not None:
                if min_fit_value is None:
                    min_fit_value = min(error_low_dict[fit_key])
                else:
                    min_fit_value = min(min_fit_value, min(error_low_dict[fit_key]))
                series_dict['1 Sigma Confidence Interval'] = {
                    'type': dp.PlotType.shade,
                    #'x_data': range(0, len(observations)),
                    'x_data': x_axis,
                    'y_data': error_low_dict[fit_key],
                    'y_shade_data': error_high_dict[fit_key],
                    'shade_colour': 'k',
                    'shade_alpha': 0.5,
                    'zorder': 2,
                    'legend': True
                }
            colour_index += 1
                
        series_dict['Solar Composition'] = {
            'type': dp.PlotType.scatter_2d,
            #'x_data': range(0, len(observations)),
            'x_data': x_axis,
            'y_data': rel_solar_abundances_for_plot,
            'line_color': 'black',
            'line_style': '--',
            'line_linewidth': 1,
            'zorder': 1,
            'legend': True
        }
        series_dict['Stellar 98\% Composition Range'] = {
            'type': dp.PlotType.shade,
            #'x_data': range(0, len(observations)),
            'x_data': x_axis,
            'y_data': lower_abundances_for_plot,
            'y_shade_data': upper_abundances_for_plot,
            'shade_colour': 'silver',
            'line_linewidth': 1,
            'zorder': 1,
            'legend': True
        }
        
        if video:
            series_dict['Pressure Text'] = {
                'type': dp.PlotType.text,
                'x_pos': 6,
                'y_pos': -1,
                'text_string': str(fit_key[1]),
                'horizontalalignment': 'left',
                'fontsize': 20,
                'legend': False
            }
            series_dict['Pressure unit Text'] = {
                'type': dp.PlotType.text,
                'x_pos': 6.7,
                'y_pos': -1,
                'text_string': 'GPa',
                'horizontalalignment': 'left',
                'fontsize': 20,
                'legend': False
            }
        extension = ['.png'] if video else ['.pdf', '.png']
        #try:
        #    y_min = -2 if min_fit_value < -2 else None
        #    y_max = 1 if video else (2 if min_fit_value < -10 else None)  # For some reason y_max starts misbehaving if the min value is very small
        #except TypeError:
        #    y_min = None
        #    y_max = None
        title_text = r'Composition of ' + obs_data_series_name + ' pollutant'
        plot_dict = {
            'complotmk2': {
                'show': False,
                'filenames': [system + '_' + file_prefix + 'compositionmk2' + e for e in extension],
                'fig_height': 5,
                'fig_width': 10,
                'dpi': 300,
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 12,
                        #'title_text': title_text,
                        'title_fontsize': 16,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Lithophiles \, \, \, \, \, \, \, \, \, \, \, \, Siderophiles \, \, \, \, \, \, \, Volatile Lithophiles \, \, \, \, Atmophiles',
                        'xlabel_fontsize': 16,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': '[X/Mg]',#'$\mathrm{log((X/Mg)/(X/Mg)}_{\mathrm{Solar}})$',
                        'ylabel_fontsize': 20,
                        'ylabel_fontweight': 'bold',
                        'x_tick_labels': x_axis,
                        'x_tick_fontsize': 20,
                        'y_tick_fontsize': 17,
                        'ylabel_pad': 0,
                        #'x_min': 0,
                        #'x_max': 100,
                        #'y_min': -1.2, # I used -1.2 for GD61
                        #'y_max': 0.9,  # I used 0.9 for GD61
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        if video:
            plotter.yield_output('<your_filepath_here>')
        else:
            plotter.yield_output(self.output_dir)

    def strip_system_suffix(self, system_name):
        obs_data_series_name = system_name
        strippable_suffixes = ['Corr', 'NoNa', 'X', 'SiO', 'Bounds']
        for suffix in strippable_suffixes:
            if obs_data_series_name.endswith(suffix):
                len_to_strip = len(suffix)
                obs_data_series_name = obs_data_series_name[:-len_to_strip]
        #if obs_data_series_name.endswith('Corr'):
        #    obs_data_series_name = obs_data_series_name[:-4]
        #if obs_data_series_name.endswith('NoNa'):
        #    obs_data_series_name = obs_data_series_name[:-4]
        #elif obs_data_series_name.endswith('X'):
        #    obs_data_series_name = obs_data_series_name[:-1]
        #elif obs_data_series_name.endswith('SiO'):
        #    obs_data_series_name = obs_data_series_name[:-3]
        #elif obs_data_series_name.endswith('Bounds'):
        #    obs_data_series_name = obs_data_series_name[:-6]
        return obs_data_series_name

    def make_composition_plot_raw(self, system, wd_type, file_prefix, all_wd_abundances, all_wd_abundance_errors, fit_dict, error_low_dict=None, error_high_dict=None, upper_limits=None, lower_limits=None, video=False, hack_legend_to_only_show_pressure=False):
        x_axis = ['Al','Ti','Ca','Ni','Fe','Cr','Mg','Si','Na','O','C','N']
        whitedwarfabundancesplot = list()
        whitedwarferrorsplot = list()
        upper_limits_for_plot = list()  # a list of bools to track which of the abundances are actually upper/lower bounds
        lower_limits_for_plot = list()
        upper_abundances_for_plot = list()
        lower_abundances_for_plot = list()
        geo_model_for_solar_abundances = gi.GeologyModel()  # We're only using this to access the solar abundances so no arguments needed
        for el, abundance in all_wd_abundances.items():
            upper_abundances_for_plot.append(geo_model_for_solar_abundances.upper_XH_ratiod_to_solar[el])
            lower_abundances_for_plot.append(geo_model_for_solar_abundances.lower_XH_ratiod_to_solar[el])
            if np.isnan(abundance) or abundance == 0:
                if upper_limits is not None and upper_limits.get(el) is not None:
                    #whitedwarfabundancesplot.append(upper_limits[el] - geo_model_for_solar_abundances.solar_abundances[el])
                    whitedwarfabundancesplot.append(upper_limits[el] - (geo_model_for_solar_abundances.solar_ratiod_to_H[el] - geo_model_for_solar_abundances.solar_ratiod_to_H[wd_type]))
                    whitedwarferrorsplot.append(0.1) # This sets arrow size
                    upper_limits_for_plot.append(True)
                    lower_limits_for_plot.append(False)
                elif lower_limits is not None and lower_limits.get(el) is not None:
                    #whitedwarfabundancesplot.append(lower_limits[el] - geo_model_for_solar_abundances.solar_abundances[el])
                    whitedwarfabundancesplot.append(lower_limits[el] - (geo_model_for_solar_abundances.solar_ratiod_to_H[el] - geo_model_for_solar_abundances.solar_ratiod_to_H[wd_type]))
                    whitedwarferrorsplot.append(0.1) # This sets arrow size
                    upper_limits_for_plot.append(False)
                    lower_limits_for_plot.append(True)
                else:
                    whitedwarfabundancesplot.append(np.nan)
                    whitedwarferrorsplot.append(0)
                    upper_limits_for_plot.append(False)
                    lower_limits_for_plot.append(False)
            else:
                #whitedwarfabundancesplot.append(all_wd_abundances[el] - geo_model_for_solar_abundances.solar_abundances[el])
                whitedwarfabundancesplot.append(all_wd_abundances[el] - (geo_model_for_solar_abundances.solar_ratiod_to_H[el] - geo_model_for_solar_abundances.solar_ratiod_to_H[wd_type]))
                whitedwarferrorsplot.append(all_wd_abundance_errors[el])
                upper_limits_for_plot.append(False)
                lower_limits_for_plot.append(False)
        print(whitedwarfabundancesplot)
        print(whitedwarferrorsplot)
        
        obs_data_series_name = self.strip_system_suffix(system)
        obs_data_series_name += ' Observational Data'
        
        series_dict = {
            obs_data_series_name: {
                'type': dp.PlotType.scatter_2d_error,
                #'x_data': range(0, len(observations)),
                'x_data': x_axis,
                'y_data': whitedwarfabundancesplot,
                'y_error_data': whitedwarferrorsplot,
                'line_type': 'ko',
                'line_linewidth': 1,
                'capsize': 2,
                'capthick': 1,
                'legend': True,
                'upper_limits': upper_limits_for_plot,
                'lower_limits': lower_limits_for_plot
            }
        }
        colours = self.colour_list
        colour_index = 0
        min_fit_value = None
        for fit_key, fit_data in fit_dict.items():
            if min_fit_value is None:
                min_fit_value = min(fit_data)
            else:
                min_fit_value = min(min_fit_value, min(fit_data))  # Potential improvement: only check values for which there is an observation
            if video or hack_legend_to_only_show_pressure:
                if fit_key == 'LP run':
                    fit_key_to_use = 'Median model (? GPa)'  # This is a hack
                else:
                    print(fit_key)
                    print(fit_key[0])
                    print(fit_key[1])
                    print(fit_key[2])
                    len_pressure = len(str(fit_key[1]))
                    padding_space = ' '*(5 - len_pressure)
                    print(fit_key[1])
                    print(len_pressure)
                    print(padding_space)
                    print(str(fit_key[1]) + padding_space + ' GPa')
                    fit_key_to_use = 'Fit at ' + str(fit_key[1]) + padding_space + ' GPa'
                print(fit_key_to_use)
            else:
                #fit_key_to_use = fit_key + ' median'
                fit_key_to_use = fit_key
            series_dict[fit_key_to_use] = {
                'type': dp.PlotType.scatter_2d,
                #'x_data': range(0, len(observations)),
                'x_data': x_axis,
                'y_data': fit_data,
                'line_color': colours[colour_index],
                'line_style': '-',
                'line_linewidth': 2,
                'zorder': 2,
                'legend': True
            }
            if error_low_dict is not None and error_low_dict.get(fit_key) is not None and error_high_dict is not None and error_high_dict.get(fit_key) is not None:
                if min_fit_value is None:
                    min_fit_value = min(error_low_dict[fit_key])
                else:
                    min_fit_value = min(min_fit_value, min(error_low_dict[fit_key]))
                series_dict['1 Sigma Confidence Interval'] = {
                    'type': dp.PlotType.shade,
                    #'x_data': range(0, len(observations)),
                    'x_data': x_axis,
                    'y_data': error_low_dict[fit_key],
                    'y_shade_data': error_high_dict[fit_key],
                    'shade_colour': 'k',
                    'shade_alpha': 0.5,
                    'zorder': 2,
                    'legend': True
                }
            colour_index += 1
                 
        solar_offset = np.nanmedian(whitedwarfabundancesplot)
        if abs(solar_offset) < 1.5:
            solar_offset = 0
        else:
            solar_offset = int(solar_offset)
        
        add_stellar_spread_point = False
        stellar_spreads = [ua - la for ua, la in zip(upper_abundances_for_plot,lower_abundances_for_plot)]
        typical_stellar_spread = np.mean(stellar_spreads)
        stellar_spread_name = 'Typical Stellar Spread'
        if wd_type != ci.Element.H:
            stellar_spread_name += ' (relative to H)'
        if add_stellar_spread_point:
            series_dict[stellar_spread_name] = {
                'type': dp.PlotType.scatter_2d_error,
                'x_data': [len(x_axis) - 1.5],
                'y_data': [solar_offset],
                'y_error_data': [typical_stellar_spread/2],
                'line_type': 'ks',
                'line_linewidth': 1,
                'capsize': 2,
                'capthick': 1,
                #'line_markersize': 0,
                'legend': True
            }
                       
        if video:
            series_dict['Pressure Text'] = {
                'type': dp.PlotType.text,
                'x_pos': 6,
                'y_pos': -1,
                'text_string': str(fit_key[1]),
                'horizontalalignment': 'left',
                'fontsize': 20,
                'legend': False
            }
            series_dict['Pressure unit Text'] = {
                'type': dp.PlotType.text,
                'x_pos': 6.7,
                'y_pos': -1,
                'text_string': 'GPa',
                'horizontalalignment': 'left',
                'fontsize': 20,
                'legend': False
            }
        extension = '.png' if video else ['.pdf', '.png']
        #try:
        #    y_min = -2 if min_fit_value < -2 else None
        #    y_max = 1 if video else (2 if min_fit_value < -10 else None)  # For some reason y_max starts misbehaving if the min value is very small
        #except TypeError:
        #    y_min = None
        #    y_max = None
        plot_dict = {
            'complotraw': {
                'show': False,
                'filenames': [system + '_' + file_prefix + 'compositionraw' + e for e in extension],
                'fig_height': 5,
                'fig_width': 10,
                'dpi': 1000,
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 12,
                        #'title_text': r'Composition of ' + obs_data_series_name + ' pollutant',
                        'title_fontsize': 16,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Lithophiles \, \, \, \, \, \, \, \, \, \, \, \, Siderophiles \, \, \, \, \, \, \, Volatile Lithophiles \, \, \, \, Atmophiles',
                        'xlabel_fontsize': 16,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': '[X/' + str(wd_type) + ']',
                        'ylabel_fontsize': 20,
                        'ylabel_fontweight': 'bold',
                        'x_tick_labels': x_axis,
                        'x_tick_fontsize': 20,
                        'y_tick_fontsize': 17,
                        'ylabel_pad': 0,
                        #'x_min': 0,
                        #'x_max': 100,
                        #'y_min': -10,
                        #'y_max': -6.5,
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        if video:
            plotter.yield_output('<your_filepath_here>')
        else:
            plotter.yield_output(self.output_dir)
        return plot_dict

    def make_earth_composition_plot(self, observed_abundances, fitted_abundances_dict, layer, tag, observed_abundance_errors=None, fitted_abundance_errors_dict=None):
        #observations, errors, x_tick_labels = self.convert_to_relative_solar_abundance(observed_abundances, observed_abundance_errors)
        fits = dict()
        ele_list_dict = dict()
        for key, fa in fitted_abundances_dict.items():
            ele_list_dict[key] = list()
            for test_ele in observed_abundances.keys():  # Make a list to guarantee consistent ordering and skip NaNs. NB: This means would ignore all observations of 0 abundance
                try:
                    if observed_abundances[test_ele][0] != 0 and fa[test_ele][layer] is not None:
                        ele_list_dict[key].append(test_ele)
                except:
                    pass
            fits[key] = list()
            for ele in ele_list_dict[key]:
                try:
                    fits[key].append(fa[ele][layer]/observed_abundances[ele][0])  # discard errors and element labels. TODO: Verify here that the fitted abundances have the same elements present!
                except IndexError:
                    #This is horrible but basically this means there's no observed abundance
                    pass
        first_key = list(ele_list_dict.keys())[0]  # This is a bit hacky - I'm assuming all the fits have the same elements...
        series_dict = {
            'McDonough 2003': {
                'type': dp.PlotType.scatter_2d_error,
                'x_data': range(0, len(ele_list_dict[first_key])),
                'y_data': [1]*len(ele_list_dict[first_key]),
                #'y_error_data': errors,
                'line_type': 'k:',
                'line_linewidth': 1,
                'legend': True,
            }
        }
        colours = ['r', 'g', 'b', 'c', 'm', 'y'] # TODO: Expand on this and loop around if necessary
        fit_index = 0
        for series_name, series in fits.items():
            series_dict[series_name] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': range(0, len(series)),
                'y_data': series,
                'line_color': colours[fit_index % len(colours)],
                'line_marker': 'x',
                'line_style': '--',
                'line_linewidth': 1,
                'legend': True,
            }
            fit_index += 1
        plot_dict = {
            'testplot': {
                'show': False,
                'filenames': ['composition_rel_' + tag + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': r'Earth ' + str(layer) + ' composition',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Element',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'Number abundance relative to McDonough 2003',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'x_tick_labels': [str(e) for e in ele_list_dict[first_key]],
                        #'x_min': 0,
                        #'x_max': 100,
                        'y_min': 0,
                        'y_max': 1.6 if str(layer) == 'mantle' else None,  # Another hack
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        
    def make_Ds_comparison_plot(self, observed_Ds_dict, predicted_Ds_dict, tag, run_name, observed_Ds_errors_dict=None, predicted_Ds_errors_dict=None):
        elements = list()
        fit = list()
        actual = list()
        errors = list()
        for el, observation in observed_Ds_dict.items():
            elements.append(el)
            actual.append(np.log10(observation))
            fit.append(np.log10(predicted_Ds_dict[el]))
            errors.append(observed_Ds_errors_dict[el])  # For some reason, I made the errors in log space but the values in linear space
        series_dict = {
            run_name: {
                'type': dp.PlotType.scatter_2d_error,
                'x_data': range(0, len(elements)),
                'y_data': actual,
                'y_error_data': errors,
                'line_marker': 'x',
                'line_markersize': 3,
                'line_color': self.colour_list[0],
                'line_style': 'None',
                #'line_type': 'k:',
                'line_linewidth': 1,
                'legend': True,
            }
        }
        series_dict['Model'] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': range(0, len(elements)),
            'y_data': fit,
            'line_color': self.colour_list[1],
            'line_marker': 'x',
            'line_style': 'None',
            'line_linewidth': 1,
            'legend': True,
        }
        #series_dict['RNText'] = {
        #    'type': dp.PlotType.text,
        #    'y_pos': 0,
        #    'x_pos': 5,
        #    'text_string': run_name,
        #    'fontsize': 12,
        #    'shade_colour': 'w',
        #    'shade_alpha': 0.5,
        #    'legend': False
        #}
        #series_dict['KEText'] = {
        #    'type': dp.PlotType.text,
        #    'y_pos': 0,
        #    'x_pos': 0,
        #    'text_string': 'Key Elements',
        #    'fontsize': 10,
        #    'shade_colour': 'k',
        #    'shade_alpha': 0.5,
        #    'legend': False
        #}
        plot_dict = {
            'testplot': {
                'show': False,
                'filenames': ['partition_coefficients_comparison_' + tag + '_' + run_name + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Element',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'log(Partition coefficient)',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'x_tick_labels': [str(e) for e in elements],
                        #'x_min': 0,
                        #'x_max': 100,
                        'y_min': -3,
                        'y_max': 3,
                        #'y_scale': 'log',
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        
    def make_Ds_comparison_multipanel_plot(self, observed_Ds_dict_list, predicted_Ds_dict_list, tag, run_name_list, observed_Ds_errors_dict_list=None, predicted_Ds_errors_dict_list=None):
        # For now assuming we only want 2 panels
        elements_lists = list()
        fit_lists = list()
        actual_lists = list()
        errors_lists = list()
        for i, observed_Ds_dict in enumerate(observed_Ds_dict_list):
            elements_lists.append(list())
            fit_lists.append(list())
            actual_lists.append(list())
            errors_lists.append(list())
            for el, observation in observed_Ds_dict.items():
                elements_lists[i].append(el)
                actual_lists[i].append(np.log10(observation))
                fit_lists[i].append(np.log10(predicted_Ds_dict_list[i][el]))
                errors_lists[i].append(observed_Ds_errors_dict_list[i][el])  # For some reason, I made the errors in log space but the values in linear space
        series_dict_panel1 = {
            'Model': {
                'type': dp.PlotType.scatter_2d,
                'x_data': range(0, len(elements_lists[0])),
                'y_data': fit_lists[0],
                'line_color': self.colour_list[1],
                'line_marker': 'None',
                'line_markersize': 5,
                'line_style': '-',
                'line_linewidth': 2,
                'legend': True,
            },
            run_name_list[0]: {
                'type': dp.PlotType.scatter_2d_error,
                'x_data': range(0, len(elements_lists[0])),
                'y_data': actual_lists[0],
                'y_error_data': errors_lists[0],
                'line_marker': 'x',
                'line_markersize': 7,
                'line_color': self.colour_list[0],
                'line_style': 'None',
                #'line_type': 'k:',
                'line_linewidth': 1,
                'legend': True,
            },
            'KeyElShade': {
                'type': dp.PlotType.shade,
                'x_data': [-1, 4.5],
                'y_data': [-3, -3],
                'y_shade_data': [3, 3],
                'shade_colour': 'g',
                'shade_alpha': 0.15,
                'zorder': 2,
                'legend': False
            },
            'KeyElText': {
                'type': dp.PlotType.text,
                'text_string': 'Key Elements',
                'x_pos': 0.1,
                'y_pos': 2.55,
                'fontsize': 24,
                'legend': False
            }
        }
        series_dict_panel2 = {
            'Model': {
                'type': dp.PlotType.scatter_2d,
                'x_data': range(0, len(elements_lists[1])),
                'y_data': fit_lists[1],
                'line_color': self.colour_list[1],
                'line_marker': 'None',
                'line_markersize': 5,
                'line_style': '-',
                'line_linewidth': 2,
                'legend': True,
            },
            run_name_list[1]: {
                'type': dp.PlotType.scatter_2d_error,
                'x_data': range(0, len(elements_lists[1])),
                'y_data': actual_lists[1],
                'y_error_data': errors_lists[1],
                'line_marker': 'o',
                'line_markersize': 5,
                'line_color': self.colour_list[0],
                'line_style': 'None',
                #'line_type': 'k:',
                'line_linewidth': 2,
                'legend': True,
            },
            'KeyElShade': {
                'type': dp.PlotType.shade,
                'x_data': [-1, 3.5],
                'y_data': [-3, -3],
                'y_shade_data': [3, 3],
                'shade_colour': 'g',
                'shade_alpha': 0.15,
                'zorder': 2,
                'legend': False
            },
            'KeyElText': {
                'type': dp.PlotType.text,
                'text_string': 'Key Elements',
                'x_pos': -0.25,
                'y_pos': -2.75,
                'fontsize': 24,
                'legend': False
            }
        }
        plot_dict = {
            'testplot': {
                'show': False,
                'filenames': ['partition_coefficients_comparison_' + tag + '_' + '_'.join(run_name_list) + '.pdf'],
                'fig_width': 15,
                'fig_height': 6,
                'gridspec_y_x': (1, 2),
                'gridspec_width_ratios': [len(elements_lists[0]), len(elements_lists[1])],
                'gridspec_wspace': 0.2,
                'subplots': {
                    'subplot1': {
                        'subplot_region': 121,
                        'legend': True,
                        #'legend_loc': 'best',
                        'legend_coord_x': 0.48,
                        'legend_coord_y': 0.25,
                        'legend_text_size': 22,
                        'title_fontsize': 28,
                        'title_fontweight': 'bold',
                        'title_text': run_name_list[0],
                        #'xlabel_text': 'Element',
                        'xlabel_fontsize': 16,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': r'$\textrm{log}_{10}(D)$',
                        'ylabel_fontsize': 30,
                        'ylabel_fontweight': 'bold',
                        'x_tick_labels': [str(e) for e in elements_lists[0]],
                        'x_tick_fontsize': 24,
                        'y_tick_fontsize': 24,
                        'x_min': -0.5,
                        #'x_max': 100,
                        'y_min': -3,
                        'y_max': 3,
                        #'y_scale': 'log',
                        'font': 'STIXGeneral',
                        'series': series_dict_panel1,
                        'gridspec_index': 0
                    },
                    'subplot2': {
                        'subplot_region': 122,
                        'legend': True,
                        'legend_loc': 'best',
                        #'legend_coord_x': 0.36,
                        #'legend_coord_y': 0.7,
                        'legend_text_size': 18,
                        'title_fontsize': 28,
                        'title_fontweight': 'bold',
                        'title_text': run_name_list[1],
                        #'xlabel_text': 'Element',
                        'xlabel_fontsize': 16,
                        'xlabel_fontweight': 'bold',
                        #'ylabel_text': 'log(Partition coefficient)',
                        'ylabel_fontsize': 16,
                        'ylabel_fontweight': 'bold',
                        'x_tick_labels': [str(e) for e in elements_lists[1]],
                        'x_tick_fontsize': 24,
                        'y_tick_fontsize': 24,
                        'x_min': -0.5,
                        #'x_max': 100,
                        #'y_min': -3,
                        #'y_max': 3,
                        #'y_scale': 'log',
                        'font': 'STIXGeneral',
                        'series': series_dict_panel2,
                        'sharey_subplot': 'subplot1',
                        'gridspec_index': 1
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
    
    def make_planet_formation_plot(self, abundances_dict, variable, variable_vals, layer_text, earth_observations, tag, ratio_against=None):
        colours = [
            '#bbbaaa',
            '#000000',
            '#FF0FD5',
            '#191970',
            '#FFA07A',
            '#9F1493',
            '#CD853F',
            '#FFD700',
            '#87CEFA',
            '#32CD32',
            '#BC8F8F',
            '#808000',
            '#FF7F50',
            '#7FFFD4',
            '#5F9EA0',
            '#FFA500',
            '#C0C0C0',
            '#800000',
            '#BA55D3',
            '#66CDAA',
            '#B0E0E6',
            '#8A2BE2',
            '#98FB98',
            '#F0E68C',
            '#006400',
            '#8B0000',
            '#DAA520',
            '#D8BFD8'
        ]
        el_index = 0
        series_dict = dict()
        for element, series in abundances_dict.items():
            line_style = '-'
            series_dict[str(element)] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': variable_vals,
                'y_data': series,
                'line_color': colours[el_index % len(colours)],
                'line_style': line_style,
                'line_linewidth': 1,
                'legend': True,
            }
            el_index += 1
        if earth_observations is not None:
            el_index = 0
            for element, value in earth_observations.items():
                if value is None or len(value) == 0:
                    el_index += 1
                    continue
                line_style = 'x'
                series_dict[str(element) + ' (Earth)'] = {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': [54 if variable == 'Pressure' else -2],
                    'y_data': value,
                    'line_color': colours[el_index % len(colours)],
                    'line_marker': line_style,
                    'line_linewidth': 1,
                    'legend': True
                }
                el_index += 1
        yaxis_extension = '' if ratio_against is None else ' relative to ' + ratio_against
        filename_extension = '' if ratio_against is None else '_rel_' + ratio_against
        plot_dict = {
            'planetplot': {
                'show': False,
                'filenames': ['planet_plot_' + layer_text + '_' + tag + '_' + variable.replace(' ', '') + filename_extension + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 4,
                        'title_text': r'Composition of Earth ' + layer_text + ' with varying ' + variable,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': variable + ' /' + self.unit_dict[variable],
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'Number abundance' + yaxis_extension,
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'x_min': min(variable_vals),
                        'x_max': max(variable_vals),
                        #'y_min': 0,
                        #'y_max': 100,
                        'y_scale': 'log',
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        focus_on_crni = True
        if focus_on_crni:
            plot_dict['planetplotcrni'] = {
                'show': False,
                'filenames': ['planet_plot_crni_' + layer_text + '_' + tag + '_' + variable.replace(' ', '') + filename_extension + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 4,
                        'title_text': r'Composition of Earth ' + layer_text + ' with varying ' + variable,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': variable + ' /' + self.unit_dict[variable],
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'Number abundance' + yaxis_extension,
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'x_min': min(variable_vals),
                        'x_max': max(variable_vals),
                        'y_min': 0,
                        #'y_max': 100,
                        'font': 'STIXGeneral',
                        'series': {'Cr': series_dict['Cr'], 'Ni': series_dict['Ni']}
                    }
                }
            }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
    
    def get_planet_formation_stacked_subplot_dict(self, abundances_dict, variable, variable_vals, file_prefix, earth_observations, colour_dict):
        first_ele_name = list(abundances_dict.keys())[0]  # Doesn't matter which - just need a key to find out how many data points we have on the next line
        series_dict = cn.OrderedDict()
        prev_series = [0]*len(abundances_dict[first_ele_name])
        for element, series in abundances_dict.items():
            line_style = '-'
            series_dict[str(element)] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': variable_vals,
                'y_data': series,
                'line_color': colour_dict[element],
                'line_style': line_style,
                'line_linewidth': 1,
                'legend': True,
            }
            where = [s is not None for s in series]  # No need to check prev_series as well: Nones will appear in the same places in series and prev_series
            series_dict[str(element) + '_shade'] = {
                'type': dp.PlotType.shade,
                'x_data': variable_vals,
                'y_data': [e if e is not None else 0 for e in prev_series],
                'y_shade_data': [e if e is not None else 0 for e in series],
                'shade_colour': colour_dict[element],
                'line_linewidth': 1,
                'shade_alpha': 0.4,
                'shade_where': where
            }
            prev_series = series
        reversed_series_dict = cn.OrderedDict()
        for series_name in reversed(series_dict.keys()):
            reversed_series_dict[series_name] = series_dict[series_name]
        #self.x_data, self.y_data, self.y_shade_data, facecolor=self.shade_colour, alpha=self.shade_alpha)
        #el_index = 0
        #for element, value in earth_observations.items():
        #    if value is None or len(value) == 0:
        #        el_index += 1
        #        continue
        #    line_style = 'x'
        #    series_dict[str(element) + ' (Earth)'] = {
        #        'type': dp.PlotType.scatter_2d,
        #        'x_data': [54 if variable == 'Pressure' else -2],
        #        'y_data': value,
        #        'line_color': colours[el_index % len(colours)],
        #        'line_marker': line_style,
        #        'line_linewidth': 1,
        #        'legend': True
        #    }
        #    el_index += 1
        unit_dict = {
            'Pressure': ' /GPa',
            'fO2':  ' (' + r'$\Delta$' + 'IW' + ')'
        }
        variable_name_to_use = r'$f_{\textrm{O}_{2}}$' if variable == 'fO2' else variable
        toret = {
            'subplot_region': 111,
            'legend': True,
            'legend_loc': 'upper left',
            'legend_text_size': 9,
            'title_text': r'Composition of Earth ' + file_prefix + ' with varying ' + variable,
            'title_fontsize': 12,
            'title_fontweight': 'bold',
            'xlabel_text': variable_name_to_use + unit_dict[variable],
            'xlabel_fontsize': 16,
            'xlabel_fontweight': 'bold',
            'ylabel_text': 'Cumulative abundance',
            'ylabel_fontsize': 16,
            'ylabel_fontweight': 'bold',
            'x_min': min(variable_vals),
            'x_max': max(variable_vals),
            'y_min': 0.001,
            'y_max': 1,
            'y_scale': 'log',
            'font': 'STIXGeneral',
            'series': reversed_series_dict
        }
        return toret
    
    def make_planet_formation_stacked_plot(self, abundances_dict, variable, variable_vals, file_prefix, earth_observations, colour_dict, run_tag):
        subplot_dict = self.get_planet_formation_stacked_subplot_dict(abundances_dict, variable, variable_vals, file_prefix, earth_observations, colour_dict)
        plot_dict = {
            'planetplot_stacked': {
                'show': False,
                'filenames': ['planet_plot_stacked_' + run_tag + '_' + file_prefix + '_' + variable + '.pdf', 'planet_plot_stacked_' + run_tag + '_' + file_prefix + '_' + variable + '.png'],
                'subplots': {
                    'subplot1': subplot_dict
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
    
    def get_planet_formation_combined_stacked_subplot_dict(self, abundances_dict, variable, variable_vals, file_prefix, earth_observations, colour_dict):
        first_ele_name = list(abundances_dict[gi.Layer.core].keys())[0]  # Doesn't matter which - just need a key to find out how many data points we have on the next line
        series_dict = cn.OrderedDict()
        prev_series = [0]*len(abundances_dict[gi.Layer.core][first_ele_name])
        for layer in [gi.Layer.core, gi.Layer.mantle]:
            for element, series in abundances_dict[layer].items():
                line_style = '-'
                el_name_to_use = str(element) if element != ci.Element.Placeholder else 'Other'
                series_name = el_name_to_use + ' (' + str(layer) + ')'
                series_dict[series_name] = {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': variable_vals,
                    'y_data': series,
                    'line_color': colour_dict[element],
                    'line_style': line_style,
                    'line_linewidth': 1,
                    #'legend': element not in elements_already_seen,
                }
                where = [s is not None for s in series]  # No need to check prev_series as well: Nones will appear in the same places in series and prev_series
                shade_series_name = el_name_to_use if layer == gi.Layer.mantle else el_name_to_use + '_shade_' + str(layer)
                series_dict[shade_series_name] = {
                    'type': dp.PlotType.shade,
                    'x_data': variable_vals,
                    'y_data': [e if e is not None else 0 for e in prev_series],
                    'y_shade_data': [e if e is not None else 0 for e in series],
                    'shade_colour': colour_dict[element],
                    'line_linewidth': 1,
                    'shade_alpha': 0.4,
                    'shade_hatch': 'x' if layer == gi.Layer.core else None,
                    'shade_where': where,
                    'legend': layer == gi.Layer.mantle
                }
                prev_series = series
        ordered_series_dict = cn.OrderedDict()
        desired_order = [ci.Element.O, ci.Element.Si, ci.Element.Mg, ci.Element.Fe, ci.Element.Al, ci.Element.Ca, ci.Element.Ni, ci.Element.Cr, ci.Element.C, ci.Element.Placeholder]
        for element in desired_order:
            el_name_to_use = str(element) if element != ci.Element.Placeholder else 'Other'
            try:
                ordered_series_dict[el_name_to_use] = series_dict[el_name_to_use]
            except KeyError:
                pass
        for series_key, series_val in series_dict.items():
            if series_key not in ordered_series_dict:
                ordered_series_dict[series_key] = series_val
        #self.x_data, self.y_data, self.y_shade_data, facecolor=self.shade_colour, alpha=self.shade_alpha)
        #el_index = 0
        #for element, value in earth_observations.items():
        #    if value is None or len(value) == 0:
        #        el_index += 1
        #        continue
        #    line_style = 'x'
        #    series_dict[str(element) + ' (Earth)'] = {
        #        'type': dp.PlotType.scatter_2d,
        #        'x_data': [54 if variable == 'Pressure' else -2],
        #        'y_data': value,
        #        'line_color': colours[el_index % len(colours)],
        #        'line_marker': line_style,
        #        'line_linewidth': 1,
        #        'legend': True
        #    }
        #    el_index += 1
        if variable == 'Pressure':
            ordered_series_dict['Core text'] = {
                'type': dp.PlotType.text,
                'x_pos': 32.5,
                'y_pos': 0.06,
                'text_string': 'Core',
                'fontsize': 28,
                'shade_colour': 'w',
                'shade_alpha': 0.75,
                'legend': False
            }
            ordered_series_dict['Mantle text'] = {
                'type': dp.PlotType.text,
                'x_pos': 30,
                'y_pos': 0.6,
                'text_string': 'Mantle',
                'fontsize': 28,
                'shade_colour': 'w',
                'shade_alpha': 0.75,
                'legend': False
            }
        unit_dict = {
            'Pressure': ' /GPa',
            'fO2':  ' (' + r'$\Delta$' + 'IW' + ')'
        }
        variable_name_to_use = r'$f_{\textrm{O}_{2}}$' if variable == 'fO2' else variable
        toret = {
            'subplot_region': 111,
            'legend': True,
            'legend_loc': 'upper left',
            'legend_text_size': 14,
            'title_text': r'Composition of Earth mantle and core with varying ' + variable,
            'title_fontsize': 12,
            'title_fontweight': 'bold',
            'xlabel_text': variable_name_to_use + unit_dict[variable],
            'xlabel_fontsize': 18,
            'xlabel_fontweight': 'bold',
            'ylabel_text': 'Cumulative abundance',
            'ylabel_fontsize': 18,
            'ylabel_fontweight': 'bold',
            'x_tick_fontsize': 14,
            'y_tick_fontsize': 14,
            'x_min': min(variable_vals),
            'x_max': max(variable_vals),
            'y_min': 0,
            'y_max': 1,
            #'y_scale': 'log',
            'font': 'Verdana',
            'series': ordered_series_dict
        }
        return toret
    
    def make_planet_formation_combined_stacked_plot(self, abundances_dict, variable, variable_vals, file_prefix, earth_observations, colour_dict, run_tag):
        subplot_dict = self.get_planet_formation_combined_stacked_subplot_dict(abundances_dict, variable, variable_vals, file_prefix, earth_observations, colour_dict)
        plot_dict = {
            'planetplot_stacked': {
                'show': False,
                'filenames': ['planet_plot_stacked_' + run_tag + '_' + file_prefix + '_' + variable + '.pdf', 'planet_plot_stacked_' + run_tag + '_' + file_prefix + '_' + variable + '.png'],
                'subplots': {
                    'subplot1': subplot_dict
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        
    def make_planet_formation_multipanel_stacked_plot(self, stacked_variable_p_vals, bulk_stacked_variable_p_vals,
        stacked_variable_fO2_vals, bulk_stacked_variable_fO2_vals, pressure_vals, fO2_vals, file_prefix,
        earth_observations, colour_dict, run_tag):
        bulk_p_subplot_dict = self.get_planet_formation_combined_stacked_subplot_dict(bulk_stacked_variable_p_vals, 'Pressure', pressure_vals, 'bulk', earth_observations, colour_dict)
        bulk_fO2_subplot_dict = self.get_planet_formation_combined_stacked_subplot_dict(bulk_stacked_variable_fO2_vals, 'fO2', fO2_vals, 'bulk', earth_observations, colour_dict)
        mantle_p_subplot_dict = self.get_planet_formation_stacked_subplot_dict(stacked_variable_p_vals[gi.Layer.mantle], 'Pressure', pressure_vals, 'mantle', earth_observations, colour_dict)
        mantle_fO2_subplot_dict = self.get_planet_formation_stacked_subplot_dict(stacked_variable_fO2_vals[gi.Layer.mantle], 'fO2', fO2_vals, 'mantle', earth_observations, colour_dict)
        core_p_subplot_dict = self.get_planet_formation_stacked_subplot_dict(stacked_variable_p_vals[gi.Layer.core], 'Pressure', pressure_vals, 'core', earth_observations, colour_dict)
        core_fO2_subplot_dict = self.get_planet_formation_stacked_subplot_dict(stacked_variable_fO2_vals[gi.Layer.core], 'fO2', fO2_vals, 'core', earth_observations, colour_dict)
        plot_dict = {
            'planetplot_multipanel_stacked': {
                'show': False,
                'filenames': ['planet_plot_stacked_' + run_tag + '_' + file_prefix + '.pdf'],
                'gridspec_y_x': (3, 2),
                'gridspec_wspace': 0.07,
                'gridspec_hspace': 0.07,
                'fig_height': 16,
                'fig_width': 13,
                'subplots': {
                    'bulk_p_subplot': bulk_p_subplot_dict,
                    'bulk_fO2_subplot': bulk_fO2_subplot_dict,
                    'mantle_p_subplot': mantle_p_subplot_dict,
                    'mantle_fO2_subplot': mantle_fO2_subplot_dict,
                    'core_p_subplot': core_p_subplot_dict,
                    'core_fO2_subplot': core_fO2_subplot_dict
                }
            }
        }
        
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_p_subplot']['title_text'] = None
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_p_subplot']['gridspec_index'] = 0
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_p_subplot']['sharex_subplot'] = 'mantle_p_subplot'
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_p_subplot']['ylabel_text'] += ' (bulk)'
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_p_subplot']['ylabel_fontsize'] = 22
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_p_subplot']['y_tick_fontsize'] = 16
        #plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_p_subplot']['y_hide_ticks'] = [0]
        
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_fO2_subplot']['title_text'] = None
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_fO2_subplot']['gridspec_index'] = 1
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_fO2_subplot']['sharey_subplot'] = 'bulk_p_subplot'
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_fO2_subplot']['sharex_subplot'] = 'mantle_fO2_subplot'
        plot_dict['planetplot_multipanel_stacked']['subplots']['bulk_fO2_subplot']['legend'] = False
        
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_p_subplot']['title_text'] = None
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_p_subplot']['y_min'] = 0.0003
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_p_subplot']['gridspec_index'] = 2
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_p_subplot']['sharex_subplot'] = 'core_p_subplot'
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_p_subplot']['ylabel_text'] += ' (mantle)'
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_p_subplot']['ylabel_fontsize'] = 22
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_p_subplot']['y_tick_fontsize'] = 16
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_p_subplot']['legend'] = False
        
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_fO2_subplot']['title_text'] = None
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_fO2_subplot']['y_min'] = 0.0003
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_fO2_subplot']['gridspec_index'] = 3
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_fO2_subplot']['sharey_subplot'] = 'mantle_p_subplot'
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_fO2_subplot']['sharex_subplot'] = 'core_fO2_subplot'
        plot_dict['planetplot_multipanel_stacked']['subplots']['mantle_fO2_subplot']['legend'] = False
        
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['title_text'] = None
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['gridspec_index'] = 4
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['ylabel_text'] += ' (core)'
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['ylabel_fontsize'] = 22
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['y_min'] = 0.003
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['legend'] = False
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['x_tick_fontsize'] = 16
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['y_tick_fontsize'] = 16
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_p_subplot']['xlabel_fontsize'] = 24
        
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_fO2_subplot']['title_text'] = None
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_fO2_subplot']['gridspec_index'] = 5
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_fO2_subplot']['sharey_subplot'] = 'core_p_subplot'
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_fO2_subplot']['y_min'] = 0.003
        #plot_dict['planetplot_multipanel_stacked']['subplots']['core_fO2_subplot']['x_hide_ticks'] = [0]
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_fO2_subplot']['legend'] = False
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_fO2_subplot']['x_tick_fontsize'] = 16
        plot_dict['planetplot_multipanel_stacked']['subplots']['core_fO2_subplot']['xlabel_fontsize'] = 24
        
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        
    def make_planet_cnf_plot(self, cnf_list, variable, variable_vals, observations_dict):
        unit_dict = {
            'Pressure': 'GPa',
            'fO2': 'ΔIW'
        }
        series_dict = {
            'cnf_series': {
                'type': dp.PlotType.scatter_2d,
                'x_data': variable_vals,
                'y_data': cnf_list,
                'line_type': 'k-',
                'line_linewidth': 1
            }
        }
        for key, value in observations_dict.items():
            line_style = 'x'
            series_dict[key] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': [54 if variable == 'Pressure' else -2],
                'y_data': [value],
                'line_type': 'k' + line_style,
                'line_linewidth': 1,
                'legend': True,
            }
        plot_dict = {
            'cnfplot': {
                'show': False,
                'filenames': ['planet_cnf_' + variable + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': r'Core number fraction with varying ' + variable,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': variable + ' /' + unit_dict[variable],
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'Number abundance',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'x_min': min(variable_vals),
                        'x_max': max(variable_vals),
                        #'y_min': 0,
                        #'y_max': 100,
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
    
    def make_convergence_plot(self, Ds_dict, run_name, colour_dict):
        # Ds_dict[element] = [D_0, D_1 ... D_iterations]
        series_dict = dict()
        i = 0
        for element, D_vals in Ds_dict.items():
            try:
                line_colour = colour_dict[element]
            except:
                line_colour = self.colour_list[i]
                i += 1
            series_dict[str(element)] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': list(range(0, len(D_vals))),
                'y_data': D_vals,
                'line_color': line_colour,
                'line_style': '-',
                'line_linewidth': 1,
                'legend': True
            }
        plot_dict = {
            'cnfplot': {
                'show': False,
                'filenames': ['convergence_' + run_name + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': 'Convergence of Partition Coefficients, D, for run $\mathrm{' + run_name.replace('_', '\_') + '}$',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Iteration',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'D',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        #'x_min': min(variable_vals),
                        #'x_max': max(variable_vals),
                        #'y_min': 0,
                        #'y_max': 100,
                        'y_scale': 'log',
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)

    def plot_log_elel_ratio(self, ratio_vals, x_vals, variable, system_name, element, ref_element, core_number_fraction,
                           observed_val, observed_error, tag, fixed_val, colour_dict, ignore_lack_of_observations=False, femg_vals=None, femg_obs=None, femg_obs_err=None):
        if observed_val is None and not ignore_lack_of_observations:
            print('No observation for ' + str(element) + '/' + str(ref_element) + ' in ' + system_name + ', skipping')
            return
        other_variable = 'Pressure' if variable == 'fO2' else 'fO2'
        if observed_val is not None:
            series_dict = {
                'log(' + str(element) + '/' + str(ref_element) + ')': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': ratio_vals,
                    'line_color': 'g', #colour_dict[element],
                    'line_style': '-',
                    'line_linewidth': 1,
                    'legend': True
                },
                'Observed log(' + str(element) + '/' + str(ref_element) + ')': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': [observed_val for x in x_vals],
                    'line_color': 'g',
                    'line_style': ':',
                    'line_linewidth': 1,
                    'legend': True
                },
                'observed + sigma': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': [observed_val + observed_error for x in x_vals],
                    'line_color': 'g',
                    'line_style': '--',
                    'line_linewidth': 1,
                    'legend': False
                },
                'observed - sigma': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': [observed_val - observed_error for x in x_vals],
                    'line_color': 'g',
                    'line_style': '--',
                    'line_linewidth': 1,
                    'legend': False
                },
                'observed + sigma shade': {
                    'type': dp.PlotType.shade,
                    'x_data': x_vals,
                    'y_data': [observed_val for x in x_vals],
                    'y_shade_data': [observed_val + observed_error for x in x_vals],
                    'shade_colour': 'g',
                    'shade_alpha': 0.3,
                    'legend': False
                },
                'observed - sigma shade': {
                    'type': dp.PlotType.shade,
                    'x_data': x_vals,
                    'y_data': [observed_val for x in x_vals],
                    'y_shade_data': [observed_val - observed_error for x in x_vals],
                    'shade_colour': 'g',
                    'shade_alpha': 0.3,
                    'legend': False
                }
            }
        else:
            series_dict = {
                'log(' + str(element) + '/' + str(ref_element) + ')': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': ratio_vals,
                    'line_color': 'g', #colour_dict[element],
                    'line_style': '-',
                    'line_linewidth': 1,
                    'legend': True
                }
            }
        if femg_obs is not None:
            series_dict_femg = {
                'log(Fe/Mg)': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': femg_vals,
                    'line_color': 'b',
                    'line_style': '-',
                    'line_linewidth': 1,
                    'legend': True
                },
                'Observed log(Fe/Mg)': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': [femg_obs for x in x_vals],
                    'line_color': 'b',
                    'line_style': ':',
                    'line_linewidth': 1,
                    'legend': True
                },
                'observed + sigma': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': [femg_obs + femg_obs_err for x in x_vals],
                    'line_color': 'b',
                    'line_style': '--',
                    'line_linewidth': 1,
                    'legend': False
                },
                'observed - sigma': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': [femg_obs - femg_obs_err for x in x_vals],
                    'line_color': 'b',
                    'line_style': '--',
                    'line_linewidth': 1,
                    'legend': False
                },
                'observed + sigma shade': {
                    'type': dp.PlotType.shade,
                    'x_data': x_vals,
                    'y_data': [femg_obs for x in x_vals],
                    'y_shade_data': [femg_obs + femg_obs_err for x in x_vals],
                    'shade_colour': 'b',
                    'shade_alpha': 0.3,
                    'legend': False
                },
                'observed - sigma shade': {
                    'type': dp.PlotType.shade,
                    'x_data': x_vals,
                    'y_data': [femg_obs for x in x_vals],
                    'y_shade_data': [femg_obs - femg_obs_err for x in x_vals],
                    'shade_colour': 'b',
                    'shade_alpha': 0.3,
                    'legend': False
                }
            }
        else:
            series_dict_femg = {
                'log(Fe/Mg)': {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals,
                    'y_data': femg_vals,
                    'line_color': 'b', #colour_dict[element],
                    'line_style': ':',
                    'line_linewidth': 1,
                    'legend': True
                }
            }
        title_var = variable if variable is not 'fcf' else 'Fragment Core Fraction'
        title = 'log(' + str(element) + '/' + str(ref_element) + ')'
        if femg_vals is not None:
            title += ' and log(Fe/Mg)'
        title += ' against ' + title_var + ' for ' + system_name
        if fixed_val is not None:
            if other_variable == 'fO2':
                title += ' (' + other_variable + ' = IW ' + str(fixed_val) + ')'
            else:
                title += ' (' + other_variable + ' = ' + str(fixed_val) + self.unit_dict[other_variable] + ')'
        plot_dict = {
            'elplot': {
                'show': False,
                'filenames': ['elel' + '_' + system_name + '_' + str(element) + '_' + str(ref_element) + '_' + variable + '_' + tag + '.pdf'],
                'subplots': {
                    'elplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': title,
                        'title_fontsize': 11,
                        'title_fontweight': 'bold',
                        'xlabel_text': variable + ' /' + self.unit_dict[variable] if variable is not 'fcf' else 'Fragment Core Fraction',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'log(' + str(element) + '/' + str(ref_element) + ')',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'x_min': min(x_vals),
                        'x_max': max(x_vals),
                        #'y_min': 0,
                        #'y_max': 100,
                        #'y_scale': 'log',
                        'font': 'STIXGeneral',
                        'series': series_dict
                    },
                    'elplot2': {
                        'subplot_region': 111,
                        'legend': True,
                        #'legend_loc': 'best',
                        #'legend_text_size': 8,
                        #'title_text': title,
                        #'title_fontsize': 11,
                        #'title_fontweight': 'bold',
                        #'xlabel_text': variable + ' /' + self.unit_dict[variable],
                        #'xlabel_fontsize': 10,
                        #'xlabel_fontweight': 'bold',
                        'ylabel_text': 'log(Fe/Mg)',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        #'x_min': min(x_vals),
                        #'x_max': max(x_vals),
                        #'y_min': 0,
                        #'y_max': 100,
                        #'y_scale': 'log',
                        'font': 'STIXGeneral',
                        'series': series_dict_femg,
                        'twin_subplot': 'elplot1'
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        
    def plot_3d_log_elel_ratio(self, ratio_vals, p_range, f_range, system_name, element, ref_element, fixed_var_val, observed_val, observed_error, tag, x_var, y_var, fixed_var, ignore_lack_of_observation=False, guidelines=list()):
        ref_element_str = str(ref_element) if ref_element != ci.Element.Placeholder else 'Hx'
        if observed_val is None and not ignore_lack_of_observation:
            print('No observation for ' + str(element) + '/' + ref_element_str + ' in ' + system_name + ', skipping')
            return
        series_dict = {
            '3delplot': {
                'type': dp.PlotType.contour_scatter,
                'x_data': p_range,
                'y_data': f_range,
                'z_data': ratio_vals,
                'fill': True,
                'cbar_label': 'log(' + str(element) + '/' + ref_element_str + ')',
                'legend': False,
                'cbar_labelrotation': 90,
                'cbar_labelfontsize': 16,
                'cbar_tickfontsize': 16
            }
        }
        if not ignore_lack_of_observation:
            series_dict['3delplot']['levels'] = [observed_val - (3*observed_error), observed_val - (2*observed_error), observed_val - observed_error, observed_val + observed_error, observed_val + (2*observed_error), observed_val + (3*observed_error)]
            series_dict['3delplot']['colours'] = ('r', 'b', 'g', 'b', 'r')
            series_dict['3delplot']['colour_below_min'] = 'k'
            series_dict['3delplot']['colour_above_max'] = 'k'
        i = 0
        for guideline in guidelines:
            name = 'log(' + str(guideline['num_el']) + '/' +  str(guideline['denom_el']) +')'
            if i > 0:
                name += ' ' + str(i+1)
            series_dict[name] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': guideline['x_data'],
                'y_data': guideline['y_data'],
                'line_color': 'r', #colour_dict[element],
                'line_style': ':',
                'line_marker': 'o',
                'line_linewidth': 1,
                'legend': True
            }
            j = 0
            assert len(guideline['x_data']) == len(guideline['y_data']) == len(guideline['z_data'])
            while j < len(guideline['x_data']):
                #try:
                #    next_x = guideline['x_data'][j+1]
                #    next_y = guideline['y_data'][j+1]
                #    delta_x = next_x - guideline['x_data'][j]
                #    delta_y = next_y - guideline['y_data'][j]
                #except IndexError:
                #    prev_x = guideline['x_data'][j-1]
                #    prev_y = guideline['y_data'][j-1]
                #    delta_x = -(prev_x - guideline['x_data'][j])
                #    delta_y = -(prev_y - guideline['y_data'][j])
                #x_offset = 5*delta_y
                #y_offset = -0.002*delta_x
                x_offset = 0.8
                y_offset = -0.04 if j != 0 else 0.02
                series_dict['guideline' + str(i) + 'text' + str(j)] = {
                    'type': dp.PlotType.text,
                    'x_pos': guideline['x_data'][j] + x_offset,
                    'y_pos': guideline['y_data'][j] + y_offset,
                    'text_string': guideline['z_data'][j],
                    'fontsize': 16,
                    'legend': False
                }
                j += 1
            i += 1
        x_var_symbol = self.symbol_dict[x_var]
        y_var_symbol = self.symbol_dict[y_var]
        x_var_unit = self.unit_dict[x_var]
        y_var_unit = self.unit_dict[y_var]
        plot_dict = {
            'elplot': {
                'show': False,
                'filenames': ['elel' + '_3d_' + system_name + '_' + str(element) + '_' + ref_element_str + '_' + tag + '.png'],
                'dpi': 300,
                'subplots': {
                    'elplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'center right',
                        'legend_text_size': 16,
                        'title_text': 'log(' + str(element) + '/' + ref_element_str + ') in ' + x_var_symbol + '/' + y_var_symbol + ' space for ' + system_name + ' (' + fixed_var + ' = ' + str(fixed_var_val) + ')' if not ignore_lack_of_observation else None,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': x_var + '/' + x_var_unit if x_var != 'fcf' else 'Fragment Core Fraction',
                        'xlabel_fontsize': 24,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': y_var + '/' + y_var_unit if y_var != 'fcf' else 'Fragment Core Fraction',
                        'ylabel_fontsize': 24,
                        'ylabel_fontweight': 'bold',
                        'x_min': min(p_range),
                        'x_max': max(p_range),
                        'y_min': -0.02,#min(f_range),
                        'y_max': 1.01,#max(f_range),
                        'x_tick_fontsize': 16,
                        'y_tick_fontsize': 16,
                        #'y_scale': 'log',
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        return plot_dict

    def plot_all_system_ratios(self, system_names, num_ref_systems, y_axis_ratios_obs, x_axis_ratios_obs, y_axis_ratios_obs_err, x_axis_ratios_obs_err, y_axis_ratios, x_axis_ratios, y_el_numerator, y_el_denominator, x_el_numerator, x_el_denominator, system_ps, system_fcfs, y_contour_vals, x_contour_vals, z_contour_vals):
        colours = [
            '#AF0FD5',
            '#191970',
            '#FFA07A',
            '#FF1493',
            '#CD853F',
            '#FFD700',
            '#87CEFA',
            '#32CD32',
            '#BC8F8F',
            '#808000',
            '#4F7F50',
            '#7FFFD4',
            '#5F9EA0',
            '#FFA500',
            '#C0C0C0',
            '#800000',
            '#BA55D3',
            '#66CDAA',
            '#B0E0E6',
            '#8A2BE2',
            '#98FB98',
            '#F0E68C',
            '#006400',
            '#8B0000',
            '#DAA520',
            '#D8BFD8'
        ]
        series_dict = dict()
        colour_index = 0
        for i, sys in enumerate(system_names):
            if i < num_ref_systems: # Then this is a reference system
                series_dict[sys] = {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': [x_axis_ratios_obs[i], x_axis_ratios[i]],
                    'y_data': [y_axis_ratios_obs[i], y_axis_ratios[i]],
                    'line_marker': 'x',
                    'line_style': 'None',
                    'line_color': 'k'
                    #'line_linewidth': 1,
                    #'line_type': 'kx--'
                }
            else:
                if x_axis_ratios_obs[i] is not None and y_axis_ratios_obs[i] is not None:
                    series_dict[sys] = {
                        'type': dp.PlotType.scatter_2d_error,
                        'x_data': [x_axis_ratios_obs[i], x_axis_ratios[i]],
                        'y_data': [y_axis_ratios_obs[i], y_axis_ratios[i]],
                        'x_error_data': [x_axis_ratios_obs_err[i], 0],
                        'y_error_data': [y_axis_ratios_obs_err[i], 0],
                        'line_color': self.colour_mixer('#FF0000', '#0000FF', system_fcfs[i-num_ref_systems]),
                        'line_style': ':',
                        #'line_linewidth': 1,
                        'line_marker': 'x',
                        'line_markersize': system_ps[i-num_ref_systems]/5
                    }
                    colour_index += 1
            if x_axis_ratios[i] is not None and y_axis_ratios[i] is not None:
                series_dict[sys + '_text'] = {
                    'type': dp.PlotType.text,
                    'x_pos': x_axis_ratios[i],
                    'y_pos': y_axis_ratios[i] + 0.035,
                    'text_string': sys,
                    'fontsize': 5,
                    'horizontalalignment': 'center'
                }
        series_dict['EXAMPLECONTOURS'] = {
            'type': dp.PlotType.tricontour_scatter,
            'x_data': x_contour_vals,
            'y_data': y_contour_vals,
            'z_data': z_contour_vals,
            'fill': True,
            'cbar_label': 'Pressure / GPa',
            'legend': False
        }
        plot_dict = {
            'bestfitratiosplot': {
                'show': False,
                'filenames': [str(y_el_numerator) + str(y_el_denominator) + '_' + str(x_el_numerator) + str(x_el_denominator) + '_ratios.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': False,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': str(y_el_numerator) + r'/' + str(y_el_denominator) + ' and ' + str(x_el_numerator) + '/' + str(x_el_denominator) + ' number ratios for selected systems',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'log(' + str(x_el_numerator) + '/' + str(x_el_denominator) + ')',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'log(' + str(y_el_numerator) + '/' + str(y_el_denominator) + ')',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        #'x_min': min(variable_vals),
                        #'x_max': max(variable_vals),
                        #'y_min': 0,
                        #'y_max': 100,
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)

    def plot_all_system_ratios_poster_version(self, system_names, num_ref_systems, y_axis_ratios_obs, x_axis_ratios_obs, y_axis_ratios_obs_err, x_axis_ratios_obs_err, y_axis_ratios, x_axis_ratios, y_el_numerator, y_el_denominator, x_el_numerator, x_el_denominator, system_ps, system_fcfs, y_contour_vals, x_contour_vals, z_contour_vals, text_offsets, example_fcf_lines, stellar_y, stellar_x, synth_stellar_y, synth_stellar_x, sinking_arrows, heating_arrows, synth_wd_y_dict, synth_wd_x_dict, systems_to_show_ellipse, y_axis_lower_bounds, x_axis_lower_bounds, y_axis_upper_bounds, x_axis_upper_bounds, system_categories):
        colours = [
            '#AF0FD5',
            '#191970',
            '#FFA07A',
            '#FF1493',
            '#CD853F',
            '#FFD700',
            '#87CEFA',
            '#32CD32',
            '#BC8F8F',
            '#808000',
            '#4F7F50',
            '#7FFFD4',
            '#5F9EA0',
            '#FFA500',
            '#C0C0C0',
            '#800000',
            '#BA55D3',
            '#66CDAA',
            '#B0E0E6',
            '#8A2BE2',
            '#98FB98',
            '#F0E68C',
            '#006400',
            '#8B0000',
            '#DAA520',
            '#D8BFD8'
        ]
        category_colours_dict = {
            'NED': 'k',
            'PF': 'k',
            'Unphysical': 'k',
            'Puncon': 'red',
            'Pdegen': 'red',
            'HPM': 'red',
            'LPC': 'red'
        }
        category_markers_dict = {
            'NED': 'x',
            'PF': 'o',
            'Unphysical': 'o',
            'Puncon': '+',
            'Pdegen': 'x',
            'HPM': 'o',
            'LPC': 's'
        }
        category_abbreviations = {
            'Unphysical': 'U',
            'Puncon': 'PU',
            'Pdegen': 'PD'
        }
        legible_fontsize = 10
        x_axis_descriptor = str(x_el_numerator) + str(x_el_denominator)
        y_axis_descriptor = str(y_el_numerator) + str(y_el_denominator)
        series_dict = dict()
        minimal = False
        #series_dict['SynthStellar'] = {
        #    'type': dp.PlotType.scatter_2d,
        #    'x_data': synth_stellar_x,
        #    'y_data': synth_stellar_y,
        #    'line_marker': 'x',
        #    'line_style': 'None',
        #    'line_color': 'b',
        #    'line_markersize': 1,
        #    #'line_linewidth': 1,
        #    #'line_type': 'kx--'
        #}
        legends_used = list()
        if not minimal:
            series_dict['SynthStellarEllipse1sig'] = {
                'type': dp.PlotType.ellipse,
                'x_data': synth_stellar_x,
                'y_data': synth_stellar_y,
                'n_std': 1,
                'face_colour': 'b',
                'face_alpha': 0.15,
                'fill': True
            }
            series_dict['SynthStellarEllipse2sig'] = {
                'type': dp.PlotType.ellipse,
                'x_data': synth_stellar_x,
                'y_data': synth_stellar_y,
                'n_std': 2,
                'face_colour': 'b',
                'face_alpha': 0.15,
                'fill': True
            }
            series_dict['SynthStellarEllipse3sig'] = {
                'type': dp.PlotType.ellipse,
                'x_data': synth_stellar_x,
                'y_data': synth_stellar_y,
                'n_std': 3,
                'face_colour': 'b',
                'face_alpha': 0.15,
                'fill': True
            }
            series_dict['Local Stars'] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': stellar_x,
                'y_data': stellar_y,
                'line_marker': 'x',
                'line_style': 'None',
                'line_color': 'b',
                'line_markersize': 1,
                #'line_linewidth': 1,
                #'line_type': 'kx--',
                'legend': True
            }
            #series_dict['StellarEllipse1sig'] = {
            #    'type': dp.PlotType.ellipse,
            #    'x_data': stellar_x,
            #    'y_data': stellar_y,
            #    'n_std': 1,
            #    'face_colour': 'g',
            #    'face_alpha': 0.15
            #}
            #series_dict['StellarEllipse2sig'] = {
            #    'type': dp.PlotType.ellipse,
            #    'x_data': stellar_x,
            #    'y_data': stellar_y,
            #    'n_std': 2,
            #    'face_colour': 'g',
            #    'face_alpha': 0.15
            #}
            #series_dict['StellarEllipse3sig'] = {
            #    'type': dp.PlotType.ellipse,
            #    'x_data': stellar_x,
            #    'y_data': stellar_y,
            #    'n_std': 3,
            #    'face_colour': 'g',
            #    'face_alpha': 0.15
            #}
        ref_colour = '#65a523'
        ref_text = 'Reference'
        for i, sys in enumerate(system_names):
            if i < num_ref_systems: # Then this is a reference system
                if not minimal and x_axis_ratios[i] is not None and y_axis_ratios[i] is not None:
                    if ref_text not in legends_used:
                        series_dict[ref_text] = {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': [x_axis_ratios_obs[i], x_axis_ratios[i]],
                            'y_data': [y_axis_ratios_obs[i], y_axis_ratios[i]],
                            'line_marker': '.',
                            'line_style': 'None',
                            'line_color': ref_colour,
                            #'line_linewidth': 1,
                            #'line_type': 'kx--',
                            'legend': True
                        }
                        legends_used.append(ref_text)
                    else:
                        series_dict[sys] = {
                            'type': dp.PlotType.scatter_2d,
                            'x_data': [x_axis_ratios_obs[i], x_axis_ratios[i]],
                            'y_data': [y_axis_ratios_obs[i], y_axis_ratios[i]],
                            'line_marker': '.',
                            'line_style': 'None',
                            'line_color': ref_colour
                            #'line_linewidth': 1,
                            #'line_type': 'kx--'
                        }
                    default_x = x_axis_ratios[i]
                    default_y = y_axis_ratios[i]
                    text_offset = text_offsets.get(sys, (0, 0))
                    series_dict[sys + '_text'] = {
                        'type': dp.PlotType.text,
                        'x_pos': default_x + text_offset[0],
                        'y_pos': default_y + text_offset[1] + 0.05,
                        'text_string': sys,
                        'fontsize': legible_fontsize,
                        'horizontalalignment': 'center',
                        'text_colour': 'k' if i < num_ref_systems else 'r'
                    }
            else:
                category = system_categories[sys]
                category_colour = category_colours_dict[category]
                category_marker = category_markers_dict[category]
                x_val_to_use = x_axis_ratios_obs[i]
                x_bound = 0
                if x_val_to_use is None:
                    x_val_to_use = x_axis_upper_bounds[i]
                    x_bound = 1
                if x_val_to_use is None:
                    x_val_to_use = x_axis_lower_bounds[i]
                    x_bound = -1
                y_val_to_use = y_axis_ratios_obs[i]
                y_bound = 0
                if y_val_to_use is None:
                    y_val_to_use = y_axis_upper_bounds[i]
                    y_bound = 1
                if y_val_to_use is None:
                    y_val_to_use = y_axis_lower_bounds[i]
                    y_bound = -1
                if x_val_to_use is not None and y_val_to_use is not None:
                    #pass
                    #series_dict[sys] = {
                    #    'type': dp.PlotType.scatter_2d_error,
                    #    'x_data': [x_axis_ratios_obs[i]],
                    #    'y_data': [y_axis_ratios_obs[i]],
                    #    'x_error_data': [x_axis_ratios_obs_err[i]],
                    #    'y_error_data': [y_axis_ratios_obs_err[i]],
                    #    'line_color': 'r',
                    #    'line_style': ':',
                    #    'line_linewidth': 1,
                    #    'line_marker': None,
                    #    'line_markersize': 0,
                    #    'error_linewidth': 0.5,
                    #    'capsize': 0.5
                    #}
                    if x_bound == 0 and y_bound == 0:
                        if category in legends_used:
                            series_dict[sys] = {
                                'type': dp.PlotType.scatter_2d,
                                'x_data': [x_val_to_use],
                                'y_data': [y_val_to_use],
                                'x_error_data': [x_axis_ratios_obs_err[i]],
                                'y_error_data': [y_axis_ratios_obs_err[i]],
                                'line_marker': category_marker,
                                'line_style': 'None',
                                'line_color': category_colour,
                                'legend': False
                            }
                        else:
                            series_dict[category_abbreviations.get(category, category)] = {
                                'type': dp.PlotType.scatter_2d,
                                'x_data': [x_val_to_use],
                                'y_data': [y_val_to_use],
                                'x_error_data': [x_axis_ratios_obs_err[i]],
                                'y_error_data': [y_axis_ratios_obs_err[i]],
                                'line_marker': category_marker,
                                'line_style': 'None',
                                'line_color': category_colour,
                                'legend': True
                            }
                            legends_used.append(category)
                    else:
                        if category in legends_used:
                            series_dict[sys] = {
                                'type': dp.PlotType.scatter_2d,
                                'x_data': [x_val_to_use],
                                'y_data': [y_val_to_use],
                                'x_error_data': [x_axis_ratios_obs_err[i]],
                                'y_error_data': [y_axis_ratios_obs_err[i]],
                                'line_marker': category_marker,
                                'line_style': 'None',
                                'line_color': category_colour,
                                'legend': False
                            }
                        else:
                            series_dict[category_abbreviations.get(category, category)] = {
                                'type': dp.PlotType.scatter_2d,
                                'x_data': [x_val_to_use],
                                'y_data': [y_val_to_use],
                                'x_error_data': [x_axis_ratios_obs_err[i]],
                                'y_error_data': [y_axis_ratios_obs_err[i]],
                                'line_marker': category_marker,
                                'line_style': 'None',
                                'line_color': category_colour,
                                'legend': True
                            }
                            legends_used.append(category)
                        series_dict[sys + '_arrow'] = {
                            'type': dp.PlotType.arrow,
                            'x_start': x_val_to_use,
                            'y_start': y_val_to_use,
                            'x_change': -0.15*x_bound,
                            'y_change': -0.15*y_bound,
                            'face_colour': category_colour,
                            'edge_colour': category_colour,
                            'face_alpha': 1,
                            'head_width': 0.03,
                            'legend': False,
                            'length_includes_head': True
                        }
                    if sys in systems_to_show_ellipse:
                        print(y_el_numerator)
                        print(sys)
                        sys_text = sys
                        if sys_text.endswith('Corr'):
                            sys_text = sys_text[:-4]
                        if sys_text.endswith('NoNa'):
                            sys_text = sys_text[:-4]
                        elif sys_text.endswith('X'):
                            sys_text = sys_text[:-1]
                        elif sys_text.endswith('SiO'):
                            sys_text = sys_text[:-3]
                        series_dict[sys + '_text'] = {
                                'type': dp.PlotType.text,
                                'x_pos': x_val_to_use + text_offsets.get(sys, (0, 0))[0],
                                'y_pos': y_val_to_use + text_offsets.get(sys, (0, 0))[1] + 0.05,
                                'text_string': sys_text,
                                'fontsize': legible_fontsize,
                                'horizontalalignment': 'center',
                                'text_colour': 'r'
                            }
                if synth_wd_x_dict.get(sys) is not None and synth_wd_y_dict.get(sys) is not None and x_axis_ratios_obs_err[i] is not None and y_axis_ratios_obs_err[i] is not None:
                    wd_ellipse_threshold = 0.05
                    if (sys in systems_to_show_ellipse) or (x_axis_ratios_obs_err[i] < wd_ellipse_threshold) or (y_axis_ratios_obs_err[i] < wd_ellipse_threshold):
                        #series_dict['Synth' + sys] = {
                        #    'type': dp.PlotType.scatter_2d,
                        #    'x_data': synth_wd_x_dict[sys],
                        #    'y_data': synth_wd_y_dict[sys],
                        #    'line_marker': 'x',
                        #    'line_style': 'None',
                        #    'line_color': 'r',
                        #    'line_markersize': 1,
                        #    #'line_linewidth': 1,
                        #    #'line_type': 'kx--'
                        #}
                        series_dict['Synth' + sys + 'Ellipse1sig'] = {
                            'type': dp.PlotType.ellipse,
                            'x_data': synth_wd_x_dict[sys],
                            'y_data': synth_wd_y_dict[sys],
                            'n_std': 1,
                            'face_colour': 'r',
                            'face_alpha': 0.3,
                            'fill': True
                        }
                        #series_dict['Synth' + sys + 'Ellipse2sig'] = {
                        #    'type': dp.PlotType.ellipse,
                        #    'x_data': synth_wd_x_dict[sys],
                        #    'y_data': synth_wd_y_dict[sys],
                        #    'n_std': 2,
                        #    'face_colour': 'r',
                        #    'face_alpha': 0.15
                        #}
                        #series_dict['Synth' + sys + 'Ellipse3sig'] = {
                        #    'type': dp.PlotType.ellipse,
                        #    'x_data': synth_wd_x_dict[sys],
                        #    'y_data': synth_wd_y_dict[sys],
                        #    'n_std': 3,
                        #    'face_colour': 'r',
                        #    'face_alpha': 0.15
                        #}
            if x_axis_ratios_obs[i] is not None and y_axis_ratios_obs[i] is not None: # This bit used to handle both the reference and non-reference systems but I messed something up
                pass
                #default_x = x_axis_ratios_obs[i]
                #default_y = y_axis_ratios_obs[i]
                #series_dict[sys + '_text'] = {
                #    'type': dp.PlotType.text,
                #    'x_pos': default_x + text_offsets[i][0],
                #    'y_pos': default_y + text_offsets[i][1] + 0.035,
                #    'text_string': sys,
                #    'fontsize': legible_fontsize,
                #    'horizontalalignment': 'center',
                #    'text_colour': 'k' if i < num_ref_systems else 'r'
                #}
        prev_3_fcf = None
        prev_3_elel_vals_dict = None
        prev_prev_fcf = None
        prev_prev_elel_vals_dict = None
        prev_fcf = None
        prev_elel_vals_dict = None
        cbar_plotted = False
        # Perhaps going in reverse, towards the falsely interpolated region, would avoid this complexity?
        for fcf, elel_vals_dict in example_fcf_lines.items():
            if prev_3_fcf is not None and prev_3_elel_vals_dict is not None:
                series_dict['EXAMPLECONTOURS' + str(fcf)] = {
                    'type': dp.PlotType.tricontour_scatter,
                    'x_data': prev_3_elel_vals_dict[x_axis_descriptor] + prev_prev_elel_vals_dict[x_axis_descriptor] + prev_elel_vals_dict[x_axis_descriptor] + elel_vals_dict[x_axis_descriptor],
                    'y_data': prev_3_elel_vals_dict[y_axis_descriptor] + prev_prev_elel_vals_dict[y_axis_descriptor] + prev_elel_vals_dict[y_axis_descriptor] + elel_vals_dict[y_axis_descriptor],
                    'z_data': prev_3_elel_vals_dict['p'] + prev_prev_elel_vals_dict['p'] + prev_elel_vals_dict['p'] + elel_vals_dict['p'],
                    'fill': True,
                    'cbar_label': None if cbar_plotted else 'Pressure / GPa',
                    'cbar_tickfontsize': 14,
                    'cbar_labelfontsize': 16,
                    'cbar_labelrotation': 90,
                    'legend': False
                }
                if not cbar_plotted:
                    print(series_dict)
                cbar_plotted = True
            prev_3_fcf = prev_prev_fcf
            prev_3_elel_vals_dict = prev_prev_elel_vals_dict
            prev_prev_fcf = prev_fcf
            prev_prev_elel_vals_dict = prev_elel_vals_dict
            prev_fcf = fcf
            prev_elel_vals_dict = elel_vals_dict
        for fcf, elel_vals_dict in example_fcf_lines.items():
            if (fcf > 0.7499 and fcf < 0.750001) or (fcf > 0.8999 and fcf < 0.90001) or (fcf > 0.01999 and fcf < 0.020001):
                series_dict['fcf test' + str(fcf)] = {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': elel_vals_dict[x_axis_descriptor],
                    'y_data': elel_vals_dict[y_axis_descriptor],
                    'line_color': 'k',
                    'line_style': ':',
                    'legend': False
                }
                text_string = str(int(fcf*100)) + '$\%$' + ' core'
                text_offset_lookup_string = str(int(fcf*100)) + ' core'
                text_offset = text_offsets.get(text_offset_lookup_string, (0, 0))
                series_dict[str(int(fcf*100)) + 'percent core'] = {
                    'type': dp.PlotType.text,
                    'x_pos': elel_vals_dict[x_axis_descriptor][0] + text_offset[0],
                    'y_pos': elel_vals_dict[y_axis_descriptor][0] + text_offset[1] + 0.05,
                    'text_string': text_string,
                    'fontsize': legible_fontsize,
                    'horizontalalignment': 'center',
                    'text_colour': 'k'
                }
        if not minimal and sinking_arrows is not None:
            text_offset = text_offsets.get('SinkingEffects', (0, 0))
            for i, arrow in enumerate(sinking_arrows):
                if arrow[2] != 0 and arrow[3] != 0:
                    series_dict['SinkingEffects' + str(i)] = {
                        'type': dp.PlotType.arrow,
                        'x_start': arrow[0],
                        'y_start': arrow[1],
                        'x_change': arrow[2],
                        'y_change': arrow[3],
                        'face_colour': 'k',
                        'face_alpha': 1,
                        'head_width': 0.03 if i == (len(sinking_arrows) - 1) else 0.0,
                        'legend': False,
                        'length_includes_head': True
                    }
            series_dict['SinkingEffects_text'] = {
                'type': dp.PlotType.text,
                'x_pos': sinking_arrows[0][0] + text_offset[0],
                'y_pos': sinking_arrows[0][1] + text_offset[1] + 0.05,
                'text_string': 'Sinking Effects',
                'fontsize': legible_fontsize,
                'horizontalalignment': 'center',
                'text_colour': 'k'
            }
        if not minimal and heating_arrows is not None:
            text_offset = text_offsets.get('HeatingEffects', (0, 0))
            for i, arrow in enumerate(heating_arrows):
                if arrow[2] != 0 and arrow[3] != 0:
                    series_dict['HeatingEffects' + str(i)] = {
                        'type': dp.PlotType.arrow,
                        'x_start': arrow[0],
                        'y_start': arrow[1],
                        'x_change': arrow[2],
                        'y_change': arrow[3],
                        'face_colour': 'k',
                        'face_alpha': 1,
                        'head_width': 0.03 if i == (len(heating_arrows) - 1) else 0.0,
                        'legend': False,
                        'length_includes_head': True
                    }
            series_dict['HeatingEffects_text'] = {
                'type': dp.PlotType.text,
                'x_pos': heating_arrows[0][0] + text_offset[0],
                'y_pos': heating_arrows[0][1] + text_offset[1] + 0.05,
                'text_string': 'Heating Effects',
                'fontsize': legible_fontsize,
                'horizontalalignment': 'center',
                'text_colour': 'k'
            }
        file_name_start = str(y_el_numerator) + str(y_el_denominator) + '_' + str(x_el_numerator) + str(x_el_denominator) + '_ratios_bowtie'
        plot_dict = {
            'bestfitratiosplot': {
                'show': False,
                'dpi': 1200,
                'filenames': [file_name_start + '_paper.pdf', file_name_start + '_paper.png'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': not minimal,
                        'legend_loc': 'upper left',
                        'legend_text_size': 12,
                        'title_text': str(y_el_numerator) + r'/' + str(y_el_denominator) + ' and ' + str(x_el_numerator) + '/' + str(x_el_denominator) + ' number ratios for selected systems',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'log(' + str(x_el_numerator) + '/' + str(x_el_denominator) + ')',
                        'xlabel_fontsize': 16,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'log(' + str(y_el_numerator) + '/' + str(y_el_denominator) + ')',
                        'ylabel_fontsize': 16,
                        'x_tick_fontsize': 14,
                        'y_tick_fontsize': 14,
                        'ylabel_fontweight': 'bold',
                        'x_min': -2,
                        'x_max': 1.55,
                        #'y_min': 0,
                        #'y_max': -0.7,
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        return plot_dict

    def multipanelise(self, list_of_plot_dicts, y_dimension, x_dimension, filenames, fig_height=15, fig_width=10, gridspec_wspace=None, gridspec_hspace=None, sharey_axes=True, sharex_axes=True):
        if not isinstance(filenames, list):
            filenames = [filenames]
        plot_dict = {
            'multipanelplot': {
                'show': False,
                'dpi': 300,
                'filenames': filenames,
                'gridspec_y_x': (y_dimension, x_dimension),
                'fig_height': fig_height,
                'fig_width': fig_width,
                #'gridspec_wspace': 0.07,
                #'gridspec_hspace': 0.07,
                #'gridspec_width_ratios': None,
                #'gridspec_height_ratios': None,
                'subplots': dict()
            }
        }
        if gridspec_wspace is not None:
            plot_dict['multipanelplot']['gridspec_wspace'] = gridspec_wspace
        if gridspec_hspace is not None:
            plot_dict['multipanelplot']['gridspec_hspace'] = gridspec_hspace
        gridspec_index = 0
        y_coord = 0
        x_coord = 0
        for pd in list_of_plot_dicts:
            for plot_name, plot in pd.items():
                for subplot_name, subplot in plot['subplots'].items():
                    new_subplot_name = 'subplot_' + str(y_coord) + '_' + str(x_coord)
                    plot_dict['multipanelplot']['subplots'][new_subplot_name] = subplot
                    if y_coord != 0:
                        plot_dict['multipanelplot']['subplots'][new_subplot_name]['title_text'] = None
                    plot_dict['multipanelplot']['subplots'][new_subplot_name]['gridspec_index'] = gridspec_index
                    plot_dict['multipanelplot']['subplots'][new_subplot_name]['sharey_subplot'] = None  # Overwriting any unhelpful metadata that may already be present...
                    plot_dict['multipanelplot']['subplots'][new_subplot_name]['sharex_subplot'] = None
                    
                    if y_coord < y_dimension - 1:
                        if sharex_axes:
                            plot_dict['multipanelplot']['subplots'][new_subplot_name]['sharex_subplot'] = 'subplot_' + str(y_coord+1) + '_' + str(x_coord)
                    if x_coord > 0:
                        if sharey_axes:
                            plot_dict['multipanelplot']['subplots'][new_subplot_name]['sharey_subplot'] = 'subplot_' + str(y_coord) + '_' + str(x_coord-1)
                    x_coord += 1
                    if x_coord >= x_dimension:
                        x_coord = 0
                        y_coord += 1
                    gridspec_index += 1
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        return plot_dict
            
    def plot_size_v_pressure(self, series_to_plot, property_string=None):
        # series_to_plot = {'Name of series': {'Pressure': [pressure vals], 'Radius': [radius vals], 'Mass': [mass vals]}}
        # property_string: either 'Mass' or 'Radius' or anything else in order to plot both
        radius_series_dict = dict()
        i = 0
        line_styles =['-', '--', ':']
        j = 0
        distinct_colours = ['k', 'r']
        for series_name, series_info in series_to_plot.items():
            series_name_to_use = series_name if property_string == 'Radius' else 'Radius (' + series_name + ')'
            radius_series_dict[series_name_to_use] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': series_info['Pressure'],
                'y_data': series_info['Radius'],
                #'line_marker': 'x',
                #'line_style': 'None',
                'line_color': distinct_colours[i],
                'line_linewidth': 2,
                'line_style': line_styles[j],
                'legend': True
            }
            j += 1
        # This bit of code is mildly horrible:
        earth_radius_series_name = 'Earth' if property_string in ['Mass', 'Radius'] else 'Object radii'
        mars_radius_series_name = 'Mars' if property_string in ['Mass', 'Radius'] else 'Mars radius'
        moon_radius_series_name = 'Moon' if property_string in ['Mass', 'Radius'] else 'Moon radius'
        radius_series_dict[earth_radius_series_name] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [54],
            'y_data': [6371],
            'line_marker': '+',
            'line_style': 'None',
            'line_color': 'green' if property_string in ['Mass', 'Radius'] else self.colour_list[i],
            'line_linewidth': 5,
            'line_markersize': 10,
            'line_type': 'k-',
            'legend': True
        }
        radius_series_dict[mars_radius_series_name] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [13],
            'y_data': [3390],
            'line_marker': '+',
            'line_style': 'None',
            'line_color': 'gold' if property_string in ['Mass', 'Radius'] else self.colour_list[i],
            'line_linewidth': 5,
            'line_markersize': 10,
            'line_type': 'k-',
            'legend': property_string in ['Mass', 'Radius']
        }
        radius_series_dict[moon_radius_series_name] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [3.5],
            'y_data': [1737.4],
            'line_marker': '+',
            'line_style': 'None',
            'line_color': 'gold' if property_string in ['Mass', 'Radius'] else self.colour_list[i],
            'line_linewidth': 5,
            'line_markersize': 10,
            'line_type': 'k-',
            'legend': property_string in ['Mass', 'Radius']
        }
        mass_series_dict = dict()
        if property_string not in ['Mass', 'Radius']:
            i = 1
        j = 0
        for series_name, series_info in series_to_plot.items():
            series_name_to_use = series_name if property_string == 'Mass' else 'Mass (' + series_name + ')'
            mass_series_dict[series_name_to_use] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': series_info['Pressure'],
                'y_data': series_info['Mass'],
                #'line_marker': 'x',
                #'line_style': 'None',
                'line_color': distinct_colours[i],
                'line_linewidth': 2,
                'line_style': line_styles[j],
                'legend': True
            }
            j += 1
        earth_mass_series_name = 'Earth' if property_string in ['Mass', 'Radius'] else 'Object masses'
        mars_mass_series_name = 'Mars' if property_string in ['Mass', 'Radius'] else 'Mars mass'
        moon_mass_series_name = 'Moon' if property_string in ['Mass', 'Radius'] else 'Moon mass'
        mass_series_dict[earth_mass_series_name] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [54],
            'y_data': [1],
            'line_marker': 'x',
            'line_style': 'None',
            'line_color': 'green' if property_string in ['Mass', 'Radius'] else distinct_colours[i],
            'line_linewidth': 5,
            'line_markersize': 10,
            'line_type': 'k-',
            'legend': True
        }
        mass_series_dict[mars_mass_series_name] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [13],
            'y_data': [0.107],
            'line_marker': 'x',
            'line_style': 'None',
            'line_color': 'gold' if property_string in ['Mass', 'Radius'] else distinct_colours[i],
            'line_linewidth': 5,
            'line_markersize': 10,
            'line_type': 'k-',
            'legend': property_string in ['Mass', 'Radius']
        }
        mass_series_dict[moon_mass_series_name] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [3.5],
            'y_data': [0.0123],
            'line_marker': 'x',
            'line_style': 'None',
            'line_color': 'gold' if property_string in ['Mass', 'Radius'] else distinct_colours[i],
            'line_linewidth': 5,
            'line_markersize': 10,
            'line_type': 'k-',
            'legend': property_string in ['Mass', 'Radius']
        }
        
        radius_subplot = {
            'subplot_region': 111,
            'legend': True,
            'legend_loc': 'lower right',
            'legend_text_size': 28 if property_string in ['Mass', 'Radius'] else 16,
            'title_fontsize': 16,
            'title_fontweight': 'bold',
            'xlabel_text': 'Mid-mantle Pressure / GPa',
            'xlabel_fontsize': 28,
            'xlabel_fontweight': 'bold',
            'ylabel_text': 'Radius / km',
            'ylabel_fontsize': 28,
            'ylabel_fontweight': 'bold',
            'x_min': 0,
            'x_max': 60,
            'y_min': 0,
            'font': 'STIXGeneral',
            'x_tick_fontsize': 20,
            'y_tick_fontsize': 20,
            'series': radius_series_dict
        }
        mass_subplot = {
            'subplot_region': 111,
            'legend': True,
            'legend_loc': 'best',
            'legend_text_size': 28,
            'title_fontsize': 16,
            'title_fontweight': 'bold',
            'xlabel_text': 'Pressure / GPa',
            'xlabel_fontsize': 28,
            'xlabel_fontweight': 'bold',
            #'ylabel_text': 'Mass / M' + r'$_{\textrm{Earth}}$',
            'ylabel_text': 'Mass / M' + r'$_{\oplus}$',
            'ylabel_fontsize': 28,
            'ylabel_fontweight': 'bold',
            'x_min': 0,
            'x_max': 60,
            'y_min': 0,
            'font': 'STIXGeneral',
            'y_label_colour': 'r',
            'y_tick_colour': 'r',
            'x_tick_fontsize': 20,
            'y_tick_fontsize': 20,
            'series': mass_series_dict
        }
        if property_string in ['Mass', 'Radius']:
            radius_subplot['title_text'] = property_string + ' against differentiation pressure'
        object_names = {
            'Mars': 13,
            'Earth': 54,
            'Moon': 3.5,
            #'Asteroids': 0
        }
        object_x_pos = {
            'Asteroids': -1,
            #'Moon': 1,
        }
        object_y_pos = {
            'Earth': 5000
        }
        if property_string == 'Mass':
            for obj_name, pressure in object_names.items():
                mass_series_dict[obj_name + '_text'] = {
                    'type': dp.PlotType.text,
                    'x_pos': pressure - object_x_pos.get(obj_name, 1.5),
                    'y_pos': object_y_pos.get(obj_name, 1),
                    'text_string': obj_name,
                    'fontsize': 24,
                    'rotation': 90
                }
        else:
            for obj_name, pressure in object_names.items():
                radius_series_dict[obj_name + '_text'] = {
                    'type': dp.PlotType.text,
                    'x_pos': pressure - object_x_pos.get(obj_name, 1.5),
                    'y_pos': object_y_pos.get(obj_name, 6000),
                    'text_string': obj_name,
                    'fontsize': 24,
                    'rotation': 90
                }
        plot_dict = {
            'pvsplot': {
                'show': False,
                'filenames': [property_string + '_v_pressure.png'] if property_string in ['Mass', 'Radius'] else ['MassRadius_v_pressure.png', 'MassRadius_v_pressure.pdf'],
                'dpi': 500,
                'fig_height': 6,
                'fig_width': 9
            }
        }
        if property_string == 'Radius':
            plot_dict['pvsplot']['subplots'] = {'radiusplot': radius_subplot}
        elif property_string == 'Mass':
            plot_dict['pvsplot']['subplots'] = {'massplot': mass_subplot}
        else:
            plot_dict['pvsplot']['subplots'] = cn.OrderedDict({'radiusplot': radius_subplot, 'massplot': mass_subplot})
            plot_dict['pvsplot']['subplots']['massplot']['twin_subplot'] = 'radiusplot'
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        
    def plot_retrieved_variable_against_fcf(self, system_name, variable, fcf_values, variable_values, variable_1s_upper_errors, variable_1s_lower_errors, variable_2s_upper_errors=None, variable_2s_lower_errors=None, tag=None):
        earth_vals = {
            'Pressure': 54,
            'Oxygen Fugacity': -2,
            'Formation Distance': 0,
            'Parent Core Fraction': 0.17
        }
        mars_vals = {
            'Pressure': 13,
            'Oxygen Fugacity': -1,
            'Formation Distance': np.log10(1.5),
            'Parent Core Fraction': 0.17  #From Table 11 in Yoshizaki+ 19, it seems its actually about the same as Earth
        }
        system_vals = {
            'Earth': earth_vals,
            'Mars': mars_vals
        }
        # Yeah this is a hack
        ylabel_dict = {
            'Pressure': 'Retrieved Pressure /GPa',
            'Oxygen Fugacity': 'Retrieved Oxygen Fugacity ' + r'$\Delta$' + 'IW',
            'Formation Distance': 'log(Retrieved Formation Distance /AU)',
            'Fragment Core Fraction': 'Retrieved Fragment Core Fraction',
            'Parent Core Fraction': 'Retrieved Parent Core Fraction'
        }
        y_limits_dict = {
            'Pressure': (0, 60),
            'Oxygen Fugacity': (-3, -1),
            'Formation Distance': (None, None),
            'Fragment Core Fraction': (0, 1),
            'Parent Core Fraction': (None, None)
        }
        text_pos_dict = {
            'Earth': {
                'Pressure': (0.27, 10),
                'Oxygen Fugacity': (0.4, -1.4),
                'Formation Distance': (0.5, -0.2),
                'Fragment Core Fraction': (0.75, 0.25),
                'Parent Core Fraction': (0.5, 0.15)
            },
            'Mars': {
                'Pressure': (0.27, 55),
                'Oxygen Fugacity': (0.7, -2.75),
                'Formation Distance': (0.5, -0.2),
                'Fragment Core Fraction': (0.75, 0.25),
                'Parent Core Fraction': (0.5, 0.09)
            }
        }
        differentiation_text_pos_dict = {
            'Earth': {
                'Pressure': (0.135, 45),
                'Oxygen Fugacity': (0.4, -1.4),
                'Formation Distance': (0.5, -0.2),
                'Fragment Core Fraction': (0.75, 0.25),
                'Parent Core Fraction': (0.5, 0.15)
            },
            'Mars': {
                'Pressure': (0.08, 55),
                'Oxygen Fugacity': (0.7, -2.75),
                'Formation Distance': (0.5, -0.2),
                'Fragment Core Fraction': (0.75, 0.25),
                'Parent Core Fraction': (0.5, 0.09)
            }
        }
        series_dict = dict()
        series_dict['System Text'] = {
            'type': dp.PlotType.text,
            'x_pos': text_pos_dict[system_name][variable][0],
            'y_pos': text_pos_dict[system_name][variable][1],
            'text_string': system_name,
            'fontsize': 32
        }
        #draw_2_sigma_bars = True
        #if draw_2_sigma_bars:
        #    if variable_2s_upper_errors is not None and variable_2s_lower_errors is not None and not np.isnan(variable_2s_upper_errors).all() and not np.isnan(variable_2s_lower_errors).all():
        #        series_dict['Retrieved value (2 sigma)'] = {
        #            'type': dp.PlotType.scatter_2d_error,
        #            'x_data': fcf_values,
        #            'y_data': variable_values,
        #            'y_error_data': [variable_2s_lower_errors, variable_2s_upper_errors],
        #            'line_marker': 'x',
        #            'line_markersize': 5,
        #            'line_style': 'None',
        #            'line_color': 'b',
        #            'line_linewidth': 1,
        #            'error_linewidth': 1,
        #            'capsize': 1,
        #            'legend': True
        #        }
        #series_dict['Retrieved value'] = {
        #    'type': dp.PlotType.scatter_2d_error,
        #    'x_data': fcf_values,
        #    'y_data': variable_values,
        #    'y_error_data': [variable_1s_lower_errors, variable_1s_upper_errors],
        #    'line_marker': 'x',
        #    'line_markersize': 5,
        #    'line_style': 'None',
        #    'line_color': 'k',
        #    'line_linewidth': 1,
        #    'error_linewidth': 1,
        #    'capsize': 1,
        #    'legend': True
        #}
        series_dict['Retrieved value'] = {
            'type': dp.PlotType.scatter_2d_error,
            'x_data': fcf_values,
            'y_data': variable_values,
            'y_error_data': [variable_2s_lower_errors, variable_2s_upper_errors],
            'line_marker': 'x',
            'line_markersize': 5,
            'line_style': 'None',
            'line_color': 'k',
            'line_linewidth': 1,
            'error_linewidth': 1,
            'capsize': 1,
            'legend': True
        }
        i = 0
        for fcf, p in zip(fcf_values, variable_values):
            if np.isnan(p):
                shade_width = 0.05 # should ideally calculate this dynamically
                series_dict['Nodiff' + str(i)] = {
                    'type': dp.PlotType.shade,
                    'x_data': [fcf - shade_width, fcf + shade_width],
                    'y_data': [y_limits_dict[variable][0], y_limits_dict[variable][0]],
                    'y_shade_data': [y_limits_dict[variable][1], y_limits_dict[variable][1]],
                    'shade_colour': 'grey',
                    'shade_alpha': 0.15,
                    'zorder': 2,
                    'legend': False
                }
                i += 1
                series_dict['No Diff Text'] = {
                    'type': dp.PlotType.text,
                    'x_pos': differentiation_text_pos_dict[system_name][variable][0],
                    'y_pos': differentiation_text_pos_dict[system_name][variable][1],
                    'text_string': 'Differentiation not invoked',
                    'fontsize': 28,
                    'rotation': 90
                }
        if system_vals[system_name].get(variable) is not None:
            series_dict['Approximate value'] = {
                'type': dp.PlotType.hline,
                'y_start': system_vals[system_name][variable],
                'x_min': -0.1,
                'x_max': 1,
                'legend': True,
                'colours': 'r',
                'line_styles': 'dashed'
            }
        if variable == 'Fragment Core Fraction':
            series_dict['Retrieved = Input'] = {
                'type': dp.PlotType.function_2d,
                'x_start': 0,
                'x_end': 1,
                'x_points': 100,
                'function': lambda x: x,
                'legend': True,
                'line_color': 'r',
                'line_style': 'dashed'
            }
        plot_dict = {
            'retrieved_variable_plot': {
                'show': False,
                'filenames': ['retrieved_' + variable.replace(' ', '_') + '_v_fcf_' + system_name + '_' + tag + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'lower right' if variable == 'Pressure' else'best',
                        'legend_text_size': 22,
                        'title_text': r'Retrieved ' + variable + ' against Fragment Core Fraction for ' + system_name + ' (' + tag + ')',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Fragment Core Fraction',
                        'xlabel_fontsize': 28,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': ylabel_dict[variable],
                        'ylabel_fontsize': 28,
                        'ylabel_fontweight': 'bold',
                        'x_min': -0.05,
                        'x_max': 1.05,
                        'y_min': y_limits_dict[variable][0],
                        'y_max': y_limits_dict[variable][1],
                        'x_tick_fontsize': 22,
                        'y_tick_fontsize': 22,
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        return plot_dict

    def plot_2D_timescales(self, timescales_to_plot, variable, x_vals, HorHe, logg, Teff, CaHe):
        x_labels = {
            'Teff': 'Temperature /K',
            'logg': 'log(g / cm s$^{-2}$)',
            'CaHe': 'log(Ca/He)',
        }
        series_dict = dict()
        i = 0
        scale_factor = 0.000000001 if HorHe == 'He' else 1
        time_unit = 'Gyr' if HorHe == 'He' else 'yr'
        elements = list()
        for element, timescales in timescales_to_plot.items():
            elements.append(element)
            series_dict[str(element)] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': x_vals,
                'y_data': [t*scale_factor for t in timescales],
                'legend': True,
                'line_color': self.colour_list[i]
            }
            i += 1
        #plot_ratio = True
        #if len(elements == 2) and plot_ratio:
        #    ratio_data = list()
        #    j = 0
        #    while j < len(timescales_to_plot[elements[0]]):
        #        ratio_data.append(timescales_to_plot[elements[0]][j] / timescales_to_plot[elements[1]][j])
        #        j += 1
        #    ratio_series_dict[str(elements[0]) + '/' + str(elements[1])] = {
        #        'type': dp.PlotType.scatter_2d,
        #        'x_data': x_vals,
        #        'y_data': ratio_data,
        #        'legend': True,
        #        'line_color': self.colour_list[i+1]
        #    }
        title_text = 'Sinking times for a ' + HorHe + ' WD at '
        vars_added = 0
        var_limit = 2 if HorHe == 'He' else 1
        if logg is not None:
            vars_added += 1
            title_text += 'log(g / cm s$^{-2}$) = ' + str(logg)
            if vars_added < var_limit:
                title_text += ', '
        if Teff is not None:
            vars_added += 1
            title_text += 'T = ' + str(Teff) + 'K'
            if vars_added < var_limit:
                title_text += ', '
        if CaHe is not None and HorHe == 'He':
            vars_added += var_limit
            title_text += 'log(Ca/He) = ' + str(CaHe)
            if vars_added < 2:
                title_text += ', '
        plot_dict = {
            '2d_timescales_plot': {
                'show': False,
                'filenames': ['timescales_v_' + variable + '_' + HorHe + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': title_text,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': x_labels[variable],
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'Sinking time /' + time_unit,
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'font': 'STIXGeneral',
                        'y_scale': 'log' if (variable == 'Teff' and HorHe == 'H') or (variable == 'logg') else None,
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        
    def plot_timescale_model_comparison(self, to_plot_new, to_plot_old, variable, x_vals, x_vals_old, HorHe, logg, Teff, CaHe, elements_to_ratio):
        x_labels = {
            'Teff': 'Temperature /K',
            'logg': 'log(g / cm s$^{-2}$)',
            'CaHe': 'log(Ca/He)',
        }
        series_dict = dict()
        i = 0
        scale_factor = 0.000000001 if HorHe == 'He' else 1
        time_unit = 'Gyr' if HorHe == 'He' else 'yr'
        for element, timescales in to_plot_new.items():
            series_dict[str(element)] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': x_vals,
                'y_data': [t*scale_factor for t in timescales],
                'legend': True,
                'line_color': self.colour_list[i]
            }
            i += 1
        if to_plot_old is not None:
            i = 0
            for element, timescales in to_plot_old.items():
                series_dict[str(element) + ' Hollands 2017'] = {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals_old,
                    'y_data': [t*scale_factor for t in timescales],
                    'legend': True,
                    'line_color': self.colour_list[i],
                    'line_style': ':'
                }
                i += 1
        print(series_dict)
        ratio_series_dict = dict()
        new_ratios = list()
        try:
            j = 0
            while j < len(x_vals):
                new_ratios.append(to_plot_new[elements_to_ratio[0]][j]/to_plot_new[elements_to_ratio[1]][j])
                j += 1
            ratio_series_dict[str(elements_to_ratio[0]) + '/' + str(elements_to_ratio[1])] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': x_vals,
                'y_data': new_ratios,
                'legend': True,
                'line_color': 'r'
            }
        except KeyError:
            print('Warning could not find one or both of elements_to_ratio in to_plot_new')
            print('Keys in elements_to_ratio:')
            print(elements_to_ratio)
            print('Keys in to_plot_new:')
            print(to_plot_new.keys())
            print('Skipping...')
        if to_plot_old is not None:
            try:
                old_ratios = list()
                j = 0
                while j < len(x_vals_old):
                    old_ratios.append(to_plot_old[elements_to_ratio[0]][j]/to_plot_old[elements_to_ratio[1]][j])
                    j += 1
                ratio_series_dict[str(elements_to_ratio[0]) + '/' + str(elements_to_ratio[1]) + ' Hollands 2017'] = {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': x_vals_old,
                    'y_data': old_ratios,
                    'legend': True,
                    'line_color': 'r',
                    'line_style': ':'
                }
            except KeyError:
                print('Warning could not find one or both of elements_to_ratio in to_plot_old')
                print('Keys in elements_to_ratio:')
                print(elements_to_ratio)
                print('Keys in to_plot_old:')
                print(to_plot_old.keys())
                print('Skipping...')
        print(ratio_series_dict)
        title_text = 'Sinking times for a ' + HorHe + ' WD at '
        vars_added = 0
        var_limit = 2 if HorHe == 'He' else 1
        if logg is not None:
            vars_added += 1
            title_text += 'log(g / cm s$^{-2}$) = ' + str(logg)
            if vars_added < var_limit:
                title_text += ', '
        if Teff is not None:
            vars_added += 1
            title_text += 'T = ' + str(Teff) + 'K'
            if vars_added < var_limit:
                title_text += ', '
        if CaHe is not None and HorHe == 'He':
            vars_added += var_limit
            title_text += 'log(Ca/He) = ' + str(CaHe)
            if vars_added < 2:
                title_text += ', '
        plot_dict = {
            '2d_timescales_plot': {
                'show': False,
                'filenames': ['timescales_v_' + variable + '_' + HorHe + '_model_comparison.pdf'],
                'subplots': cn.OrderedDict({
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': title_text,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': x_labels[variable],
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'Sinking time /' + time_unit,
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'font': 'STIXGeneral',
                        'y_scale': 'log' if (variable == 'Teff' and HorHe == 'H') or (variable == 'logg') else None,
                        'series': series_dict
                    },
                    'ratio_subplot': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        #'title_text': title_text,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        #'xlabel_text': x_labels[variable],
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': str(elements_to_ratio[0]) + '/' + str(elements_to_ratio[1]),
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'font': 'STIXGeneral',
                        'series': ratio_series_dict,
                        'twin_subplot': 'subplot1'
                    }
                })
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)

    def plot_el_v_Teff(self, element, el_vals, teff_vals, el_upper_bound_vals, teff_upper_bound_vals):
        series_dict = dict()
        He_key = str(element) + '/He'
        H_key = str(element) + '/H'
        series_dict[He_key] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': teff_vals['He'],
            'y_data': el_vals['He'],
            'legend': True,
            'line_color': 'b',
            'line_style': 'None',
            'line_marker': 'o',
        }
        series_dict[H_key] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': teff_vals['H'],
            'y_data': el_vals['H'],
            'legend': True,
            'line_color': 'r',
            'line_style': 'None',
            'line_marker': 'o',
        }
        series_dict[He_key + ' upper bounds'] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': teff_upper_bound_vals['He'],
            'y_data': el_upper_bound_vals['He'],
            'legend': True,
            'line_color': 'b',
            'line_style': 'None',
            'line_marker': 'x',
        }
        series_dict[H_key + ' upper bounds'] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': teff_upper_bound_vals['H'],
            'y_data': el_upper_bound_vals['H'],
            'legend': True,
            'line_color': 'r',
            'line_style': 'None',
            'line_marker': 'x',
        }
        plot_dict = {
            'el_v_Teff': {
                'show': False,
                'filenames': [str(element) + '_v_Teff.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Teff',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': str(element) + '/Hx',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'font': 'STIXGeneral',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)

    def make_fcrit_plot(self, variable_fcrit_dict, pressure_vals, fO2):
        series_dict = dict()
        i = 0
        for element, f_crit_list in variable_fcrit_dict.items():
            if element in [ci.Element.Si, ci.Element.Cr, ci.Element.Ti, ci.Element.Fe, ci.Element.Ni, ci.Element.C]:
                series_dict[str(element)] = {
                    'type': dp.PlotType.scatter_2d,
                    'x_data': pressure_vals,
                    'y_data': f_crit_list,
                    'legend': True,
                    'line_color': self.colour_list[i]
                }
                i += 1
        plot_dict = {
            'f_crit_plot': {
                'show': False,
                'filenames': ['fcrit_v_pressure.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': 'Critical Fragment Core Fraction against Pressure at fO2 = ' + str(fO2),
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Pressure / GPa',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'Critical Fragment Core Fraction',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'font': 'STIXGeneral',
                        'y_min': 0,
                        'y_max': 1,
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)

    def make_dlogX_dP_plot(self, variable_dlogX_dP_dict, fcf_vals, fO2):
        series_dict = dict()
        linestyles = ['-', '--', ':']
        j = 0
        for element, variable_dlogX_dP_el_dict in variable_dlogX_dP_dict.items():
            if element not in [ci.Element.Si, ci.Element.Cr, ci.Element.Ni]:
                continue
            i = 0
            for p, dlogX_dP_list in variable_dlogX_dP_el_dict.items():
                if element in [ci.Element.Si, ci.Element.Cr, ci.Element.Ni]:
                    print(fcf_vals)
                    print(dlogX_dP_list)
                    series_dict[str(element)+', ' + str(p) + ' GPa'] = {
                        'type': dp.PlotType.scatter_2d,
                        'x_data': fcf_vals,
                        'y_data': dlogX_dP_list,
                        'legend': True,
                        'line_color': self.colour_list[j],
                        'line_style': linestyles[i]
                    }
                    i += 1
            j += 1
        resolvability_limit = 0.01
        series_dict['hline0'] = {
            'type': dp.PlotType.hline,
            'y_start': 0,
            'x_min': 0,
            'x_max': 1,
            'line_styles': ':'
        }
        series_dict['hlineresolvability_limitupper'] = {
            'type': dp.PlotType.hline,
            'y_start': resolvability_limit,
            'x_min': 0,
            'x_max': 1,
            'line_styles': ':'
        }
        series_dict['hlineresolvability_limitlower'] = {
            'type': dp.PlotType.hline,
            'y_start': -resolvability_limit,
            'x_min': 0,
            'x_max': 1,
            'line_styles': ':'
        }
        series_dict['Resolvable'] = {
            'type': dp.PlotType.shade,
            'x_data': [0, 1],
            'y_data': [resolvability_limit, resolvability_limit],
            'y_shade_data': [100*resolvability_limit, 100*resolvability_limit],
            'shade_colour': 'g',
            'shade_alpha': 0.15,
            'zorder': 2,
            'legend': False
        }
        series_dict['ResolvableText1'] = {
            'type': dp.PlotType.text,
            'x_pos': 0.1,
            'y_pos': 0.012,
            'text_string': 'Resolvable',
            'horizontalalignment': 'center',
            'fontsize': 16,
            'legend': False
        }
        series_dict['ResolvableLower'] = {
            'type': dp.PlotType.shade,
            'x_data': [0, 1],
            'y_data': [-resolvability_limit, -resolvability_limit],
            'y_shade_data': [-100*resolvability_limit, -100*resolvability_limit],
            'shade_colour': 'g',
            'shade_alpha': 0.15,
            'zorder': 2,
            'legend': False
        }
        series_dict['ResolvableText2'] = {
            'type': dp.PlotType.text,
            'x_pos': 0.1,
            'y_pos': -0.013,
            'text_string': 'Resolvable',
            'horizontalalignment': 'center',
            'fontsize': 16,
            'legend': False
        }
        plot_dict = {
            'f_crit_plot': {
                'show': False,
                'filenames': ['dlogX_dP_v_fcf.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_coord_x': 0.25,
                        'legend_coord_y': 0.5,
                        'legend_text_size': 9,
                        #'title_text': 'Pressure Sensitivity against Fragment Core Fraction  (Bulk Earth at fO2 = IW' + str(fO2) + ')',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Fragment Core Fraction',
                        'xlabel_fontsize': 20,
                        'xlabel_fontweight': 'bold',
                        'x_tick_fontsize': 16,
                        'y_tick_fontsize': 16,
                        'ylabel_text': 'Pressure Sensitivity / GPa$^{-1}$',
                        'ylabel_fontsize': 20,
                        'ylabel_fontweight': 'bold',
                        'font': 'STIXGeneral',
                        'x_min': 0.001,
                        'x_max': 1,                       
                        'y_min': -0.015,
                        'y_max': 0.015,
                        'x_scale': 'log',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)

    def plot_timesince_v_accretiontime(self, x_data, y_data, system_name, file_prefix, t_Mg):
        series_dict = dict()
        xy_min = 0
        xy_max = 8
        buss_border_mg_timescales = 5
        series_dict['dec_border'] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [xy_min, xy_max],
            'y_data': [xy_min, xy_max],
            'legend': False,
            'line_color': 'grey',
            'line_styles': '--'
        }
        series_dict['buss_border'] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [np.log10(buss_border_mg_timescales*t_Mg), np.log10(buss_border_mg_timescales*t_Mg)],
            'y_data': [xy_max, np.log10(buss_border_mg_timescales*t_Mg)],
            'legend': False,
            'line_color': 'grey',
            'line_styles': '--'
        }

        if (0.3 <= t_Mg < 2.5) or t_Mg > 200000:
            if 0.3 <= t_Mg < 2.5:
                x_start = 1.2
                y_start = 7.4
                x_change = -0.9
                y_change = 0
            elif 200000 < t_Mg <= 1500000:
                x_start = 7
                y_start = 6.2
                x_change = 0
                y_change = 1
            else:
                x_start = 7.6
                y_start = 6.2
                x_change = 0
                y_change = 1.5
            series_dict['arrow'] = {
                'type': dp.PlotType.arrow,
                'x_start': x_start,
                'y_start': y_start,
                'x_change': x_change,
                'y_change': y_change,
                'head_width': 0.12,
                'line_linewidth': 1.25,
                'face_colour': 'k'
            }
        series_dict['dec_phase_text'] = {
            'type': dp.PlotType.text,
            'x_pos': 4,
            'y_pos': 0.5,
            'text_string': 'Declining Phase',
            'horizontalalignment': 'center'
        }
        series_dict['ss_text'] = {
            'type': dp.PlotType.text,
            'x_pos': (7.5+np.log10(5*t_Mg))/2 if t_Mg <= 200000 else 7,
            'y_pos': 7.5 if t_Mg <= 200000 else 5.8,
            'text_string': 'Steady State',
            'horizontalalignment': 'center'
        }
        series_dict['ss_phase_text'] = {
            'type': dp.PlotType.text,
            'x_pos': (7.5+np.log10(5*t_Mg))/2 if t_Mg <= 200000 else 7,
            'y_pos': 7.1 if t_Mg <= 200000 else 5.4,
            'text_string': 'Phase',
            'horizontalalignment': 'center'
        }
        if t_Mg >= 0.3:
            series_dict['bu_text'] = {
                'type': dp.PlotType.text,
                'x_pos': 1.7 if t_Mg < 2.5 else np.log10(5*t_Mg)/2,
                'y_pos': 7.5,
                'text_string': 'Build-Up',
                'horizontalalignment': 'center'
            }
            series_dict['bu_phase_text'] = {
                'type': dp.PlotType.text,
                'x_pos': 1.7 if t_Mg < 2.5 else np.log10(5*t_Mg)/2,
                'y_pos': 7.1,
                'text_string': 'Phase',
                'horizontalalignment': 'center'
            }
        series_dict['sys_name'] = {
            'type': dp.PlotType.text,
            'x_pos': 7,
            'y_pos': 1,
            'text_string': self.strip_system_suffix(system_name),
            'horizontalalignment': 'center'
        }
        bin1d = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,6.6,6.8,7,7.2,7.4,7.6,7.8,8]
        bins = (bin1d, bin1d)
        series_dict['hist'] = {
            'type': dp.PlotType.hist2d,
            'x_data': x_data,
            'y_data': y_data,
            'bins': bins,
            'normed': True,
            'cmap': self.discrete_cmap(8, 'Greys'),
            'cbar_label': 'Probability density',
            'cbar_ticks': [0.00,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.00,1.125,1.25,1.375,1.5,1.625,1.75],
            'cbar_ticklabels': [0.00/25,0.125/25,0.25/25,0.375/25,0.5/25,0.625/25,0.75/25,0.875/25,1.00/25,1.125/25,1.25/25,1.375/25,1.5/25,1.625/25,1.75/25],
            'cbar_labelfontsize': 14,
            'cbar_labelrotation': 270,
            'cbar_labelpad': 18
        }
        plot_dict = {
            'time_since_plot': {
                'show': False,
                'filenames': [system_name + '_' + file_prefix + 'timesince_v_accretiontime.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,         #plt.legend(loc='lower right', frameon = False, handletextpad=0.5, fontsize=9)
                        'legend_loc': 'lower right',
                        'legend_text_size': 9,
                        'xlabel_text': 'log(Time since Accretion Started/Yrs)',
                        'xlabel_fontsize': 14,
                        'ylabel_text': 'log(Accretion Event Lifetime/Yrs)',
                        'ylabel_fontsize': 14,
                        'font': 'STIXGeneral',
                        'y_min': xy_min,
                        'y_max': xy_max,
                        'x_min': xy_min,
                        'x_max': xy_max,
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)                                              

    def plot_pressure_v_oxygen_fugacity(self, x_data, y_data, system_name, file_prefix):
        series_dict = dict()
        #x_bins = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60]
        x_bins = [0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60]
        y_bins = [-3,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1]
        bins = (x_bins, y_bins)
        series_dict['hist'] = {
            'type': dp.PlotType.hist2d,
            'x_data': x_data,
            'y_data': y_data,
            'bins': bins,
            'normed': True,
            'cmap': self.discrete_cmap(5, 'Greys'),
            'cbar_label': 'Probability density',
            'cbar_ticks': [0.0,0.005,0.01,0.015,0.02],
            'cbar_ticklabels': [0.0,0.005,0.01,0.015,0.02],
            #'cbar_ticklabels': [0.00/25,0.125/25,0.25/25,0.375/25,0.5/25,0.625/25,0.75/25,0.875/25,1.00/25,1.125/25,1.25/25,1.375/25,1.5/25,1.625/25,1.75/25],
            'cbar_labelfontsize': 20,
            'cbar_tickfontsize': 16,
            'cbar_labelrotation': 90,
            'cbar_labelpad': 18
        }
        plot_dict = {
            'p_fO2_plot': {
                'show': False,
                'filenames': [system_name + '_' + file_prefix + 'pressure_v_fO2.pdf', system_name + '_' + file_prefix + 'pressure_v_fO2.png'],
                'dpi': 300,
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': False,         #plt.legend(loc='lower right', frameon = False, handletextpad=0.5, fontsize=9)
                        'legend_loc': 'lower right',
                        'legend_text_size': 9,
                        'xlabel_text': 'Pressure /GPa',
                        'xlabel_fontsize': 24,
                        'ylabel_text': 'Oxygen Fugacity ' + r'$\Delta$' + 'IW',
                        'ylabel_fontsize': 24,
                        'x_tick_fontsize': 16,
                        'y_tick_fontsize': 16,
                        'y_tick_locations': [-3, -2.5, -2, -1.5, -1],
                        'font': 'STIXGeneral',
                        'y_min': -3,
                        'y_max': -1,
                        'x_min': 0,
                        'x_max': 60,
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)

    def make_histogram(self, xbar, all_heights, wd_names, file_prefix, bar_width, relative_text_height, xlabel, file_suffix, additional_text_dict, additional_line_dict, additional_x_axis_dict=None, relative_y_max=1.25):
        series_dict = dict()
        single_hist = len(all_heights) < 2
        max_height = 0
        use_hack = xlabel == 'Temperature /K'
        i = 0
        colours_to_use = self.rainbow if (use_hack and not single_hist) else self.colour_list
        for heights in all_heights:
            if np.max(heights) > max_height:
                max_height = np.max(heights)
            heights = list(heights)
            colour_index = i % len(colours_to_use)
            series_name = self.strip_system_suffix(wd_names[i])
            series_dict[series_name] = {
                'type': dp.PlotType.bar,
                'step': not single_hist,
                'x_data': xbar if single_hist else [xbar[0]-bar_width] + xbar + [xbar[-1]+bar_width],
                'y_data': heights if single_hist else [0] + heights + [0],
                'bar_width': bar_width,
                'colours': colours_to_use[colour_index],
                'edge_colours': colours_to_use[colour_index],
                #'face_alpha': 1 if single_hist else 0,
                'legend': True
            }
            i += 1
        text_height = max_height*relative_text_height
        if additional_text_dict is not None:
            for text_name, text_dict in additional_text_dict.items():
                if isinstance(text_dict, dict):
                    series_dict[text_name] = {
                        'type': dp.PlotType.text,
                        'x_pos': text_dict['x_pos'],
                        'y_pos': text_height,
                        'text_string': text_dict['text_string'],
                        'horizontalalignment': text_dict.get('horizontalalignment'),
                        'verticalalignment': text_dict.get('verticalalignment'),
                        'rotation': text_dict.get('rotation'),
                        'fontsize': text_dict.get('fontsize')
                    }
        else:
            additional_text_dict = {'title_text': None}
        if additional_line_dict is not None:
            for line_name, line_dict in additional_line_dict.items():
                series_dict[line_name] = {
                    'type': dp.PlotType.vline,
                    'x_start': line_dict['x_start'],
                    'y_min': 0,
                    'y_max': max_height*relative_y_max,
                    'line_styles': 'dotted',
                    'colours': 'grey',
                    'line_linewidth': 2
                }
        hack = 10 if use_hack else 0.5
        plot_dict = {
            'hist_plot': {
                'show': False,
                'filenames': [file_prefix + file_suffix + '.pdf', file_prefix + file_suffix + '.png'],
                'dpi': 300,
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,         #frameon = False
                        'legend_loc': 'center right' if use_hack else 'best',
                        'legend_text_size': 18,
                        'xlabel_text': xlabel,
                        'xlabel_fontsize': 28,
                        'ylabel_text': 'Probability density',
                        'ylabel_fontsize': 28,
                        'x_tick_fontsize': 20,
                        'y_tick_fontsize': 20,
                        'font': 'STIXGeneral',
                        'x_min': xbar[0] - (hack*bar_width),
                        'x_max': xbar[-1] + (0.5*bar_width),
                        'y_min': 0,
                        'y_max': max_height*relative_y_max,
                        'series': series_dict,
                        'title_text': additional_text_dict.get('title_text'),
                        'title_fontsize': 36,
                    }
                }
            }
        }
        if additional_x_axis_dict is not None:
            plot_dict['hist_plot']['subplots']['subplot2'] = {
                'subplot_region': 111,
                'xlabel_text': additional_x_axis_dict['xlabel_text'],
                'x_tick_locations': additional_x_axis_dict['x_tick_locations'],#[0, 10, 20, 30, 45, 60],
                'x_tick_labels': additional_x_axis_dict['x_tick_labels'],#[5, 6, 8, 11, 124],
                'x_max': xbar[-1] + (0.5*bar_width),
                #'x_max': additional_x_axis_dict['x_max'],
                'xlabel_fontsize': 22,
                'x_tick_fontsize': 18,
                #'series': None,
                'twin_subplot': 'subplot1',
                'twin_on_x': True
            }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)     
        return plot_dict
    
    def plot_mix_pie_chart(self, fcf, mix_dict):
        series_dict = dict()
        x_data = list()
        labels = list()
        colour_list = list()
        els_we_care_about = [ci.Element.O, ci.Element.Mg, ci.Element.Si, ci.Element.Fe, ci.Element.Ni, ci.Element.Placeholder]
        for element in els_we_care_about:
            toappend = mix_dict.get(element, None)
            if toappend is not None:
                x_data.append(toappend)
                labels.append(str(element) if element != ci.Element.Placeholder else 'Other')
                colour_list.append(self.colour_dict[element])
        series_dict[str(fcf)] = {
            'type': dp.PlotType.pie,
            'x_data': x_data,
            'labels': None,
            'colours': colour_list,
            'legend': True
        }
        plot_dict = {
            'hist_plot': {
                'show': False,
                'filenames': ['pie_mix_' + str(fcf) + '.pdf'],
                'subplots': {
                    'pie1': {
                        'subplot_region': 111,
                        'legend': False,         #frameon = False
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
    
    def plot_observability(self, data_dict, x_element, y_element, lithophile_element, siderophile_element):
        series_dict = dict()
        linestyles = ['-', '--', ':']
        j = 0
        x_key = str(x_element) + '/Hx'
        y_key = str(y_element) + '/Hx'
        known_fcfs = list()
        known_Ps = list()
        for fcf_P_tuple in data_dict:
            fcf = fcf_P_tuple[0]
            P = fcf_P_tuple[1]
            if fcf not in known_fcfs:
                known_fcfs.append(fcf)
            if P not in known_Ps:
                known_Ps.append(P)
        known_fcfs.sort()
        known_Ps.sort()
        for fcf_P_tuple, data in data_dict.items():
            run_name = 'fcf = ' + str(fcf_P_tuple[0]) + ', P = ' + str(fcf_P_tuple[1]) + ' GPa'
            fcf_index = known_fcfs.index(fcf_P_tuple[0])
            P_index = known_Ps.index(fcf_P_tuple[1])
            series_dict[run_name] = {
                'type': dp.PlotType.scatter_2d,
                'x_data': data[x_key],
                'y_data': data[y_key],
                'legend': True,
                'line_color': 'g' if fcf_index == 0 else self.colour_list[fcf_index],
                'line_style': linestyles[P_index]
            }
        interpolate_between_fcfs = len(known_fcfs) > 1
        interpolate_between_pressure = len(known_Ps) > 1 and not interpolate_between_fcfs
        # Do some maths to calculate the fcf indicator. Firstly identify the beginning and end points:
        if interpolate_between_fcfs:
            start_x = data_dict[(known_fcfs[0], known_Ps[0])][x_key][0] - 3.5
            start_y = data_dict[(known_fcfs[0], known_Ps[0])][y_key][0] - 3.5
            end_x = data_dict[(known_fcfs[-1], known_Ps[0])][x_key][0] - 3.5
            end_y = data_dict[(known_fcfs[-1], known_Ps[0])][y_key][0] - 3.5
            anchor_values = known_fcfs
        elif interpolate_between_pressure:
            start_x = data_dict[(known_fcfs[0], known_Ps[0])][x_key][0] - 3
            start_y = data_dict[(known_fcfs[0], known_Ps[0])][y_key][0] - 3
            end_x = data_dict[(known_fcfs[0], known_Ps[-1])][x_key][0] - 3
            end_y = data_dict[(known_fcfs[0], known_Ps[-1])][y_key][0] - 3
            anchor_values = known_Ps
        perpendicular_vector = [end_y - start_y, start_x - end_x]
        sum_sq = 0
        for pv in perpendicular_vector:
            sum_sq += pv**2
        perpendicular_vector = perpendicular_vector/np.sqrt(sum_sq)
        # Now generate a bunch of points along the line:
        spacing = 0.25
        if interpolate_between_pressure:
            spacing = 1
        current_position = 0
        rel_tol_on_comparison = 0.0000000001
        while current_position <= 1:
            point_x = start_x + current_position*(end_x - start_x)
            point_y = start_y + current_position*(end_y - start_y)
            # Identify current x_el/y_el value:
            xelyel = point_x - point_y
            # Now find the nearest neighbouring actual xelyel values, and how far between them we are: (assume no ties)
            lower_neighbour = None
            upper_neighbour = None
            lower_av = None
            upper_av = None
            for av in anchor_values:
                if interpolate_between_fcfs:
                    test_xelyel = data_dict[(av, known_Ps[0])][x_key][0] - data_dict[(av, known_Ps[0])][y_key][0]
                elif interpolate_between_pressure:
                    test_xelyel = data_dict[(known_fcfs[0], av)][x_key][0] - data_dict[(known_fcfs[0], av)][y_key][0]
                if lower_av is None:
                    if test_xelyel < xelyel or math.isclose(test_xelyel, xelyel, rel_tol=rel_tol_on_comparison):  # Then this is a starting point
                        lower_av = av
                        lower_neighbour = test_xelyel
                else:
                    if (test_xelyel < xelyel or math.isclose(test_xelyel, xelyel, rel_tol=rel_tol_on_comparison)) and test_xelyel > lower_neighbour: # Then this is a tighter bound
                        lower_av = av
                        lower_neighbour = test_xelyel
                if upper_av is None:
                    if test_xelyel > xelyel or math.isclose(test_xelyel, xelyel, rel_tol=rel_tol_on_comparison):  # Then this is a starting point
                        upper_av = av
                        upper_neighbour = test_xelyel
                else:
                    if (test_xelyel > xelyel or math.isclose(test_xelyel, xelyel, rel_tol=rel_tol_on_comparison)) and test_xelyel < upper_neighbour: # Then this is a tighter bound
                        upper_av = av
                        upper_neighbour = test_xelyel
            if math.isclose(upper_neighbour, lower_neighbour, rel_tol=rel_tol_on_comparison):
                distance_to_interpolate = 0
            else:
                distance_to_interpolate = (xelyel - lower_neighbour) / (upper_neighbour - lower_neighbour)
            # Now interpolate!
            lith_key = str(lithophile_element) + '/Hx'
            sid_key = str(siderophile_element) + '/Hx'
            if interpolate_between_fcfs:
                lower_LithSid = data_dict[(lower_av, known_Ps[0])][lith_key][0] - data_dict[(lower_av, known_Ps[0])][sid_key][0]
                upper_LithSid = data_dict[(upper_av, known_Ps[0])][lith_key][0] - data_dict[(upper_av, known_Ps[0])][sid_key][0]
            elif interpolate_between_pressure:
                lower_LithSid = data_dict[(known_fcfs[0], lower_av)][lith_key][0] - data_dict[(known_fcfs[0], lower_av)][sid_key][0]
                upper_LithSid = data_dict[(known_fcfs[0], upper_av)][lith_key][0] - data_dict[(known_fcfs[0], upper_av)][sid_key][0]
            interpolated_LithSid = lower_LithSid + distance_to_interpolate*(upper_LithSid - lower_LithSid)
            series_dict['fcf text ' +  str(current_position)] = {
                'type': dp.PlotType.text,
                'x_pos': point_x + (0.2*perpendicular_vector[0]),
                'y_pos': point_y + (0.2*perpendicular_vector[1]),
                'text_string': str(round(interpolated_LithSid, 2)),
                'legend': False
            }
            current_position += spacing
        series_dict['log(' + str(lithophile_element) + '/' + str(siderophile_element) + ')'] = {
            'type': dp.PlotType.scatter_2d,
            'x_data': [start_x, end_x],
            'y_data': [start_y, end_y],
            'line_color': 'k',
            'line_style': '-',
            'line_width': 10,
            'legend': True
        }
        plot_dict = {
            'observability_plot': {
                'show': False,
                'filenames': ['observability_' + str(x_element)  + '_v_' + str(y_element) + '.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        #'title_text': 'Observability Test Plot',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'log(' + x_key + ')',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'log(' + y_key + ')',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'font': 'STIXGeneral',
                        'x_min': -8,
                        'x_max': -2,
                        #'y_min': -0.015,
                        'y_max': -2,
                        #'x_scale': 'log',
                        'series': series_dict
                    }
                }
            }
        }
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
        return plot_dict
        
    def plot_fcf_v_delta_time(self, fcfs, fcfs_upper_errors, fcfs_lower_errors, delta_times, delta_times_upper_errors, delta_times_lower_errors, t_Mgs, wd_types, scaled=False):
        H_delta_times = list()
        H_delta_times_lower_errors = list()
        H_delta_times_upper_errors = list()
        H_fcfs = list()
        H_fcfs_lower_errors = list()
        H_fcfs_upper_errors = list()
        He_delta_times = list()
        He_delta_times_lower_errors = list()
        He_delta_times_upper_errors = list()
        He_fcfs = list()
        He_fcfs_lower_errors = list()
        He_fcfs_upper_errors = list()
        for i, dt in enumerate(delta_times):
            if scaled:
                scaling_factor = t_Mgs[i]
            else:
                scaling_factor = 1
            if wd_types[i] == ci.Element.H:
                H_delta_times.append(dt/scaling_factor)
                H_delta_times_lower_errors.append(delta_times_lower_errors[i]/scaling_factor)
                H_delta_times_upper_errors.append(delta_times_upper_errors[i]/scaling_factor)
                H_fcfs.append(fcfs[i])
                H_fcfs_lower_errors.append(fcfs_lower_errors[i])
                H_fcfs_upper_errors.append(fcfs_upper_errors[i])
            elif wd_types[i] == ci.Element.He:
                He_delta_times.append(dt/scaling_factor)
                He_delta_times_lower_errors.append(delta_times_lower_errors[i]/scaling_factor)
                He_delta_times_upper_errors.append(delta_times_upper_errors[i]/scaling_factor)
                He_fcfs.append(fcfs[i])
                He_fcfs_lower_errors.append(fcfs_lower_errors[i])
                He_fcfs_upper_errors.append(fcfs_upper_errors[i])
            else:
                raise
        series_dict = dict()
        if not scaled:
            series_dict['H'] = {
                'type': dp.PlotType.scatter_2d_error,
                'x_data': H_delta_times,
                'y_data': H_fcfs,
                'x_error_data': [H_delta_times_lower_errors, H_delta_times_upper_errors],
                'y_error_data': [H_fcfs_lower_errors, H_fcfs_upper_errors],
                'line_color': 'r',
                'line_style': 'None',
                'line_width': 10,
                'legend': True,
                'line_marker': '+',
                'line_markersize': 5
            }
        series_dict['He'] = {
            'type': dp.PlotType.scatter_2d_error,
            'x_data': He_delta_times,
            'y_data': He_fcfs,
            'x_error_data': [He_delta_times_lower_errors, He_delta_times_upper_errors],
            'y_error_data': [He_fcfs_lower_errors, He_fcfs_upper_errors],
            'line_color': 'k',
            'line_style': 'None',
            'line_width': 10,
            'legend': True,
            'line_marker': '+',
            'line_markersize': 5
        }
        plot_dict = {
            'fcf_v_delta_time_plot': {
                'show': False,
                'filenames': ['fcf_v_delta_time.pdf'],
                'subplots': {
                    'subplot1': {
                        'subplot_region': 111,
                        'legend': True,
                        'legend_loc': 'best',
                        'legend_text_size': 8,
                        'title_text': 'Fragment Core Fraction v Delta Time',
                        'title_fontsize': 12,
                        'title_fontweight': 'bold',
                        'xlabel_text': 'Delta Time/$t_{Mg}$' if scaled else 'Delta Time',
                        'xlabel_fontsize': 10,
                        'xlabel_fontweight': 'bold',
                        'ylabel_text': 'Fragment Core Fraction',
                        'ylabel_fontsize': 10,
                        'ylabel_fontweight': 'bold',
                        'font': 'STIXGeneral',
                        #'x_min': -8,
                        #'x_max': 0.05,
                        'y_min': 0,
                        'y_max': 1,
                        #'x_scale': 'log',
                        'series': series_dict
                    }
                }
            }
        }
        if scaled:
            plot_dict['fcf_v_delta_time_plot']['subplots']['subplot1']['x_scale'] = 'symlog'
            plot_dict['fcf_v_delta_time_plot']['subplots']['subplot1']['x_max'] = 0.00002
            plot_dict['fcf_v_delta_time_plot']['subplots']['subplot1']['x_min'] = -0.0000025
        plotter = dp.DictPlotter(plot_dict)
        plotter.draw()
        plotter.yield_output(self.output_dir)
