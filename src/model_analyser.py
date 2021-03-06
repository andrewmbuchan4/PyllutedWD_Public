#!/usr/bin/env python
# -*- coding: utf-8 -*-

import corner
import csv
import matplotlib.pyplot as plt
import numpy as np
import pymultinest as pn
from scipy.special import erfcinv
from scipy.special import lambertw as W
from scipy import stats as stat
import os

import abundance_model as am
import chemistry_info as ci
import complete_model as cm
import disc_model as dm
import enhancement_model as em
import geology_info as gi
import graph_factory as gf
import live_data as ld
import model_parameters as mp
import pollution_model as pm
import pwd_utils as pu

class ModelAnalyser:
    
    def __init__(self, graph_dir):
        self.graph_dir = graph_dir
        self.graph_fac = gf.GraphFactory(self.graph_dir)
    
    def confidence_intervals(self, sample_draws, array, length):
        prob = np.empty(sample_draws)
        prob.fill(1.0/sample_draws)   # Each sample assumed to be equally likely
            
        sig_1 = 0.5 + 0.6826/2.0  # Cummulative propabiltiy of +1 sigma level
        sig_2 = 0.5 + 0.954/2.0   # Cummulative propabiltiy of +2 sigma level
        sig_3 = 0.5 + 0.997/2.0   # Cummulative propabiltiy of +3 sigma level
            
        if (length > 0):  # If array is a vector (e.g. a line / curve where we want contours for each x value)
                
            arr_low3 = np.zeros(shape=(length))
            arr_low2 = np.zeros(shape=(length))
            arr_low1 = np.zeros(shape=(length))
            arr_median = np.zeros(shape=(length))
            arr_high1 = np.zeros(shape=(length))
            arr_high2 = np.zeros(shape=(length))
            arr_high3 = np.zeros(shape=(length))
                
            for i in range(length):
            
                arr_ordered = list(zip(prob[:], array[:, i]))
                arr_ordered.sort(key=lambda x: x[1])
                arr_ordered = np.array(arr_ordered)
        
                arr_ordered[:,0] = arr_ordered[:,0].cumsum()
        
                arr_ordered_interp = lambda x: np.interp(x, arr_ordered[:,0], arr_ordered[:,1],
                                                         left=arr_ordered[0,1], right=arr_ordered[-1,1])
        
                arr_low3[i] = arr_ordered_interp(1-sig_3)
                arr_low2[i] = arr_ordered_interp(1-sig_2)
                arr_low1[i] = arr_ordered_interp(1-sig_1)
                arr_median[i] = arr_ordered_interp(0.5)
                arr_high1[i] = arr_ordered_interp(sig_1)
                arr_high2[i] = arr_ordered_interp(sig_2) 
                arr_high3[i] = arr_ordered_interp(sig_3) 
                
            return arr_low3, arr_low2, arr_low1, arr_median, arr_high1, arr_high2, arr_high3
                
        if (length == 0):  # If quantity is just a float
                
            arr_ordered = list(zip(prob[:], array[:]))
            arr_ordered.sort(key=lambda x: x[1])
            arr_ordered = np.array(arr_ordered)
        
            arr_ordered[:,0] = arr_ordered[:,0].cumsum()
        
            arr_ordered_interp = lambda x: np.interp(x, arr_ordered[:,0], arr_ordered[:,1],
                                                     left=arr_ordered[0,1], right=arr_ordered[-1,1])
        
            arr_low3 = arr_ordered_interp(1-sig_3)
            arr_low2 = arr_ordered_interp(1-sig_2)
            arr_low1 = arr_ordered_interp(1-sig_1)
            arr_median = arr_ordered_interp(0.5)
            arr_high1 = arr_ordered_interp(sig_1)
            arr_high2 = arr_ordered_interp(sig_2) 
            arr_high3 = arr_ordered_interp(sig_3) 
                    
            return arr_low3, arr_low2, arr_low1, arr_median, arr_high1, arr_high2, arr_high3   

    def z_to_sigma(self, ln_Z1, ln_Z2):  
        np.set_printoptions(precision=50)
        B = np.exp(ln_Z1 - ln_Z2) 
        p = np.real(np.exp(W((-1.0/(B*np.exp(1))),-1)))
        sigma = np.sqrt(2)*erfcinv(p)
        return B, sigma  
    
    def detection_significance(self, additional_parameter, err_data, n_params, n_removed, basename_M1, basename_M2):
    
        norm_log = (-0.5*np.log(2.0*np.pi*err_data*err_data)).sum()
        
        retrieved_M1 = pn.Analyzer(n_params = n_params, outputfiles_basename = basename_M1)
        retrieved_M2 = pn.Analyzer(n_params = (n_params-n_removed), outputfiles_basename = basename_M2)
        s_M1 = retrieved_M1.get_stats()
        s_M2 = retrieved_M2.get_stats()
        
        best_fit_M1 = retrieved_M1.get_best_fit()                       # Best-fitting parameter combination from model 1
        max_likelihood_M1 = best_fit_M1['log_likelihood']               # Maximum log-likelihood for model 1
        best_chi_square_M1 = -2.0 * (max_likelihood_M1 - norm_log)      # Best chi-square for model 1
        
        ln_Z_M1 = s_M1['global evidence']     # Natural log of Bayesian evidence from model 1 
        
        print("Evidence with " + additional_parameter + " = " + str(ln_Z_M1))
                    
        best_fit_M2 = retrieved_M2.get_best_fit()                                      # Best-fitting parameter combination from model 2
        max_likelihood_M2 = best_fit_M2['log_likelihood']                              # Maximum log-likelihood for model 2
        best_chi_square_M2 = -2.0 * (max_likelihood_M2 - norm_log)                     # Best chi-square for model 2
        
        ln_Z_M2 = s_M2['global evidence']     # Natural log of Bayesian evidence from model 1 
        
        print("Evidence without " + additional_parameter + " = " + str(ln_Z_M2))
                    
        Bayes_factor, n_sigma = self.z_to_sigma(ln_Z_M1, ln_Z_M2)
        
        print("Bayes_factor of " + additional_parameter + " = " + str(Bayes_factor))
        print("Sigma significance of " + additional_parameter + " = " + str(n_sigma))
        
        return ln_Z_M1, ln_Z_M2, Bayes_factor, n_sigma, best_chi_square_M1, best_chi_square_M2   

    def compare(self, model_dict, system, number_of_data_points):
        # model_dict should be an OrderedDict with the base model in first position. Comparison of base to itself is to get ln_Z
        base_model_name = None
        best_model_so_far_name = None
        for model_name, model in model_dict.items():
            if base_model_name is None:
                base_model_name = model_name
            else:
                base_model_name = best_model_so_far_name
            base_model = model_dict[base_model_name]
            
            self.compare_two_models(model, base_model, system, number_of_data_points)
            
            if best_model_so_far_name is None:
                best_model_so_far_name = base_model_name
            else:
                ln_Z_model = model.comparison[base_model_name]['ln_Z_model']
                ln_Z_base = model.comparison[base_model_name]['ln_Z_base']
                if ln_Z_model > ln_Z_base:
                    best_model_so_far_name = model_name
    
    def compare_two_models(self, model, base_model, system, number_of_data_points):
        print('Comparing result for model ' + str(model.basename) + ' with base_model ' + str(base_model.basename))
        good_fit_threshold = 2
        ln_Z_model, ln_Z_base, Bayes_factor_model_base, n_sigma_model_base, chi_model, chi_base = self.detection_significance(
            'Additional parameter name goes here for readability',
            ld._live_non_zero_wd_errors,
            model.get_n_dims(),
            (model.get_n_dims() - base_model.get_n_dims()),
            model.get_full_prefix(system),
            base_model.get_full_prefix(system)
        )
        model.comparison[base_model.basename] = {
            'ln_Z_model': ln_Z_model,
            'ln_Z_base': ln_Z_base,
            'Bayes_factor_model_base': Bayes_factor_model_base,
            'n_sigma_model_base': n_sigma_model_base,
            'chi_model': chi_model,
            'chi_base': chi_base,
            'chi_model_per_data_point': chi_model/number_of_data_points,
            'chi_base_per_data_point': chi_base/number_of_data_points,
            'good_fit': chi_model/number_of_data_points < good_fit_threshold
        }
    
    def find_diff_sigma(self, models, system, number_of_data_points, stats_file):
        # First need to find best non-differentiated system
        fcf_string = mp.model_parameter_strings[mp.ModelParameter.fragment_core_frac]
        best_non_diff_model_name = None
        best_non_diff_model_lnZ = None
        for model_name, model in models.items():
            if model.prior_name == 'Default' and fcf_string not in model.get_model_params():
                if best_non_diff_model_name is None:
                    best_non_diff_model_name = model_name
                    if len(model.comparison.items()) == 0:
                        self.compare_two_models(model, model, system, number_of_data_points)  # Compare model to itself - we just need the ln Z so actual comparison is arbitrary
                    assert len(model.comparison.items()) > 0
                    for dummy_name, arbitrary_comparison in model.comparison.items():
                        best_non_diff_model_lnZ = arbitrary_comparison['ln_Z_model']
                else:
                    if len(model.comparison.items()) == 0:
                        self.compare_two_models(model, model, system, number_of_data_points)  # Compare model to itself - we just need the ln Z so actual comparison is arbitrary
                    assert len(model.comparison.items()) > 0
                    test_model_lnZ = None
                    for dummy_name, arbitrary_comparison in model.comparison.items():
                        test_model_lnZ = arbitrary_comparison['ln_Z_model']
                    if test_model_lnZ > best_non_diff_model_lnZ:
                        best_non_diff_model_lnZ = test_model_lnZ
                        best_non_diff_model_name = model_name
        if best_non_diff_model_name is None:
            # Then there was no model which actually omitted fcf, so no point continuing
            print('Could not find a non-differentiated model, no diff_sigma could be calculated')
            return
        else:
            for model_name, model in models.items():
                if model_name == best_non_diff_model_name:
                    model.best_nondiff_model = True
        # Now find best model:
        best_model_name = None
        for model_name, model in models.items():
            if model.prior_name == 'Default' and model.best_model:
                best_model_name = model_name
        # Now compare them if necessary
        if best_model_name != best_non_diff_model_name:
            if best_non_diff_model_name not in models[best_model_name].comparison:
                self.compare_two_models(models[best_model_name], models[best_non_diff_model_name], system, number_of_data_points)
                with open(stats_file, 'a', newline='', encoding='utf-8') as f:
                    to_write = csv.writer(f)
                    to_write.writerow([
                        best_model_name,
                        best_non_diff_model_name,
                        models[best_model_name].comparison[best_non_diff_model_name]['ln_Z_model'],
                        models[best_model_name].comparison[best_non_diff_model_name]['ln_Z_base'],
                        models[best_model_name].comparison[best_non_diff_model_name]['Bayes_factor_model_base'],
                        models[best_model_name].comparison[best_non_diff_model_name]['n_sigma_model_base'],
                        models[best_model_name].comparison[best_non_diff_model_name]['chi_model'],
                        models[best_model_name].comparison[best_non_diff_model_name]['chi_base'],
                        models[best_model_name].comparison[best_non_diff_model_name]['chi_model_per_data_point'],
                        models[best_model_name].comparison[best_non_diff_model_name]['chi_base_per_data_point'],
                        models[best_model_name].get_model_params(),
                        'N/A' if models[best_model_name].best_model is None else models[best_model_name].best_model,
                        models[best_model_name].comparison[best_non_diff_model_name]['good_fit']
                    ])
            diff_sigma_str = str(models[best_model_name].comparison[best_non_diff_model_name]['n_sigma_model_base'])
            bayes_factor_str = str(models[best_model_name].comparison[best_non_diff_model_name]['Bayes_factor_model_base'])
            chi_sq_nondiff_str = str(models[best_model_name].comparison[best_non_diff_model_name]['chi_base_per_data_point'])
        else:
            diff_sigma_str = 'N/A'
            bayes_factor_str = 'N/A'
            chi_sq_nondiff_str = 'N/A'
        with open(stats_file, 'a', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow([
                'Best model name:',
                best_model_name
            ])
            to_write.writerow([
                'Best non-differentiated model name:',
                best_non_diff_model_name
            ])
            to_write.writerow([
                'Diff Sigma:',
                diff_sigma_str
            ])
            to_write.writerow([
                'Corresponding Bayes Factor:',
                bayes_factor_str
            ])
            to_write.writerow([
                'Chi squared per data point (best non-diff):',
                chi_sq_nondiff_str
            ])
    
    def dump_model_stats(self, model_name, model, stats_file):
        print(model.comparison)
        with open(stats_file, 'a', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            for base_model_name, comp_vals in model.comparison.items():
                to_write.writerow([
                    model_name,
                    base_model_name,
                    comp_vals['ln_Z_model'],
                    comp_vals['ln_Z_base'],
                    comp_vals['Bayes_factor_model_base'],
                    comp_vals['n_sigma_model_base'],
                    comp_vals['chi_model'],
                    comp_vals['chi_base'],
                    comp_vals['chi_model_per_data_point'],
                    comp_vals['chi_base_per_data_point'],
                    model.get_model_params(),
                    'N/A' if model.best_model is None else model.best_model,
                    comp_vals['good_fit']
                ])
    
    def make_all_plots(self, model_name, model, observation_number, all_wd_timescales, wd_name, all_wd_abundances, all_wd_upper_bounds, all_wd_lower_bounds, all_wd_abundance_errors, enhancement_model, wd_type):
        self.make_corner_plot(model_name, model, observation_number, wd_name)
        self.make_composition_plot_mk2(model_name, model, observation_number, all_wd_timescales, wd_name, all_wd_abundances, all_wd_upper_bounds, all_wd_lower_bounds, all_wd_abundance_errors, enhancement_model, wd_type)
        self.make_time_plot(model, model_name, wd_name, observation_number)
        self.make_p_v_fO2_plot(model, model_name, wd_name, observation_number)
        temp_stats = self.make_temperature_plot(model, model_name, wd_name, observation_number)
        mass_stats = self.make_mass_plot(model, model_name, wd_name, observation_number, enhancement_model)
        self.make_pressure_plot(model, model_name, wd_name, observation_number)
        return temp_stats, mass_stats
        
    def dump_model_fit(self, model, N_wd, stats_file, temp_stats, mass_stats):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, N_wd)
        parameter_indices = mp.parameter_indices(ld._live_model)
        t_sinceaccretion_index = parameter_indices[mp.ModelParameter.t_sinceaccretion]
        accretion_timescale_index = parameter_indices[mp.ModelParameter.accretion_timescale]
        medians = list()
        percentile_5 = list()
        percentile_16 = list()
        percentile_84 = list()
        percentile_95 = list()
        modes = list()
        for index, param in enumerate(model.get_model_params()):
            medians.append(np.percentile(weightpost[:,index], 50))
            percentile_5.append(np.percentile(weightpost[:,index], 5))
            percentile_16.append(np.percentile(weightpost[:,index], 16))
            percentile_84.append(np.percentile(weightpost[:,index], 84))
            percentile_95.append(np.percentile(weightpost[:,index], 95))
            modes.append(np.asscalar(stat.mode(np.around(weightpost[:,index], decimals=2))[0]))
        print('medians')
        print(medians)
        input_values = list()
        for param in mp.get_model_params_in_order():
            if mp.model_parameter_strings[param] in model.get_model_params():
                index = model.get_model_params().index(mp.model_parameter_strings[param])
                to_append = medians[index]
            else:
                to_append = mp.default_values[ld._enhancement_model][param]
            input_values.append(to_append)
        print('ma221 believes the final medians are: ')
        print(input_values)
        print('Also assuming t_formation=1.5 and normalise_abundances=True...')
        model_result, diagnostics = cm.complete_model_calculation(
            input_values[0],
            input_values[1],
            input_values[2],
            input_values[3],
            input_values[4],
            input_values[5],
            input_values[6],
            input_values[7],
            input_values[8],
            10**input_values[9],
            input_values[10],
            input_values[11],
            model.enhancement_model_name
        )
        print('output')
        print(model_result)
        print(diagnostics)
        try:
            elements_list = list()
            for el in model_result.keys():
                elements_list.append(el)
        except AttributeError:
            # Model could potentially be None
            print('Warning! Median model gave None result!')
            elements_list = ci.usual_elements
        elements_results = list()
        disc_comp = list()
        bulk_comp = list()
        core_comp = list()
        mantle_comp = list()
        ds = list()
        for el in elements_list:
            try:
                elements_results.append(model_result[el])
            except TypeError:
                pass
            try:
                disc_comp.append(diagnostics['DiscAbundances'][el])
            except (TypeError, KeyError):
                pass
            try:
                bulk_comp.append(diagnostics['Enhancements']['Abundances'][el][gi.Layer.bulk])
            except (TypeError, KeyError):
                pass
            try:
                core_comp.append(diagnostics['Enhancements']['Abundances'][el][gi.Layer.core])
            except (TypeError, KeyError):
                pass
            try:
                mantle_comp.append(diagnostics['Enhancements']['Abundances'][el][gi.Layer.mantle])
            except (TypeError, KeyError):
                pass
            try:
                ds.append(diagnostics['Enhancements']['Ds'][el])
            except (TypeError, KeyError):
                pass
        radius = None
        mass = None
        radius_lower_limit = None
        mass_lower_limit = None
        if mp.model_parameter_strings[mp.ModelParameter.pressure] in model.get_model_params():
            geo_model = gi.GeologyModel()  # Argument here irrelevent, we will supply the needed abundances
            radius_mass = geo_model.calculate_radius_and_mass(input_values[10], diagnostics['Enhancements']['Abundances'])
            radius = radius_mass[0]
            mass = radius_mass[1]
            radius_mass_lower_limit = geo_model.calculate_radius_and_mass(input_values[10], diagnostics['Enhancements']['Abundances'], True)
            radius_lower_limit = radius_mass_lower_limit[0]
            mass_lower_limit = radius_mass_lower_limit[1]
            
        accretion_timescale_extract = 10**(weightpost[:,accretion_timescale_index] - 6)
        t_sinceaccretion_extract = weightpost[:,t_sinceaccretion_index]

        tmedian, terrorplus, terrorminus = self.extract_median_and_error(accretion_timescale_extract - t_sinceaccretion_extract)
        
        Tmedian, Terrorplus, Terrorminus = self.extract_median_and_error(temp_stats)

        Mmedian, Merrorplus, Merrorminus = self.extract_median_and_error(mass_stats)
        
        pressure_stats = self.get_untransformed_parameter_stats(model, N_wd, mp.ModelParameter.pressure)
        Pmedian, Perrorplus, Perrorminus = self.extract_median_and_error(pressure_stats)
        P_10_percentile = self.extract_xth_percentile(pressure_stats, 10)
        P_90_percentile = self.extract_xth_percentile(pressure_stats, 90)
        P_percentile_of_5 = self.extract_percentile_of_x(pressure_stats, 5)
        P_percentile_of_10 = self.extract_percentile_of_x(pressure_stats, 10)
        P_percentile_of_15 = self.extract_percentile_of_x(pressure_stats, 15)
        P_percentile_of_45 = self.extract_percentile_of_x(pressure_stats, 45)
        P_percentile_of_50 = self.extract_percentile_of_x(pressure_stats, 50)
        P_percentile_of_55 = self.extract_percentile_of_x(pressure_stats, 55)
        
        with open(stats_file, 'a', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow([])
            to_write.writerow(['Parameter:'] + model.get_model_params())
            to_write.writerow(['5th percentile:'] + percentile_5)
            to_write.writerow(['16th percentile:'] + percentile_16)
            to_write.writerow(['Median:'] + medians)
            to_write.writerow(['84th percentile:'] + percentile_84)
            to_write.writerow(['95th percentile:'] + percentile_95)
            to_write.writerow(['Mode:'] + modes)
            to_write.writerow([])
            to_write.writerow(['Output from median run:'])
            to_write.writerow(['Elements:'] + elements_list)
            to_write.writerow(['Disc Composition:'] + disc_comp)
            #to_write.writerow(['Partition Coefficients (X_c for Oxygen):'] + ds) # For old Oxygen logic
            to_write.writerow(['Partition Coefficients:'] + ds)
            to_write.writerow(['Bulk Composition:'] + bulk_comp)
            to_write.writerow(['Core Composition:'] + core_comp)
            to_write.writerow(['Mantle Composition:'] + mantle_comp)
            to_write.writerow(['Final Result:'] + elements_results)
            try:
                to_write.writerow(['Parent Core Number Fraction:', diagnostics['Enhancements']['ParentCoreNumberFraction']])
            except TypeError:
                to_write.writerow(['Parent Core Number Fraction:'])
            to_write.writerow(['Radius /km:', radius])
            to_write.writerow(['Mass /M_Earth:', mass])
            to_write.writerow(['Radius lower limit /km:', radius_lower_limit])
            to_write.writerow(['Mass lower limit /M_Earth:', mass_lower_limit])
            to_write.writerow(['delta time/Myrs, +error, -error:', tmedian, terrorplus, terrorminus])
            to_write.writerow(['Temperature /K, +error, -error:', Tmedian, Terrorplus, Terrorminus])
            to_write.writerow(['log(Mass /kg), +error, -error:', Mmedian, Merrorplus, Merrorminus])
            to_write.writerow(['Pressure /GPa, +error, -error:', Pmedian, Perrorplus, Perrorminus])
            to_write.writerow(['90% chance of Pressure above:', P_10_percentile])
            to_write.writerow(['90% chance of Pressure below:', P_90_percentile])
            to_write.writerow(['Percentile of 5 GPa:', P_percentile_of_5])
            to_write.writerow(['Percentile of 10 GPa:', P_percentile_of_10])
            to_write.writerow(['Percentile of 15 GPa:', P_percentile_of_15])
            to_write.writerow(['Percentile of 45 GPa:', P_percentile_of_45])
            to_write.writerow(['Percentile of 50 GPa:', P_percentile_of_50])
            to_write.writerow(['Percentile of 55 GPa:', P_percentile_of_55])
    
    def extract_median_and_error(self, stats_object):
        if stats_object is None:
            return None, None, None
        median = np.percentile(stats_object, 50)
        errorplus = np.percentile(stats_object, 84) - median
        errorminus = median - np.percentile(stats_object, 16)
        return median, errorplus, errorminus
        
    def extract_xth_percentile(self, stats_object, percentile):
        if stats_object is None:
            return None
        return np.percentile(stats_object, percentile)
        
    def extract_percentile_of_x(self, stats_object, value):
        if stats_object is None:
            return None
        count = 0
        for num in stats_object:
            if num < value:
                count += 1
        return count / len(stats_object)
    
    def get_stats_weightpost_and_best_fit(self, model, observation_number):
        a = pn.Analyzer(n_params = model.get_n_dims(), outputfiles_basename = model.get_full_prefix(observation_number))
        stats = a.get_stats()
        weightpost = a.get_equal_weighted_posterior()[:, 0:model.get_n_dims()]  # This is basically just excluding the final column of the ...post_equal_weights.dat file
        best_fit = a.get_best_fit()
        return stats, weightpost, best_fit['log_likelihood'], best_fit['parameters']

    def make_corner_plot(self, model_name, model, observation_number, wd_name):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number)
        delta_free_params = list()
        for dfp in model.get_model_params():
            if '??' in dfp:
                # TODO: See if replacing it with r'$\Delta$' solves it
                delta_free_params.append(dfp.replace('??', ''))  # The corner function doesn't handle this character when plt.rcParams['text.usetex'] = 'True' in dict_plotter
            else:
                delta_free_params.append(dfp)
        corner.corner(weightpost, labels=delta_free_params, label_kwargs=dict(fontsize=12), levels = (0.39346934,0.86466472,0.988891), smooth=True)
        prefix = self.graph_dir + '/' + wd_name + '_' + model.get_prefix(observation_number)
        plt.savefig(prefix + 'corner.pdf')  # TODO: do this within DictPlotter

    def make_composition_plot_mk2(self, model_name, model, observation_number, all_wd_timescales, wd_name, all_wd_abundances, all_wd_upper_bounds, all_wd_lower_bounds, all_wd_abundance_errors, enhancement_model, wd_type):
        fixed_pressure_version = False # Set this to true to produce a version of this plot with an extra line corresponding to a certain fixed pressure
        fixed_pressure = 60
        fixed_pressure_name = 'High Pressure'
        
        include_best_fit = False
        
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number)

        number_of_possible_samples = len(weightpost[:,0])
        max_number_of_samples = 10000
        number_of_samples = min(max_number_of_samples, number_of_possible_samples)
        samples = np.random.choice(number_of_possible_samples, number_of_samples, replace=False)
        len_x_axis = len(ci.usual_elements) - 1
        model_output = np.zeros(shape=(number_of_samples, len_x_axis))
        model_raw_output = np.zeros(shape=(number_of_samples, len_x_axis+1))
        input_values = list()
        parameter_indices = mp.parameter_indices(ld._live_model)
        geo_model_for_solar_abundances = gi.GeologyModel()  # We're only using this to access the solar abundances so no arguments needed
        
        if fixed_pressure_version:
            fp_model_output = np.zeros(shape=(number_of_samples, len_x_axis))
            fp_model_raw_output = np.zeros(shape=(number_of_samples, len_x_axis + 1))
        print('Sampling for Composition plot')
        number_of_samples_used = 0
        for i in range(number_of_samples):
            if (i % 250) == 0:
                print('Composition sample ' + str(i))
            fe_star = weightpost[samples[i], parameter_indices[mp.ModelParameter.metallicity]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.metallicity) else mp.default_values[enhancement_model][mp.ModelParameter.metallicity]
            t_sinceaccretion = weightpost[samples[i], parameter_indices[mp.ModelParameter.t_sinceaccretion]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.t_sinceaccretion) else mp.default_values[enhancement_model][mp.ModelParameter.t_sinceaccretion]
            d_formation = weightpost[samples[i], parameter_indices[mp.ModelParameter.formation_distance]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.formation_distance) else mp.default_values[enhancement_model][mp.ModelParameter.formation_distance]
            t_formation = 1.5
            z_formation = weightpost[samples[i], parameter_indices[mp.ModelParameter.feeding_zone_size]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.feeding_zone_size) else mp.default_values[enhancement_model][mp.ModelParameter.feeding_zone_size]
            N_c = weightpost[samples[i], parameter_indices[mp.ModelParameter.parent_core_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.parent_core_frac) else mp.default_values[enhancement_model][mp.ModelParameter.parent_core_frac]
            N_o = weightpost[samples[i], parameter_indices[mp.ModelParameter.parent_crust_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.parent_crust_frac) else mp.default_values[enhancement_model][mp.ModelParameter.parent_crust_frac]
            f_c = weightpost[samples[i], parameter_indices[mp.ModelParameter.fragment_core_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.fragment_core_frac) else mp.default_values[enhancement_model][mp.ModelParameter.fragment_core_frac]
            f_o = weightpost[samples[i], parameter_indices[mp.ModelParameter.fragment_crust_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.fragment_crust_frac) else mp.default_values[enhancement_model][mp.ModelParameter.fragment_crust_frac]
            pollutionfraction = weightpost[samples[i], parameter_indices[mp.ModelParameter.pollution_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.pollution_frac) else mp.default_values[enhancement_model][mp.ModelParameter.pollution_frac]
            log_t_disc = weightpost[samples[i], parameter_indices[mp.ModelParameter.accretion_timescale]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.accretion_timescale) else mp.default_values[enhancement_model][mp.ModelParameter.accretion_timescale]
            pressure = weightpost[samples[i], parameter_indices[mp.ModelParameter.pressure]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.pressure) else mp.default_values[enhancement_model][mp.ModelParameter.pressure]
            fO2 = weightpost[samples[i], parameter_indices[mp.ModelParameter.oxygen_fugacity]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.oxygen_fugacity) else mp.default_values[enhancement_model][mp.ModelParameter.oxygen_fugacity]
            expressions, ignore = cm.complete_model_calculation(
                fe_star,
                t_sinceaccretion,
                d_formation,
                z_formation,
                N_c,
                N_o,
                f_c,
                f_o,
                pollutionfraction,
                10**log_t_disc,
                pressure,
                fO2,
                enhancement_model,
                t_formation
            )
            if fixed_pressure_version:
                fp_expressions, fp_ignore = cm.complete_model_calculation(
                    fe_star,
                    t_sinceaccretion,
                    d_formation,
                    z_formation,
                    N_c,
                    N_o,
                    f_c,
                    f_o,
                    pollutionfraction,
                    10**log_t_disc,
                    fixed_pressure,
                    fO2,
                    enhancement_model,
                    t_formation
                )
            if expressions is None or (fixed_pressure_version and fp_expressions is None):
                print('Warning! The sampled model does not converge! (Or otherwise returns None). Skipping this sample. Printing attempted params:')
                print('fe_star: ' + str(fe_star))
                print('t_sinceaccretion: ' + str(t_sinceaccretion))
                print('d_formation: ' + str(d_formation))
                print('z_formation: ' + str(z_formation))
                print('N_c: ' + str(N_c))
                print('N_o: ' + str(N_o))
                print('f_c: ' + str(f_c))
                print('f_o: ' + str(f_o))
                print('pollutionfraction: ' + str(pollutionfraction))
                print('log_t_disc: ' + str(log_t_disc))
                print('pressure: ' + str(pressure))
                print('fO2: ' + str(fO2))
                print('t_formation: ' + str(t_formation))
                continue
            number_of_samples_used += 1
            output = list()
            raw_output = list()
            for element in ci.usual_elements:
                if element is not ci.Element.Mg:
                    output.append(expressions[element] - expressions[ci.Element.Mg] - geo_model_for_solar_abundances.solar_abundances[element])
                raw_output.append(expressions[element] - (geo_model_for_solar_abundances.solar_ratiod_to_H[element] - geo_model_for_solar_abundances.solar_ratiod_to_H[wd_type]))
            model_output[i,:] = output
            model_raw_output[i,:] = raw_output
            if fixed_pressure_version:
                fp_output = list()
                fp_raw_output = list()
                for element in ci.usual_elements:
                    if element is not ci.Element.Mg:
                        fp_output.append(fp_expressions[element] - fp_expressions[ci.Element.Mg] - geo_model_for_solar_abundances.solar_abundances[element])
                    fp_raw_output.append(fp_expressions[element] - (geo_model_for_solar_abundances.solar_ratiod_to_H[element] - geo_model_for_solar_abundances.solar_ratiod_to_H[wd_type]))
                fp_model_output[i,:] = fp_output
                fp_model_raw_output[i,:] = fp_raw_output
        model_low3, model_low2, model_low1, model_median, model_high1, model_high2, model_high3 = self.confidence_intervals(number_of_samples_used, model_output, len_x_axis)
        fit_dict = {'Model median': model_median}
        error_low_dict = {'Model median': model_low1}
        error_high_dict = {'Model median': model_high1}
        raw_model_low3, raw_model_low2, raw_model_low1, raw_model_median, raw_model_high1, raw_model_high2, raw_model_high3 = self.confidence_intervals(number_of_samples_used, model_raw_output, len_x_axis + 1)
        raw_fit_dict = {'Model median': raw_model_median}
        raw_error_low_dict = {'Model median': raw_model_low1}
        raw_error_high_dict = {'Model median': raw_model_high1}
        if fixed_pressure_version:
            fp_model_low3, fp_model_low2, fp_model_low1, fp_model_median, fp_model_high1, fp_model_high2, fp_model_high3 = self.confidence_intervals(number_of_samples_used, fp_model_output, len_x_axis)
            fit_dict[fixed_pressure_name] = fp_model_median
            fp_raw_model_low3, fp_raw_model_low2, fp_raw_model_low1, fp_raw_model_median, fp_raw_model_high1, fp_raw_model_high2, fp_raw_model_high3 = self.confidence_intervals(number_of_samples_used, fp_model_raw_output, len_x_axis+1)
            raw_fit_dict[fixed_pressure_name] = fp_raw_model_median
        if include_best_fit:
            fe_star = best_fit_params[parameter_indices[mp.ModelParameter.metallicity]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.metallicity) else mp.default_values[enhancement_model][mp.ModelParameter.metallicity]
            t_sinceaccretion = best_fit_params[parameter_indices[mp.ModelParameter.t_sinceaccretion]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.t_sinceaccretion) else mp.default_values[enhancement_model][mp.ModelParameter.t_sinceaccretion]
            d_formation = best_fit_params[parameter_indices[mp.ModelParameter.formation_distance]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.formation_distance) else mp.default_values[enhancement_model][mp.ModelParameter.formation_distance]
            t_formation = 1.5
            z_formation = best_fit_params[parameter_indices[mp.ModelParameter.feeding_zone_size]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.feeding_zone_size) else mp.default_values[enhancement_model][mp.ModelParameter.feeding_zone_size]
            N_c = best_fit_params[parameter_indices[mp.ModelParameter.parent_core_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.parent_core_frac) else mp.default_values[enhancement_model][mp.ModelParameter.parent_core_frac]
            N_o = best_fit_params[parameter_indices[mp.ModelParameter.parent_crust_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.parent_crust_frac) else mp.default_values[enhancement_model][mp.ModelParameter.parent_crust_frac]
            f_c = best_fit_params[parameter_indices[mp.ModelParameter.fragment_core_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.fragment_core_frac) else mp.default_values[enhancement_model][mp.ModelParameter.fragment_core_frac]
            f_o = best_fit_params[parameter_indices[mp.ModelParameter.fragment_crust_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.fragment_crust_frac) else mp.default_values[enhancement_model][mp.ModelParameter.fragment_crust_frac]
            pollutionfraction = best_fit_params[parameter_indices[mp.ModelParameter.pollution_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.pollution_frac) else mp.default_values[enhancement_model][mp.ModelParameter.pollution_frac]
            log_t_disc = best_fit_params[parameter_indices[mp.ModelParameter.accretion_timescale]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.accretion_timescale) else mp.default_values[enhancement_model][mp.ModelParameter.accretion_timescale]
            pressure = best_fit_params[parameter_indices[mp.ModelParameter.pressure]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.pressure) else mp.default_values[enhancement_model][mp.ModelParameter.pressure]
            fO2 = best_fit_params[parameter_indices[mp.ModelParameter.oxygen_fugacity]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.oxygen_fugacity) else mp.default_values[enhancement_model][mp.ModelParameter.oxygen_fugacity]
            expressions, ignore = cm.complete_model_calculation(
                fe_star,
                t_sinceaccretion,
                d_formation,
                z_formation,
                N_c,
                N_o,
                f_c,
                f_o,
                pollutionfraction,
                10**log_t_disc,
                pressure,
                fO2,
                enhancement_model,
                t_formation
            )
            output = list()
            raw_output = list()
            for element in ci.usual_elements:
                if element is not ci.Element.Mg:
                    output.append(expressions[element] - expressions[ci.Element.Mg] - geo_model_for_solar_abundances.solar_abundances[element])
                raw_output.append(expressions[element] - (geo_model_for_solar_abundances.solar_ratiod_to_H[element] - geo_model_for_solar_abundances.solar_ratiod_to_H[wd_type]))
            fit_dict['Best fit'] = output
            raw_fit_dict['Best fit'] = raw_output
            if fixed_pressure_version:
                fp_expressions, fp_ignore = cm.complete_model_calculation(
                    fe_star,
                    t_sinceaccretion,
                    d_formation,
                    z_formation,
                    N_c,
                    N_o,
                    f_c,
                    f_o,
                    pollutionfraction,
                    10**log_t_disc,
                    fixed_pressure,
                    fO2,
                    enhancement_model,
                    t_formation
                )
                output = list()
                raw_output = list()
                for element in ci.usual_elements:
                    if element is not ci.Element.Mg:
                        output.append(fp_expressions[element] - fp_expressions[ci.Element.Mg] - geo_model_for_solar_abundances.solar_abundances[element])
                    raw_output.append(fp_expressions[element] - (geo_model_for_solar_abundances.solar_ratiod_to_H[element] - geo_model_for_solar_abundances.solar_ratiod_to_H[wd_type]))
                fit_dict[fixed_pressure_name + ' fit'] = output
                raw_fit_dict[fixed_pressure_name + ' fit'] = raw_output
        self.graph_fac.make_composition_plot_mk2(
            wd_name,
            model.get_prefix(observation_number),
            all_wd_abundances,
            all_wd_abundance_errors,
            fit_dict,
            error_low_dict,
            error_high_dict,
            all_wd_upper_bounds,
            all_wd_lower_bounds
        )
        self.graph_fac.make_composition_plot_raw(
            wd_name,
            wd_type,
            model.get_prefix(observation_number),
            all_wd_abundances,
            all_wd_abundance_errors,
            raw_fit_dict,
            raw_error_low_dict,
            raw_error_high_dict,
            all_wd_upper_bounds,
            all_wd_lower_bounds
        )

    def get_untransformed_parameter_stats(self, model, observation_number, model_param):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number)
        parameter_indices = mp.parameter_indices(ld._live_model)
        try:
            param_index = parameter_indices[model_param]
        except KeyError:
            # If the parameter was not present, this calculation is irrelevant
            return None
        toret = weightpost[:, param_index]
        return toret

    def make_pressure_plot(self, model, model_name, wd_name, observation_number):
        p_stats = self.get_untransformed_parameter_stats(model, observation_number, mp.ModelParameter.pressure)
        if p_stats is None:
            print('Abandoning pressure plot (no pressure stats available)')
            return
        
        half_bin_size = 0.5
        bins, x_bar_centres = self.generate_bins_and_bar_centres(0, 60, half_bin_size)
        
        heights, bins2 = np.histogram(
            p_stats,
            bins,
            density=True
        )
        geo_model = gi.GeologyModel()
        pressure_vals = [0, 18, 29.5, 39.5, 48.7, 57.5]
        mass_vals = list()
        for p in pressure_vals:
            if p == 0:
                mass_vals.append(0)
            else:
                rm_tuple = geo_model.calculate_radius_and_mass(p, geo_model.element_info)
                mass_vals.append(str(round(rm_tuple[1],2)).rstrip('0').rstrip('.'))
        additional_x_axis_dict = {
            'xlabel_text': 'Mass / M' + r'$_{\oplus}$',
            'x_tick_locations': pressure_vals,
            'x_tick_labels': mass_vals
            #'x_max': 60
        }
        self.graph_fac.make_histogram(x_bar_centres, [heights], [wd_name], wd_name + '_' + model.get_prefix(observation_number), half_bin_size*2, 1.1, 'Pressure /GPa', 'pressure_dist', None, None, additional_x_axis_dict)
        return p_stats

    def build_model(self, hierarchy_name, hierarchy_level_list):
        # Based on the equivalent function in manager.py
        # Make a pollution model consisting of the parameters in the relevent hierarchy levels of the parameter hierarchy
        model_name = hierarchy_name + '_levels_'
        parameters_to_use = list()
        for hl in hierarchy_level_list:
            parameters_to_use += mp.hierarchy_definitions_dict[hierarchy_name][hl]
            model_name += str(hl)
        mp.model_definitions_dict[model_name] = dict()
        for potential_param in mp.ModelParameter:
            mp.model_definitions_dict[model_name][potential_param] = potential_param in parameters_to_use
        #TODO: Add a warning when some incompatible parameter combination is entered. e.g. fragment_crust without fragment_core
        ld._live_model = model_name
        return model_name

    def recreate_model(self, observation_number):
        # This is hacky! Sorry
        candidate_files = list()
        for filename in os.listdir(self.graph_dir):
            if filename.endswith('pressure_dist.pdf'):
                try:
                    test_obs_number = int(filename.split('obs')[1].split('_')[0])
                    if test_obs_number == observation_number:
                        candidate_files.append(filename)
                except:
                    pass
        if len(candidate_files) > 1:
            raise # Now we're stuck
        elif len(candidate_files) < 1:
            return None
        else:
            # Parse file name
            filename = candidate_files[0]
            system_name = filename.split('_')[0]
            hierarchy_name = filename.split(system_name + '_')[1].split('_levels')[0]
            hierarchy_level_list_str = filename.split('levels_')[1].split('_')[0]
            hierarchy_level_list = [int(c) for c in hierarchy_level_list_str]
            model_name = hierarchy_name + '_levels_' + hierarchy_level_list_str
            self.build_model(hierarchy_name, hierarchy_level_list)
            enhancement_model_abbreviation = filename.split('obs')[1].split('_')[1]
            prior_abbreviation = filename.split('obs')[1].split('_')[2]
            live_points = int(filename.split('lp')[1].split('_')[0])
            enhancement_model_name = None
            prior_name = None
            for name, abbrev in pu.abbreviations.items():
                if abbrev == prior_abbreviation:
                    prior_name = name
                if abbrev == enhancement_model_abbreviation:
                    enhancement_model_name = name
            path_to_output = self.graph_dir.split('graphs')[0] + 'output'
        return pm.PollutionModel(model_name, path_to_output, enhancement_model_name, prior_name, live_points)
    
    def make_multisystem_pressure_plot(self, observation_numbers, labels=None, additional_text_dict=None):
        p_stats_list = list()
        used_observation_numbers_list = list()
        for observation_number in observation_numbers:
            recreated_model = self.recreate_model(observation_number)
            p_stats = self.get_untransformed_parameter_stats(recreated_model, observation_number, mp.ModelParameter.pressure)
            if p_stats is None:
                print('Abandoning pressure plot (no pressure stats available)')
                continue
            p_stats_list.append(p_stats)
            used_observation_numbers_list.append(str(observation_number))
        
        half_bin_size = 0.5
        bins, x_bar_centres = self.generate_bins_and_bar_centres(0, 60, half_bin_size)
        
        all_heights = list()
        for p_stats in p_stats_list:
            heights, bins2 = np.histogram(
                p_stats,
                bins,
                density=True
            )
            all_heights.append(heights)
        
        prefix = 'multisys'
        for uon in used_observation_numbers_list:
            prefix += '_' + uon

        hist_plot = self.graph_fac.make_histogram(x_bar_centres, all_heights, used_observation_numbers_list if labels is None else labels, prefix, half_bin_size*2, 0.8, 'Pressure /GPa', 'pressure_dist', additional_text_dict, None, None, 1.1)
        return hist_plot

    def get_temperature_stats(self, model, observation_number):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number)
        parameter_indices = mp.parameter_indices(ld._live_model)
        try:
            distance_index = parameter_indices[mp.ModelParameter.formation_distance]
        except KeyError:
            # If there was no formation_distance parameter, this calculation is irrelevant
            return None
        T = np.zeros(len(weightpost))
        for i in range(0,len(weightpost)):
            T[i] = dm.T_disc(10**(weightpost[i, distance_index]), 1.5)  # TODO: Consider having a variable like default_t_formation = 1.5 in some high level module
        return T

    def generate_bins_and_bar_centres(self, min_bin_edge, max_bin_edge, half_bin_size):
        bin_size = 2 * half_bin_size
        bins = np.arange(min_bin_edge, max_bin_edge + bin_size, bin_size)
        bar_centres = [x + half_bin_size for x in bins]
        bar_centres.pop()
        return bins, bar_centres

    def make_temperature_plot(self, model, model_name, wd_name, observation_number):
        temp_stats = self.get_temperature_stats(model, observation_number)
        if temp_stats is None:
            print('Abandoning temperature plot (no temperature stats available)')
            return
        half_bin_size = 12.5
        bins, x_bar_centres = self.generate_bins_and_bar_centres(0, 3000, half_bin_size)
        heights, bins2 = np.histogram(
            temp_stats,
            bins,
            density=True
        )
        text_dict = {
            'icy_text': {
                'x_pos': 130, # Change to 0 for plotting separate dists...,
                'text_string': 'Icy',
                'horizontalalignment': 'center'
            },
            'dry_text': {
                'x_pos': 700,
                'text_string': 'Dry',
                'horizontalalignment': 'center'
            },
            'eh_text': {
                'x_pos': 2175,
                'text_string': 'Extreme Heating',
                'horizontalalignment': 'center'
            }
        }
        line_dict = {
            'vline1': {
                'x_start': 220
            },
            'vline2': {
                'x_start': 1250
            }
        }
        self.graph_fac.make_histogram(x_bar_centres, [heights], [wd_name], wd_name + '_' + model.get_prefix(observation_number), half_bin_size*2, 1.1, 'Formation Temperature /K', 'temp_dist', text_dict, line_dict)
        return temp_stats
        
    def get_mass_stats(self, model, observation_number, enhancement_model):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number)
        number_of_possible_samples = len(weightpost[:,0])
        max_number_of_samples = 10000
        number_of_samples = min(max_number_of_samples, number_of_possible_samples)
        samples = np.random.choice(number_of_possible_samples, number_of_samples, replace=False)
        massofpollutant = np.zeros(number_of_samples)

        geo_model_for_solar_mass = gi.GeologyModel()
        solar_mass_in_kg = geo_model_for_solar_mass.M_Sun
        HorHe = ld._live_type
        q = ld._live_q
        M_wd = ld._live_mass
        parameter_indices = mp.parameter_indices(ld._live_model)
        print('Sampling for mass plot')
        for j in range(number_of_samples):
            if (j % 250) == 0:
                print('Mass plot sample ' + str(j))
            fe_star = weightpost[samples[j], parameter_indices[mp.ModelParameter.metallicity]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.metallicity) else mp.default_values[enhancement_model][mp.ModelParameter.metallicity]
            t_sinceaccretion = weightpost[samples[j], parameter_indices[mp.ModelParameter.t_sinceaccretion]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.t_sinceaccretion) else mp.default_values[enhancement_model][mp.ModelParameter.t_sinceaccretion]
            d_formation = weightpost[samples[j], parameter_indices[mp.ModelParameter.formation_distance]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.formation_distance) else mp.default_values[enhancement_model][mp.ModelParameter.formation_distance]
            z_formation = weightpost[samples[j], parameter_indices[mp.ModelParameter.feeding_zone_size]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.feeding_zone_size) else mp.default_values[enhancement_model][mp.ModelParameter.feeding_zone_size]
            N_c = weightpost[samples[j], parameter_indices[mp.ModelParameter.parent_core_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.parent_core_frac) else mp.default_values[enhancement_model][mp.ModelParameter.parent_core_frac]
            N_o = weightpost[samples[j], parameter_indices[mp.ModelParameter.parent_crust_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.parent_crust_frac) else mp.default_values[enhancement_model][mp.ModelParameter.parent_crust_frac]
            f_c = weightpost[samples[j], parameter_indices[mp.ModelParameter.fragment_core_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.fragment_core_frac) else mp.default_values[enhancement_model][mp.ModelParameter.fragment_core_frac]
            f_o = weightpost[samples[j], parameter_indices[mp.ModelParameter.fragment_crust_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.fragment_crust_frac) else mp.default_values[enhancement_model][mp.ModelParameter.fragment_crust_frac]
            pollutionfraction = weightpost[samples[j], parameter_indices[mp.ModelParameter.pollution_frac]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.pollution_frac) else mp.default_values[enhancement_model][mp.ModelParameter.pollution_frac]
            log_t_disc = weightpost[samples[j], parameter_indices[mp.ModelParameter.accretion_timescale]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.accretion_timescale) else mp.default_values[enhancement_model][mp.ModelParameter.accretion_timescale]
            pressure = weightpost[samples[j], parameter_indices[mp.ModelParameter.pressure]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.pressure) else mp.default_values[enhancement_model][mp.ModelParameter.pressure]
            fO2 = weightpost[samples[j], parameter_indices[mp.ModelParameter.oxygen_fugacity]] if mp.model_uses_parameter(ld._live_model, mp.ModelParameter.oxygen_fugacity) else mp.default_values[enhancement_model][mp.ModelParameter.oxygen_fugacity]
            t_disc = 10**log_t_disc
            t_formation = 1.5
            
            
            expressions, ignore = cm.complete_model_calculation(
                fe_star,
                t_sinceaccretion,
                d_formation,
                z_formation,
                N_c,
                N_o,
                f_c,
                f_o,
                pollutionfraction,
                t_disc,
                pressure,
                fO2,
                enhancement_model,
                t_formation,
                True,
                False
            )
            
            if expressions is None:
                print('Warning! The sampled model does not converge! (Or otherwise returns None). Skipping this mass plot sample. Printing attempted params:')
                print('fe_star: ' + str(fe_star))
                print('t_sinceaccretion: ' + str(t_sinceaccretion))
                print('d_formation: ' + str(d_formation))
                print('z_formation: ' + str(z_formation))
                print('N_c: ' + str(N_c))
                print('N_o: ' + str(N_o))
                print('f_c: ' + str(f_c))
                print('f_o: ' + str(f_o))
                print('pollutionfraction: ' + str(pollutionfraction))
                print('log_t_disc: ' + str(log_t_disc))
                print('pressure: ' + str(pressure))
                print('fO2: ' + str(fO2))
                print('t_formation: ' + str(t_formation))
                continue
            
            element_masses = dict() # in kg. At the moment, storing this in a dict is overly complicated but we may want to break it down by element at some point
            for element in ci.usual_elements:
                # This calculation assumes that M_cvz = M_HorHe, which should be a pretty safe assumption
                log_el_HorHe = expressions[element] # log of number ratio of lifetime element to present H (or He)
                el_to_HorHe_mass_ratio = ci.get_element_mass(element)/ci.get_element_mass(HorHe)
                element_mass = (10**(q + log_el_HorHe)) * el_to_HorHe_mass_ratio * M_wd * solar_mass_in_kg
                element_masses[element] = element_mass
            total_mass = 0
            for el, mass in element_masses.items():
                total_mass += mass
            massofpollutant[j] = total_mass
        print(massofpollutant)
        logmassofpollutant = np.log10(massofpollutant)
        return logmassofpollutant
                
    def make_mass_plot(self, model, model_name, wd_name, observation_number, enhancement_model):
        logmassofpollutant = self.get_mass_stats(model, observation_number, enhancement_model)
        
        half_bin_size = 0.05
        bins, x_bar_centres = self.generate_bins_and_bar_centres(8, 25, half_bin_size)
        
        heights, bins2 = np.histogram(
            logmassofpollutant,
            bins,
            density=True
        )
        
        text_dict = {
            'Earth': {
                'x_pos': 24.78,
                'text_string': 'Earth',
                'horizontalalignment': 'right',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Mars': {
                'x_pos': 23.81,
                'text_string': 'Mars',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Moon': {
                'x_pos': 22.87,
                'text_string': 'Moon',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Pluto': {
                'x_pos': 22.11,
                'text_string': 'Pluto',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Ceres': {
                'x_pos': 20.97,
                'text_string': 'Ceres',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Vesta': {
                'x_pos': 20.41,
                'text_string': 'Vesta',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Hygiea': {
                'x_pos': 19.94,
                'text_string': 'Hygiea',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Psyche': {
                'x_pos': 19.38,
                'text_string': 'Psyche',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Flora': {
                'x_pos': 18.93,
                'text_string': 'Flora',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Epimetheus': {
                'x_pos': 17.72,
                'text_string': 'Epimetheus',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Ophelia': {
                'x_pos': 16.72,
                'text_string': 'Ophelia',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Phobos': {
                'x_pos': 16.02,
                'text_string': 'Phobos',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Deimos': {
                'x_pos': 15.17,
                'text_string': 'Deimos',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Comet Halley': {
                'x_pos': 14.34,
                'text_string': 'Comet Halley',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Toutatis': {
                'x_pos': 13.70,
                'text_string': 'Toutatis',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Comet 67P': {
                'x_pos': 13.00,
                'text_string': 'Comet 67P',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Comet SL9': {
                'x_pos': 12.18,
                'text_string': 'Comet SL9',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Ryugu': {
                'x_pos': 11.65,
                'text_string': 'Ryugu',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Bennu': {
                'x_pos': 11.15,
                'text_string': 'Bennu',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            'Itokawa': {
                'x_pos': 10.55,
                'text_string': 'Itokawa',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            },
            '1994 WR12': {
                'x_pos': 9.46,
                'text_string': '1994 WR12',
                'verticalalignment': 'top',
                'rotation': 'vertical'
            }
        }
        
        self.graph_fac.make_histogram(x_bar_centres, [heights], [wd_name], wd_name + '_' + model.get_prefix(observation_number), half_bin_size*2, 1.2, 'Log(Mass of Pollutant/kg)', 'mass_dist', text_dict, None)
        return logmassofpollutant # This is just so that the later functions don't need to recalculate it!
           
    def make_time_plot(self, model, model_name, wd_name, observation_number):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number)
        parameter_indices = mp.parameter_indices(ld._live_model)
        t_sinceaccretion_index = parameter_indices[mp.ModelParameter.t_sinceaccretion]
        accretion_timescale_index = parameter_indices[mp.ModelParameter.accretion_timescale]
        x_data = np.log10(weightpost[:,t_sinceaccretion_index]*1000000)
        y_data = weightpost[:,accretion_timescale_index]
        self.graph_fac.plot_timesince_v_accretiontime(x_data, y_data, wd_name, model.get_prefix(observation_number), ld._live_t_mg)
        
    def make_p_v_fO2_plot(self, model, model_name, wd_name, observation_number):
        x_data = self.get_untransformed_parameter_stats(model, observation_number, mp.ModelParameter.pressure)
        y_data = self.get_untransformed_parameter_stats(model, observation_number, mp.ModelParameter.oxygen_fugacity)
        if x_data is None or y_data is None:
            print('Abandoning pressure v fO2 plot (no pressure and/or fO2 stats available)')
            return
        self.graph_fac.plot_pressure_v_oxygen_fugacity(x_data, y_data, wd_name, model.get_prefix(observation_number))
        
