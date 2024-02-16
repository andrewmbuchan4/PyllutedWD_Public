#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections as cn
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
import excess_oxygen_calculator as eoc
import geology_info as gi
import graph_factory as gf
import live_data as ld
import model_parameters as mp
import pollution_model as pm
import pwd_utils as pu
import solar_abundances as sa
import white_dwarf as wd

class ModelAnalyser:

    def __init__(self, graph_dir):
        self.graph_dir = graph_dir
        self.graph_fac = gf.GraphFactory(self.graph_dir)
        self.excess_oxygen_calculator = eoc.ExcessOxygenCalculator()
        self.absolute_string = 'absolute'
        self.fractional_string = 'fractional'

    def update_graph_dir(self, graph_dir):
        self.graph_dir = graph_dir
        self.graph_fac.output_dir = graph_dir

    def confidence_intervals(self, input_array):
        sig_1 = 50 + 68.26/2  # Cumulative probability of +1 sigma level (0.8413) (0.1587)
        sig_2 = 50 + 95.4/2   # Cumulative probability of +2 sigma level (0.977) (0.023)
        sig_3 = 50 + 99.7/2   # Cumulative probability of +3 sigma level (0.9985) (0.0015)

        percentile_vals = [100-sig_3, 100-sig_2, 100-sig_1, 50, sig_1, sig_2, sig_3]

        if not isinstance(input_array, np.ndarray):
            input_array = np.array(input_array)

        length = input_array.ndim
        if length > 1:
            toret = np.zeros((len(input_array[0,:]), len(percentile_vals)))
            for i, sample in enumerate(input_array.T):
                sample_filter = [a for a in sample if a is not None and not np.isnan(a)]
                toret[i,:] = [np.percentile(sample_filter, pv) for pv in percentile_vals]
            return toret[:,0], toret[:,1], toret[:,2], toret[:,3], toret[:,4], toret[:,5], toret[:,6]
        else:
            sample_filter = [a for a in input_array if a is not None and not np.isnan(a)]
            toret = [np.percentile(sample_filter, pv) for pv in percentile_vals]
            return toret[0], toret[1], toret[2], toret[3], toret[4], toret[5], toret[6]

    def z_to_sigma(self, ln_Z1, ln_Z2):
        np.set_printoptions(precision=50)
        B = np.exp(ln_Z1 - ln_Z2)
        p = np.real(np.exp(W((-1.0/(B*np.exp(1))),-1)))
        sigma = np.sqrt(2)*erfcinv(p)
        return B, sigma

    def detection_significance(self, err_data, n_params, n_removed, basename_M1, basename_M2):

        norm_log = (-0.5*np.log(2.0*np.pi*err_data*err_data)).sum()

        retrieved_M1 = pn.Analyzer(n_params = n_params, outputfiles_basename = basename_M1)
        retrieved_M2 = pn.Analyzer(n_params = (n_params-n_removed), outputfiles_basename = basename_M2)
        s_M1 = retrieved_M1.get_stats()
        s_M2 = retrieved_M2.get_stats()

        best_fit_M1 = retrieved_M1.get_best_fit()                       # Best-fitting parameter combination from model 1
        max_likelihood_M1 = best_fit_M1['log_likelihood']               # Maximum log-likelihood for model 1
        best_chi_square_M1 = -2.0 * (max_likelihood_M1 - norm_log)      # Best chi-square for model 1

        ln_Z_M1 = s_M1['global evidence']     # Natural log of Bayesian evidence from model 1


        best_fit_M2 = retrieved_M2.get_best_fit()                                      # Best-fitting parameter combination from model 2
        max_likelihood_M2 = best_fit_M2['log_likelihood']                              # Maximum log-likelihood for model 2
        best_chi_square_M2 = -2.0 * (max_likelihood_M2 - norm_log)                     # Best chi-square for model 2

        ln_Z_M2 = s_M2['global evidence']     # Natural log of Bayesian evidence from model 1


        Bayes_factor, n_sigma = self.z_to_sigma(ln_Z_M1, ln_Z_M2)
        print('Comparing ' + basename_M1 + ' with ' + basename_M2)
        print('Bayes Factor = ' + str(Bayes_factor))

        return ln_Z_M1, ln_Z_M2, Bayes_factor, n_sigma, best_chi_square_M1, best_chi_square_M2

    def compare(self, model_dict, output_dir, observation_number, number_of_data_points):
        # model_dict should be an OrderedDict with the base model in first position. Comparison of base to itself is to get ln_Z
        base_model_name = None
        best_model_so_far_name = None
        for model_name, model in model_dict.items():
            if base_model_name is None:
                base_model_name = model_name
            else:
                base_model_name = best_model_so_far_name
            base_model = model_dict[base_model_name]

            self.compare_two_models(model, base_model, output_dir, observation_number, number_of_data_points)

            if best_model_so_far_name is None:
                best_model_so_far_name = base_model_name
            else:
                ln_Z_model = model.comparison[base_model_name]['ln_Z_model']
                ln_Z_base = model.comparison[base_model_name]['ln_Z_base']
                if ln_Z_model > ln_Z_base:
                    best_model_so_far_name = model_name

    def compare_two_models(self, model, base_model, output_dir, observation_number, number_of_data_points):
        print('Comparing result for model ' + str(model.basename) + ' with base_model ' + str(base_model.basename))
        good_fit_threshold = 2
        ln_Z_model, ln_Z_base, Bayes_factor_model_base, n_sigma_model_base, chi_model, chi_base = self.detection_significance(
            ld._live_white_dwarf.get_errors_for_present_elements_as_array(),
            model.get_n_dims(),
            (model.get_n_dims() - base_model.get_n_dims()),
            model.get_full_prefix(output_dir, observation_number),
            base_model.get_full_prefix(output_dir, observation_number)
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

    def find_best_heated_model(self, models, system, number_of_data_points, stats_file, model_output_dir):
        heating_string = mp.model_parameter_strings[mp.ModelParameter.formation_distance]
        best_heated_model_name = None
        best_heated_model_lnZ = None
        for model_name, model in models.items():
            if model.prior_name == 'Default' and heating_string in model.get_model_params():
                if best_heated_model_name is None:
                    best_heated_model_name = model_name
                    if len(model.comparison.items()) == 0:
                        self.compare_two_models(model, model, model_output_dir, system, number_of_data_points)  # Compare model to itself - we just need the ln Z so actual comparison is arbitrary
                    assert len(model.comparison.items()) > 0
                    for dummy_name, arbitrary_comparison in model.comparison.items():
                        best_heated_model_lnZ = arbitrary_comparison['ln_Z_model']
                else:
                    if len(model.comparison.items()) == 0:
                        self.compare_two_models(model, model, model_output_dir, system, number_of_data_points)  # Compare model to itself - we just need the ln Z so actual comparison is arbitrary
                    assert len(model.comparison.items()) > 0
                    test_model_lnZ = None
                    for dummy_name, arbitrary_comparison in model.comparison.items():
                        test_model_lnZ = arbitrary_comparison['ln_Z_model']
                    if test_model_lnZ > best_heated_model_lnZ:
                        best_heated_model_lnZ = test_model_lnZ
                        best_heated_model_name = model_name
        if best_heated_model_name is None:
            # Then there was no model which actually included heating, so no point continuing
            print('Could not find a heated model')
            return
        else:
            for model_name, model in models.items():
                if model_name == best_heated_model_name:
                    model.best_heated_model = True
            with open(stats_file, 'a', newline='', encoding='utf-8') as f:
                to_write = csv.writer(f)
                to_write.writerow([])
                to_write.writerow([
                    'Best heated model name:',
                    best_heated_model_name
                ])

    def find_parameter_sigma(self, process_name, parameter, models, system, number_of_data_points, stats_file, model_output_dir):
        # Find sigma significance of a certain parameter that describes a certain process
        # Do this by identifying the best model without that parameter included,
        # and seeing how much better the best overall model is
        parameter_string = mp.model_parameter_strings[parameter]
        best_non_param_model_name = None
        best_non_param_model_lnZ = None
        for model_name, model in models.items():
            if model.prior_name == 'Default' and parameter_string not in model.get_model_params():
                if best_non_param_model_name is None:
                    best_non_param_model_name = model_name
                    if len(model.comparison.items()) == 0:
                        self.compare_two_models(model, model, model_output_dir, system, number_of_data_points)  # Compare model to itself - we just need the ln Z so actual comparison is arbitrary
                    assert len(model.comparison.items()) > 0
                    for dummy_name, arbitrary_comparison in model.comparison.items():
                        best_non_param_model_lnZ = arbitrary_comparison['ln_Z_model']
                else:
                    if len(model.comparison.items()) == 0:
                        self.compare_two_models(model, model, model_output_dir, system, number_of_data_points)  # Compare model to itself - we just need the ln Z so actual comparison is arbitrary
                    assert len(model.comparison.items()) > 0
                    test_model_lnZ = None
                    for dummy_name, arbitrary_comparison in model.comparison.items():
                        test_model_lnZ = arbitrary_comparison['ln_Z_model']
                    if test_model_lnZ > best_non_param_model_lnZ:
                        best_non_param_model_lnZ = test_model_lnZ
                        best_non_param_model_name = model_name
        if best_non_param_model_name is None:
            # Then there was no model which actually omitted the parameter, so no point continuing
            print('Could not find a model without parameter ' + str(parameter) + ', no sigma significance could be calculated')
            return
        else:
            for model_name, model in models.items():
                if model_name == best_non_param_model_name:
                    model.best_model_without_parameter[parameter] = True
        # Now find best model:
        best_model_name = None
        for model_name, model in models.items():
            if model.prior_name == 'Default' and model.best_model:
                best_model_name = model_name
        # Now compare them if necessary
        if best_model_name != best_non_param_model_name:
            if best_non_param_model_name not in models[best_model_name].comparison:
                self.compare_two_models(models[best_model_name], models[best_non_param_model_name], model_output_dir, system, number_of_data_points)
                with open(stats_file, 'a', newline='', encoding='utf-8') as f:
                    to_write = csv.writer(f)
                    to_write.writerow([])
                    to_write.writerow([
                        best_model_name,
                        best_non_param_model_name,
                        models[best_model_name].comparison[best_non_param_model_name]['ln_Z_model'],
                        models[best_model_name].comparison[best_non_param_model_name]['ln_Z_base'],
                        models[best_model_name].comparison[best_non_param_model_name]['Bayes_factor_model_base'],
                        models[best_model_name].comparison[best_non_param_model_name]['n_sigma_model_base'],
                        models[best_model_name].comparison[best_non_param_model_name]['chi_model'],
                        models[best_model_name].comparison[best_non_param_model_name]['chi_base'],
                        models[best_model_name].comparison[best_non_param_model_name]['chi_model_per_data_point'],
                        models[best_model_name].comparison[best_non_param_model_name]['chi_base_per_data_point'],
                        models[best_model_name].get_model_params(),
                        'N/A' if models[best_model_name].best_model is None else models[best_model_name].best_model,
                        models[best_model_name].comparison[best_non_param_model_name]['good_fit']
                    ])
            param_sigma_str = str(models[best_model_name].comparison[best_non_param_model_name]['n_sigma_model_base'])
            bayes_factor_str = str(models[best_model_name].comparison[best_non_param_model_name]['Bayes_factor_model_base'])
            chi_sq_nondiff_str = str(models[best_model_name].comparison[best_non_param_model_name]['chi_base_per_data_point'])
        else:
            param_sigma_str = 'N/A'
            bayes_factor_str = 'N/A'
            chi_sq_nondiff_str = 'N/A'
        with open(stats_file, 'a', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow([
                'Best model name:',
                best_model_name
            ])
            to_write.writerow([
                'Best model name without ' + process_name + ':',
                best_non_param_model_name
            ])
            to_write.writerow([
                process_name + ' Sigma:',
                param_sigma_str
            ])
            to_write.writerow([
                process_name + ' Bayes Factor:',
                bayes_factor_str
            ])
            to_write.writerow([
                'Chi squared per data point (best model without ' + process_name + '):',
                chi_sq_nondiff_str
            ])

    def dump_model_stats(self, model_name, model, stats_file):
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

    def make_plots_and_dump_fit(self, white_dwarf, N_wd, chains_dir, stats_file, model_name, model, enhancement_model, bonus_fits, bonus_error_lows, bonus_error_highs, suppress_graphical_output=False):
        temp_stats, mass_stats, comp_fits_and_errors, eo_samples, semisampled_eo_dict, disc_abundance_medians, bulk_composition_medians, core_composition_medians, mantle_composition_medians, partition_coefficient_medians, pcnf_median = self.make_all_plots(
            white_dwarf,
            chains_dir,
            model_name,
            model,
            N_wd,
            enhancement_model,
            bonus_fits,
            bonus_error_lows,
            bonus_error_highs,
            suppress_graphical_output
        )
        #if bonus_plots != []:
        #    for key, fit in mk2_fits_and_errors[0].items():
        #        bonus_fits[bonus_plot_names[model_name] + ' ' + key] = fit
        #    for key, error_lows in mk2_fits_and_errors[1].items():
        #        bonus_error_lows[bonus_plot_names[model_name] + ' ' + key] = error_lows
        #    for key, error_highs in mk2_fits_and_errors[2].items():
        #        bonus_error_highs[bonus_plot_names[model_name] + ' ' + key] = error_highs
        self.dump_model_fit(
            white_dwarf,
            chains_dir,
            model,
            N_wd,
            stats_file,
            comp_fits_and_errors,
            temp_stats,
            mass_stats,
            model_name,
            eo_samples,
            semisampled_eo_dict,
            disc_abundance_medians,
            bulk_composition_medians,
            core_composition_medians,
            mantle_composition_medians,
            partition_coefficient_medians,
            pcnf_median,
            suppress_graphical_output
        )

    #def make_all_plots(self, chains_dir, model_name, model, observation_number, all_wd_timescales, wd_name, all_wd_abundances, all_wd_upper_bounds, all_wd_lower_bounds, all_wd_abundance_errors, excluded_wd_abundances, excluded_wd_upper_bounds, excluded_wd_lower_bounds, excluded_wd_abundance_errors, enhancement_model, wd_type, bonus_fits=None, bonus_error_lows=None, bonus_error_highs=None, suppress_graphical_output=False):
    def make_all_plots(self, white_dwarf, chains_dir, model_name, model, observation_number, enhancement_model, bonus_fits=None, bonus_error_lows=None, bonus_error_highs=None, suppress_graphical_output=False):
        #TODO put results from each model in a separate subdirectory
        wd_name = white_dwarf.name
        #mk2_fits_and_errors = self.make_composition_plot_mk2(chains_dir, model_name, model, observation_number, all_wd_timescales, wd_name, all_wd_abundances, all_wd_upper_bounds, all_wd_lower_bounds, all_wd_abundance_errors, excluded_wd_abundances, excluded_wd_upper_bounds, excluded_wd_lower_bounds, excluded_wd_abundance_errors, enhancement_model, wd_type, bonus_fits, bonus_error_lows, bonus_error_highs, suppress_graphical_output)
        comp_fits_and_errors, eo_samples, disc_abundance_medians, bulk_composition_medians, core_composition_medians, mantle_composition_medians, partition_coefficient_medians, pcnf_median = self.make_composition_plot(chains_dir, model_name, model, observation_number, white_dwarf, enhancement_model, bonus_fits, bonus_error_lows, bonus_error_highs, suppress_graphical_output)
        temp_stats = self.make_temperature_plot(chains_dir, model, model_name, wd_name, observation_number, suppress_graphical_output)
        mass_stats = self.make_mass_plot(chains_dir, model, model_name, wd_name, observation_number, enhancement_model, suppress_graphical_output)
        #semisampled_eo_dict = self.make_semisampled_eo_distibution_plot(chains_dir, model_name, model, observation_number, all_wd_timescales, wd_name, all_wd_abundances, all_wd_upper_bounds, all_wd_lower_bounds, all_wd_abundance_errors, suppress_graphical_output)
        semisampled_eo_dict = self.make_semisampled_eo_distibution_plot(chains_dir, model_name, model, observation_number, white_dwarf, comp_fits_and_errors, suppress_graphical_output)
        if not suppress_graphical_output:
            self.make_corner_plot(model_name, model, observation_number, wd_name, chains_dir)
            self.make_time_plot(chains_dir, model, model_name, wd_name, observation_number)
            self.make_scaled_tevent_plot(chains_dir, model, model_name, white_dwarf, observation_number)
            self.make_p_v_fO2_plot(chains_dir, model, model_name, wd_name, observation_number)
            self.make_pressure_plot(chains_dir, model, model_name, wd_name, observation_number)
            self.make_eo_distribution_plot(model, model_name, wd_name, observation_number, eo_samples)
        return temp_stats, mass_stats, comp_fits_and_errors, eo_samples, semisampled_eo_dict, disc_abundance_medians, bulk_composition_medians, core_composition_medians, mantle_composition_medians, partition_coefficient_medians, pcnf_median

    def get_eo_sigma(self, eo_samples):
        if eo_samples is None:
            return None, None, None, None, None, None
        eo_sample_stats_no_nans = [eo for eo in eo_samples if eo is not None and not np.isnan(eo)]
        median_excess_oxygen = np.percentile(eo_sample_stats_no_nans, 50)
        lower_error_excess_oxygen = np.percentile(eo_sample_stats_no_nans, 16)
        upper_error_excess_oxygen = np.percentile(eo_sample_stats_no_nans, 84)
        negative_eos = sum(1 for i in eo_sample_stats_no_nans if i < 0)
        positive_eos = sum(1 for i in eo_sample_stats_no_nans if i > 0)
        if median_excess_oxygen >= 0:
            probability_of_result = positive_eos/len(eo_sample_stats_no_nans)
        else:
            probability_of_result = negative_eos/len(eo_sample_stats_no_nans)
        eo_sigma_2 = stat.norm.ppf(probability_of_result)
        return median_excess_oxygen, upper_error_excess_oxygen - median_excess_oxygen, median_excess_oxygen - lower_error_excess_oxygen, eo_sigma_2, len(eo_sample_stats_no_nans), len(eo_samples), probability_of_result

    def get_predicted_abundances(self, model, parameter_values):
        raise # I don't think we should do this any more, we should use the abundance fits found in make_composition_plot
        input_values = list()
        for param in mp.get_model_params_in_order():
            if mp.model_parameter_strings[param] in model.get_model_params():
                index = model.get_model_params().index(mp.model_parameter_strings[param])
                to_append = parameter_values[index]
            else:
                to_append = mp.default_values[ld._enhancement_model][param]
            input_values.append(to_append)
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
        return model_result, diagnostics, input_values

    def get_param_vals_by_percentile(self, model, weightpost, percentile):
        toret = list()
        for index, param in enumerate(model.get_model_params()):
            toret.append(np.percentile(weightpost[:,index], percentile))
        return toret

    def reconstruct_composition_dict(self, bulk_composition, core_composition, mantle_composition):
        toret = dict()
        for element in ci.usual_elements:
            toret[element] = {
                gi.Layer.bulk: bulk_composition[element],
                gi.Layer.core: core_composition[element],
                gi.Layer.mantle: mantle_composition[element]
            }
        return toret

    def dump_model_fit(self, white_dwarf, chains_dir, model, N_wd, stats_file, fit_data, temp_stats, mass_stats, model_name, eo_sample_stats, semisampled_eo_dict, disc_abundance_medians, bulk_composition_medians, core_composition_medians, mantle_composition_medians, partition_coefficient_medians, pcnf_median, suppress_graphical_output=False):
        print_all_eo_sample_vals = False # Switch these to False by default!
        print_all_time_sample_vals = False
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, N_wd, chains_dir)
        parameter_indices = mp.parameter_indices(ld._live_model)
        t_sinceaccretion_index = parameter_indices[mp.ModelParameter.t_sinceaccretion]
        accretion_timescale_index = parameter_indices[mp.ModelParameter.accretion_timescale]
        modes = list()
        for index, param in enumerate(model.get_model_params()):
            modes.append(np.asscalar(stat.mode(np.around(weightpost[:,index], decimals=2))[0]))
        medians = self.get_param_vals_by_percentile(model, weightpost, 50)
        percentile_5 = self.get_param_vals_by_percentile(model, weightpost, 5)
        percentile_16 = self.get_param_vals_by_percentile(model, weightpost, 16)
        percentile_84 = self.get_param_vals_by_percentile(model, weightpost, 84)
        percentile_95 = self.get_param_vals_by_percentile(model, weightpost, 95)
        #model_result, diagnostics, input_values = self.get_predicted_abundances(model, medians)
        try:
            elements_list = list()
            for el in fit_data[white_dwarf.get_atmospheric_type().value]['Model median'].keys():
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
        pollutant_comp_list = list()
        ds = list()
        for el in elements_list:
            try:
                elements_results.append(fit_data[white_dwarf.get_atmospheric_type().value]['Model median'][el])
            except TypeError:
                pass
            try:
                disc_comp.append(disc_abundance_medians[el])
            except (TypeError, KeyError):
                pass
            try:
                bulk_comp.append(bulk_composition_medians[el])
            except (TypeError, KeyError):
                pass
            try:
                core_comp.append(core_composition_medians[el])
            except (TypeError, KeyError):
                pass
            try:
                mantle_comp.append(mantle_composition_medians[el])
            except (TypeError, KeyError):
                pass
            try:
                ds.append(partition_coefficient_medians[el])
            except (TypeError, KeyError):
                pass
            #try:
            #    wd_abundances_list.append(wd_abundances[el])
            #except (TypeError, KeyError):
            #    pass
            #try:
            #    wd_abundance_upper_bounds_list.append(wd_abundance_upper_bounds[el])
            #except (TypeError, KeyError):
            #    pass
            #try:
            #    wd_abundance_lower_bounds_list.append(wd_abundance_lower_bounds[el])
            #except (TypeError, KeyError):
            #    pass
            #try:
            #    wd_errors_list.append(wd_errors[el])
            #except (TypeError, KeyError):
            #    pass
            #wd_excluded_abundances_list.append(excluded_wd_abundances.get(el, ''))
            #wd_excluded_abundance_upper_bounds_list.append(excluded_wd_upper_bounds.get(el, ''))
            #wd_excluded_abundance_lower_bounds_list.append(excluded_wd_lower_bounds.get(el, ''))
            #wd_excluded_errors_list.append(excluded_wd_abundance_errors.get(el, ''))

            wd_abundance_els_list, wd_abundances_list, wd_errors_list, wd_lower_errors_list = white_dwarf.get_measurements_as_lists(
                True,
                False,
                ci.usual_elements,
                True
            )
            wd_excluded_abundance_els_list, wd_excluded_abundances_list, wd_excluded_errors_list, wd_excluded_lower_errors_list = white_dwarf.get_measurements_as_lists(
                False,
                True,
                ci.usual_elements,
                True
            )
            wd_upper_bound_els_list, wd_abundance_upper_bounds_list = white_dwarf.get_upper_bounds_as_lists(
                True,
                False,
                ci.usual_elements,
                True
            )
            wd_excluded_upper_bound_els_list, wd_excluded_abundance_upper_bounds_list = white_dwarf.get_upper_bounds_as_lists(
                False,
                True,
                ci.usual_elements,
                True
            )
            wd_lower_bound_els_list, wd_abundance_lower_bounds_list = white_dwarf.get_lower_bounds_as_lists(
                True,
                False,
                ci.usual_elements,
                True
            )
            wd_excluded_lower_bound_els_list, wd_excluded_abundance_lower_bounds_list = white_dwarf.get_lower_bounds_as_lists(
                False,
                True,
                ci.usual_elements,
                True
            )
        # Check that the outputs above will be unambiguous:
        assert len(elements_results) in [0, len(elements_list)]
        assert len(disc_comp) in [0, len(elements_list)]
        assert len(bulk_comp) in [0, len(elements_list)]
        assert len(core_comp) in [0, len(elements_list)]
        assert len(mantle_comp) in [0, len(elements_list)]
        assert len(ds) in [0, len(elements_list)]
        assert len(wd_abundances_list) in [0, len(elements_list)]
        assert len(wd_abundance_upper_bounds_list) in [0, len(elements_list)]
        assert len(wd_abundance_lower_bounds_list) in [0, len(elements_list)]
        assert len(wd_errors_list) in [0, len(elements_list)]
        radius = None
        mass = None
        radius_lower_limit = None
        mass_lower_limit = None
        eoc_stats = None
        model_prediction_string = 'Model prediction'
        cmf_median = None
        cmf_percentile_16 = None
        cmf_percentile_84 = None
        median_composition_dict = self.reconstruct_composition_dict(bulk_composition_medians, core_composition_medians, mantle_composition_medians)
        if mp.model_parameter_strings[mp.ModelParameter.pressure] in model.get_model_params():
            fcf_index = parameter_indices[mp.ModelParameter.fragment_core_frac]
            pressure_index = parameter_indices[mp.ModelParameter.pressure]
            median_fcf = medians[fcf_index]
            median_pressure = medians[pressure_index]
            geo_model = gi.GeologyModel()  # Argument here irrelevent, we will supply the needed abundances
            cmf_median = geo_model.convert_core_number_fraction_to_core_mass_fraction(
                median_fcf,
                median_composition_dict # Caveat: using this composition for all fcf values is an approximation - there may be systematic changes in composition with fcf
            )
            cmf_percentile_16 = geo_model.convert_core_number_fraction_to_core_mass_fraction(
                percentile_16[fcf_index],
                median_composition_dict # Caveat: using this composition for all fcf values is an approximation - there may be systematic changes in composition with fcf
            )
            cmf_percentile_84 = geo_model.convert_core_number_fraction_to_core_mass_fraction(
                percentile_84[fcf_index],
                median_composition_dict # Caveat: using this composition for all fcf values is an approximation - there may be systematic changes in composition with fcf
            )
            pollutant_composition = geo_model.get_earth_mix(median_fcf, median_composition_dict)
            for el in elements_list:
                try:
                    pollutant_comp_list.append(pollutant_composition[el])
                except KeyError:
                    pollutant_comp_list.append(None)
            radius_mass = geo_model.calculate_radius_and_mass(median_pressure, median_composition_dict)
            radius = radius_mass[0]
            mass = radius_mass[1]
            radius_mass_lower_limit = geo_model.calculate_radius_and_mass(median_pressure, median_composition_dict, True)
            radius_lower_limit = radius_mass_lower_limit[0]
            mass_lower_limit = radius_mass_lower_limit[1]
            mass_abundance_dict = geo_model.convert_number_abundances_to_mass(median_composition_dict)
            mantle_mass_fraction = 1 - geo_model.get_core_fraction_from_abundance_dict(mass_abundance_dict)
            mantle_mass = mantle_mass_fraction*mass
            # Temporary hack for Marc's thing! But this should be incorporated properly, it's quite useful. And should really take a white dwarf object as an argument rather than all the little bits
            input_dict = {model_prediction_string: pollutant_composition}
            eoc_stats = self.extract_and_plot_excess_oxygen(model, N_wd, input_dict, white_dwarf.name, white_dwarf.get_abundance_values_dict(), white_dwarf.get_error_values_dict(), white_dwarf.get_timescale_values_dict(), 1, suppress_graphical_output)
        else:
            for el in elements_list:
                try:
                    pollutant_comp_list.append(disc_abundance_medians[el])
                except KeyError:
                    pollutant_comp_list.append(None)
            input_dict = {model_prediction_string: disc_abundance_medians}
            eoc_stats = self.extract_and_plot_excess_oxygen(model, N_wd, input_dict, white_dwarf.name, white_dwarf.get_abundance_values_dict(), white_dwarf.get_error_values_dict(), white_dwarf.get_timescale_values_dict(), 1, suppress_graphical_output)

        accretion_timescale_extract = 10**(weightpost[:,accretion_timescale_index] - 6)
        t_sinceaccretion_extract = weightpost[:,t_sinceaccretion_index]
        delta_time = t_sinceaccretion_extract - accretion_timescale_extract

        tmedian, terrorplus, terrorminus = self.extract_median_and_error(delta_time)

        bu_count = 0
        ss_count = 0
        dec_count = 0
        bu_ss_cutoff = (5*white_dwarf.timescale_dict[ci.Element.Mg])/1000000 #Important to convert to Myr! Physically we're saying it takes 5 Mg sinking times to reach steady state
        for t_acc_t_obs_pair in zip(accretion_timescale_extract, t_sinceaccretion_extract):
            accretion_timescale = t_acc_t_obs_pair[0]
            time_since_accretion = t_acc_t_obs_pair[1]
            if time_since_accretion > accretion_timescale:
                dec_count += 1
            elif time_since_accretion < bu_ss_cutoff:
                bu_count += 1
            else:
                ss_count += 1

        bu_percent = (100*bu_count)/(bu_count + ss_count + dec_count)
        ss_percent = (100*ss_count)/(bu_count + ss_count + dec_count)
        dec_percent = (100*dec_count)/(bu_count + ss_count + dec_count)

        Tmedian, Terrorplus, Terrorminus = self.extract_median_and_error(temp_stats)

        Mmedian, Merrorplus, Merrorminus = self.extract_median_and_error(mass_stats)

        pressure_stats = self.get_untransformed_parameter_stats(chains_dir, model, N_wd, mp.ModelParameter.pressure)
        Pmedian, Perrorplus, Perrorminus = self.extract_median_and_error(pressure_stats)
        P_10_percentile = self.extract_xth_percentile(pressure_stats, 10)
        P_90_percentile = self.extract_xth_percentile(pressure_stats, 90)
        P_percentile_of_5 = self.extract_percentile_of_x(pressure_stats, 5)
        P_percentile_of_10 = self.extract_percentile_of_x(pressure_stats, 10)
        P_percentile_of_15 = self.extract_percentile_of_x(pressure_stats, 15)
        P_percentile_of_45 = self.extract_percentile_of_x(pressure_stats, 45)
        P_percentile_of_50 = self.extract_percentile_of_x(pressure_stats, 50)
        P_percentile_of_55 = self.extract_percentile_of_x(pressure_stats, 55)
        oa_elements = dict()
        oa_assignations = dict()
        for oa_strat, oa_stats in eoc_stats['Model prediction'].items():
            oa_elements[oa_strat] = list()
            oa_assignations[oa_strat] = list()
            for element, value in oa_stats[2].items():
                oa_elements[oa_strat].append(element)
                oa_assignations[oa_strat].append(value)

        try:
            cmf_error_plus = cmf_percentile_84 - cmf_median
        except TypeError:
            cmf_error_plus = None
        try:
            cmf_error_minus = cmf_median - cmf_percentile_16
        except TypeError:
            cmf_error_minus = None

        sampled_eo_stats = dict()

        for ox_strat, stats in eo_sample_stats.items():
            sampled_eo_stats[ox_strat] = dict()
            absolute_sample = stats.get(self.absolute_string)
            sampled_eo_stats[ox_strat][self.absolute_string] = self.get_eo_sigma(absolute_sample)
            fractional_sample = stats.get(self.fractional_string)
            sampled_eo_stats[ox_strat][self.fractional_string] = self.get_eo_sigma(fractional_sample)

        with open(stats_file, 'a', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow([])
            to_write.writerow(['Results from model:', model_name])
            to_write.writerow([])
            to_write.writerow(['Parameter:'] + model.get_model_params())
            to_write.writerow(['5th percentile:'] + percentile_5)
            to_write.writerow(['16th percentile:'] + percentile_16)
            to_write.writerow(['Median:'] + medians)
            to_write.writerow(['84th percentile:'] + percentile_84)
            to_write.writerow(['95th percentile:'] + percentile_95)
            to_write.writerow(['Mode:'] + modes)
            to_write.writerow([])
            to_write.writerow(['Elements:'] + elements_list)
            to_write.writerow(['Input Abundances:'] + wd_abundances_list)
            to_write.writerow(['Input Errors:'] + wd_errors_list)
            to_write.writerow(['Input Upper Bounds:'] + wd_abundance_upper_bounds_list)
            to_write.writerow(['Input Lower Bounds:'] + wd_abundance_lower_bounds_list)
            to_write.writerow(['Excluded Abundances:'] + wd_excluded_abundances_list)
            to_write.writerow(['Excluded Errors:'] + wd_excluded_errors_list)
            to_write.writerow(['Excluded Upper Bounds:'] + wd_excluded_abundance_upper_bounds_list)
            to_write.writerow(['Excluded Lower Bounds:'] + wd_excluded_abundance_lower_bounds_list)
            to_write.writerow(['Output from median run:'])
            to_write.writerow(['Disc Composition:'] + disc_comp)
            #to_write.writerow(['Partition Coefficients (X_c for Oxygen):'] + ds) # For old Oxygen logic
            to_write.writerow(['Partition Coefficients:'] + ds)
            to_write.writerow(['Bulk Composition:'] + bulk_comp)
            to_write.writerow(['Core Composition:'] + core_comp)
            to_write.writerow(['Mantle Composition:'] + mantle_comp)
            to_write.writerow(['Pollutant Composition:'] + pollutant_comp_list)
            to_write.writerow(['Final Result:'] + elements_results)
            #try:
            to_write.writerow(['Parent Core Number Fraction:', pcnf_median])
            #except TypeError:
            #    to_write.writerow(['Parent Core Number Fraction:'])
            to_write.writerow(['Fragment core mass fraction, +error, -error:', cmf_median, cmf_error_plus, cmf_error_minus])
            to_write.writerow(['Radius /km:', radius])
            to_write.writerow(['Mass /M_Earth:', mass])
            to_write.writerow(['Radius lower limit /km:', radius_lower_limit])
            to_write.writerow(['Mass lower limit /M_Earth:', mass_lower_limit])
            to_write.writerow(['delta time/Myrs, +error, -error:', tmedian, terrorplus, terrorminus])
            to_write.writerow(['Build Up % (sampled):', bu_percent])
            to_write.writerow(['Steady State % (sampled):', ss_percent])
            to_write.writerow(['Declining % (sampled):', dec_percent])
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
            for oa_strat, oa_stats in eoc_stats['Model prediction'].items():
                to_write.writerow(['With Oxidation Strategy:', oa_strat])
                to_write.writerow(['Using composition of:', oa_stats[9]])
                to_write.writerow(['Excess Oxygen:', oa_stats[0]])
                to_write.writerow(['Fractional Excess Oxygen:', oa_stats[1]])
                to_write.writerow(['Excess O Sigma Significance:', oa_stats[8]])
                to_write.writerow(['Water Mass:', oa_stats[3]])
                to_write.writerow(['Water Mass Fraction:', oa_stats[4]])
                to_write.writerow(['Excess Oxygen Mass:', oa_stats[5]])
                to_write.writerow(['Oxygen Assignation Elements:'] + oa_elements[oa_strat])
                to_write.writerow(['Oxygen Assignations:'] + oa_assignations[oa_strat])
            if eoc_stats.get('Model with O data point') is not None:
                for oa_strat, oa_stats in eoc_stats['Model prediction, data O'].items():
                    to_write.writerow(['Excess Oxygen using O data point:', oa_stats[0]])
            if eoc_stats.get('Model prediction, data (SS)') is not None:
                for oa_strat, oa_stats in eoc_stats['Model prediction, data (SS)'].items():
                    list_of_els = list()
                    list_of_el_abundances = list()
                    for el, a in oa_stats[6].items():
                        list_of_els.append(el)
                        list_of_el_abundances.append(a[gi.Layer.bulk])
                    to_write.writerow(['Elements using SS corrected data points (oxidation ' + str(oa_strat) + '):'] + list_of_els)
                    to_write.writerow(['Composition using SS corrected data points (oxidation ' + str(oa_strat) + '):'] + list_of_el_abundances)
                    to_write.writerow(['Excess Oxygen using SS corrected data points (oxidation ' + str(oa_strat) + '):', oa_stats[0]])
                    to_write.writerow(['Fractional Excess Oxygen using SS corrected data points (oxidation ' + str(oa_strat) + '):', oa_stats[1]])
                    to_write.writerow(['Excess O Sigma Significance using SS corrected data points (oxidation ' + str(oa_strat) + '):', oa_stats[8]])
                    to_write.writerow(['Water Mass using SS corrected data points (oxidation ' + str(oa_strat) + '):', oa_stats[3]])
                    to_write.writerow(['Water Mass Fraction using SS corrected data points (oxidation ' + str(oa_strat) + '):', oa_stats[4]])
                    to_write.writerow(['Excess Oxygen Mass using SS corrected data points (oxidation ' + str(oa_strat) + '):', oa_stats[5]])
            to_write.writerow([])
            to_write.writerow(['Excess Oxygen Sampling results:'])
            for ox_strat, stats in sampled_eo_stats.items():
                to_write.writerow([])
                to_write.writerow(['For Oxidation Strat: ' + str(ox_strat)])
                to_write.writerow(['Sampled Excess Oxygen, +error, -error:', stats[self.absolute_string][0], stats[self.absolute_string][1], stats[self.absolute_string][2]])
                to_write.writerow(['Sampled Excess Oxygen Sigma Significance:', stats[self.absolute_string][3]])
                to_write.writerow(['Sampled Oxygen Excess or Deficit Rate (whichever greater):', stats[self.absolute_string][6]])
                to_write.writerow(['Samples used:', str(stats[self.absolute_string][4]) + '/' + str(stats[self.absolute_string][5])])
                to_write.writerow(['Sampled Fractional Excess Oxygen, +error, -error:', stats[self.fractional_string][0], stats[self.fractional_string][1], stats[self.fractional_string][2]])
                to_write.writerow(['Sampled Fractional Excess Oxygen Sigma Significance:', stats[self.fractional_string][3]])
                to_write.writerow(['Sampled Fractional Oxygen Excess or Deficit Rate (whichever greater):', stats[self.fractional_string][6]])
                to_write.writerow(['Samples used (fractional):', str(stats[self.fractional_string][4]) + '/' + str(stats[self.fractional_string][5])])
            to_write.writerow([])
            to_write.writerow(['Excess Oxygen Semi-Sampling results:'])
            for key, value in semisampled_eo_dict.items():
                to_write.writerow([key, str(value)])
        if print_all_eo_sample_vals:
            sample_dump = stats_file.split('.')[0] + '_eo_samples_' + model_name + '.csv'
            with open(sample_dump, 'w', newline='', encoding='utf-8') as f:
                to_write = csv.writer(f)
                lists_to_zip = list()
                names_of_zipped = list()
                for ox_strat, stats in eo_sample_stats.items():
                    for key, list_to_zip in stats.items():
                        names_of_zipped.append(str(ox_strat) + ', ' + key)
                        lists_to_zip.append(list_to_zip)
                to_write.writerow(names_of_zipped)
                zipped = zip(*lists_to_zip)
                for z in zipped:
                    to_write.writerow(z)
        if print_all_time_sample_vals:
            sample_dump = stats_file.split('.')[0] + '_time_samples_' + model_name + '.csv'
            with open(sample_dump, 'w', newline='', encoding='utf-8') as f:
                to_write = csv.writer(f)
                to_write.writerow(['Accretion Timescale /Myr', 'Time Since Accretion /Myr', 'Delta /Myr'])
                zipped = zip(accretion_timescale_extract, t_sinceaccretion_extract, delta_time)
                for z in zipped:
                    to_write.writerow(z)

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

    def get_stats_weightpost_and_best_fit(self, model, observation_number, chains_dir):
        a = pn.Analyzer(n_params = model.get_n_dims(), outputfiles_basename = model.get_full_prefix(chains_dir, observation_number))
        stats = a.get_stats()
        weightpost = a.get_equal_weighted_posterior()[:, 0:model.get_n_dims()]  # This is basically just excluding the final column of the ...post_equal_weights.dat file
        best_fit = a.get_best_fit()
        return stats, weightpost, best_fit['log_likelihood'], best_fit['parameters']

    def extract_and_plot_excess_oxygen(self, model, observation_number, input_dict, system_name, observations, errors, wd_timescales=None, mass=1, suppress_graphical_output=False):
        eoc_stats = self.excess_oxygen_calculator.run_full_calculation(input_dict, mass, True, True, observations, errors, wd_timescales)
        if not suppress_graphical_output:
            self.graph_fac.make_excess_oxygen_plot(
                eoc_stats,
                system_name + '_' + model.get_prefix(observation_number)
            )
        return eoc_stats

    def make_corner_plot(self, model_name, model, observation_number, wd_name, chains_dir):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number, chains_dir)
        delta_free_params = list()
        for dfp in model.get_model_params():
            if 'Δ' in dfp:
                # TODO: See if replacing it with r'$\Delta$' solves it
                delta_free_params.append(dfp.replace('Δ', ''))  # The corner function doesn't handle this character when plt.rcParams['text.usetex'] = 'True' in dict_plotter
            else:
                delta_free_params.append(dfp)
        corner.corner(weightpost, labels=delta_free_params, label_kwargs=dict(fontsize=12), levels = (0.39346934,0.86466472,0.988891), smooth=True)
        prefix = self.graph_dir + wd_name + '_' + model.get_prefix(observation_number)
        plt.savefig(prefix + 'corner.pdf')  # TODO: do this within DictPlotter

    def make_composition_plot(self, chains_dir, model_name, model, observation_number, white_dwarf, enhancement_model, bonus_mk2_fits=None, bonus_mk2_error_lows=None, bonus_mk2_error_highs=None, suppress_graphical_output=False):
        reference_elements = [ci.Element.Mg, white_dwarf.get_atmospheric_type().value]

        fixed_pressure_version = True # Set this to true to produce a version of this plot with an extra line corresponding to a certain fixed pressure
        fixed_pressure = 0
        fixed_pressure_name = 'Low Pressure'

        include_best_fit = True

        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number, chains_dir)

        number_of_possible_samples = len(weightpost[:,0])
        max_number_of_samples = 10000
        number_of_samples = min(max_number_of_samples, number_of_possible_samples)
        samples = np.random.choice(number_of_possible_samples, number_of_samples, replace=False)
        len_x_axis = len(ci.usual_elements)
        fit_data = dict()
        model_output = dict()
        for element in reference_elements:
            fit_data[element] = dict()
            model_output[element] = np.zeros(shape=(number_of_samples, len_x_axis))
        #model_raw_output = np.zeros(shape=(number_of_samples, len_x_axis+1))
        sampled_disc_abundances = np.zeros(shape=(number_of_samples, len_x_axis))
        sampled_bulk_compositions = np.zeros(shape=(number_of_samples, len_x_axis))
        sampled_core_compositions = np.zeros(shape=(number_of_samples, len_x_axis))
        sampled_mantle_compositions = np.zeros(shape=(number_of_samples, len_x_axis))
        sampled_partition_compositions = np.zeros(shape=(number_of_samples, len_x_axis))
        sampled_pcnfs = list()

        input_values = list()
        parameter_indices = mp.parameter_indices(ld._live_model)
        geo_model_for_solar_abundances = gi.GeologyModel()  # We're only using this to access the solar abundances so no arguments needed
        eo_samples = {
            eoc.OxidationStrategy.default: {
                self.absolute_string: list(),
                self.fractional_string: list()
            },
            eoc.OxidationStrategy.conservative: {
                self.absolute_string: list(),
                self.fractional_string: list()
            }
        }

        if fixed_pressure_version:
            fp_model_output = dict()
            for element in reference_elements:
                fp_model_output[element] = np.zeros(shape=(number_of_samples, len_x_axis))
        print('Sampling for Composition plot')
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
            expressions, diagnostics = cm.complete_model_calculation(
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
                eo_samples[eoc.OxidationStrategy.default][self.absolute_string].append(None)
                eo_samples[eoc.OxidationStrategy.default][self.fractional_string].append(None)
                eo_samples[eoc.OxidationStrategy.conservative][self.absolute_string].append(None)
                eo_samples[eoc.OxidationStrategy.conservative][self.fractional_string].append(None)
                continue
            for ref_element in reference_elements:
                scaled_abundances = sa.scale_abundances_to_solar(expressions, ref_element)
                model_output[ref_element][i,:] = [scaled_abundances[element] for element in ci.usual_elements]
            try:
                eoc_comp = geo_model_for_solar_abundances.get_earth_mix(f_c, diagnostics['Enhancements']['Abundances'])
            except (TypeError, KeyError):
                eoc_comp = diagnostics['DiscAbundances']
            eoc_stats = self.excess_oxygen_calculator.run_full_calculation({'Sample ' + str(i): eoc_comp})
            try:
                eo_samples[eoc.OxidationStrategy.default][self.absolute_string].append(eoc_stats['Sample ' + str(i)][eoc.OxidationStrategy.default][0])
                eo_samples[eoc.OxidationStrategy.default][self.fractional_string].append(eoc_stats['Sample ' + str(i)][eoc.OxidationStrategy.default][1])
                eo_samples[eoc.OxidationStrategy.conservative][self.absolute_string].append(eoc_stats['Sample ' + str(i)][eoc.OxidationStrategy.conservative][0])
                eo_samples[eoc.OxidationStrategy.conservative][self.fractional_string].append(eoc_stats['Sample ' + str(i)][eoc.OxidationStrategy.conservative][1])
            except TypeError:
                # It was None
                eo_samples[eoc.OxidationStrategy.default][self.absolute_string].append(None)
                eo_samples[eoc.OxidationStrategy.default][self.fractional_string].append(None)
                eo_samples[eoc.OxidationStrategy.conservative][self.absolute_string].append(None)
                eo_samples[eoc.OxidationStrategy.conservative][self.fractional_string].append(None)
            sampled_disc_abundances[i,:] = [diagnostics['DiscAbundances'][element] for element in ci.usual_elements]
            if 'ParentCoreNumberFraction' in diagnostics['Enhancements']:
                sampled_pcnfs.append(diagnostics['Enhancements']['ParentCoreNumberFraction'])
            if 'Abundances' in diagnostics['Enhancements']:
                sampled_bulk_compositions[i,:] = [diagnostics['Enhancements']['Abundances'][element][gi.Layer.bulk] for element in ci.usual_elements]
                sampled_core_compositions[i,:] = [diagnostics['Enhancements']['Abundances'][element][gi.Layer.core] for element in ci.usual_elements]
                sampled_mantle_compositions[i,:] = [diagnostics['Enhancements']['Abundances'][element][gi.Layer.mantle] for element in ci.usual_elements]
            if 'Ds' in diagnostics['Enhancements']:
                sampled_partition_compositions[i,:] = [diagnostics['Enhancements']['Ds'][element] for element in ci.usual_elements]
            if fixed_pressure_version:
                for ref_element in reference_elements:
                    scaled_abundances = sa.scale_abundances_to_solar(fp_expressions, ref_element)
                    fp_model_output[ref_element][i,:] = [scaled_abundances[element] for element in ci.usual_elements]
        disc_abundance_medians = dict(zip(ci.usual_elements, self.confidence_intervals(sampled_disc_abundances)[3]))
        bulk_composition_medians = dict(zip(ci.usual_elements, self.confidence_intervals(sampled_bulk_compositions)[3]))
        core_composition_medians = dict(zip(ci.usual_elements, self.confidence_intervals(sampled_core_compositions)[3]))
        mantle_composition_medians = dict(zip(ci.usual_elements, self.confidence_intervals(sampled_mantle_compositions)[3]))
        partition_coefficient_medians = dict(zip(ci.usual_elements, self.confidence_intervals(sampled_partition_compositions)[3]))
        pcnf_median = np.percentile(sampled_pcnfs, 50)
        for ref_element in reference_elements:
            model_low3, model_low2, model_low1, model_median, model_high1, model_high2, model_high3 = self.confidence_intervals(model_output[ref_element])
            fit_data[ref_element]['Model median'] = (dict(zip(ci.usual_elements, model_median)), dict(zip(ci.usual_elements, model_low1)), dict(zip(ci.usual_elements, model_high1)))

            if fixed_pressure_version:
                fp_model_low3, fp_model_low2, fp_model_low1, fp_model_median, fp_model_high1, fp_model_high2, fp_model_high3 = self.confidence_intervals(fp_model_output[ref_element])
                fit_data[ref_element][fixed_pressure_name] = (dict(zip(ci.usual_elements, fp_model_median)), dict(zip(ci.usual_elements, fp_model_low1)), dict(zip(ci.usual_elements, fp_model_high1)))
                #if fixed_pressure_version:
                #    fp_output, fp_ignore = cm.complete_model_calculation(
                #        fe_star,
                #        t_sinceaccretion,
                #        d_formation,
                #        z_formation,
                #        N_c,
                #        N_o,
                #        f_c,
                #        f_o,
                #        pollutionfraction,
                #        10**log_t_disc,
                #        fixed_pressure,
                #        fO2,
                #        enhancement_model,
                #        t_formation
                #    )
                #    fit_data[ref_element][fixed_pressure_name + ' fit'] = fp_output
            for key, fit in bonus_mk2_fits.items():
                fit_dict[ref_element][key] = fit
            for key, error_lows in bonus_mk2_error_lows.items():
                error_low_dict[ref_element][key] = error_lows
            for key, error_highs in bonus_mk2_error_highs.items():
                error_high_dict[ref_element][key] = error_highs
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
            output, ignore = cm.complete_model_calculation(
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
            for ref_element in reference_elements:
                scaled_abundances = sa.scale_abundances_to_solar(output, ref_element)
                fit_data[ref_element]['Best fit'] = (scaled_abundances, None, None)
        if not suppress_graphical_output:
            for ref_element in reference_elements:
                self.graph_fac.make_composition_plot_mk3(
                    white_dwarf,
                    ci.usual_elements,
                    fit_data[ref_element],
                    ref_element,
                    model.get_prefix(observation_number)
                )
        return fit_data, eo_samples, disc_abundance_medians, bulk_composition_medians, core_composition_medians, mantle_composition_medians, partition_coefficient_medians, pcnf_median

    def get_untransformed_parameter_stats(self, chains_dir, model, observation_number, model_param):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number, chains_dir)
        parameter_indices = mp.parameter_indices(ld._live_model)
        try:
            param_index = parameter_indices[model_param]
        except KeyError:
            # If the parameter was not present, this calculation is irrelevant
            return None
        toret = weightpost[:, param_index]
        return toret

    def make_pressure_plot(self, chains_dir, model, model_name, wd_name, observation_number):
        p_stats = self.get_untransformed_parameter_stats(chains_dir, model, observation_number, mp.ModelParameter.pressure)
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

    def make_eo_distribution_plot(self, model, model_name, wd_name, observation_number, eo_stats):
        half_bin_size = 0.025
        for ox_strat, stats in eo_stats.items():
            for abs_or_frac in [self.absolute_string, self.fractional_string]:
                eo_stats_to_use = stats.get(abs_or_frac)
                if eo_stats_to_use is not None:
                    max_bin = -3
                    eo_stats_no_nans = [eo for eo in eo_stats_to_use if eo is not None and not np.isnan(eo)]
                    while max_bin < max(eo_stats_no_nans):
                        max_bin += 0.1
                    min_bin = 3
                    while min_bin > min(eo_stats_no_nans):
                        min_bin -= 0.1
                    bins, x_bar_centres = self.generate_bins_and_bar_centres(min_bin, max_bin, half_bin_size)

                    heights, bins2 = np.histogram(
                        eo_stats_no_nans,
                        bins,
                        density=True
                    )

                    nan_count = len(eo_stats_to_use) - len(eo_stats_no_nans)
                    text_to_write = str(ox_strat).capitalize() + ' oxidation'
                    if nan_count > 0:
                        text_to_write += ', ' + str(nan_count) + ' runs excluded'
                    text_dict = {
                        'ox_strat_string': {
                            'x_pos': min_bin + ((max_bin-min_bin)/10),
                            'text_string': text_to_write
                        }
                    }
                    file_suffix = 'eo_dist_' + str(ox_strat) + '_' + abs_or_frac
                    x_label = 'Excess Oxygen' if abs_or_frac == self.absolute_string else 'Fractional Excess Oxygen'
                    self.graph_fac.make_histogram(x_bar_centres, [heights], [wd_name], wd_name + '_' + model.get_prefix(observation_number), half_bin_size*2, 1.1, x_label, file_suffix, text_dict, None, None)
        return None

    def make_semisampled_eo_distibution_plot(self, chains_dir, model_name, model, observation_number, white_dwarf, resampled_fit_data, suppress_graphical_output=False):
        # Adapted from Marc Brouwers' code
        N_O = None
        half_bin_size = 0.025

        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number, chains_dir)

        #medians = self.get_param_vals_by_percentile(model, weightpost, 50)
        #predicted_median_abundances, diagnostics, input_values_used = self.get_predicted_abundances(model, medians)

        # I think really the underlying problem here is that the median abundances, and the abundances calculated by forward modelling with the medians might not be the same!
        # We should really just reuse the median abundances calculated from the composition resampling, which should hopefully never be None

        parameter_indices = mp.parameter_indices(ld._live_model)
        t_sinceaccretion_index = parameter_indices[mp.ModelParameter.t_sinceaccretion]
        accretion_timescale_index = parameter_indices[mp.ModelParameter.accretion_timescale]
        t_obs_data = weightpost[:,t_sinceaccretion_index]*1000000 # convert from Myr to yr
        t_event_data = 10**weightpost[:,accretion_timescale_index] # convert from log yr to yr

        assert len(t_obs_data) == len(t_event_data)
        N_samples = len(t_event_data)

        O_acc_def = np.zeros(N_samples)
        O_acc_cons = np.zeros(N_samples)
        O_acc_def_ss = np.zeros(N_samples)
        els_to_iterate = list(self.excess_oxygen_calculator.oxidation_dict[eoc.OxidationStrategy.default].keys())
        if ci.Element.O not in els_to_iterate:
            els_to_iterate.append(ci.Element.O)
        abundances_to_use = dict()
        errors_to_use = dict()
        for el in ci.usual_elements:
            relevant_data = white_dwarf.get_abundance(el)
            if relevant_data is None or relevant_data.data_point_type != wd.WhiteDwarfDataPointType.measurement:
                abundances_to_use[el] = resampled_fit_data[white_dwarf.get_atmospheric_type().value]['Model median'][0][el] # Really the model name should not be hardcoded as a magic value like this! It basically just needs to be the same name that we saved the resampled median fits to in make_composition_plot
                #abundances_to_use[el] = predicted_median_abundances[el]
                errors_to_use[el] = 0.4
            else:
                abundances_to_use[el] = relevant_data.value
                errors_to_use[el] = relevant_data.upper_error
            #if white_dwarf.get_abundanceall_wd_abundances.get(el, 0) != 0:
            #    abundances_to_use[el] = all_wd_abundances[el]
            #else:
            #    # If that element didn't have an abundance, or was only an upper bound, we'll use the abundance predicted by the model for the median parameter values
            #    abundances_to_use[el] = predicted_median_abundances[el]
            #if all_wd_abundance_errors.get(el, 0) != 0:
            #    errors_to_use[el] = all_wd_abundance_errors[el]
            #else:
            #    # If that element didn't have an abundance, or was only an upper bound, we'll assume a (large) error
            #    errors_to_use[el] = 0.4
        for el in els_to_iterate:
            mu = abundances_to_use[el]
            sigma = errors_to_use[el]
            log_el_hx =  np.random.normal(loc=mu, scale=sigma, size=N_samples)
            t_el = white_dwarf.timescale_dict[el]
            N_el = 10**log_el_hx / (t_el*(np.exp((np.minimum(t_event_data, t_obs_data)-t_obs_data)/t_el) - np.exp(-t_obs_data/t_el)))
            N_el_ss = 10**log_el_hx / t_el
            if el == ci.Element.O:
                N_O = N_el
                N_O_ss = N_el_ss
            else:
                O_acc_def += self.excess_oxygen_calculator.oxidation_dict[eoc.OxidationStrategy.default][el]*N_el
                O_acc_cons += self.excess_oxygen_calculator.oxidation_dict[eoc.OxidationStrategy.conservative][el]*N_el
                O_acc_def_ss += self.excess_oxygen_calculator.oxidation_dict[eoc.OxidationStrategy.default][el]*N_el_ss

        O_excess_def = (N_O-O_acc_def)/N_O
        O_excess_cons = (N_O-O_acc_cons)/N_O
        median_def = np.median(O_excess_def)
        median_cons = np.median(O_excess_cons)
        O_excess_def_ss = (N_O_ss - O_acc_def_ss)/N_O_ss

        bins_many, x_bar_centres = self.generate_bins_and_bar_centres(
            np.floor(min(min(O_excess_def), min(O_excess_cons), min(O_excess_def_ss))), # Possible TODO: This is sometimes unnecessarily (and stupidly) low, but truncating it will change the later maths
            np.ceil(max(max(O_excess_def), max(O_excess_cons), max(O_excess_def_ss))),
            half_bin_size
        )
        plot_scaling_factor = 100 # We want to plot percent rather than the raw fraction
        x_bar_centres_percent = [plot_scaling_factor*x for x in x_bar_centres]

        heights_default, bins2 = np.histogram(
            O_excess_def,
            bins_many,
            density=True
        )

        n_cons, edges_def = np.histogram(O_excess_cons, bins=bins_many)

        lower_sigma_def = np.median(O_excess_def) - bins_many[np.where(np.cumsum(heights_default)>0.158*np.sum(heights_default))][0]
        upper_sigma_def = bins_many[np.where(np.cumsum(heights_default)>0.841*np.sum(heights_default))][0] - np.median(O_excess_def)
        # stat.norm.ppf(x) is the area under the curve of Gaussian(mu=0, sigma=1) between -inf and x
        # Satisfies e.g. stat.norm.ppf(0.5) = 0, i.e. the sigma significance of a 50:50 event is 0
        sig_excess_def = abs(stat.norm.ppf(len(O_excess_def[np.where(O_excess_def>0)])/len(O_excess_def)))
        p_excess_def = len(O_excess_def[np.where(O_excess_def>0)])/len(O_excess_def)
        p_deficit_def = len(O_excess_def[np.where(O_excess_def<0)])/len(O_excess_def)
        sig_max_def = abs(stat.norm.ppf(1/len(O_excess_def)))
        if np.isinf(sig_excess_def):
            sig_excess_def = sig_max_def

        lower_sigma_cons = np.median(O_excess_cons) - bins_many[np.where(np.cumsum(n_cons)>0.158*np.sum(n_cons))][0]
        upper_sigma_cons = bins_many[np.where(np.cumsum(n_cons)>0.841*np.sum(n_cons))][0] - np.median(O_excess_cons)

        median_def_ss = np.median(O_excess_def_ss)
        n_def_ss, edges_ss = np.histogram(O_excess_def_ss, bins=bins_many)
        lower_sigma_def_ss = np.median(O_excess_def_ss) - bins_many[np.where(np.cumsum(n_def_ss)>0.158*np.sum(n_def_ss))][0]
        upper_sigma_def_ss = bins_many[np.where(np.cumsum(n_def_ss)>0.841*np.sum(n_def_ss))][0] - np.median(O_excess_def_ss)

        text_dict = {
            'sigma_display_text': {
                'x_pos': 50 if p_excess_def > p_deficit_def else -100,
                'text_string': '$' + str(np.round(sig_excess_def,1)) + '\sigma$',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 16,
                'y_pos': 0.2*max(heights_default)
            }
        }
        text_dict_thesis = {
            'sigma_display_text': {
                'x_pos': 0.8 if p_excess_def > p_deficit_def else 0.2,
                'text_string': '$' + str(np.round(sig_excess_def,1)) + '\sigma$',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 16,
                'position_text_relative_to_plot': True,
                'y_pos': 0.2
            },
            'name_text': {
                'x_pos': 0.02,
                'text_string': white_dwarf.name,
                'horizontalalignment': 'left',
                'verticalalignment': 'center',
                'fontsize': 18,
                'position_text_relative_to_plot': True,
                'y_pos': 0.65
            }
        }

        line_dict = {
            'vline1': {
                'x_start': 0
            }
        }
        if not suppress_graphical_output:
            self.graph_fac.make_histogram(
                x_bar_centres_percent,
                [heights_default],
                [white_dwarf.name],
                white_dwarf.name + '_' + model.get_prefix(observation_number),
                half_bin_size*2*plot_scaling_factor,
                1.2,
                'Excess Oxygen (\%)',
                'semisampled_oxygen_excess',
                text_dict,
                line_dict,
                {'supress_second_axis': True, 'x_min': -150, 'x_max': 100},
                1.25,
                dict(),
                False,
                ['#ADD3E6'] if sig_excess_def > 2 else ['#979797']
            )

            self.graph_fac.make_histogram(
                x_bar_centres_percent,
                [heights_default],
                [white_dwarf.name],
                white_dwarf.name + '_' + model.get_prefix(observation_number),
                half_bin_size*2*plot_scaling_factor,
                1.2,
                'Excess Oxygen (\%)',
                'semisampled_oxygen_excess_thesis',
                text_dict_thesis,
                line_dict,
                {'supress_second_axis': True, 'x_min': -150, 'x_max': 100},
                1.25,
                dict(),
                False,
                ['#ADD3E6'] if sig_excess_def > 2 else ['#979797'],
                True
            )

        semisampled_eo_dict = cn.OrderedDict({
            'Sigma excess (default)': sig_excess_def,
            'P(excess) (default)': p_excess_def,
            'P(deficit) (default)': p_deficit_def,
            'Sigma max (default)': sig_max_def,
            'Median fractional excess (default)': median_def,
            'Upper error fractional excess (default)': upper_sigma_def,
            'Lower error fractional excess (default)': lower_sigma_def,

            #'O_excess_def': O_excess_def,
            #'O_excess_cons': O_excess_cons,

            'Median fractional excess (conservative)': median_cons,
            'Upper error fractional excess (conservative)': upper_sigma_cons,
            'Lower error fractional excess (conservative)': lower_sigma_cons,

            'Median fractional excess (default, ss)': median_def_ss,
            'Upper error fractional excess (default, ss)': upper_sigma_def_ss,
            'Lower error fractional excess (default, ss)': lower_sigma_def_ss
        })

        return semisampled_eo_dict

    def build_model(self, hierarchy_name, hierarchy_level_list):
        # Based on the equivalent function in manager.py
        # Make a pollution model consisting of the parameters in the relevent hierarchy levels of the parameter hierarchy
        model_name = pu.hierarchy_abbreviations[hierarchy_name]
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
            hierarchy_abbreviation = filename.split(system_name + '_')[1].split('_lv')[0]
            hierarchy_level_list_str = filename.split('lv_')[1].split('_')[0]
            hierarchy_level_list = [int(c) for c in hierarchy_level_list_str]
            model_name = hierarchy_abbreviation + '_lv_' + hierarchy_level_list_str
            for name, abbrev in pu.hierarchy_abbreviations.items():
                if abbrev == hierarchy_abbreviation:
                    hierarchy_name = name
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
        return pm.PollutionModel(model_name, enhancement_model_name, prior_name, live_points)

    def make_variable_distribution_plot(self, variable_name, variable_values_dict, x_min, x_max, half_bin_size, text_dict=None, file_prefix=None, cumulative=False, text_size_dict=dict(), weights_dict=dict(), comparison_dist_dict=dict()):
        bins, x_bar_centres = self.generate_bins_and_bar_centres(x_min, x_max, half_bin_size)

        all_heights = list()
        names = list()
        for config_name, variable_values in variable_values_dict.items():
            heights, bins2 = np.histogram(
                variable_values,
                bins,
                weights=weights_dict.get(config_name),
                density=True
            )
            names.append(config_name)
            all_heights.append(heights)

        file_suffix = '_' + str(variable_name).replace('/', '').replace(' ', '')
        xlabel = str(variable_name)
        if isinstance(variable_name, ci.Element):
            xlabel = 'log(' + str(variable_name) + '/Hx)'
        if isinstance(variable_name, tuple):
            # Then these should be 2 elements
            if isinstance(variable_name[0], ci.Element) and isinstance(variable_name[1], ci.Element):
                xlabel = 'log(' + str(variable_name[0]) + '/' + str(variable_name[1]) + ')'
                file_suffix = '_' + str(variable_name[0]) + '_' + str(variable_name[1])
                #if variable_name[0] == ci.Element.Ca and variable_name[1] == ci.Element.Fe:
                #    # Add Hollands et al 2017 data for comparison
                #    heights, bins2 = np.histogram(
                #    ha.get_hollands_ca_fe_values(),
                #    bins,
                #    density=True
                #)
                #names.append('Hollands')
                #all_heights.append(heights)
        for comparison_dist_name, comparison_dist in comparison_dist_dict.items():
            names.append(comparison_dist_name)
            heights, bins2 = np.histogram(
                comparison_dist,
                bins,
                density=True
            )
            all_heights.append(heights)
        hist_plot = self.graph_fac.make_histogram(x_bar_centres, all_heights, names, file_prefix, half_bin_size*2, 0.8, xlabel, file_suffix, text_dict, None, None, 1.1, text_size_dict, cumulative)
        return hist_plot

    def make_multisystem_pressure_plot(self, chains_dir, observation_numbers, labels=None, additional_text_dict=None):
        p_stats_list = list()
        used_observation_numbers_list = list()
        for observation_number in observation_numbers:
            recreated_model = self.recreate_model(observation_number)
            p_stats = self.get_untransformed_parameter_stats(chains_dir, recreated_model, observation_number, mp.ModelParameter.pressure)
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

    def get_temperature_stats(self, chains_dir, model, observation_number):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number, chains_dir)
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

    def make_temperature_plot(self, chains_dir, model, model_name, wd_name, observation_number, suppress_graphical_output=False):
        temp_stats = self.get_temperature_stats(chains_dir, model, observation_number)
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
        if not suppress_graphical_output:
            self.graph_fac.make_histogram(
                x_bar_centres,
                [heights],
                [wd_name],
                wd_name + '_' + model.get_prefix(observation_number),
                half_bin_size*2,
                1.1,
                'Formation Temperature /K',
                'temp_dist',
                text_dict,
                line_dict
            )
        return temp_stats

    def get_mass_stats(self, chains_dir, model, observation_number, enhancement_model):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number, chains_dir)
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
        logmassofpollutant = np.log10(massofpollutant)
        return logmassofpollutant

    def make_mass_plot(self, chains_dir, model, model_name, wd_name, observation_number, enhancement_model, suppress_graphical_output=False):
        logmassofpollutant = self.get_mass_stats(chains_dir, model, observation_number, enhancement_model)

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
        if not suppress_graphical_output:
            self.graph_fac.make_histogram(
                x_bar_centres,
                [heights],
                [wd_name],
                wd_name + '_' + model.get_prefix(observation_number),
                half_bin_size*2,
                1.2,
                'Log(Mass of Pollutant/kg)',
                'mass_dist',
                text_dict,
                None
            )
        return logmassofpollutant # This is just so that the later functions don't need to recalculate it!

    def make_time_plot(self, chains_dir, model, model_name, wd_name, observation_number):
        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number, chains_dir)
        parameter_indices = mp.parameter_indices(ld._live_model)
        t_sinceaccretion_index = parameter_indices[mp.ModelParameter.t_sinceaccretion]
        accretion_timescale_index = parameter_indices[mp.ModelParameter.accretion_timescale]
        x_data = np.log10(weightpost[:,t_sinceaccretion_index]*1000000)
        y_data = weightpost[:,accretion_timescale_index]
        self.graph_fac.plot_timesince_v_accretiontime(x_data, y_data, wd_name, model.get_prefix(observation_number), ld._live_t_mg)

    def make_scaled_tevent_plot(self, chains_dir, model, model_name, white_dwarf, observation_number):
        # Based on Marc Brouwers' code, used for Brouwers et al 2022b
        wd_name = white_dwarf.name
        t_ss_scaling_factor = 5 # Marc used 3
        steady_state_scaling_factor = 1 # The steady state panel takes up steady_state_scaling_factor*t_ss width on the plot. Marc used 1
        half_bin_size = 0.1
        critical_element = ci.Element.Mg # Marc used O

        stats, weightpost, best_fit_lnz, best_fit_params = self.get_stats_weightpost_and_best_fit(model, observation_number, chains_dir)
        parameter_indices = mp.parameter_indices(ld._live_model)
        t_sinceaccretion_index = parameter_indices[mp.ModelParameter.t_sinceaccretion]
        accretion_timescale_index = parameter_indices[mp.ModelParameter.accretion_timescale]
        t_obs_data = weightpost[:,t_sinceaccretion_index]*1000000 # convert from Myr to yr
        t_event_data = 10**weightpost[:,accretion_timescale_index] # convert from log yr to yr

        t_element = white_dwarf.timescale_dict[critical_element]

        t_ss = t_ss_scaling_factor*t_element # time to reach steady state

        # This bit calculates a scaled time coordinate for all possible phases of accretion (apart from build-up because that's just t_obs_data). We will pick the right one to use later
        scaled_ss_times = ((steady_state_scaling_factor*t_ss)*((t_obs_data-t_ss)/(t_event_data-t_ss))) + t_ss # i.e., what fraction of the way through steady state are we? Then rescale by arbitrary multiplicative factor so it looks ok on the plot, then add a non-arbitrary offset so that it goes after the build-up phase
        scaled_declining_times = (t_obs_data-t_event_data) + t_ss + (steady_state_scaling_factor*t_ss) # i.e., how far into declining phase are you? Now add the size of the build-up and steady state panels

        # This bit picks out the indices to which each possible accretion phase applies
        i_buildup = np.where((t_obs_data <= t_ss) & (t_obs_data <= t_event_data)) # build-up defined as t_ss after start.
        i_ss = np.where((t_obs_data > t_ss) & (t_obs_data <= t_event_data))
        i_declining = np.where(t_obs_data > t_event_data)

        scaled_times = np.concatenate((t_obs_data[i_buildup]/t_element, scaled_ss_times[i_ss]/t_element, scaled_declining_times[i_declining]/t_element), axis=None) # This picks out the correct phase for each sample and puts them all together in one array

        try:
            max_bin = max(np.amax(scaled_declining_times[i_declining]/t_element), (2 + steady_state_scaling_factor)*t_ss_scaling_factor)
        except ValueError:
            # Probably means there were no declining phase samples
            max_bin = (2 + steady_state_scaling_factor)*t_ss_scaling_factor
        bins, x_bar_centres = self.generate_bins_and_bar_centres(0, max_bin, half_bin_size)
        heights, bins2 = np.histogram(
            scaled_times,
            bins,
            density=True
        )

        names = [wd_name]
        text_dict = {
            'bu_text': {
                'x_pos': t_ss_scaling_factor/2,
                'text_string': str(np.round(100*len(t_obs_data[i_buildup])/len(t_obs_data),1)) + '\%',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 20,
                'y_pos': max(heights)*0.4
            },
            'ss_text': {
                'x_pos': t_ss_scaling_factor + steady_state_scaling_factor*t_ss_scaling_factor*0.5,
                'text_string': str(np.round(100*len(t_obs_data[i_ss])/len(t_obs_data),1)) + '\%',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 20,
                'y_pos': max(heights)*0.4
            },
            'dec_text': {
                'x_pos': (1.5 + steady_state_scaling_factor)*t_ss_scaling_factor,
                'text_string': str(np.round(100*len(t_obs_data[i_declining])/len(t_obs_data),1)) + '\%',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 20,
                'y_pos': max(heights)*0.4
            },
            'bu_display_text': {
                'x_pos': t_ss_scaling_factor/2,
                'text_string': 'Build-up',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 16,
                'y_pos': 0.8*max(heights)
            },
            'ss_display_text': {
                'x_pos': t_ss_scaling_factor + steady_state_scaling_factor*t_ss_scaling_factor*0.5,
                'text_string': 'Steady State',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 16,
                'y_pos': 0.8*max(heights)
            },
            'dec_display_text': {
                'x_pos': (1.5 + steady_state_scaling_factor)*t_ss_scaling_factor,
                'text_string': 'Declining',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 16,
                'y_pos': 0.8*max(heights)
            }
        }
        text_dict_thesis = {
            'name_text': {
                'x_pos': 0.02,
                'text_string': wd_name,
                'horizontalalignment': 'left',
                'verticalalignment': 'center',
                'fontsize': 18,
                'position_text_relative_to_plot': True,
                'y_pos': 0.65
            },
            'bu_text': {
                'x_pos': 0.17,
                'text_string': str(np.round(100*len(t_obs_data[i_buildup])/len(t_obs_data),1)) + '\%',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 20,
                'position_text_relative_to_plot': True,
                'y_pos': 0.2
            },
            'ss_text': {
                'x_pos': 0.5,
                'text_string': str(np.round(100*len(t_obs_data[i_ss])/len(t_obs_data),1)) + '\%',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 20,
                'position_text_relative_to_plot': True,
                'y_pos': 0.2
            },
            'dec_text': {
                'x_pos': 0.83,
                'text_string': str(np.round(100*len(t_obs_data[i_declining])/len(t_obs_data),1)) + '\%',
                'horizontalalignment': 'center',
                'verticalalignment': 'center',
                'fontsize': 20,
                'position_text_relative_to_plot': True,
                'y_pos': 0.2
            }
        }
        line_dict = {
            'vline1': {
                'x_start': t_ss_scaling_factor
            },
            'vline2': {
                'x_start': t_ss_scaling_factor*(1 + steady_state_scaling_factor)
            }
        }

        x_axis_dict = {
            'supress_second_axis': True,
            'x_tick_locations_override': [0, t_ss_scaling_factor, t_ss_scaling_factor*(1 + steady_state_scaling_factor), t_ss_scaling_factor*(2 + steady_state_scaling_factor)],
            'x_tick_labels_override': ['0', str(t_ss_scaling_factor) + r'$\tau_{' + str(critical_element) +'}$', r'$t_{event}$', r'$t_{event}$ + ' + str(t_ss_scaling_factor) + r'$\tau_{' + str(critical_element) +'}$'],
        }

        x_axis_dict_thesis = {
            'supress_second_axis': True,
            'x_tick_locations_override': [0, t_ss_scaling_factor, t_ss_scaling_factor*(1 + steady_state_scaling_factor), t_ss_scaling_factor*(2 + steady_state_scaling_factor)],
            'x_tick_labels_override': ['0', str(t_ss_scaling_factor) + r'$\tau_{' + str(critical_element) +'}$', r'$t_{event}$', r'$t_{event}$ + ' + str(t_ss_scaling_factor) + r'$\tau_{' + str(critical_element) +'}$'],
            'x_min': 0,
            'x_max': 15
        }

        self.graph_fac.make_histogram(
            x_bar_centres,
            [heights],
            names,
            wd_name + '_' + model.get_prefix(observation_number),
            half_bin_size*2,
            1.2,
            'Time, $t$',
            'scaled_times',
            text_dict,
            line_dict,
            x_axis_dict,
            1.25,
            dict(),
            False,
            ['#979797']
        )

        self.graph_fac.make_histogram(
            x_bar_centres,
            [heights],
            names,
            wd_name + '_' + model.get_prefix(observation_number),
            half_bin_size*2,
            1.2,
            'Time, $t$',
            'scaled_times_thesis',
            text_dict_thesis,
            line_dict,
            x_axis_dict_thesis,
            1.25,
            dict(),
            False,
            ['#979797'],
            True
        )

    def make_p_v_fO2_plot(self, chains_dir, model, model_name, wd_name, observation_number):
        x_data = self.get_untransformed_parameter_stats(chains_dir, model, observation_number, mp.ModelParameter.pressure)
        y_data = self.get_untransformed_parameter_stats(chains_dir, model, observation_number, mp.ModelParameter.oxygen_fugacity)
        if x_data is None or y_data is None:
            print('Abandoning pressure v fO2 plot (no pressure and/or fO2 stats available)')
            return
        self.graph_fac.plot_pressure_v_oxygen_fugacity(x_data, y_data, wd_name, model.get_prefix(observation_number))

