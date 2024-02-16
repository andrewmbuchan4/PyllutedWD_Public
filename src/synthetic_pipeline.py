#!/usr/bin/env python
# -*- coding: utf-8 -*-

from argparse import Namespace
from enum import Enum
import numpy as np
import os
import scipy.stats as st
import sys

import chemistry_info as ci
import forty_parsec_abundances as fp
import graph_factory as gf
import hollands_abundances as ha
import manager as mn
import model_analyser as ma
import model_parameters as mp
import multivariate_tests as mv
import mwdd_abundances as mwdd
import pwd_utils as pu
import synthetic_modeller as sm
import synthetic_population as sp
import synthetic_observer as so
import timescale_interpolator as ti

sys.path.append(pu.get_path_to_utils())

import dict_plotter as dp

run_mv_tests = mv.imported_correctly

class Pipeline:

    def __init__(self, name, population_parameter_dict, observer_dict, modeller_dict):
        self.name = name
        self.population_parameter_dict = population_parameter_dict
        self.population_dict = dict()
        self.observer_dict = observer_dict
        self.modeller_dict = modeller_dict
        self.base_dir = pu.get_path_to_pipeline_base_dir()
        self.results_dict = dict()
        self.p_threshold = 0.05
        self.timescale_interpolator = ti.TimescaleInterpolator()
        self.manager = mn.Manager(  # Load manager once here, and send it to hollands_abundances, rather than load a new one every time
            Namespace(
                wd_data_filename='WDInputData.csv',
                stellar_compositions_filename='StellarCompositionsSortFE.csv',
                n_live_points = 0,
                pollution_model_names=['Model_24'],
                enhancement_model='Earthlike',
                base_dir=pu.get_path_to_data()
            )
        )
        self.ternary_scaling_factors = {
            ci.Element.Ca: 10,
            ci.Element.Fe: 1,
            ci.Element.Mg: 1
        }

    def get_pipeline_dir(self, pop_name, obs_name, mod_name):
        path_name = self.get_base_pipeline_dir() + self.get_pipeline_name(pop_name, obs_name, mod_name) + '/'
        os.makedirs(path_name, exist_ok=True)
        return path_name

    def get_base_pipeline_dir(self):
        path_name = self.base_dir + self.name + '/'
        os.makedirs(path_name, exist_ok=True)
        return path_name

    def get_pipeline_name(self, pop_name, obs_name, mod_name):
        return pop_name + '_' + obs_name + '_' + mod_name

    def get_popdump_dir(self):
        return self.base_dir + 'popdumps/'

    def create_mv_plots(self, base_pop1_name, base_pop2_name, errors, threshold_offsets, N_values, modeller, base_obs_name, pool=False, resample_count=50, comparison_et_value=None, most_polluted=None):
        if not run_mv_tests:
            return
        # dprob = probability of distinguishing
        graph_fac = gf.GraphFactory(self.get_base_pipeline_dir())

        # For this to be remotely practical, assume that the modeller name is always the same, and we're only interested in 'Modelled' output
        # (Unless populations are the same!)
        # Let's require that the obs names we are interested in should all end in errxpyz which is parsed as error = x.yz dex

        mv_p_values_dict = dict()
        mv_dprob_values_dict = dict()

        for test_type in mv.get_mv_test_method_dict():
            mv_p_values_dict[test_type] = dict()
            mv_dprob_values_dict[test_type] = dict()

        # Also assume that we change error size OR threshold offsets, not both. Need to detect which mode we're using:
        errors_or_thresholds = list()
        error_mode = None
        if len(errors) == 1:
            if len(threshold_offsets) == 1:
                # Give a warning!
                print('WARNING! It\'s not clear whether we\'re varying error at constant threshold offset, or varying threshold offset at constant error, because only one value was supplied for each')
                print('Assuming it\'s the former.')
                errors_or_thresholds = errors
                error_mode = True
            else:
                errors_or_thresholds = threshold_offsets
                error_mode = False
        else:
            assert len(threshold_offsets) == 1
            errors_or_thresholds = errors
            error_mode = True
        pop1_output_name = base_pop1_name + base_obs_name
        pop2_output_name = base_pop2_name + base_obs_name
        for test_type in mv.get_mv_test_method_dict():
            for N in N_values:
                mv_p_values_dict[test_type][N] = list()
                mv_dprob_values_dict[test_type][N] = list()
        for N in N_values:
            for et in errors_or_thresholds:
                if error_mode:
                    obs_name = base_obs_name + 'err' + str(et).replace('.', 'p')
                    comparison_obs_name = obs_name if comparison_et_value is None else base_obs_name + 'err' + str(comparison_et_value).replace('.', 'p')
                else:
                    obs_name = base_obs_name + 'off' + str(et).replace('-', 'm').replace('.', 'p')
                    comparison_obs_name = obs_name if comparison_et_value is None else base_obs_name + 'off' + str(comparison_et_value).replace('-', 'm').replace('.', 'p')
                if not pool:
                    raise NotImplementedError
                else:
                    # Then assume we are being directed to resample a pop of size N from a larger pool
                    # Assuming this becomes the standard protocol, can assume that the populations will no longer be named ...Nxyz
                    pop1_name = base_pop1_name
                    pop2_name = base_pop2_name
                    if pop2_name.startswith('Hollands'):
                        el_list_to_compare = list()
                        for hname in pop2_name.split('_'):
                            try:
                                el_list_to_compare.append(ci.Element[hname])
                            except KeyError:
                                pass
                        pop2_output_name = base_pop2_name
                        mv_combination1 = (modeller, obs_name, pop1_name, 'Observed', N)
                        hollands_data = ha.get_hollands_common_el_values(el_list_to_compare)
                        mv_results_dict, mv_distinguishable_dict = self.run_mv_tests_against_comparison(mv_combination1, 'Hollands', hollands_data, resample_count, most_polluted, el_list_to_compare)
                    elif pop2_name.startswith('40pc'):
                        raise NotImplementedError
                    elif comparison_et_value is not None:
                        raise NotImplementedError
                    else:
                        raise NotImplementedError
                for test_type in mv.get_mv_test_method_dict():
                    mv_dprob_values_dict[test_type][N].append(mv_distinguishable_dict[test_type][0]/mv_distinguishable_dict[test_type][1])
                    mv_p_values_dict[test_type][N].append(mv_results_dict[test_type])
        y_label = 'P(Input != Output)' if comparison_et_value is None else 'P(Output != No Noise)'
        graph_fac.make_dprob_v_error_plot(errors_or_thresholds, mv_dprob_values_dict, None, pop1_output_name, pop2_output_name, error_mode, y_label, True)
        graph_fac.make_ks_v_error_plot(errors_or_thresholds, mv_p_values_dict, self.p_threshold, None, pop1_output_name, pop2_output_name, error_mode, None, True)

    def run_mv_tests_against_comparison(self, combination, comparison_name, comparison_values, resample_count=1, most_polluted=None, filter_on_elements=None):
        # combinations should be a tuple containing:
        #      0         1         2          3                                                              4
        # (mod_name, obs_name, pop_name, variable_type ('Input', 'Pollution', 'Observed', 'Modelled'), population size)
        # (population_size can be omitted)

        normalise_abundance_values = True

        pop_including_unobserved = self.results_dict[combination[0]][combination[1]][combination[2]]
        if filter_on_elements is not None:
            # Take all systems that are have detections of the requested elements
            pop = pop_including_unobserved.get_subset_with_detected_elements(filter_on_elements)
        else:
            # Take all systems that are actually observed:
            pop = pop_including_unobserved.get_observed_subset()
        try:
            pop_size = combination[4]
        except IndexError:
            pop_size = None
        number_of_times_resampled = 0
        results_dict = dict()
        for test_type in mv.get_mv_test_method_dict():
            results_dict[test_type] = list()


        if normalise_abundance_values:
            comparison_values_normalised = self.normalise_abundance_lists(comparison_values)
            data_set1_r = mv.convert_data_to_r_format(comparison_values_normalised)
        else:
            data_set1_r = mv.convert_data_to_r_format(comparison_values)
        max_no_ternary_series = 2
        ternary_abundance_dict = {
            'Real': dict()
        }
        for i, el in enumerate(filter_on_elements):
            if normalise_abundance_values:
                ternary_abundance_dict['Real'][el] = comparison_values_normalised[i]
            else:
                ternary_abundance_dict['Real'][el] = comparison_values[i]
        while number_of_times_resampled < resample_count:
            if most_polluted is not None:
                # Then we will use (as our 1 sample) the most polluted white dwarfs - this roughly mimics Hollands et al 2017
                pop_size = most_polluted[0]
                pol_element = most_polluted[1]
                print('Sampling the ' + str(pop_size) + ' most polluted WDs (by ' + str(pol_element) + ')')
                pop_sample = pop.get_most_polluted_subset(pop_size, pol_element)
                resample_count = 1
            else:
                pop_sample = pop.get_random_subset(pop_size, True)
            if pop_sample is not None:
                if normalise_abundance_values:
                    sample_values_normalised = self.normalise_abundance_lists(pop_sample.get_common_el_values(filter_on_elements))
                    assert len(comparison_values_normalised[0]) == len(sample_values_normalised[0]) # Just a quick sanity check in case I forget to set the synthetic sample size to match the real sample size
                    data_set2_r = mv.convert_data_to_r_format(sample_values_normalised)
                else:
                    sample_values = pop_sample.get_common_el_values(filter_on_elements)
                    assert len(comparison_values[0]) == len(sample_values[0]) # Just a quick sanity check in case I forget to set the synthetic sample size to match the real sample size
                    data_set2_r = mv.convert_data_to_r_format(sample_values)
                if len(ternary_abundance_dict.keys()) < max_no_ternary_series:
                    ternary_abundance_dict['Synthetic'] = dict()
                    for i, el in enumerate(filter_on_elements):
                        if normalise_abundance_values:
                            ternary_abundance_dict['Synthetic'][el] = sample_values_normalised[i]
                        else:
                            ternary_abundance_dict['Synthetic'][el] = sample_values[i]
                for test_type, test_function in mv.get_mv_test_method_dict().items():
                    p_value = test_function(
                        data_set1_r,
                        data_set2_r
                    )
                    results_dict[test_type].append(p_value)
                number_of_times_resampled += 1
            else:
                break
        distinguishable_dict = self.report_mv_test_results(combination, comparison_name, pop_size, len(comparison_values[0]), results_dict, True)
        if len(filter_on_elements) == 3:
            graph_fac = gf.GraphFactory(self.get_pipeline_dir(combination[2], combination[1], combination[0]))
            graph_fac.make_ternary_plot(
                filter_on_elements,
                ternary_abundance_dict,
                'ternary_plot_sampled_' + combination[2] + '_' + combination[1],
                self.ternary_scaling_factors
                #additional_text_dict
            )
        return results_dict, distinguishable_dict

    def normalise_abundance_lists(self, abundance_lists):
        number_of_dimensions = len(abundance_lists)
        number_of_data_vectors = len(abundance_lists[0]) # Hopefully this is the same for all the lists

        toret = list()
        j = 0
        while j < number_of_dimensions:
            toret.append(list())
            j += 1

        i = 0
        while i < number_of_data_vectors:
            total_abundance = 0
            for abundance_list in abundance_lists:
                total_abundance += 10**(abundance_list[i])
            k = 0
            for k, abundance_list in enumerate(abundance_lists):
                toret[k].append(10**(abundance_list[i])/total_abundance)
            i += 1
        return toret

    def create_dprob_ks_v_error_plot(self, base_pop1_name, base_pop2_name, errors, threshold_offsets, N_values, modeller, base_obs_name, variable, plot_descriptions, pool=False, resample_count=50, comparison_et_value=None, most_polluted=None, filter_on_elements=None):
        # dprob = probability of distinguishing
        graph_fac = gf.GraphFactory(self.get_base_pipeline_dir())

        # For this to be remotely practical, assume that the modeller name is always the same, and we're only interested in 'Modelled' output
        # (Unless populations are the same!)
        # Let's require that the obs names we are interested in should all end in errxpyz which is parsed as error = x.yz dex

        dprob_values_dict = dict()
        ks_values_dict = dict()
        # Also assume that we change error size OR threshold offsets, not both. Need to detect which mode we're using:
        errors_or_thresholds = list()
        error_mode = None
        if len(errors) == 1:
            if len(threshold_offsets) == 1:
                # Give a warning!
                print('WARNING! It\'s not clear whether we\'re varying error at constant threshold offset, or varying threshold offset at constant error, because only one value was supplied for each')
                print('Assuming it\'s the former.')
                errors_or_thresholds = errors
                error_mode = True
            else:
                errors_or_thresholds = threshold_offsets
                error_mode = False
        else:
            assert len(threshold_offsets) == 1
            errors_or_thresholds = errors
            error_mode = True
        pop1_output_name = base_pop1_name + base_obs_name
        pop2_output_name = base_pop2_name + base_obs_name
        for N in N_values:
            dprob_values_dict[N] = list()
            ks_values_dict[N] = list()
            for et in errors_or_thresholds:
                if error_mode:
                    obs_name = base_obs_name + 'err' + str(et).replace('.', 'p')
                    comparison_obs_name = obs_name if comparison_et_value is None else base_obs_name + 'err' + str(comparison_et_value).replace('.', 'p')
                else:
                    obs_name = base_obs_name + 'off' + str(et).replace('-', 'm').replace('.', 'p')
                    comparison_obs_name = obs_name if comparison_et_value is None else base_obs_name + 'off' + str(comparison_et_value).replace('-', 'm').replace('.', 'p')
                if not pool:
                    pop1_name = base_pop1_name if len(N_values) == 1 else base_pop1_name + 'N' + str(N)
                    pop2_name = base_pop2_name if len(N_values) == 1 else base_pop2_name + 'N' + str(N)

                    if pop2_name.startswith('Hollands'):
                        hollands_el1_el2 = pop2_name.split('_')
                        el1 = ci.Element[hollands_el1_el2[1]]
                        el2 = ci.Element[hollands_el1_el2[2]]
                        pop2_output_name = base_pop2_name
                        ks_combination1 = (modeller, obs_name, pop1_name, 'Observed')
                        ks_results, distinguishable_count, total_count = self.run_ks_test_against_comparison(ks_combination1, 'Hollands', ha.get_hollands_el_el_values(el1, el2, self.manager), variable, plot_descriptions[variable], most_polluted, filter_on_elements)
                    elif pop2_name.startswith('40pc'):
                        el1_el2 = pop2_name.split('_')
                        el1 = ci.Element[el1_el2[1]]
                        el2 = ci.Element[el1_el2[2]]
                        pop2_output_name = base_pop2_name
                        ks_combination1 = (modeller, obs_name, pop1_name, 'Observed')
                        ks_results, distinguishable_count, total_count = self.run_ks_test_against_comparison(ks_combination1, '40pc', fp.get_40pc_el_el_values(el1, el2, self.manager), variable, plot_descriptions[variable], most_polluted, filter_on_elements)
                    elif comparison_et_value is not None:
                        # Then we're being told to compare all outputs to the output of one specific run
                        ks_combination1 = (modeller, obs_name, pop1_name, 'Modelled')
                        ks_combination2 = (modeller, comparison_obs_name, pop2_name, 'Modelled')
                        ks_results, distinguishable_count, total_count = self.run_ks_test_on_results(ks_combination1, ks_combination2, variable, plot_descriptions[variable])
                    else:
                        ks_combination1 = (modeller, obs_name, pop1_name, 'Modelled') if pop1_name != pop2_name else (modeller, obs_name, pop1_name, 'Input')
                        ks_combination2 = (modeller, obs_name, pop2_name, 'Modelled')
                        ks_results, distinguishable_count, total_count = self.run_ks_test_on_results(ks_combination1, ks_combination2, variable, plot_descriptions[variable])
                else:
                    # Then assume we are being directed to resample a pop of size N from a larger pool
                    # Assuming this becomes the standard protocol, can assume that the populations will no longer be named ...Nxyz
                    pop1_name = base_pop1_name
                    pop2_name = base_pop2_name
                    if pop2_name.startswith('Hollands'):
                        hollands_el1_el2 = pop2_name.split('_')
                        el1 = ci.Element[hollands_el1_el2[1]]
                        el2 = ci.Element[hollands_el1_el2[2]]
                        pop2_output_name = base_pop2_name
                        ks_combination1 = (modeller, obs_name, pop1_name, 'Observed', N)
                        ks_results, distinguishable_count, total_count = self.run_ks_test_against_comparison(ks_combination1, 'Hollands', ha.get_hollands_el_el_values(el1, el2, self.manager), variable, plot_descriptions[variable], resample_count, most_polluted, filter_on_elements)
                    elif pop2_name.startswith('40pc'):
                        el1_el2 = pop2_name.split('_')
                        el1 = ci.Element[el1_el2[1]]
                        el2 = ci.Element[el1_el2[2]]
                        pop2_output_name = base_pop2_name
                        ks_combination1 = (modeller, obs_name, pop1_name, 'Observed', N)
                        ks_results, distinguishable_count, total_count = self.run_ks_test_against_comparison(ks_combination1, '40pc', fp.get_40pc_el_el_values(el1, el2), variable, plot_descriptions[variable], resample_count, most_polluted, filter_on_elements)
                    elif comparison_et_value is not None:
                        # Then we're being told to compare all outputs to the output of one specific run
                        ks_combination1 = (modeller, obs_name, pop1_name, 'Modelled', N)
                        ks_combination2 = (modeller, comparison_obs_name, pop2_name, 'Modelled', N)
                        ks_results, distinguishable_count, total_count = self.run_ks_test_on_results(ks_combination1, ks_combination2, variable, plot_descriptions[variable], resample_count)
                    else:
                        ks_combination1 = (modeller, obs_name, pop1_name, 'Modelled', N) if pop1_name != pop2_name else (modeller, obs_name, pop1_name, 'Input', N)
                        ks_combination2 = (modeller, obs_name, pop2_name, 'Modelled', N)
                        ks_results, distinguishable_count, total_count = self.run_ks_test_on_results(ks_combination1, ks_combination2, variable, plot_descriptions[variable], resample_count)
                dprob_values_dict[N].append(distinguishable_count/total_count)
                ks_values_dict[N].append(ks_results)
        y_label = 'P(Input != Output)' if comparison_et_value is None else 'P(Output != No Noise)'
        graph_fac.make_dprob_v_error_plot(errors_or_thresholds, dprob_values_dict, plot_descriptions[variable][0][0], pop1_output_name, pop2_output_name, error_mode, y_label)
        graph_fac.make_ks_v_error_plot(errors_or_thresholds, ks_values_dict, self.p_threshold, plot_descriptions[variable][0][0], pop1_output_name, pop2_output_name, error_mode)

    def create_dprob_ks_v_N_plot(self, base_pop1_name, base_pop2_name, errors, threshold_offsets, N_values, modeller, base_obs_name, variable, plot_descriptions, pool=False, resample_count=50):
        # dprob = probability of distinguishing
        graph_fac = gf.GraphFactory(self.get_base_pipeline_dir())

        # For this to be remotely practical, assume that the modeller name is always the same, and we're only interested in 'Modelled' output
        # (Unless populations are the same!)
        # Let's require that the obs names we are interested in should all end in errxpyz which is parsed as error = x.yz dex

        dprob_values_dict = dict()
        ks_values_dict = dict()
        dprob_grid = np.empty((len(errors), len(N_values)))

        for error in errors: # can generalise to errors or thresholds later!
            obs_name = base_obs_name + 'err' + str(error).replace('.', 'p')
            dprob_values_dict[error] = list()
            ks_values_dict[error] = list()
            for N in N_values:
                # Assume that we are always resampling a pop of size N from a larger pool
                ks_combination1 = (modeller, obs_name, base_pop1_name, 'Modelled', N) if base_pop1_name != base_pop2_name else (modeller, obs_name, base_pop1_name, 'Input', N)
                ks_combination2 = (modeller, obs_name, base_pop2_name, 'Modelled', N)
                ks_results, distinguishable_count, total_count = self.run_ks_test_on_results(ks_combination1, ks_combination2, variable, plot_descriptions[variable], resample_count)
                dprob = distinguishable_count/total_count
                dprob_values_dict[error].append(dprob)
                ks_values_dict[error].append(ks_results)
                dprob_grid[errors.index(error), N_values.index(N)] = dprob
        graph_fac.make_dprob_v_N_plot(N_values, dprob_values_dict, plot_descriptions[variable][0][0], base_pop1_name, base_pop2_name, 'P(Input != Output)')
        graph_fac.make_ks_v_N_plot(N_values, ks_values_dict, self.p_threshold, plot_descriptions[variable][0][0], base_pop1_name, base_pop2_name, False)
        graph_fac.make_dprob_N_error_heatmap(N_values, errors, dprob_grid, plot_descriptions[variable][0][0], base_pop1_name, base_pop2_name)

    def create_dprob_ks_v_N_multipop_plot(self, base_pop1_name, base_pop2_name, base_obs1_name, base_obs2_name, modeller1, modeller2, errors, N_values, N_scaling_factor, variable, plot_descriptions, resample_count=50):
        # dprob = probability of distinguishing
        graph_fac = gf.GraphFactory(self.get_base_pipeline_dir())

        # For this to be remotely practical, assume that the modeller name is always the same, and we're only interested in 'Modelled' output
        # (Unless populations are the same!)
        # Let's require that the obs names we are interested in should all end in errxpyz which is parsed as error = x.yz dex

        pop1_key = base_pop1_name + '_'  + base_obs1_name + '_' + modeller1
        pop2_key = base_pop2_name + '_'  + base_obs2_name + '_' + modeller2
        dprob_values_dict = {pop1_key: dict(), pop2_key: dict()}
        ks_values_dict = {pop1_key: dict(), pop2_key: dict()}
        N_values_dict = {pop1_key: dict(), pop2_key: dict()}

        for error in errors: # can generalise to errors or thresholds later!
            obs1_name = base_obs1_name + 'err' + str(error).replace('.', 'p')
            obs2_name = base_obs2_name + 'err' + str(error).replace('.', 'p')
            dprob_values_dict[pop1_key][error] = list()
            dprob_values_dict[pop2_key][error] = list()
            ks_values_dict[pop1_key][error] = list()
            ks_values_dict[pop2_key][error] = list()
            N_values_dict[pop1_key][error] = list()
            N_values_dict[pop2_key][error] = list()
            for N in N_values:
                # Assume that we are always resampling a pop of size N from a larger pool
                ks_combination_pop1_1 = (modeller1, obs1_name, base_pop1_name, 'Input', N)
                ks_combination_pop1_2 = (modeller1, obs1_name, base_pop1_name, 'Modelled', N)
                ks_results_pop1, distinguishable_count_pop1, total_count_pop1 = self.run_ks_test_on_results(ks_combination_pop1_1, ks_combination_pop1_2, variable, plot_descriptions[variable], resample_count)
                dprob_values_dict[pop1_key][error].append(distinguishable_count_pop1/total_count_pop1)
                ks_values_dict[pop1_key][error].append(ks_results_pop1)
                N_values_dict[pop1_key][error].append(N)
                pop2_sample_size = N*N_scaling_factor
                ks_combination_pop2_1 = (modeller2, obs2_name, base_pop2_name, 'Input', pop2_sample_size)
                ks_combination_pop2_2 = (modeller2, obs2_name, base_pop2_name, 'Modelled', pop2_sample_size)
                ks_results_pop2, distinguishable_count_pop2, total_count_pop2 = self.run_ks_test_on_results(ks_combination_pop2_1, ks_combination_pop2_2, variable, plot_descriptions[variable], resample_count)
                dprob_values_dict[pop2_key][error].append(distinguishable_count_pop2/total_count_pop2)
                ks_values_dict[pop2_key][error].append(ks_results_pop2)
                N_values_dict[pop2_key][error].append(pop2_sample_size)
        for error in errors:
            graph_fac.make_dprob_v_N_multipop_plot(N_values_dict, dprob_values_dict, plot_descriptions[variable][0][0], pop1_key, pop2_key, error, 'P(Input != Output)')
            graph_fac.make_ks_v_N_multipop_plot(N_values_dict, ks_values_dict, self.p_threshold, plot_descriptions[variable][0][0], pop1_key, pop2_key, error, False, None, False)

    # This could be a @staticmethod    ?
    def run_ks_test_on_samples(self, data1, data2, min_bin_edge, max_bin_edge, half_bin_size, precalculated_weights1=None, precalculated_weights2=None):
        n1 = len(data1)
        n2 = len(data2)
        bin_size = 2*half_bin_size
        bins = np.arange(min_bin_edge, max_bin_edge + bin_size, bin_size)
        cdf1 = self.convert_sample_to_cdf(data1, bins, precalculated_weights1)
        cdf2 = self.convert_sample_to_cdf(data2, bins, precalculated_weights2)
        return self.run_ks_test_on_cdf(cdf1, cdf2, n1, n2)

    def convert_sample_to_cdf(self, data, bins, precalculated_weights=None):
        weights_to_use = list()
        vals = list()
        if precalculated_weights is not None:
            # Assume data and precalculated_weights are already the vals and weights_to_use we need
            vals = data
            weights_to_use = precalculated_weights
        else:
            for val in data:
                # Expect each item to be a list. Need to collapse this list of lists into a single list, but keep track of the weight each should be assigned
                # (eg If a system has 3 possible output values, they need to each by weighted by 1/3 to avoid distorting the histogram)
                if isinstance(val, list):
                    weight = 1/len(val)
                    for entry in val:
                        vals.append(entry)
                        weights_to_use.append(weight)
                else:
                    # Assume it's a scalar
                    vals.append(val)
                    weights_to_use.append(1)
        heights, bins_ignore = np.histogram(
            vals,
            bins,
            weights=weights_to_use,
            density=False
        )
        cdf = np.cumsum(heights)  # Ideally we would avoid this by setting density to True, but when I tested that, it sometimes weighted things incorrectly
        scaling_factor = cdf[-1]
        cdf = np.divide(cdf, scaling_factor)
        return cdf

    def run_ks_test_on_cdf(self, cdf1, cdf2, n1, n2):
        # NB: The cdfs must have been binned the same way!!

        # This function has caused me some Bother
        # Basically, the scipy implementation of ks_2samp changed from v1.3.0 onwards
        # For previous versions, it would always calculate an approximation
        # My implementation here in run_ks_test_on_cdf is based on this approximation
        # From 1.3.0 onwards, it only uses the approximation for large samples ( > 10,000)
        # and calculates an exact value for smaller samples
        # One of the tests (test_ks_tests) checks the output of run_ks_test_on_cdf against the
        # output of ks_2samp (on a small sample)
        # Therefore, if the installed version of scipy is >= 1.3.0, there will be a discrepancy between run_ks_test_on_cdf and ks_2samp, because run_ks_test_on_cdf will calculate an approximate value and ks_2samp will calculate the exact value
        # The release notes (https://docs.scipy.org/doc/scipy/release/1.3.0-notes.html) claim that calling ks_2samp with alternative="two-sided", method="asymp" should replicate the old, approximate method
        # But when I tried this, there was still a discrepancy. I imagine probably caused by different internal cdf conversion logic?
        # (also note there is a separate bug in ks_2samp for versions between 1.3.0 and 1.5.0 which results in wildly incorrect p-values sometimes - avoid these versions!!!)

        # Further comments, following a bit more testing (see plot_trial_ks_tests.py)
        # The discrepency isn't too bad
        # My implementation seems to calculate slightly higher p-values than ks_2samp
        # But the difference is small (certainly not an order-of-magnitude effect)
        # So the far bigger caveat to these tests is the arbitrary significance threshold
        # If anything, overestimating the p-values is conservative
        # And when the samples are actually from the same dist, ks_samp does seem to sometimes produce p-values that look too low by eye

        # For future reference: the reason I didn't just use ks_2samp in the first place is that it takes, as input, 2 samples of data, and it is not expecting any of these data points to be multi-valued
        # But various outputs from the synthetic modeller are lists, e.g. if there are multiple best fit nebular compositions, the inferred Fe/H index is a list of all the best fit indices, assumed to be equally likely
        # So I dug into the source code and realised that it firstly converts to cdf, then performs the ks test on the cdf
        # So my solution was to write some logic to convert my outputs to a cdf, then duplicate the cdf-based ks test code from scipy
        if n1 == 0 or n2 == 0:
            return np.nan, np.nan
        d = np.max(np.absolute(cdf1 - cdf2))
        en = np.sqrt(n1 * n2 / (n1 + n2))
        prob = st.distributions.kstwobign.sf((en + 0.12 + 0.11 / en) * d)
        # The above line looks odd. But see https://stats.stackexchange.com/questions/264445/constants-in-scipy-implementation-of-kolmogorov-smirnov-two-sample-test
        return d, prob # This is the ks statistic and the p value

    def run_all_combinations(self):
        for mod_name, mod_description in self.modeller_dict.items():
            self.results_dict[mod_name] = dict()
            modeller = sm.Modeller(mod_description[0], mod_description[1])
            for obs_name, obs_description in self.observer_dict.items():
                self.results_dict[mod_name][obs_name] = dict()
                error_dict = dict()
                for el in ci.usual_elements:
                    if obs_description[1] > 0:
                        error_dict[el] = obs_description[1]
                try:
                    threshold_type = obs_description[2]
                except IndexError:
                    threshold_type = 'Default'
                observer = so.Observer(obs_description[0], error_dict, threshold_type)
                print('About to generate populations in pipeline ' + self.name)
                generated_pops = self.generate_populations(self.population_parameter_dict, obs_name, mod_name)
                print('About to observe in pipeline ' + self.name)
                observer.observe_populations(generated_pops)
                print('About to model in pipeline ' + self.name)
                modeller.model_populations(generated_pops)
                for pop_name, pop in generated_pops.items():
                    self.results_dict[mod_name][obs_name][pop_name] = pop
                    pop.assess_core_mantle_trustworthiness()

    def find_sample_size_for_given_distinguishability_probability(self, combination1, combination2, variable, plot_description, p_threshold=None):
        # combinations should be a tuple containing:
        #      0         1         2          3                                                              4
        # (mod_name, obs_name, pop_name, variable_type ('Input', 'Pollution', 'Observed', 'Modelled'), population size)
        # population_size is optional
        if p_threshold is None:
            p_threshold = self.p_threshold
        # variable, plot_description should be a key, value pair from plot_descriptions
        pop1_including_unobserved = self.results_dict[combination1[0]][combination1[1]][combination1[2]]
        pop2_including_unobserved = self.results_dict[combination2[0]][combination2[1]][combination2[2]]
        # Now take only the systems that are actually observed:
        pop1 = pop1_including_unobserved.get_subset_with_modelled_properties([mp.ModelParameter.fragment_core_frac]) # we only care about systems which have an fcf estimate
        pop2 = pop2_including_unobserved.get_subset_with_modelled_properties([mp.ModelParameter.fragment_core_frac])

        sample_size = 5 # For now, to keep simple, I'll just increase this until P goes above 90%
        distinguishability_probability_threshold = 0.9
        max_number_of_times_to_resample = 100
        sample_sizes_tried = list()
        while True:
            results_list = list()
            number_of_times_sampled = 0
            pop1.reset_unsampled_indices()
            pop2.reset_unsampled_indices()
            while number_of_times_sampled < max_number_of_times_to_resample:
                sample_sizes_tried.append(sample_size)
                pop1_sample = pop1.get_random_subset(sample_size, True)
                if combination1[0] == combination2[0] and combination1[1] == combination2[1] and combination1[2] == combination2[2]:
                    pop2_sample = pop1_sample  # Then they're actually the same population and we're comparing input and output, so the sample should be the same
                else:
                    pop2_sample = pop2.get_random_subset(sample_size, True)
                if pop1_sample is not None and pop2_sample is not None:
                    variable_dict1, weights_dict1, variable_below1, variable_above1 = self.extract_plottables(combination1[2], pop1_sample, {variable: plot_description})
                    variable_dict2, weights_dict2, variable_below2, variable_above2 = self.extract_plottables(combination2[2], pop2_sample, {variable: plot_description})
                    sample1 = variable_dict1[combination1[3]][variable][combination1[2]]
                    sample2 = variable_dict2[combination2[3]][variable][combination2[2]]

                    weights1 = weights_dict1[combination1[3]][variable][combination1[2]]
                    weights2 = weights_dict2[combination2[3]][variable][combination2[2]]

                    assert abs(sum(weights1) - sample_size) < 0.00001 # NB sample1 and sample 2 could be longer than sample_size if a system has multiple multiple values, should probably rather check the sum of weights
                    assert abs(sum(weights2) - sample_size) < 0.00001

                    if combination1[3] in self.get_io_type_variables():
                        assert combination2[3] in self.get_io_type_variables()
                        min_bin_edge = plot_description[0][1]
                        max_bin_edge = plot_description[0][2]
                        half_bin_size = plot_description[0][3]
                    else:
                        assert combination2[3] not in self.get_io_type_variables()
                        min_bin_edge = plot_description[1][1]
                        max_bin_edge = plot_description[1][2]
                        half_bin_size = plot_description[1][3]
                    result = self.run_ks_test_on_samples(sample1, sample2, min_bin_edge, max_bin_edge, half_bin_size, weights1, weights2)
                    results_list.append(result)
                    number_of_times_sampled += 1
                else:
                    break   # This will occur when we run out of WDs to sample (and the get_random_subset function returns None). We now want to break out of the loop and report the results
            distinguishable_count, total_count = self.report_ks_test_results(combination1, combination2, variable, sample_size, sample_size, results_list, True, p_threshold)
            if distinguishable_count/total_count > distinguishability_probability_threshold:
                break
            else:
                # Some logic to increase the sample size more if we're really far away from making the threshold
                #sample_size /= 0.48*(np.tanh(((distinguishable_count/total_count)*6)-3) + 1.05)
                #sample_size /= (distinguishable_count/total_count)**1.2
                sample_size *= 1.1
                if distinguishable_count/total_count < 0.3:
                    sample_size *= 2
                #if distinguishable_count/total_count < 0.3:
                #    sample_size *= 6
                sample_size = int(sample_size)
        pop1.reset_unsampled_indices()
        pop2.reset_unsampled_indices()
        return sample_size #Note this will always be an overestimate (by up to ~10%), but it's non-deterministic (and boundary can be inherently fuzzy) so we're only in the game of making v rough estimates anyway. I just round down to a roughly round number

    def run_ks_test_against_comparison(self, combination, comparison_name, comparison_values, variable, plot_description, resample_count=1, most_polluted=None, filter_on_elements=None):
        # combinations should be a tuple containing:
        #      0         1         2          3                                                              4
        # (mod_name, obs_name, pop_name, variable_type ('Input', 'Pollution', 'Observed', 'Modelled'), population size)
        # (population_size can be omitted)

        # variable, plot_description should be a key, value pair from plot_descriptions

        #mv_tests affects the functionality massively and should maybe not have a default value? If True, it doesn't actually run a ks test but instead pulls out "observed" Ca, Mg and Fe data and runs multivariate statistical tests on those
        pop_including_unobserved = self.results_dict[combination[0]][combination[1]][combination[2]]
        if filter_on_elements is not None:
            # Take all systems that are have detections of the requested elements
            pop = pop_including_unobserved.get_subset_with_detected_elements(filter_on_elements)
        else:
            # Take all systems that are actually observed:
            pop = pop_including_unobserved.get_observed_subset()
        try:
            pop_size = combination[4]
        except IndexError:
            pop_size = None
        number_of_times_resampled = 0
        results_list = list()
        while number_of_times_resampled < resample_count:
            if most_polluted is not None:
                # Then we will use (as our 1 sample) the most polluted white dwarfs - this roughly mimics Hollands et al 2017
                pop_size = most_polluted[0]
                pol_element = most_polluted[1]
                print('Sampling the ' + str(pop_size) + ' most polluted WDs (by ' + str(pol_element) + ')')
                pop_sample = pop.get_most_polluted_subset(pop_size, pol_element)
                resample_count = 1
            else:
                pop_sample = pop.get_random_subset(pop_size, True)
            if pop_sample is not None:
                variable_dict, weights_dict, variable_below, variable_above = self.extract_plottables(combination[2], pop_sample, {variable: plot_description})

                sample = variable_dict[combination[3]][variable][combination[2]]

                weights = weights_dict[combination[3]][variable][combination[2]]

                if combination[3] in self.get_io_type_variables():
                    min_bin_edge = plot_description[0][1]
                    max_bin_edge = plot_description[0][2]
                    half_bin_size = plot_description[0][3]
                else:
                    min_bin_edge = plot_description[1][1]
                    max_bin_edge = plot_description[1][2]
                    half_bin_size = plot_description[1][3]
                result = self.run_ks_test_on_samples(sample, comparison_values, min_bin_edge, max_bin_edge, half_bin_size, weights, [1]*len(comparison_values))
                results_list.append(result)
                number_of_times_resampled += 1
            else:
                break
        distinguishable_count, total_count = self.report_ks_test_results(combination, comparison_name, variable, pop_size, len(comparison_values), results_list, True)
        return results_list, distinguishable_count, total_count

    def run_ks_test_on_results(self, combination1, combination2, variable, plot_description, resample_count=1):
        # combinations should be a tuple containing:
        #      0         1         2          3                                                              4
        # (mod_name, obs_name, pop_name, variable_type ('Input', 'Pollution', 'Observed', 'Modelled'), population size)

        # variable, plot_description should be a key, value pair from plot_descriptions
        pop1_including_unobserved = self.results_dict[combination1[0]][combination1[1]][combination1[2]]
        pop2_including_unobserved = self.results_dict[combination2[0]][combination2[1]][combination2[2]]

        #Need to pick out the right subset - bear in mind that extract_plottables can return a smaller subset of values than expected if it receives systems that it shouldn't. This is bad (it throws off the KS test)
        if combination1[3] == 'Modelled':
            pop1 = pop1_including_unobserved.get_subset_with_modelled_properties([plot_description[0][0]])
        else:
            raise NotImplementedError('For other variable types, need to pick the appropriate subset. e.g. for Observed, should pick out systems with detections of relevant elements')

        if combination2[3] == 'Modelled':
            pop2 = pop2_including_unobserved.get_subset_with_modelled_properties([plot_description[0][0]])
        else:
            raise NotImplementedError('For other variable types, need to pick the appropriate subset. e.g. for Observed, should pick out systems with detections of relevant elements')

        try:
            pop1_size = combination1[4]
        except IndexError:
            pop1_size = None
        try:
            pop2_size = combination2[4]
        except IndexError:
            pop2_size = None

        number_of_times_resampled = 0
        results_list = list()
        while number_of_times_resampled < resample_count:
            pop1_sample = pop1.get_random_subset(pop1_size, True)
            if combination1[0] == combination2[0] and combination1[1] == combination2[1] and combination1[2] == combination2[2]:
                pop2_sample = pop1_sample  # Then they're actually the same population and we're comparing input and output, so the sample should be the same
            else:
                pop2_sample = pop2.get_random_subset(pop2_size, True)
            if pop1_sample is not None and pop2_sample is not None:
                variable_dict1, weights_dict1, variable_below1, variable_above1 = self.extract_plottables(combination1[2], pop1_sample, {variable: plot_description})
                variable_dict2, weights_dict2, variable_below2, variable_above2 = self.extract_plottables(combination2[2], pop2_sample, {variable: plot_description})

                sample1 = variable_dict1[combination1[3]][variable][combination1[2]]
                sample2 = variable_dict2[combination2[3]][variable][combination2[2]]

                weights1 = weights_dict1[combination1[3]][variable][combination1[2]]
                weights2 = weights_dict2[combination2[3]][variable][combination2[2]]

                if combination1[3] in self.get_io_type_variables():
                    assert combination2[3] in self.get_io_type_variables()
                    min_bin_edge = plot_description[0][1]
                    max_bin_edge = plot_description[0][2]
                    half_bin_size = plot_description[0][3]
                else:
                    assert combination2[3] not in self.get_io_type_variables()
                    min_bin_edge = plot_description[1][1]
                    max_bin_edge = plot_description[1][2]
                    half_bin_size = plot_description[1][3]
                result = self.run_ks_test_on_samples(sample1, sample2, min_bin_edge, max_bin_edge, half_bin_size, weights1, weights2)
                results_list.append(result)
                number_of_times_resampled += 1
            else:
                break   # This will occur when we run out of WDs to sample (and the get_random_subset function returns None). We now want to break out of the loop and report the results
        distinguishable_count, total_count = self.report_ks_test_results(combination1, combination2, variable, pop1_size, pop2_size, results_list)
        return results_list, distinguishable_count, total_count

    def report_ks_test_results(self, combination1, combination2, variable, pop1_size, pop2_size, results_list, append=False, p_threshold=None):
        if p_threshold is None:
            p_threshold = self.p_threshold
        to_write = list()
        if append:
            write_mode = 'a'
        else:
            write_mode = 'w'
        to_write.append('Comparison between: ')
        to_write.append('Modeller = ' + combination1[0] + ', Observer = ' + combination1[1] + ', Pop = ' + combination1[2] + ' (N = ' + str(pop1_size) + '), Variable = ' + combination1[3] + ' ' + variable)
        if isinstance(combination2, str):
            to_write.append('Comparison sample: ' + combination2)
        else:
            to_write.append('Modeller = ' + combination2[0] + ', Observer = ' + combination2[1] + ', Pop = ' + combination2[2] + ' (N = ' + str(pop2_size) + '), Variable = ' + combination2[3] + ' ' + variable)
        distinguishable_count = 0
        total_count = 0
        for result in results_list:
            to_write.append('KS statistic = ' + str(result[0]))
            to_write.append('KS p-value = ' + str(result[1]))
            if result[1] is not None and not np.isnan(result[1]):
                if result[1] > p_threshold:
                    to_write.append('Could not distinguish distributions')
                else:
                    to_write.append('Samples appear to be from different distributions')
                    distinguishable_count += 1
                total_count += 1
            else:
                to_write.append('(skipping)')
            to_write.append('')
        to_write.append('Rate of distinguishing samples: ' + str(distinguishable_count) + '/' + str(total_count))
        to_write.append('')
        print('Writing to ' + self.get_base_pipeline_dir() + 'ks_results.txt')
        with open(self.get_base_pipeline_dir() + 'ks_results.txt', write_mode, newline='', encoding='utf-8') as f:
            for line in to_write:
                f.write(line)
                f.write('\n')
        for line in to_write:
            print(line)
        return distinguishable_count, total_count

    def report_mv_test_results(self, combination1, combination2, pop1_size, pop2_size, results_dict, append=False):
        to_write = list()
        if append:
            write_mode = 'a'
        else:
            write_mode = 'w'
        to_write.append('Comparison between: ')
        to_write.append('Modeller = ' + combination1[0] + ', Observer = ' + combination1[1] + ', Pop = ' + combination1[2] + ' (N = ' + str(pop1_size) + '), ' + combination1[3])
        if isinstance(combination2, str):
            to_write.append('Comparison sample: ' + combination2)
        else:
            to_write.append('Modeller = ' + combination2[0] + ', Observer = ' + combination2[1] + ', Pop = ' + combination2[2] + ' (N = ' + str(pop2_size) + '), ' + combination2[3] + ' ' + variable)
        toret = dict()
        for test_type in mv.get_mv_test_method_dict():
            distinguishable_count = 0
            total_count = 0
            for result in results_dict[test_type]:
                to_write.append(str(test_type) + ' p-value = ' + str(result))
                if result is not None and not np.isnan(result):
                    if result > self.p_threshold:
                        to_write.append('Could not distinguish distributions')
                    else:
                        to_write.append('Samples appear to be from different distributions')
                        distinguishable_count += 1
                    total_count += 1
                else:
                    to_write.append('(skipping)')
                to_write.append('')
            to_write.append('Rate of distinguishing samples: ' + str(distinguishable_count) + '/' + str(total_count))
            to_write.append('')
            toret[test_type] = (distinguishable_count, total_count)
        print('Writing to ' + self.get_base_pipeline_dir() + 'mv_results.txt')
        with open(self.get_base_pipeline_dir() + 'mv_results.txt', write_mode, newline='', encoding='utf-8') as f:
            for line in to_write:
                f.write(line)
                f.write('\n')
        for line in to_write:
            print(line)
        return toret

    def generate_populations(self, population_parameter_dict, obs_name, mod_name):
        toret = dict()
        for pop_name, params_dict in population_parameter_dict.items():
            # 3 cases: Never seen this combo before, seen pop_name before but not one/both of obs and mod, seen exact combo before
            # In the 1st and 3rd cases, load up with no premade pop - in the 1st case it should realise it's seen this before and load the right pop
            # In the 2nd case, need to be careful. We have got pop data saved but must only load the raw population data, not any outputs

            file_to_look_for = self.get_popdump_dir() + 'popdump_' + self.get_pipeline_name(pop_name, obs_name, mod_name) + '.csv'

            premade_pop = None
            if not os.path.exists(file_to_look_for) and pop_name in self.results_dict:
                dummy_obs = self.results_dict[pop_name].keys().next()
                dummy_mod = self.results_dict[pop_name][dummy_obs].keys().next()
                premade_pop = self.results_dict[pop_name][dummy_obs][dummy_mod].extract_raw_population()
            toret[pop_name] = sp.SyntheticPopulation(
                params_dict[sp.PopulationParameter.size],
                params_dict[sp.PopulationParameter.wd_config],
                params_dict[sp.PopulationParameter.pollution_config],
                file_to_look_for,
                premade_pop,
                self.timescale_interpolator
            )
            #toret[pop_name].check_sinking_timescales() # This is an optional diagnostic (time consuming!)
        return toret

    def plot_results(self, plot_descriptions, mwdd_hx=None, most_polluted=None, filter_on_elements=None):
        #mwdd_hx refers to the reference element, i.e. either ci.Element.H or ci.Element.He
        # Warning: if this is supplied, some of the plots will then implicitly assume that all the white dwarfs in the synthetic pops are DAs/DBs respectively
        # I anticipate we probably don't want to do every combo, so comment stuff out!

        #most_polluted indicates whether to only pull out the most heavily polluted WDs for analysis (in which case most_polluted should be an (int, element) tuple indicating how many to pull out and what element to use to determine pollution level. Can be None if using all elements)
        # filter_on_elements is a list of elements (or None). If not None, we only pick out WDs for analysis with detections of the specified elements.
        add_excluded_systems_to_scatter_plot = True
        io_comparison_plots_dict = dict()
        for mod_name, mod_results in self.results_dict.items():
            io_comparison_plots_dict[mod_name] = dict()
            for obs_name, results in mod_results.items():
                io_comparison_plots_dict[mod_name][obs_name] = dict()
                for pop_name, pop in results.items():
                    variable_dict, weights_dict, variable_below_threshold_dict, variable_above_threshold_dict = self.extract_plottables(pop_name, pop, plot_descriptions, most_polluted, filter_on_elements)
                    model_analyser = ma.ModelAnalyser(self.get_pipeline_dir(pop_name, obs_name, mod_name))
                    individual_plots_for_multipanel = self.plot_individual_population_variables(model_analyser, plot_descriptions, variable_dict, weights_dict, variable_below_threshold_dict, variable_above_threshold_dict, mwdd_hx)
                    self.plot_multipanels(model_analyser, individual_plots_for_multipanel, plot_descriptions)
                    io_comparison_plots_dict[mod_name][obs_name][pop_name] = self.plot_individual_population_io_comparison(model_analyser, pop_name, plot_descriptions, variable_dict, weights_dict, variable_below_threshold_dict, variable_above_threshold_dict)
                    self.plot_individual_population_io_scatter(model_analyser, pop, pop_name, plot_descriptions, add_excluded_systems_to_scatter_plot, mwdd_hx)
        self.plot_specific_io_comparisons(io_comparison_plots_dict)

    def perform_knn_comparisons(self):
        hollands_pop = sp.generate_hollands_cool_dz_sample()
        for mod_name, mod_results in self.results_dict.items():
            for obs_name, results in mod_results.items():
                for pop_name, pop in results.items():
                    relevant_subset = pop.get_subset_with_detected_elements([ci.Element.Ca, ci.Element.Mg, ci.Element.Fe])
                    random_subset = relevant_subset.get_random_subset(len(hollands_pop), True)
                    correct, incorrect, p_value = random_subset.knn_comparison_against_other_pop(hollands_pop)
                    print()
                    print(pop_name)
                    print(obs_name)
                    print(correct)
                    print(incorrect)
                    print(p_value)

    def plot_ternary_comparisons(self, most_polluted_tuple=None):
        manager_for_hollands_abundances = ha.load_manager()
        el_list = [ci.Element.Ca, ci.Element.Fe, ci.Element.Mg]
        for mod_name, mod_results in self.results_dict.items():
            for obs_name, results in mod_results.items():
                for pop_name, pop in results.items():
                    graph_fac = gf.GraphFactory(self.get_pipeline_dir(pop_name, obs_name, mod_name))
                    el_abundance_dict = dict()
                    el_abundance_dict[pop_name] = dict()
                    for element in el_list:
                        el_abundance_dict[pop_name][element] = list()
                    if most_polluted_tuple is not None:
                        pop_size = most_polluted_tuple[0]
                        element_to_use = most_polluted_tuple[1]
                        pop_to_use = pop.get_most_polluted_subset(pop_size, element_to_use)
                    else:
                        pop_to_use = pop
                    for wd in pop_to_use:
                        if not wd.observed:
                            continue
                        all_elements_observed = True
                        for element in el_list:
                            if wd.observed_abundances.get(element) is None:
                                all_elements_observed = False
                        if all_elements_observed:
                            for element in el_list:
                                el_abundance_dict[pop_name][element].append(10**wd.observed_abundances[element])
                    #Now add in Hollands abundances
                    el_abundance_dict['Hollands'] = dict()
                    for element in el_list:
                        el_abundance_dict['Hollands'][element] = [10**hav for hav in ha.get_hollands_abundance_values(element, manager_for_hollands_abundances)]
                    if run_mv_tests:
                        synthetic_data_for_mv_test = [
                            el_abundance_dict[pop_name][ci.Element.Ca], # Things to check: does it matter if I just use 2 of these (as long as its normalised), does the scaling matter, does the test type matter
                            el_abundance_dict[pop_name][ci.Element.Mg], # Also, these need to be normalised anyway: we only care about the relative quantities surely? Maybe this is another thing to test
                            el_abundance_dict[pop_name][ci.Element.Fe]
                        ]
                        hollands_data_for_mv_test = [
                            el_abundance_dict['Hollands'][ci.Element.Ca],
                            el_abundance_dict['Hollands'][ci.Element.Mg],
                            el_abundance_dict['Hollands'][ci.Element.Fe]
                        ]
                        synthetic_converted = mv.convert_data_to_r_format(synthetic_data_for_mv_test)
                        hollands_converted = mv.convert_data_to_r_format(hollands_data_for_mv_test)
                        print('Performing Cramér test...')
                        p_value_cramer = mv.cramer_test_p_value(
                            synthetic_converted,
                            hollands_converted
                        )
                        print('Done')
                        #print('Performing Fasano Franceschini test...')
                        #p_value_fasano_franceschini = mv.fasano_franceschini_test_p_value(
                        #    synthetic_converted,
                        #    hollands_converted
                        #)
                        #print('Done')
                        additional_text_dict = {
                            'Cramér': {
                                'x_pos': -0.08,
                                'y_pos': 0.87,
                                'text_string': 'Cramér p-value: ' + str(round(p_value_cramer, 7))
                            }
                            #'FF': {
                            #    'x_pos': -0.08,
                            #    'y_pos': 0.82,
                            #    'text_string': 'FF p-value: ' + str(round(p_value_fasano_franceschini, 7))
                            #}
                        }
                    graph_fac.make_ternary_plot(el_list, el_abundance_dict, 'ternary_plot_' + pop_name + '_' + obs_name, self.ternary_scaling_factors, additional_text_dict)

    def plot_specific_io_comparisons(self, io_comparison_plots_dict):
        # The following is to be activated on an as/when basis:
        activate_fcf_comparison = False
        activate_strategy_comparison = False
        if activate_fcf_comparison:
            graph_fac = gf.GraphFactory(self.get_base_pipeline_dir())
            plots_to_multipanelise = list()
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0']['DATidal']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['series']['test'] = {
                    'type': dp.PlotType.text,
                    'x_pos': 0.5,
                    'y_pos': 8,
                    'text_string': 'Tidal, 0 dex error',
                    'fontsize': 24
                }
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0']['DACollisional']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['series']['test'] = {
                'type': dp.PlotType.text,
                'x_pos': 0.5,
                'y_pos': 5,
                'text_string': 'Collisional, 0 dex error',
                'fontsize': 24
            }
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0p2']['DATidal']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['series']['test'] = {
                'type': dp.PlotType.text,
                'x_pos': 0.5,
                'y_pos': 8,
                'text_string': 'Tidal, 0.2 dex error',
                'fontsize': 24
            }
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0p2']['DACollisional']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['series']['test'] = {
                'type': dp.PlotType.text,
                'x_pos': 0.5,
                'y_pos': 5,
                'text_string': 'Collisional, 0.2 dex error',
                'fontsize': 24
            }
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0p2']['DACollisional']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['x_hide_ticks'] = [0]
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0']['DATidal']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['legend_text_size'] = 26
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0p2']['DATidal']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['legend'] = False
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0']['DACollisional']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['legend'] = False
            io_comparison_plots_dict['StandardModeller']['RealisticObservererr0p2']['DACollisional']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['legend'] = False
            plots_to_multipanelise.append(io_comparison_plots_dict['StandardModeller']['RealisticObservererr0']['DATidal']['fcf']['pdf'])
            plots_to_multipanelise.append(io_comparison_plots_dict['StandardModeller']['RealisticObservererr0p2']['DATidal']['fcf']['pdf'])
            plots_to_multipanelise.append(io_comparison_plots_dict['StandardModeller']['RealisticObservererr0']['DACollisional']['fcf']['pdf'])
            plots_to_multipanelise.append(io_comparison_plots_dict['StandardModeller']['RealisticObservererr0p2']['DACollisional']['fcf']['pdf'])
            graph_fac.multipanelise(plots_to_multipanelise, 2, 2, ['tidalvcollisional_multipanel.pdf'], 15, 18, 0, 0, True, True)
        if activate_strategy_comparison:
            combined_plot = io_comparison_plots_dict['StandardModeller']['ShortObservation']['DAControlPopLarge']['fcf']['pdf']
            combined_plot['hist_plot']['filenames'] = ['strategy_io_comparisons.pdf']
            combined_plot['hist_plot']['subplots']['subplot1']['y_max'] = max(combined_plot['hist_plot']['subplots']['subplot1']['y_max'], io_comparison_plots_dict['StandardModeller']['LongObservation']['DAControlPopSmall']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['y_max'])
            combined_plot['hist_plot']['subplots']['subplot1']['legend_text_size'] = 11
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Large Sample Initial'] = combined_plot['hist_plot']['subplots']['subplot1']['series'].pop('Initial')
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Large Sample Final'] = combined_plot['hist_plot']['subplots']['subplot1']['series'].pop('Final')
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Large Sample Final']['step'] = True
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Small Sample Initial'] = io_comparison_plots_dict['StandardModeller']['LongObservation']['DAControlPopSmall']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['series']['Initial']
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Small Sample Initial']['colours'] = '#DF4443'
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Small Sample Initial']['edge_colours'] = '#DF4443'
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Small Sample Final'] = io_comparison_plots_dict['StandardModeller']['LongObservation']['DAControlPopSmall']['fcf']['pdf']['hist_plot']['subplots']['subplot1']['series']['Final']
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Small Sample Final']['colours'] = '#CD853F'
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Small Sample Final']['edge_colours'] = '#CD853F'
            combined_plot['hist_plot']['subplots']['subplot1']['series']['Small Sample Final']['step'] = True
            plotter = dp.DictPlotter(combined_plot)
            plotter.draw()
            plotter.yield_output(self.get_base_pipeline_dir())

    def extract_plottables(self, pop_name, pop, plot_descriptions, most_polluted=None, filter_on_elements=None):
        variable_list = self.get_variable_type_list()
        variable_dict = dict()
        weights_dict = dict()
        variable_below_threshold_dict = dict()
        variable_above_threshold_dict = dict()
        for vl in variable_list:
            variable_dict[vl] = dict()
            weights_dict[vl] = dict()
            variable_below_threshold_dict[vl] = dict()
            variable_above_threshold_dict[vl] = dict()
            for plot_name, description in plot_descriptions.items():
                variable_dict[vl][plot_name] = dict()
                weights_dict[vl][plot_name] = dict()
                variable_below_threshold_dict[vl][plot_name] = dict()
                variable_above_threshold_dict[vl][plot_name] = dict()
                vfilter = self.get_vfilter(vl, description)
                if vl in ['Input', 'Pollution']:
                    population = pop
                else:
                    if filter_on_elements is None:
                        population = pop.get_observed_subset()
                    else:
                        # filter_on_elements should be a list of elements. We only want the subset where those elements are detected
                        population = pop.get_subset_with_detected_elements(filter_on_elements)
                    if most_polluted is not None:
                        pop_size = most_polluted[0]
                        element_to_use = most_polluted[1]
                        population = population.get_most_polluted_subset(pop_size, element_to_use)
                if vl == 'Input':
                    prefilter = population.input_values(vfilter[0])
                    weights = [1]*len(prefilter)
                elif vl == 'Pollution':
                    prefilter = population.pollution_abundance_values(vfilter[0])
                    weights = [1]*len(prefilter)
                elif vl == 'Observed':
                    if len(vfilter[0]) == 2:
                        # Then we're dealing with a tuple of 2 elements, and we want to return the ratio
                        prefilter = population.observed_abundance_log_ratios(vfilter[0][0], vfilter[0][1])
                    else:
                        prefilter = population.observed_abundances(vfilter[0])
                    weights = [1]*len(prefilter)
                else:
                    raw_modelled_values = population.modelled_values(vfilter[0])
                    weights = list()
                    prefilter = list()
                    for i, rmv in enumerate(raw_modelled_values):
                        # Expect each item to be a list. Need to collapse this list of lists into a single list, but keep track of the weight each should be assigned
                        # (eg If a system has 3 possible output values, they need to each by weighted by 1/3 to avoid distorting the histogram)
                        weight = 1/len(rmv)
                        added_weight = list()
                        for entry in rmv:
                            prefilter.append(entry)
                            weights.append(weight)
                            added_weight.append(weight)
                        assert abs(sum(added_weight) - 1) < 0.0000001
                variable_dict[vl][plot_name][pop_name] = list()
                weights_dict[vl][plot_name][pop_name] = list()
                variable_below_threshold_dict[vl][plot_name][pop_name] = 0
                variable_above_threshold_dict[vl][plot_name][pop_name] = 0
                for i, val in enumerate(prefilter):
                    try:
                        if val is not None and ~np.isnan(val):
                            variable_dict[vl][plot_name][pop_name].append(val)
                            weights_dict[vl][plot_name][pop_name].append(weights[i])
                            if val < vfilter[1]:
                                variable_below_threshold_dict[vl][plot_name][pop_name] += weights[i]
                                print(val)
                                print(vfilter[1])
                                raise ValueError('For now, this should be forbidden! (Just increase the bounds in the plot_description if it\'s an issue)')
                            if val > vfilter[2]:
                                variable_above_threshold_dict[vl][plot_name][pop_name] += weights[i]
                                print(val)
                                print(vfilter[2])
                                raise ValueError('For now, this should be forbidden! (Just increase the bounds in the plot_description if it\'s an issue)')
                    except TypeError:
                        pass
        return variable_dict, weights_dict, variable_below_threshold_dict, variable_above_threshold_dict

    def plot_individual_population_io_comparison(self, model_analyser, pop_name, plot_descriptions, variable_dict, weights_dict, variable_below_threshold_dict, variable_above_threshold_dict):
        variables_to_plot = ['fcf']
        hist_modes = ['pdf', 'cdf']
        paper_version = False
        io_comparison_plots_dict = {vtp: dict() for vtp in variables_to_plot}
        for plot_name in variables_to_plot:
            io_variable_dict = dict()
            io_weights_dict = dict()
            io_text_dict = dict()
            vfilter = self.get_vfilter('Input', plot_descriptions[plot_name])
            initial_text_height = 1
            for pop_name, values in variable_dict['Input'][plot_name].items():
                series_name = 'Initial' if paper_version else pop_name + ' Initial'
                io_variable_dict[series_name] = values
                io_weights_dict[series_name] = weights_dict['Input'][plot_name][pop_name]
                if not paper_version:
                    io_text_dict[pop_name + ' Initial excluded'] = {
                        'text_string': pop_name + ' Initial' + ': ' + str(sum(weights_dict['Input'][plot_name][pop_name])) + ' total, ' + str(variable_below_threshold_dict['Input'][plot_name][pop_name]) + ' below, ' + str(variable_above_threshold_dict['Input'][plot_name][pop_name]) + ' above',
                        'x_pos': vfilter[1] + ((vfilter[2] - vfilter[1])/5),
                        'y_pos': initial_text_height,
                        #'horizontalalignment'
                        #'verticalalignment'
                        #'rotation'
                        'fontsize': 12
                    }
                    initial_text_height -= 0.1
            for pop_name, values in variable_dict['Modelled'][plot_name].items():
                series_name = 'Final' if paper_version else pop_name + ' Final'
                io_variable_dict[series_name] = values
                io_weights_dict[series_name] = weights_dict['Modelled'][plot_name][pop_name]
                if not paper_version:
                    io_text_dict[pop_name + ' Final excluded'] = {
                        'text_string': pop_name + ' Final' + ': ' + str(sum(weights_dict['Modelled'][plot_name][pop_name])) + ' total, ' + str(variable_below_threshold_dict['Modelled'][plot_name][pop_name]) + ' below, ' + str(variable_above_threshold_dict['Modelled'][plot_name][pop_name]) + ' above',
                        'x_pos': vfilter[1] + ((vfilter[2] - vfilter[1])/5),
                        'y_pos': initial_text_height,
                        #'horizontalalignment'
                        #'verticalalignment'
                        #'rotation'
                        'fontsize': 12
                    }
                    initial_text_height -= 0.1
            for mode in hist_modes:
                io_comparison_plot = model_analyser.make_variable_distribution_plot(
                    vfilter[0],
                    io_variable_dict,
                    vfilter[1],
                    vfilter[2],
                    vfilter[3],
                    io_text_dict,
                    plot_name + '_io_dist_' + pop_name + '_' + mode,
                    mode == 'cdf',
                    dict(),
                    io_weights_dict
                )
                io_comparison_plots_dict[plot_name][mode] = io_comparison_plot
        return io_comparison_plots_dict

    def plot_individual_population_io_scatter(self, model_analyser, population, pop_name, plot_descriptions, add_excluded_systems_to_scatter_plot, hx):
        # Possible TODO: should this have a mode for picking out the most polluted systems?
        for plot_name, plot_description in plot_descriptions.items():
            vfilter = self.get_vfilter('Input', plot_description)
            observed_pop = population.get_observed_subset()
            modelled_input_vals, modelled_output_vals, unmodelled_input_vals = observed_pop.modelled_io_pairs(vfilter[0])
            bins, x_bar_centres = model_analyser.generate_bins_and_bar_centres(vfilter[1], vfilter[2], vfilter[3])
            heights, bins2 = np.histogram(
                unmodelled_input_vals,
                bins,
                density=False
            )
            io_scatter_plot = model_analyser.graph_fac.make_io_scatter_plot(modelled_input_vals, modelled_output_vals, heights, x_bar_centres, plot_name + '_io_scatter_' + pop_name + '.pdf', vfilter[0], vfilter[1], vfilter[2], vfilter[3], add_excluded_systems_to_scatter_plot)
        plot_elel_scatter = False
        if plot_elel_scatter:
            el1 = ci.Element.Ca
            el2 = ci.Element.Fe
            el1_values = population.pollution_abundance_values(el1)
            el2_values = population.pollution_abundance_values(el2)
            el1_el2_observed_bools = population.observed_abundance_bools(el1, el2)
            #dummy_observer = so.Observer(so.ObservationType.TeffIndividualElementCutoff)
            #spectral_type = 'DA' if hx == ci.Element.H else 'DB'
            #teff = 5000
            #el1_threshold = dummy_observer.calculate_element_cutoff_threshold(el1, spectral_type, teff)
            #el2_threshold = dummy_observer.calculate_element_cutoff_threshold(el2, spectral_type, teff)
            population_limit = 1000
            if population_limit is not None:
                el1_values = el1_values[0:population_limit]
                el2_values = el2_values[0:population_limit]
                el1_el2_observed_bools = el1_el2_observed_bools[0:population_limit]
            io_scatter_plot = model_analyser.graph_fac.make_elel_scatter_plot(el1_values, el2_values, el1_el2_observed_bools, el1, el2, None, None, str(el1) + '_' + str(el2) + '_scatter_' + pop_name + '.pdf', hx)

    def plot_individual_population_variables(self, model_analyser, plot_descriptions, variable_dict, weights_dict, variable_below_threshold_dict, variable_above_threshold_dict, mwdd_hx=None):
        variable_list = self.get_variable_type_list_for_multipanel()
        vls_to_make_standalone = ['Observed']
        plots_for_multipanel = dict()
        hist_modes = ['pdf', 'cdf']
        for vl in variable_list:
            plots_for_multipanel[vl] = dict()
            for plot_name, description in plot_descriptions.items():
                vfilter = self.get_vfilter(vl, description)
                excluded_text_dict = dict()
                initial_text_height = 1
                for pop_name in variable_dict[vl][plot_name].keys():
                    excluded_text_dict[pop_name + '_excluded'] = {
                        'text_string': pop_name + ': ' + str(sum(weights_dict[vl][plot_name][pop_name])) + ' total, ' + str(variable_below_threshold_dict[vl][plot_name][pop_name]) + ' below, ' + str(variable_above_threshold_dict[vl][plot_name][pop_name]) + ' above',
                        'x_pos': vfilter[1] + ((vfilter[2] - vfilter[1])/5),
                        'y_pos': initial_text_height,
                        #'horizontalalignment'
                        #'verticalalignment'
                        #'rotation'
                        'fontsize': 12
                    }
                    initial_text_height -= 0.1
                plots_for_multipanel[vl][plot_name] = dict()
                comparison_dist_dict = dict()
                if isinstance(vfilter[0], tuple) and len(vfilter[0]) == 2:
                    # Assume that it's a 2-tuple of 2 elements: comment one of these out as appropriate!
                    comparison_dist_dict['Hollands'] = ha.get_hollands_el_el_values(vfilter[0][0], vfilter[0][1], self.manager)
                    #comparison_dist_dict['40pc'] = fp.get_40pc_el_el_values(vfilter[0][0], vfilter[0][1])
                variable_name = vfilter[0]
                if isinstance(vfilter[0], ci.Element) and mwdd_hx is not None:
                    #comment one of these out as appropriate!
                    comparison_dist_dict['MWDD ' + str(vfilter[0]) + '/' + str(mwdd_hx)] = mwdd.get_mwdd_abundances(vfilter[0], mwdd_hx)
                    comparison_dist_dict['Hollands+ 2017 ' + str(vfilter[0]) + '/' + str(mwdd_hx)] = ha.get_hollands_abundance_values(vfilter[0], self.manager)
                    variable_name = (vfilter[0], mwdd_hx)
                for mode in hist_modes:
                    vl_plot = model_analyser.make_variable_distribution_plot(
                        variable_name, #variable_name
                        variable_dict[vl][plot_name], #variable_values_dict,
                        vfilter[1], #x_min
                        vfilter[2], #x_max
                        vfilter[3], #half_bin_size
                        excluded_text_dict, #text_dict
                        plot_name + '_' + vl + '_' + mode, #file_prefix
                        mode == 'cdf', #cumulative (True or False)
                        dict(), #text_size_dict
                        weights_dict[vl][plot_name], #weights_dict
                        comparison_dist_dict
                    )
                    plots_for_multipanel[vl][plot_name][mode] = vl_plot
                    if vl in vls_to_make_standalone:
                        standalone_plot = model_analyser.make_variable_distribution_plot(
                            variable_name,
                            variable_dict[vl][plot_name],
                            vfilter[1],
                            vfilter[2],
                            vfilter[3],
                            None,
                            plot_name + '_' + vl + '_' + mode + '_standalone',
                            mode == 'cdf',
                            {'legend_text_size': 14, 'xlabel_fontsize': 22, 'ylabel_fontsize': 22, 'x_tick_fontsize': 14, 'y_tick_fontsize': 14},
                            weights_dict[vl][plot_name],
                            comparison_dist_dict
                        )
        return plots_for_multipanel

    def plot_multipanels(self, model_analyser, individual_pop_plots_dict, plot_descriptions):
        #for plot_name in plot_descriptions.keys():
        hist_modes = ['pdf', 'cdf']
        for plot_name in ['fcf']:
            for mode in hist_modes:
                plots_to_multipanelise = list()
                for vl in self.get_variable_type_list_for_multipanel():
                    plots_to_multipanelise.append(individual_pop_plots_dict[vl][plot_name][mode])
                model_analyser.graph_fac.multipanelise(plots_to_multipanelise, 2, 2, [plot_name + '_FullPipeline_' + mode + '.pdf'], 15, 18, 0.2, 0.2, False, False)

    def get_variable_type_list(self):
        return ['Input', 'Pollution', 'Observed', 'Modelled']

    def get_variable_type_list_for_multipanel(self):
        # This is in a different order to arrange the panels in a certain way
        return ['Input', 'Pollution', 'Modelled', 'Observed']

    def get_io_type_variables(self):
        return ['Input', 'Modelled']

    def get_vfilter(self, variable_type, plot_description):
        if variable_type in self.get_io_type_variables():
            return plot_description[0]
        else:
            return plot_description[1]

def run_control_pipeline():
    # KEYS IN THESE DICTS SHOULD IDEALLY CONTAIN NO UNDERSCORES OR WHITESPACE
    population_parameter_dict = {
        'DAControlPop': {
            sp.PopulationParameter.size: 2000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DAControl'
        }
    }
    observer_dict = {
        'IdealObservererr0': [
            so.ObservationType.NoCut,
            0
        ]
    }
    modeller_dict = {
        'StandardModeller': [
            sm.ModellerType.AnalyticApproximation,
            [False, 'synthetic_grid_dummy.csv']
        ]
    }
    pipeline = Pipeline('control', population_parameter_dict, observer_dict, modeller_dict)
    pipeline.run_all_combinations()

    plot_descriptions = {
         # input description, output description (each is variable, min bin value, max bin value, half bin size)
        'fcf': ((mp.ModelParameter.fragment_core_frac, 0, 1, 0.025), (ci.Element.Fe, -20, -2, 0.1)),
        'distance': ((mp.ModelParameter.formation_distance, -2, 1, 0.0125), (ci.Element.Na, -20, -4, 0.1)),
        'time': ((mp.ModelParameter.t_sinceaccretion, 0, 100000, 1), (ci.Element.Ca, -20, -4, 0.1)),
        'metallicity': ((mp.ModelParameter.metallicity, 0, 958, 24), (ci.Element.Ca, -20, -4, 0.1))
    }

    pipeline.plot_results(plot_descriptions)

def run_control_pipeline_with_realistic_observer():
    # KEYS IN THESE DICTS SHOULD IDEALLY CONTAIN NO UNDERSCORES OR WHITESPACE
    population_parameter_dict = {
        'DAControlPop': {
            sp.PopulationParameter.size: 100000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DAControl'
        }
    }
    observer_dict = {
        'RealisticObservererr0': [
            so.ObservationType.TeffIndividualElementCutoff,
            0
        ],
        'RealisticObservererr0p025': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.025
        ],
        'RealisticObservererr0p05': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.05
        ],
        'RealisticObservererr0p075': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.075
        ],
        'RealisticObservererr0p1': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.1
        ],
        'RealisticObservererr0p125': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.125
        ],
        'RealisticObservererr0p15': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.15
        ],
        'RealisticObservererr0p175': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.175
        ],
        'RealisticObservererr0p2': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.2
        ],
        'RealisticObservererr0p225': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.225
        ],
        'RealisticObservererr0p25': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.25
        ],
        'RealisticObservererr0p275': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.275
        ],
        'RealisticObservererr0p3': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.3
        ],
        'RealisticObservererr0p325': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.325
        ],
        'RealisticObservererr0p35': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.35
        ],
        'RealisticObservererr0p375': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.375
        ],
        'RealisticObservererr0p4': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.4
        ]
    }
    modeller_dict = {
        'StandardModeller': [
            sm.ModellerType.AnalyticApproximation,
            [False, 'synthetic_grid_dummy.csv']
        ]
    }
    pipeline = Pipeline('controlrealistic', population_parameter_dict, observer_dict, modeller_dict)
    pipeline.run_all_combinations()

    plot_descriptions = {
         # input description, output description (each is variable, min bin value, max bin value, half bin size)
        'fcf': ((mp.ModelParameter.fragment_core_frac, 0, 1, 0.025), (ci.Element.Fe, -20, -2, 0.1)),
        'distance': ((mp.ModelParameter.formation_distance, -2, 1, 0.0125), (ci.Element.Na, -20, -4, 0.1)),
        'time': ((mp.ModelParameter.t_sinceaccretion, 0, 100000, 1), (ci.Element.Ca, -20, -4, 0.1)),
        'metallicity': ((mp.ModelParameter.metallicity, 0, 958, 24), (ci.Element.Ca, -20, -4, 0.1))
    }

    pipeline.plot_results(plot_descriptions, ci.Element.H)

def run_control_pipeline_with_variable_observer_offset():
    # Currently unused, but keeping it here to demonstrate this functionality
    population_parameter_dict = {
        'DAControlPop': {
            sp.PopulationParameter.size: 100000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DAControl'
        }
    }
    observer_dict = {
        'DefaultObserveroff2': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            2
        ],
        'DefaultObserveroff1p5': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            1.5
        ],
        'DefaultObserveroff1': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            1
        ],
        'DefaultObserveroff0p5': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            0.5
        ],
        'DefaultObserveroff0': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            0
        ],
        'DefaultObserveroffm0p5': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            -0.5
        ],
        'DefaultObserveroffm1': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            -1
        ],
        'DefaultObserveroffm1p5': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            -1.5
        ],
        'DefaultObserveroffm2': [
            so.ObservationType.TeffIndividualElementCutoff,
            0,
            -2
        ]
    }
    modeller_dict = {
        'StandardModeller': [
            sm.ModellerType.AnalyticApproximation,
            [False, 'synthetic_grid_dummy.csv']
        ]
    }
    pipeline = Pipeline('controlvariable', population_parameter_dict, observer_dict, modeller_dict)
    pipeline.run_all_combinations()

    plot_descriptions = {
         # input description, output description (each is variable, min bin value, max bin value, half bin size)
        'fcf': ((mp.ModelParameter.fragment_core_frac, 0, 1, 0.025), (ci.Element.Fe, -20, -2, 0.1)),
        'distance': ((mp.ModelParameter.formation_distance, -2, 1, 0.0125), (ci.Element.Na, -20, -4, 0.1)),
        'time': ((mp.ModelParameter.t_sinceaccretion, 0, 100000, 1), (ci.Element.Ca, -20, -4, 0.1)),
        'metallicity': ((mp.ModelParameter.metallicity, 0, 958, 24), (ci.Element.Ca, -20, -4, 0.1))
    }

    pipeline.plot_results(plot_descriptions, ci.Element.H)

def run_fcf_comparison_pipeline():
    # KEYS IN THESE DICTS SHOULD IDEALLY CONTAIN NO UNDERSCORES OR WHITESPACE
    population_parameter_dict = {
        'DATidal': {
            sp.PopulationParameter.size: 100000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DATidal'
        },
        'DACollisional': {
            sp.PopulationParameter.size: 100000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DACollisional'
        },
        'DADeltaFcfPop': {
            sp.PopulationParameter.size: 100000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DADeltaFcf'
        }
    }
    observer_dict = {
        'RealisticObservererr0': [
            so.ObservationType.TeffIndividualElementCutoff,
            0
        ],
        'RealisticObservererr0p025': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.025
        ],
        'RealisticObservererr0p05': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.05
        ],
        'RealisticObservererr0p075': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.075
        ],
        'RealisticObservererr0p1': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.1
        ],
        'RealisticObservererr0p125': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.125
        ],
        'RealisticObservererr0p15': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.15
        ],
        'RealisticObservererr0p175': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.175
        ],
        'RealisticObservererr0p2': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.2
        ],
        'RealisticObservererr0p225': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.225
        ],
        'RealisticObservererr0p25': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.25
        ],
        'RealisticObservererr0p275': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.275
        ],
        'RealisticObservererr0p3': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.3
        ],
        'RealisticObservererr0p325': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.325
        ],
        'RealisticObservererr0p35': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.35
        ],
        'RealisticObservererr0p375': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.375
        ],
        'RealisticObservererr0p4': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.4
        ]
    }
    modeller_dict = {
        'StandardModeller': [
            sm.ModellerType.AnalyticApproximation,
            [False, 'synthetic_grid_dummy.csv']
        ]
    }
    pipeline = Pipeline('fcfcomparisons', population_parameter_dict, observer_dict, modeller_dict)
    pipeline.run_all_combinations()

    plot_descriptions = {
         # input description, output description (each is variable, min bin value, max bin value, half bin size)
        'fcf': ((mp.ModelParameter.fragment_core_frac, 0, 1, 0.025), (ci.Element.Fe, -20, -2, 0.1))
    }

    make_p_threshold_table = False
    if make_p_threshold_table:
        p_thresholds = [0.005, 0.01, 0.05, 0.1]
        sample_sizes_cd = list()
        sample_sizes_co = list()
        sample_sizes_od = list()
        for p_threshold in p_thresholds:
            sample_size_cd = pipeline.find_sample_size_for_given_distinguishability_probability(
                ['StandardModeller', 'RealisticObservererr0p2', 'DACollisional', 'Modelled'],
                ['StandardModeller', 'RealisticObservererr0p2', 'DADeltaFcfPop', 'Modelled'],
                'fcf',
                plot_descriptions['fcf'],
                p_threshold
            )
            sample_size_od = pipeline.find_sample_size_for_given_distinguishability_probability(
                ['StandardModeller', 'RealisticObservererr0p2', 'DATidal', 'Modelled'],
                ['StandardModeller', 'RealisticObservererr0p2', 'DADeltaFcfPop', 'Modelled'],
                'fcf',
                plot_descriptions['fcf'],
                p_threshold
            )
            sample_size_co = pipeline.find_sample_size_for_given_distinguishability_probability(
                ['StandardModeller', 'RealisticObservererr0p2', 'DACollisional', 'Modelled'],
                ['StandardModeller', 'RealisticObservererr0p2', 'DATidal', 'Modelled'],
                'fcf',
                plot_descriptions['fcf'],
                p_threshold
            )
            sample_sizes_cd.append(sample_size_cd)
            sample_sizes_co.append(sample_size_co)
            sample_sizes_od.append(sample_size_od)

    plot_results = False
    if plot_results:
        pipeline.plot_results(plot_descriptions, ci.Element.H)

    pipeline.create_dprob_ks_v_error_plot(
        'DATidal',
        'DACollisional',
        [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4],
        [0],
        [20, 100, 200, 500, 1000, 2000],
        'StandardModeller',
        'RealisticObserver',
        'fcf',
        plot_descriptions,
        True,
        5000,
        None
    )

    pipeline.create_dprob_ks_v_error_plot(
        'DADeltaFcfPop',
        'DATidal',
        [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4],
        [0],
        [20, 100, 200, 500, 1000, 2000],
        'StandardModeller',
        'RealisticObserver',
        'fcf',
        plot_descriptions,
        True,
        5000,
        None
    )

    pipeline.create_dprob_ks_v_error_plot(
        'DADeltaFcfPop',
        'DACollisional',
        [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4],
        [0],
        [20, 100, 200, 500, 1000, 2000],
        'StandardModeller',
        'RealisticObserver',
        'fcf',
        plot_descriptions,
        True,
        5000,
        None
    )

def run_delta_fcf_variable_error_pipeline():
    # For the version where Cr is included as a fourth element, uncomment/swap over the lines commented and marked with 'CR'
    # Although technically there's no real reason you couldn't just do it all together in one go
    population_parameter_dict = {
        'SyntheticHollandsDeltaPop': {
            sp.PopulationParameter.size: 50000, # Should be enough to give a decent number of resamples during mv testing
            sp.PopulationParameter.wd_config: 'HollandsDBs',
            sp.PopulationParameter.pollution_config: 'DBDeltaFcf'
        },
        'SyntheticHollandsTidal': {
            sp.PopulationParameter.size: 50000,
            sp.PopulationParameter.wd_config: 'HollandsDBs',
            sp.PopulationParameter.pollution_config: 'DBTidal'
        },
        'SyntheticHollandsCollisional': {
            sp.PopulationParameter.size: 50000,
            sp.PopulationParameter.wd_config: 'HollandsDBs',
            sp.PopulationParameter.pollution_config: 'DBCollisional'
        }
    }
    observer_dict = dict()
    list_of_error_tuples = [ # Doing this explicitly rather than by iteration to avoid floating point errors
        (0, '0'),
        (0.025, '0p025'),
        (0.05, '0p05'),
        (0.075, '0p075'),
        (0.1, '0p1'),
        (0.125, '0p125'),
        (0.15, '0p15'),
        (0.175, '0p175'),
        (0.2, '0p2'),
        (0.225, '0p225'),
        (0.25, '0p25'),
        (0.275, '0p275'),
        (0.3, '0p3'),
        (0.325, '0p325'),
        (0.35, '0p35'),
        (0.375, '0p375'),
        (0.4, '0p4')
    ]
    for error_tuple in list_of_error_tuples:
        observer_dict['HollandsObservererr' + error_tuple[1]] = [
            so.ObservationType.TeffIndividualElementCutoff,
            error_tuple[0],
            'Hollands'
        ]
        # CR
        #observer_dict['HollandsCrObservererr' + error_tuple[1]] = [
        #    so.ObservationType.TeffIndividualElementCutoff,
        #    error_tuple[0],
        #    'HollandsCr'
        #]

    modeller_dict = {
        'NullModeller': [
            sm.ModellerType.Null,
            None
        ]
    }
    pipeline = Pipeline('dzcomparison', population_parameter_dict, observer_dict, modeller_dict)
    pipeline.run_all_combinations()

    sample_size = 202
    sample_size_cr = 22

    element_for_most_polluted = None
    #element_for_most_polluted = ci.Element.Ca
    most_polluted_tuple = None
    #most_polluted_tuple = (sample_size, element_for_most_polluted)

    filter_on_elements = None
    #filter_on_elements = [ci.Element.Ca, ci.Element.Fe, ci.Element.Mg]

    #pipeline.perform_knn_comparisons() # This didn't work

    pipeline.plot_ternary_comparisons(most_polluted_tuple)

    pops_to_iterate = ['SyntheticHollandsDeltaPop', 'SyntheticHollandsTidal', 'SyntheticHollandsCollisional']
    comparison_tuples_to_iterate = [] # Not comparing 2-element tuples any more...
    mv_comparisons_to_iterate = ['Hollands_Ca_Mg_Fe'] #...instead compare 3 elements using a multivariate test
    #CR
    #mv_comparisons_to_iterate = ['Hollands_Ca_Mg_Fe_Cr'] # Now try comparing 4!

    for pop_to_iterate in pops_to_iterate:
        for comparison in mv_comparisons_to_iterate:
            pipeline.create_mv_plots(
                    pop_to_iterate,
                    comparison,
                    [et[0] for et in list_of_error_tuples],
                    [0],
                    [sample_size],
                    'NullModeller',
                    'HollandsObserver',
                    True,
                    248, #Theoretical maximum (= 50000/202)
                    None,
                    most_polluted_tuple
                )
            #CR
            #pipeline.create_mv_plots(
            #        pop_to_iterate,
            #        comparison,
            #        [et[0] for et in list_of_error_tuples],
            #        [0],
            #        [sample_size_cr],
            #        'NullModeller',
            #        'HollandsCrObserver',
            #        True,
            #        248, #Theoretical maximum (= 50000/202)
            #        None,
            #        most_polluted_tuple
            #    )
        for comparison_tuple_to_iterate in comparison_tuples_to_iterate:
                pipeline.create_dprob_ks_v_error_plot(
                    pop_to_iterate,
                    comparison_tuple_to_iterate[0],
                    [et[0] for et in list_of_error_tuples],
                    [0],
                    [sample_size],
                    'NullModeller',
                    'HollandsObserver',
                    'fcf',
                    comparison_tuple_to_iterate[1],
                    True,
                    248,
                    None,
                    most_polluted_tuple,
                    filter_on_elements
                )

def run_strategy_comparison_pipelines():
# KEYS IN THESE DICTS SHOULD IDEALLY CONTAIN NO UNDERSCORES OR WHITESPACE
    population_parameter_dict = {
        'DACollisional': {
            sp.PopulationParameter.size: 100000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DACollisional'
        },
        'DADeltaFcfPop': {
            sp.PopulationParameter.size: 100000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DADeltaFcf'
        }
    }
    observer_dict = {
        'RealisticObservererr0': [
            so.ObservationType.TeffIndividualElementCutoff,
            0
        ],
        'RealisticObservererr0p025': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.025
        ],
        'RealisticObservererr0p05': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.05
        ],
        'RealisticObservererr0p075': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.075
        ],
        'RealisticObservererr0p1': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.1
        ],
        'RealisticObservererr0p125': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.125
        ],
        'RealisticObservererr0p15': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.15
        ],
        'RealisticObservererr0p175': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.175
        ],
        'RealisticObservererr0p2': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.2
        ],
        'RealisticObservererr0p225': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.225
        ],
        'RealisticObservererr0p25': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.25
        ],
        'RealisticObservererr0p275': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.275
        ],
        'RealisticObservererr0p3': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.3
        ],
        'RealisticObservererr0p325': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.325
        ],
        'RealisticObservererr0p35': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.35
        ],
        'RealisticObservererr0p375': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.375
        ],
        'RealisticObservererr0p4': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.4
        ]
    }
    modeller_dict = {
        'StandardModeller': [
            sm.ModellerType.AnalyticApproximation,
            [False, 'synthetic_grid_dummy.csv']
        ]
    }

    pipeline = Pipeline('strategy_comp', population_parameter_dict, observer_dict, modeller_dict)
    pipeline.run_all_combinations()

    plot_descriptions = {
         # input description, output description (each is variable, min bin value, max bin value, half bin size)
        'fcf': ((mp.ModelParameter.fragment_core_frac, 0, 1, 0.025), (ci.Element.Fe, -20, -2, 0.1))
    }

    pipeline.create_dprob_ks_v_N_plot(
        'DADeltaFcfPop',
        'DACollisional',
        [0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4],
        [0],
        [5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 400, 500, 1000, 2000],
        'StandardModeller',
        'RealisticObserver',
        'fcf',
        plot_descriptions,
        True,
        500 #5000
    )

def run_profiling_pipeline():
    # This is a test pipeline for checking performance!
    # KEYS IN THESE DICTS SHOULD IDEALLY CONTAIN NO UNDERSCORES OR WHITESPACE
    population_parameter_dict = {
        'DAProfile': {
            sp.PopulationParameter.size: 1000,
            sp.PopulationParameter.wd_config: 'RealisticDAs',
            sp.PopulationParameter.pollution_config: 'DAControl'
        }
    }
    observer_dict = {
        'ProfileObservererr0': [
            so.ObservationType.TeffIndividualElementCutoff,
            0
        ],
        'ProfileObservererr0p05': [
            so.ObservationType.TeffIndividualElementCutoff,
            0.05
        ]
    }
    modeller_dict = {
        'StandardModeller': [
            sm.ModellerType.AnalyticApproximation,
            [False, 'synthetic_grid_dummy.csv']
        ]
    }
    pipeline = Pipeline('profiling', population_parameter_dict, observer_dict, modeller_dict)
    pipeline.run_all_combinations()

    plot_descriptions = {
         # input description, output description (each is variable, min bin value, max bin value, half bin size)
        'fcf': ((mp.ModelParameter.fragment_core_frac, 0, 1, 0.025), (ci.Element.Fe, -20, -2, 0.1)),
        'distance': ((mp.ModelParameter.formation_distance, -2, 1, 0.0125), (ci.Element.Na, -20, -4, 0.1)),
        'time': ((mp.ModelParameter.t_sinceaccretion, 0, 100000, 1), (ci.Element.Ca, -20, -4, 0.1)),
        'metallicity': ((mp.ModelParameter.metallicity, 0, 958, 24), (ci.Element.Ca, -20, -4, 0.1))
    }
    pipeline.plot_results(plot_descriptions, ci.Element.H)


def main():
    run_profiling_pipeline()
    run_control_pipeline()
    run_control_pipeline_with_realistic_observer()
    #run_control_pipeline_with_variable_observer_offset() # Unused
    run_fcf_comparison_pipeline()
    run_delta_fcf_variable_error_pipeline()
    run_strategy_comparison_pipelines()

if __name__ == '__main__':
    main()
