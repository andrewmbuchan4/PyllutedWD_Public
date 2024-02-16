#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats as st

import model_analyser as ma
import pwd_utils as pu
import synthetic_pipeline as spi

np.random.seed(809)

def trial_ks_tests():
    pipeline = spi.Pipeline('TestPipeline', dict(), dict(), dict())
    mu1 = 7
    sigma1 = 1
    mu2 = 6.9
    sigma2 = 1.1
    min_bin_edge = 2
    max_bin_edge = 12


    configs = {
        #'small': {
        #    'N': 10,
        #    'half_bin_size': 0.5
        #},
        'med': {
            'N': 100,
            'half_bin_size': 0.25
        },
        'large': {
            'N': 1000,
            'half_bin_size': 0.2
        },
        'vlarge': {
            'N': 3000,
            'half_bin_size': 0.15
        },
        'vvlarge': {
            'N': 10001, # 10001 is just over the asymp activation threshold. 10000 is too low
            'half_bin_size': 0.1
        }
    }

    model_analyser = ma.ModelAnalyser(pu.get_path_to_default_graphs())
    list_of_plot_dicts = list()

    for config_name, config in configs.items():
        data1 = np.random.normal(mu1, sigma1, config['N'])
        data2 = np.random.normal(mu2, sigma2, config['N'])
        ks_statistic_original, p_value_original = pipeline.run_ks_test_on_samples(data1, data2, min_bin_edge, max_bin_edge, config['half_bin_size'])

        new_result = st.ks_2samp(data1, data2, "two-sided")

        equiv_result = st.ks_2samp(data1, data2, "two-sided", "asymp") # Always use the approx. version (usually deployed only for N > 10000) this is supposedly equivalent to the original one but appears to be different - maybe because of intermediate step I perform of converting to a cdf

        plot1 = model_analyser.make_variable_distribution_plot(
            'X',
            {
                'Gaussian 1 ($\mu_X$ = ' + str(mu1) + ', $\sigma_X$ = ' + str(sigma1) + ')': data1,
                'Gaussian 2 ($\mu_X$ = ' + str(mu2) + ', $\sigma_X$ = ' + str(sigma2) + ')': data2
            },
            min_bin_edge,
            max_bin_edge,
            config['half_bin_size'],
            {
                'Ntext': {
                    'x_pos': 8.5,
                    'y_pos': 0.7,
                    'text_string': 'N = ' + str(config['N'])
                },
                'ptext': {
                    'x_pos': 8.5,
                    'y_pos': 0.6,
                    'text_string': 'p-value 1 = ' + str(p_value_original)
                },
                'ptext2': {
                    'x_pos': 8.5,
                    'y_pos': 0.55,
                    'text_string': 'p-value 2 = ' + str(new_result[1])
                },
                'ptext3': {
                    'x_pos': 8.5,
                    'y_pos': 0.5,
                    'text_string': 'p-value 3 = ' + str(equiv_result[1])
                }
            },
            config_name
        )

        list_of_plot_dicts.append(plot1)

    model_analyser.graph_fac.multipanelise(
        list_of_plot_dicts,
        2,
        2,
        ['ks_trials.pdf', 'ks_trials.png'],
        15,
        20,
        None,
        None,
        True,
        True
    )

def main():
    trial_ks_tests()

if __name__ == '__main__':
    main()
