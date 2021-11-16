#!/usr/bin/env python
# -*- coding: utf-8 -*-

import kd as kd
import constants as constants

experiments = {
    'Corgne': {
        'PL-169': {
            'P': 3.6,
            'T': 2123,
            'fO2': -0.8,
            'nbot': 2.1
        },
        'PL-185': {
            'P': 3.6,
            'T': 2273,
            'fO2': -0.9,
            'nbot': 2.1
        },
        'PR-369': {
            'P': 3.6,
            'T': 2473,
            'fO2': -1.1,
            'nbot': 2.1
        },
        'PR-373': {
            'P': 7.7,
            'T': 2273,
            'fO2': -0.6,
            'nbot': 2.2
        },
        'PR-376': {
            'P': 3.6,
            'T': 2273,
            'fO2': -3.2,
            'nbot': 1.8
        },
        'PR-383': {
            'P': 3.6,
            'T': 2273,
            'fO2': -3.1,
            'nbot': 1.7
        },
        'PR-375': {
            'P': 3.6,
            'T': 2273,
            'fO2': -3.0,
            'nbot': 1.7
        },
        'PR-368': {
            'P': 3.6,
            'T': 2473,
            'fO2': -3.1,
            'nbot': 1.8
        },
        'PR-365': {
            'P': 7.7,
            'T': 2273,
            'fO2': -3.2,
            'nbot': 1.7
        }
    },
    'Fischer': {
        'DAC1': {
            'P': 31,
            'T': 4250,
            'fO2': -1.0
        },
        'DAC2': {
            'P': 39,
            'T': 3600,
            'fO2': -1.2
        },
        'DAC3': {
            'P': 56,
            'T': 4360,
            'fO2': -1.1
        },
        'DAC4': {
            'P': 57,
            'T': 4440,
            'fO2': -0.9
        },
        'DAC5': {
            'P': 100,
            'T': 5700,
            'fO2': -1.1
        },
        'MAP1': {
            'P': 25,
            'T': 2740,
            'fO2': -2.2
        },
        'MAP2': {
            'P': 25,
            'T': 2740,
            'fO2': -2.7
        },
        'MAP3': {
            'P': 25,
            'T': 2850,
            'fO2': -2.4
        },
        'MAP4': {
            'P': 25,
            'T': 2850,
            'fO2': -2.8
        },
        'MAP5': {
            'P': 25,
            'T': 2860,
            'fO2': -2.3
        },
        'MAP6': {
            'P': 25,
            'T': 2950,
            'fO2': -2.1
        },
        'MAP7': {
            'P': 25,
            'T': 2950,
            'fO2': -2.5
        },
        'MAP8': {
            'P': 25,
            'T': 2950,
            'fO2': -1.9
        },
        'MAP9': {
            'P': 25,
            'T': 2950,
            'fO2': -3.1
        },
        'MAP10': {
            'P': 25,
            'T': 2950,
            'fO2': -2.6
        },
        'MAP11': {
            'P': 15,
            'T': 2807,
            'fO2': -3.7
        },
        'MAP12': {
            'P': 15,
            'T': 2807,
            'fO2': -2.5
        },
        'MAP13': {
            'P': 15,
            'T': 2802,
            'fO2': -2.0
        },
        'MAP14': {
            'P': 15,
            'T': 2720,
            'fO2': -2.8
        },
        'MAP15': {
            'P': 20.5,
            'T': 2727,
            'fO2': -2.8
        }
    },
    'Test cases': {
        'Earth': {
            'P': 54,
            'T': constants.Tpdliq(54),
            'fO2': -2
        }
    }
}

results = {
    'Corgne': {
        'PL-169': {
            'fe': 9.6,
            'mn': 0.013,
            'ni': 750,
            'cr': 0.51,
            'ga': 4.0,
            'si': 0.0014,
            'nb': 0.03,
            'ta': 0.0092,
            'ti': None,
            'cu': 20,
            'zn': 1.14
        },
        'PL-185': {
            'fe': 9.9,
            'mn': 0.03,
            'ni': 583,
            'cr': 0.86,
            'ga': 3.6,
            'si': 0.0014,
            'nb': 0.24,
            'ta': 0.02,
            'ti': None,
            'cu': 19,
            'zn': 0.39
        },
        'PR-369': {
            'fe': 12.5,
            'mn': 0.069,
            'ni': 301,
            'cr': 1.92,
            'ga': 4.7,
            'si': 0.0017,
            'nb': 0.72,
            'ta': 0.03,
            'ti': None,
            'cu': 11,
            'zn': 0.36
        },
        'PR-373': {
            'fe': 8.8,
            'mn': 0.03,
            'ni': 229,
            'cr': 0.85,
            'ga': 2.0,
            'si': 0.001,
            'nb': 0.11,
            'ta': 0.017,
            'ti': None,
            'cu': 13,
            'zn': 0.76
        },
        'PR-376': {
            'fe': 150,
            'mn': 0.84,
            'ni': 5051,
            'cr': 17.0,
            'ga': 64,
            'si': 0.33,
            'nb': 211,
            'ta': 10.9,
            'ti': 0.19,
            'cu': 86,
            'zn': 11.9
        },
        'PR-383': {
            'fe': 144,
            'mn': 0.86,
            'ni': 6109,
            'cr': 18.0,
            'ga': 51,
            'si': 0.34,
            'nb': 197,
            'ta': 10.8,
            'ti': 0.2,
            'cu': 92,
            'zn': 10.2
        },
        'PR-375': {
            'fe': 120,
            'mn': 1.1,
            'ni': 1861,
            'cr': 16.4,
            'ga': 40,
            'si': 0.36,
            'nb': 90,
            'ta': 5.9,
            'ti': 0.17,
            'cu': 49,
            'zn': 7.8
        },
        'PR-368': {
            'fe': 139,
            'mn': 1.12,
            'ni': 2619,
            'cr': 17.5,
            'ga': 45,
            'si': 0.35,
            'nb': 207,
            'ta': 11.1,
            'ti': 0.16,
            'cu': 50,
            'zn': 4.7
        },
        'PR-365': {
            'fe': 148,
            'mn': 1.05,
            'ni': 4665,
            'cr': 17.0,
            'ga': 72,
            'si': 0.33,
            'nb': 128,
            'ta': 9.9,
            'ti': 0.11,
            'cu': 92,
            'zn': 6.2
        }
    },
    'Fischer': {
        'DAC1': {
            'ni': 12,
            'co': 6.4,
            'v': 0.43,
            'cr': None,
            'si': 0.129,
            'fe': 3.08
        },
        'DAC2': {
            'ni': 13,
            'co': 10,
            'v': None,
            'cr': None,
            'si': 0.064,
            'fe': 3.8
        },
        'DAC3': {
            'ni': 7,
            'co': 6,
            'v': None,
            'cr': 2,
            'si': 0.098,
            'fe': 3.6
        },
        'DAC4': {
            'ni': 5.1,
            'co': 4.3,
            'v': 0.5,
            'cr': 1.7,
            'si': 0.053,
            'fe': 2.8
        },
        'DAC5': {
            'ni': 4.9,
            'co': 5.2,
            'v': 0.86,
            'cr': None,
            'si': 0.352,
            'fe': 3.6
        },
        'MAP1': {
            'ni': 150,
            'co': 100,
            'v': 0.216,
            'cr': None,
            'si': 0.08,
            'fe': 13
        },
        'MAP2': {
            'ni': 190,
            'co': 120,
            'v': None,
            'cr': None,
            'si': None,
            'fe': 23
        },
        'MAP3': {
            'ni': 140,
            'co': 90,
            'v': 0.3,
            'cr': None,
            'si': 0.101,
            'fe': 16.4
        },
        'MAP4': {
            'ni': 180,
            'co': 120,
            'v': 0.42,
            'cr': 1.6,
            'si': 0.14,
            'fe': 24
        },
        'MAP5': {
            'ni': 170,
            'co': 100,
            'v': 0.21,
            'cr': None,
            'si': 0.06,
            'fe': 14
        },
        'MAP6': {
            'ni': 99,
            'co': 64,
            'v': None,
            'cr': None,
            'si': 0.08,
            'fe': 11.6
        },
        'MAP7': {
            'ni': 130,
            'co': 110,
            'v': 0.34,
            'cr': 1.4,
            'si': 0.099,
            'fe': 18.5
        },
        'MAP8': {
            'ni': 71,
            'co': 65,
            'v': 0.16,
            'cr': 0.7,
            'si': 0.02,
            'fe': 8.6
        },
        'MAP9': {
            'ni': 210,
            'co': 160,
            'v': 0.83,
            'cr': 2.98,
            'si': 0.307,
            'fe': 37
        },
        'MAP10': {
            'ni': 140,
            'co': None,
            'v': 0.36,
            'cr': 1.48,
            'si': 0.109,
            'fe': 19.2
        },
        'MAP11': {
            'ni': 370,
            'co': 220,
            'v': 1.9,
            'cr': 5.1,
            'si': 0.348,
            'fe': 67
        },
        'MAP12': {
            'ni': 130,
            'co': 90,
            'v': 0.64,
            'cr': 2.0,
            'si': 0.058,
            'fe': 17.7
        },
        'MAP13': {
            'ni': 110,
            'co': 60,
            'v': 0.27,
            'cr': 1.0,
            'si': 0.009,
            'fe': 9.9
        },
        'MAP14': {
            'ni': 240,
            'co': 150,
            'v': 0.53,
            'cr': 1.8,
            'si': 0.068,
            'fe': 25
        },
        'MAP15': {
            'ni': 170,
            'co': 120,
            'v': 0.48,
            'cr': 1.7,
            'si': 0.07,
            'fe': 23.7
        }
    },
    'Test cases': {
        'Earth': { #Values from Rudge 2010 supplementary table 5
            'w': 10**1.513,
            'ni': 10**1.418,
            'p': 10**1.398,
            'co': 10**1.381,
            'pb': 10**1.159,
            'fe': 10**1.136,
            'cu': 10**0.801,
            'v': 10**0.262,
            'cr': 10**0.195,
            'mn': 10**-0.155,
            'nb': 10**-0.276,
            'ta': 10**-0.611,
            'si': 10**-0.728,
            'zn': 10**-0.824,
            'ga': 10**-1.0,
            'ti': 0
        }
    }
}

def main():
    ele_set = ['hf','u','ta','pb','nb','si','mn','zn','ga','v','cr','cu','ti','fe','w','p','co','ni','o']
    D = kd.D(
        #'data/part_param.dat',
        'data/part_param_fischer_update.dat',
        #'data/part_param_fischer_epsilon_update.dat',
        #'data/int_param.dat',
        'data/int_param_fischer_blanchard_update.dat',
        'data/composition.dat',
        True
    )
    for ex_name, ex_runs in experiments.items():
        print('Running experiments from ' + ex_name)
        for run_name, run_params in ex_runs.items():
            print('Running experiment ' + run_name)
            P = run_params['P']
            T = run_params['T']
            fO2 = run_params['fO2']
            nbot = run_params.get('nbot', constants.nbot)
            print('P = ' + str(P) + ', T = ' + str(T) + ', fO2 = ' + str(fO2) + ', nbot = ' + str(nbot))
            Ds, Ds_sd = D.mkd(P, T, fO2, constants.nbot, constants.gammaFe_sil, ele_set)
            result = results[ex_name][run_name]
            for element, element_result in result.items():
                prediction = Ds[element]
                try:
                    error = (prediction - element_result) / element_result
                except (TypeError, ZeroDivisionError):
                    error = None
                print('For ' + element + ', actual vs predicted = ' + str(element_result) + ', ' + str(prediction) + '  [error: ' + str(error) + ']')
    
if __name__ == '__main__':
    main()
