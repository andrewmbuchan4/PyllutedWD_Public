#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import chemistry_info as ci
import geology_info as gi
import graph_factory as gf

# Yoshizaki 2020, ignoring H
mars_comp = {
    ci.Element.O: 0.53,
    ci.Element.Mg: 0.15,
    ci.Element.Si: 0.14,
    ci.Element.Fe: 0.1,
    ci.Element.Al: 0.013,
    ci.Element.S: 0.009,
    ci.Element.Ca: 0.01,
    ci.Element.Ni: 0.005
}

geo_model_earth = gi.GeologyModel()
geo_model_mars = gi.GeologyModel(mars_comp)

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
        #'MAP11': {
        #    'P': 15,
        #    'T': 2807,
        #    'fO2': -3.7
        #}, # This one fails to converge
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
        'Earth (assumed)': {
            'P': geo_model_earth.get_earth_differentiation_pressure(),
            'T': None,
            'fO2': geo_model_earth.get_earth_oxygen_fugacity()
        },
        'Mars (assumed)': {
            'P': geo_model_earth.get_mars_differentiation_pressure(),
            'T': None,
            'fO2': geo_model_earth.get_mars_oxygen_fugacity()
        },
        'Earth': {
            'P': 45,
            'T': None,
            'fO2': -1.3
        },
        'Mars': {
            'P': 5,
            'T': None,
            'fO2': -1.1
        }
    }
}

results = {
    'Corgne': {
        'PL-169': {
            ci.Element.Fe: 9.6,
            ci.Element.Mn: 0.013,
            ci.Element.Ni: 750,
            ci.Element.Cr: 0.51,
            ci.Element.Ga: 4.0,
            ci.Element.Si: 0.0014,
            ci.Element.Nb: 0.03,
            ci.Element.Ta: 0.0092,
            ci.Element.Ti: None,
            ci.Element.Cu: 20,
            ci.Element.Zn: 1.14
        },
        'PL-185': {
            ci.Element.Fe: 9.9,
            ci.Element.Mn: 0.03,
            ci.Element.Ni: 583,
            ci.Element.Cr: 0.86,
            ci.Element.Ga: 3.6,
            ci.Element.Si: 0.0014,
            ci.Element.Nb: 0.24,
            ci.Element.Ta: 0.02,
            ci.Element.Ti: None,
            ci.Element.Cu: 19,
            ci.Element.Zn: 0.39
        },
        'PR-369': {
            ci.Element.Fe: 12.5,
            ci.Element.Mn: 0.069,
            ci.Element.Ni: 301,
            ci.Element.Cr: 1.92,
            ci.Element.Ga: 4.7,
            ci.Element.Si: 0.0017,
            ci.Element.Nb: 0.72,
            ci.Element.Ta: 0.03,
            ci.Element.Ti: None,
            ci.Element.Cu: 11,
            ci.Element.Zn: 0.36
        },
        'PR-373': {
            ci.Element.Fe: 8.8,
            ci.Element.Mn: 0.03,
            ci.Element.Ni: 229,
            ci.Element.Cr: 0.85,
            ci.Element.Ga: 2.0,
            ci.Element.Si: 0.001,
            ci.Element.Nb: 0.11,
            ci.Element.Ta: 0.017,
            ci.Element.Ti: None,
            ci.Element.Cu: 13,
            ci.Element.Zn: 0.76
        },
        'PR-376': {
            ci.Element.Fe: 150,
            ci.Element.Mn: 0.84,
            ci.Element.Ni: 5051,
            ci.Element.Cr: 17.0,
            ci.Element.Ga: 64,
            ci.Element.Si: 0.33,
            ci.Element.Nb: 211,
            ci.Element.Ta: 10.9,
            ci.Element.Ti: 0.19,
            ci.Element.Cu: 86,
            ci.Element.Zn: 11.9
        },
        'PR-383': {
            ci.Element.Fe: 144,
            ci.Element.Mn: 0.86,
            ci.Element.Ni: 6109,
            ci.Element.Cr: 18.0,
            ci.Element.Ga: 51,
            ci.Element.Si: 0.34,
            ci.Element.Nb: 197,
            ci.Element.Ta: 10.8,
            ci.Element.Ti: 0.2,
            ci.Element.Cu: 92,
            ci.Element.Zn: 10.2
        },
        'PR-375': {
            ci.Element.Fe: 120,
            ci.Element.Mn: 1.1,
            ci.Element.Ni: 1861,
            ci.Element.Cr: 16.4,
            ci.Element.Ga: 40,
            ci.Element.Si: 0.36,
            ci.Element.Nb: 90,
            ci.Element.Ta: 5.9,
            ci.Element.Ti: 0.17,
            ci.Element.Cu: 49,
            ci.Element.Zn: 7.8
        },
        'PR-368': {
            ci.Element.Fe: 139,
            ci.Element.Mn: 1.12,
            ci.Element.Ni: 2619,
            ci.Element.Cr: 17.5,
            ci.Element.Ga: 45,
            ci.Element.Si: 0.35,
            ci.Element.Nb: 207,
            ci.Element.Ta: 11.1,
            ci.Element.Ti: 0.16,
            ci.Element.Cu: 50,
            ci.Element.Zn: 4.7
        },
        'PR-365': {
            ci.Element.Fe: 148,
            ci.Element.Mn: 1.05,
            ci.Element.Ni: 4665,
            ci.Element.Cr: 17.0,
            ci.Element.Ga: 72,
            ci.Element.Si: 0.33,
            ci.Element.Nb: 128,
            ci.Element.Ta: 9.9,
            ci.Element.Ti: 0.11,
            ci.Element.Cu: 92,
            ci.Element.Zn: 6.2
        }
    },
    'Fischer': {
        'DAC1': {
            ci.Element.Ni: 12,
            ci.Element.Co: 6.4,
            ci.Element.V: 0.43,
            ci.Element.Cr: None,
            ci.Element.Si: 0.129,
            ci.Element.Fe: 3.08
        },
        'DAC2': {
            ci.Element.Ni: 13,
            ci.Element.Co: 10,
            ci.Element.V: None,
            ci.Element.Cr: None,
            ci.Element.Si: 0.064,
            ci.Element.Fe: 3.8
        },
        'DAC3': {
            ci.Element.Ni: 7,
            ci.Element.Co: 6,
            ci.Element.V: None,
            ci.Element.Cr: 2,
            ci.Element.Si: 0.098,
            ci.Element.Fe: 3.6
        },
        'DAC4': {
            ci.Element.Ni: 5.1,
            ci.Element.Co: 4.3,
            ci.Element.V: 0.5,
            ci.Element.Cr: 1.7,
            ci.Element.Si: 0.053,
            ci.Element.Fe: 2.8
        },
        'DAC5': {
            ci.Element.Ni: 4.9,
            ci.Element.Co: 5.2,
            ci.Element.V: 0.86,
            ci.Element.Cr: None,
            ci.Element.Si: 0.352,
            ci.Element.Fe: 3.6
        },
        'MAP1': {
            ci.Element.Ni: 150,
            ci.Element.Co: 100,
            ci.Element.V: 0.216,
            ci.Element.Cr: None,
            ci.Element.Si: 0.08,
            ci.Element.Fe: 13
        },
        'MAP2': {
            ci.Element.Ni: 190,
            ci.Element.Co: 120,
            ci.Element.V: None,
            ci.Element.Cr: None,
            ci.Element.Si: None,
            ci.Element.Fe: 23
        },
        'MAP3': {
            ci.Element.Ni: 140,
            ci.Element.Co: 90,
            ci.Element.V: 0.3,
            ci.Element.Cr: None,
            ci.Element.Si: 0.101,
            ci.Element.Fe: 16.4
        },
        'MAP4': {
            ci.Element.Ni: 180,
            ci.Element.Co: 120,
            ci.Element.V: 0.42,
            ci.Element.Cr: 1.6,
            ci.Element.Si: 0.14,
            ci.Element.Fe: 24
        },
        'MAP5': {
            ci.Element.Ni: 170,
            ci.Element.Co: 100,
            ci.Element.V: 0.21,
            ci.Element.Cr: None,
            ci.Element.Si: 0.06,
            ci.Element.Fe: 14
        },
        'MAP6': {
            ci.Element.Ni: 99,
            ci.Element.Co: 64,
            ci.Element.V: None,
            ci.Element.Cr: None,
            ci.Element.Si: 0.08,
            ci.Element.Fe: 11.6
        },
        'MAP7': {
            ci.Element.Ni: 130,
            ci.Element.Co: 110,
            ci.Element.V: 0.34,
            ci.Element.Cr: 1.4,
            ci.Element.Si: 0.099,
            ci.Element.Fe: 18.5
        },
        'MAP8': {
            ci.Element.Ni: 71,
            ci.Element.Co: 65,
            ci.Element.V: 0.16,
            ci.Element.Cr: 0.7,
            ci.Element.Si: 0.02,
            ci.Element.Fe: 8.6
        },
        'MAP9': {
            ci.Element.Ni: 210,
            ci.Element.Co: 160,
            ci.Element.V: 0.83,
            ci.Element.Cr: 2.98,
            ci.Element.Si: 0.307,
            ci.Element.Fe: 37
        },
        'MAP10': {
            ci.Element.Ni: 140,
            ci.Element.Co: None,
            ci.Element.V: 0.36,
            ci.Element.Cr: 1.48,
            ci.Element.Si: 0.109,
            ci.Element.Fe: 19.2
        },
        'MAP11': {
            ci.Element.Ni: 370,
            ci.Element.Co: 220,
            ci.Element.V: 1.9,
            ci.Element.Cr: 5.1,
            ci.Element.Si: 0.348,
            ci.Element.Fe: 67
        },
        'MAP12': {
            ci.Element.Ni: 130,
            ci.Element.Co: 90,
            ci.Element.V: 0.64,
            ci.Element.Cr: 2.0,
            ci.Element.Si: 0.058,
            ci.Element.Fe: 17.7
        },
        'MAP13': {
            ci.Element.Ni: 110,
            ci.Element.Co: 60,
            ci.Element.V: 0.27,
            ci.Element.Cr: 1.0,
            ci.Element.Si: 0.009,
            ci.Element.Fe: 9.9
        },
        'MAP14': {
            ci.Element.Ni: 240,
            ci.Element.Co: 150,
            ci.Element.V: 0.53,
            ci.Element.Cr: 1.8,
            ci.Element.Si: 0.068,
            ci.Element.Fe: 25
        },
        'MAP15': {
            ci.Element.Ni: 170,
            ci.Element.Co: 120,
            ci.Element.V: 0.48,
            ci.Element.Cr: 1.7,
            ci.Element.Si: 0.07,
            ci.Element.Fe: 23.7
        }
    },
    'Test cases': {
        'Earth (assumed)': { #Values from Rudge 2010 supplementary table 5
            ci.Element.Ni: 10**1.418,
            ci.Element.Fe: 10**1.136,
            ci.Element.Cr: 10**0.195,
            #ci.Element.O: 0.1344,  # for use with Fischer O
            ci.Element.Si: 10**-0.728, # McDonough 2003 implies this should be somewhat higher
            ci.Element.O: 10**-0.6505,  # For use with Badro O Calculation: Badro 2015 gives core O as 2.7-5 wt%., i.e. ~ 9.42 - 17.45 mol %, i.e. 0.1344 +- 0.0402. Take mantle O as 0.6011. D = core concentration/ mantle conc. = 0.2236 +- 0.0669, so logD = -0.6505 += 0.1544 ish (taking lower bound to give largest error, to be conservative)
            ci.Element.W: 10**1.513,
            ci.Element.P: 10**1.398,
            ci.Element.Co: 10**1.381,
            #ci.Element.Pb: 10**1.159,
            ci.Element.Cu: 10**0.801,
            ci.Element.V: 10**0.262,
            ci.Element.Mn: 10**-0.155,
            ci.Element.Nb: 10**-0.276,
            ci.Element.Ta: 10**-0.611,
            ci.Element.Zn: 10**-0.824,
            ci.Element.Ga: 10**-1.0
            #ci.Element.Ti: 0
        },
        'Mars (assumed)': { # Values taken from table 10 in Yoshizaki 2020. D = Core/BSM
            ci.Element.Ni: 4.2/0.01,
            ci.Element.Fe: 48/4.4,
            ci.Element.O: 11/59,
            ci.Element.Si: 0,
            #ci.Element.Mg: 0,
            #ci.Element.Al: 0,
            ci.Element.S: 6.9/0.02,
            #ci.Element.Ca: 0,
            #ci.Element.H: 30/0.03
        },
        'Earth': { #Values from Rudge 2010 supplementary table 5
            ci.Element.Ni: 10**1.418,
            ci.Element.Fe: 10**1.136,
            ci.Element.Cr: 10**0.195,
            #ci.Element.O: 0.1344,  # for use with Fischer O
            ci.Element.O: 10**-0.6505,  # For use with Badro O Calculation: Badro 2015 gives core O as 2.7-5 wt%., i.e. ~ 9.42 - 17.45 mol %, i.e. 0.1344 +- 0.0402. Take mantle O as 0.6011. D = core concentration/ mantle conc. = 0.2236 +- 0.0669, so logD = -0.6505 += 0.1544 ish (taking lower bound to give largest error, to be conservative)
            ci.Element.Si: 10**-0.728, # McDonough 2003 implies this should be somewhat higher
            ci.Element.W: 10**1.513,
            ci.Element.P: 10**1.398,
            ci.Element.Co: 10**1.381,
            #ci.Element.Pb: 10**1.159,
            ci.Element.Cu: 10**0.801,
            ci.Element.V: 10**0.262,
            ci.Element.Mn: 10**-0.155,
            ci.Element.Nb: 10**-0.276,
            ci.Element.Ta: 10**-0.611,
            ci.Element.Zn: 10**-0.824,
            ci.Element.Ga: 10**-1.0
            #ci.Element.Ti: 0
        },
        'Mars': { # Values taken from table 10 in Yoshizaki 2020. D = Core/BSM
            ci.Element.Ni: 4.2/0.01,
            ci.Element.Fe: 48/4.4,
            ci.Element.Si: 0,
            ci.Element.O: 11/59,
            #ci.Element.Mg: 0,
            #ci.Element.Al: 0,
            ci.Element.S: 6.9/0.02,
            #ci.Element.Ca: 0,
            #ci.Element.H: 30/0.03
        }
    }
}

errors_on_logD = {  #NB these are errors on log(D), not on D!
    'Test cases': {
        'Earth (assumed)': {
            ci.Element.W: 0.077,
            ci.Element.Ni: 0.017,
            ci.Element.P: 0.14,
            ci.Element.Co: 0.013,
            ci.Element.Pb: 0.118,
            ci.Element.Fe: 0,  # This is listed without an error...
            ci.Element.Cu: 0.099,
            ci.Element.V: 0.042,
            ci.Element.Cr: 0.175,
            ci.Element.Mn: 0.274,
            ci.Element.Nb: 0.211,
            ci.Element.Ta: 0.195,
            ci.Element.Si: 0.136,
            ci.Element.Zn: 0.301,
            ci.Element.Ga: 0.301,
            #ci.Element.Ti: 0,  # This is listed without an error...
            ci.Element.O: 0.1544  # This is calculated by me! see above
        },
        'Mars (assumed)': { # Values not filled in yet
            ci.Element.O: 0,
            ci.Element.Mg: 0,
            ci.Element.Si: 0,
            ci.Element.Fe: 0,
            ci.Element.Al: 0,
            ci.Element.S: 0,
            ci.Element.Ca: 0,
            ci.Element.Ni: 0
            #ci.Element.H: 0
        },
        'Earth': {
            ci.Element.W: 0.077,
            ci.Element.Ni: 0.017,
            ci.Element.P: 0.14,
            ci.Element.Co: 0.013,
            ci.Element.Pb: 0.118,
            ci.Element.Fe: 0,  # This is listed without an error...
            ci.Element.Cu: 0.099,
            ci.Element.V: 0.042,
            ci.Element.Cr: 0.175,
            ci.Element.Mn: 0.274,
            ci.Element.Nb: 0.211,
            ci.Element.Ta: 0.195,
            ci.Element.Si: 0.136,
            ci.Element.Zn: 0.301,
            ci.Element.Ga: 0.301,
            #ci.Element.Ti: 0,  # This is listed without an error...
            ci.Element.O: 0.1544  # This is calculated by me! see above
        },
        'Mars': { # Values not filled in yet
            ci.Element.O: 0,
            ci.Element.Mg: 0,
            ci.Element.Si: 0,
            ci.Element.Fe: 0,
            ci.Element.Al: 0,
            ci.Element.S: 0,
            ci.Element.Ca: 0,
            ci.Element.Ni: 0
            #ci.Element.H: 0
        }
    }
}

observed_mantle_abundances = {
    'Earth (assumed)': { # By weight. McDonough 2003
        ci.Element.Ni: 0.001960,
        ci.Element.Fe: 0.0626,
        ci.Element.Cr: 0.002625,
        ci.Element.Si: 0.21,
        ci.Element.O: 0.44
    },
    'Mars (assumed)': { # This is by number. Yoshizaki 2020
        ci.Element.Ni: 0.0001,
        ci.Element.Fe: 0.044,
        ci.Element.Si: 0.16,
        ci.Element.O: 0.5959
    },
    'Earth': { # By weight. McDonough 2003
        ci.Element.Ni: 0.001960,
        ci.Element.Fe: 0.0626,
        ci.Element.Cr: 0.002625,
        ci.Element.Si: 0.21,
        ci.Element.O: 0.44
    },
    'Mars': { # This is by number. Yoshizaki 2020
        ci.Element.Ni: 0.0001,
        ci.Element.Fe: 0.044,
        ci.Element.Si: 0.16,
        ci.Element.O: 0.5959
    }
}

def main():
    graph_fac = gf.GraphFactory()
    Ds_dict = dict()
    for ex_name, ex_runs in experiments.items():
        if ex_name != 'Test cases':
            continue
        Ds_dict[ex_name] = dict()
        print('Running experiments from ' + ex_name)
        for run_name, run_params in ex_runs.items():
            print('Running experiment ' + run_name)
            P = run_params['P']
            T = run_params['T']
            fO2 = run_params['fO2']
            nbot = run_params.get('nbot', None)
            print('P = ' + str(P) + ', T = ' + str(T) + ', fO2 = ' + str(fO2) + ', nbot = ' + str(nbot))
            if run_name == 'Mars':
                a, w, Ds, all_Ds = geo_model_mars.form_a_planet_iteratively(P, fO2, T, nbot)
                geo_model_mars.tabulate_output(a, w, Ds)
                a_by_weight = geo_model_mars.convert_number_abundances_to_mass(a)
            else:
                a, w, Ds, all_Ds = geo_model_earth.form_a_planet_iteratively(P, fO2, T, nbot)
                geo_model_earth.tabulate_output(a, w, Ds)
                a_by_weight = geo_model_earth.convert_number_abundances_to_mass(a)
            Ds_dict[ex_name][run_name] = Ds
            result = results[ex_name][run_name]
            for element, element_result in result.items():
                prediction = Ds[element]
                try:
                    error = (prediction - element_result) / element_result
                except (TypeError, ZeroDivisionError):
                    error = None
                print('For ' + str(element) + ', actual vs predicted = ' + str(element_result) + ', ' + str(prediction) + '  [error: ' + str(error) + ']')
                log_prediction = np.log10(prediction)
                log_element_result = np.log10(element_result)
                try:
                    log_error = (log_prediction - log_element_result) / log_element_result
                except (TypeError, ZeroDivisionError):
                    log_error = None
                print('For ' + str(element) + ', log(actual) vs log(predicted) = ' + str(log_element_result) + ', ' + str(log_prediction) + '  [error: ' + str(log_error) + ']')
                if run_name == 'Mars':
                    predicted_mantle_abundance = a[element][gi.Layer.mantle]
                else:
                    #use weight
                    predicted_mantle_abundance = a_by_weight[element][gi.Layer.mantle]
                try:
                    observed_mantle_abundance = observed_mantle_abundances[run_name][element]
                except KeyError:
                    observed_mantle_abundance = None
                try:
                    error_ma = (predicted_mantle_abundance - observed_mantle_abundance) / observed_mantle_abundance
                except (TypeError, ZeroDivisionError):
                    error_ma = None
                print('For ' + str(element) + ', predicted_mantle_abundance vs observed_mantle_abundances = ' + str(predicted_mantle_abundance) + ', ' + str(observed_mantle_abundance) + '  [error: ' + str(error_ma) + ']')
            if ex_name == 'Test cases':
                print(ex_name)
                print(run_name)
                print(P)
                print(T)
                print(fO2)
                print(nbot)
                print(Ds)
                print(w)
                graph_fac.make_Ds_comparison_plot(results[ex_name][run_name], Ds_dict[ex_name][run_name], 'sisi', run_name, errors_on_logD[ex_name][run_name])
    print(results['Test cases']['Earth'])
    print(Ds_dict['Test cases']['Earth'])
    graph_fac.make_Ds_comparison_multipanel_plot(
        [results['Test cases']['Earth'], results['Test cases']['Mars']],
        [Ds_dict['Test cases']['Earth'], Ds_dict['Test cases']['Mars']],
        'sisi',
        ['Earth', 'Mars'],
        [errors_on_logD['Test cases']['Earth'], errors_on_logD['Test cases']['Mars']]
    )
    
if __name__ == '__main__':
    main()
