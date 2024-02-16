#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import chemistry_info as ci
import geology_info as gi


target_logds_dict = {
    'Earth': {
        #Values from Rudge 2010 supplementary table 5
        ci.Element.Ni: 1.418,
        ci.Element.Fe: 1.136,
        ci.Element.Cr: 0.195,
        ci.Element.O: -0.6505,  # For use with Badro O Calculation: Badro 2015 gives core O as 2.7-5 wt%., i.e. ~ 9.42 - 17.45 mol %, i.e. 0.1344 +- 0.0402. Take mantle O as 0.6011. D = core concentration/ mantle conc. = 0.2236 +- 0.0669, so logD = -0.6505 += 0.1544 ish (taking lower bound to give largest error, to be conservative)
        ci.Element.Si: -0.728, # McDonough 2003 implies this should be somewhat higher
        ci.Element.W: 1.513,
        ci.Element.P: 1.398,
        ci.Element.Co: 1.381,
        ci.Element.Cu: 0.801,
        ci.Element.V: 0.262,
        ci.Element.Mn: -0.155,
        ci.Element.Nb: -0.276,
        ci.Element.Ta: -0.611,
        ci.Element.Zn: -0.824,
        ci.Element.Ga: -1.0
    },
    'Mars': { # Values taken from table 10 in Yoshizaki 2020. D = Core/BSM
        ci.Element.Ni: np.log10(4.2/0.01),
        ci.Element.Fe: np.log10(48/4.4),
        ci.Element.O: np.log10(11/59),
        ci.Element.Si: np.log10(0.001),
        ci.Element.S: np.log10(6.9/0.02)
    }
}

target_logd_errors_dict = {
    'Earth': {
        ci.Element.W: 0.077,
        ci.Element.Ni: 0.017,
        ci.Element.P: 0.14,
        ci.Element.Co: 0.013,
        ci.Element.Pb: 0.118,
        ci.Element.Fe: 0.013, # This is arbitrary - no error was given because it's well constrained! Just made it equal to Ni's error (should maybe be Co?)
        ci.Element.Cu: 0.099,
        ci.Element.V: 0.042,
        ci.Element.Cr: 0.175,
        ci.Element.Mn: 0.274,
        ci.Element.Nb: 0.211,
        ci.Element.Ta: 0.195,
        ci.Element.Si: 0.136,
        ci.Element.Zn: 0.301,
        ci.Element.Ga: 0.301,
        ci.Element.O: 0.1544  # This is calculated by me! see above
    },
    'Mars': { # Unknown! Assuming all equal for simplicity (although O might be larger?)
        ci.Element.Ni: 0.1,
        ci.Element.Fe: 0.1,
        ci.Element.O: 0.1,
        ci.Element.Si: 0.1,
        ci.Element.S: 0.1
    }
}

# For if we care more about some elements that others
penalty_weights = {
    ci.Element.W: 1,
    ci.Element.Ni: 1,
    ci.Element.P: 1,
    ci.Element.Co: 1,
    ci.Element.Pb: 1,
    ci.Element.Fe: 1,
    ci.Element.Cu: 1,
    ci.Element.V: 1,
    ci.Element.Cr: 1,
    ci.Element.Mn: 1,
    ci.Element.Nb: 1,
    ci.Element.Ta: 1,
    ci.Element.Si: 1,
    ci.Element.Zn: 1,
    ci.Element.Ga: 1,
    ci.Element.O: 1,
    ci.Element.S: 1
}

def chi_squared(ds, target_name):
    target_logds = target_logds_dict[target_name]
    target_logd_errors = target_logd_errors_dict[target_name]
    total_penalty = 0
    degrees_of_freedom = 0
    for element, target_value in target_logds.items():
        error = target_logd_errors[element]
        calculated_value = np.log10(ds.get(element))
        if calculated_value is None:
            print('Warning! No calculated value for ' + str(element))
            print(ds)
            raise
        weight = penalty_weights[element]
        penalty = (((target_value - calculated_value)**2)/(error**2))*weight # should probably divide by degrees of freedom to make this a reduced chi squared
        total_penalty += penalty
    return total_penalty

def run_search(target_name):
    geo_model = gi.GeologyModel()
    pressure_values = np.arange(0, 60, 0.5)
    fo2_values = np.arange(-3, -1, 0.01)
    min_penalty = np.inf
    best_combo = None
    for p in pressure_values:
        for fo2 in fo2_values:
            composition, pcnf, ds, diagnostics = geo_model.form_a_planet_iteratively(p, fo2)
            if pcnf is None or np.isnan(pcnf) or ds is None:
                print('Warning! Model did not converge')
            else:
                chi_sq = chi_squared(ds, target_name)
                print()
                print(' --- Trying case:')
                print(' --- Pressure = ' + str(p) + ' GPa')
                print(' --- fo2 = IW' + str(fo2))
                print(' --- Penalty = ' + str(chi_sq))
                if chi_sq < min_penalty:
                    min_penalty = chi_sq
                    best_combo = (p, fo2)
    print()
    print()
    print()
    print('The Best Combo Is...')
    print(' --- Pressure = ' + str(best_combo[0]) + ' GPa')
    print(' --- fo2 = IW' + str(best_combo[1]))
    print(' --- Penalty = ' + str(min_penalty))
    composition, pcnf, ds, diagnostics = geo_model.form_a_planet_iteratively(best_combo[0], best_combo[1])
    geo_model.tabulate_output(composition, pcnf, ds)
    
def main():
    run_search('Earth')
    #run_search('Mars')
    
if __name__ == '__main__':
    main()
