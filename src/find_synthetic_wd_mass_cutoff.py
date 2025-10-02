#!/usr/bin/env python

import numpy as np

import complete_model as cm
import manager as mn
import model_parameters as mp
import physical_constants as pc
import synthetic_observer as so
import synthetic_population as sp
import timescale_interpolator as ti

def find_mdot_cutoff_as_function_of_teff(spectral_type, threshold_type, teff): # Returns log10(kg per Myr)

    tolerance = 1 # The nearest gram ser second should be adequate!
    converged = False
    lower_mdot = 1 # grams per second
    upper_mdot = 1000000000000000000 # grams per second
    while not converged:
        mdot = (lower_mdot + upper_mdot)/2
        print()
        print(mdot)
        accretion_timescale = 100 # Myr  -> w.l.o.g just fix this to a large number which can accommodate any realistic situation (need to allow time to settle into steady state)

        mass = np.log10((mdot/1000) * accretion_timescale * 1000000 * pc.seconds_per_year) # in log(kg)

        print(mass)

        input_dict = {
            mp.ModelParameter.metallicity: 478,
            mp.ModelParameter.t_sinceaccretion: 10, #Myr -> should be long enough to reach steady state, but less than accretion_timescale
            mp.ModelParameter.formation_distance: 2,
            mp.ModelParameter.feeding_zone_size: 0.05,
            mp.ModelParameter.parent_core_frac: None,
            mp.ModelParameter.parent_crust_frac: None,
            mp.ModelParameter.fragment_core_frac: 0.17,
            mp.ModelParameter.fragment_crust_frac: 0,
            mp.ModelParameter.fragment_mass: mass,
            mp.ModelParameter.accretion_timescale: accretion_timescale,
            mp.ModelParameter.pressure: 45,
            mp.ModelParameter.oxygen_fugacity: -2
        }

        pollution_abundances, diagnostics = cm.complete_model_calculation(
            input_dict[mp.ModelParameter.metallicity],
            input_dict[mp.ModelParameter.t_sinceaccretion],
            input_dict[mp.ModelParameter.formation_distance],
            input_dict[mp.ModelParameter.feeding_zone_size],
            input_dict[mp.ModelParameter.parent_core_frac],
            input_dict[mp.ModelParameter.parent_crust_frac],
            input_dict[mp.ModelParameter.fragment_core_frac],
            input_dict[mp.ModelParameter.fragment_crust_frac],
            input_dict[mp.ModelParameter.fragment_mass],
            input_dict[mp.ModelParameter.accretion_timescale],
            input_dict[mp.ModelParameter.pressure],
            input_dict[mp.ModelParameter.oxygen_fugacity],
            'NonEarthlike',
            False
        )
        wd_properties = {
            mp.WDParameter.spectral_type: spectral_type,
            mp.WDParameter.temperature: teff
        }
        pollution_properties = dict()
        test_wd = sp.SyntheticSystem(wd_properties, pollution_properties, pollution_abundances)

        error_dict = dict()

        observer = so.Observer(so.ObservationType.TeffIndividualElementCutoff, error_dict, threshold_type)
        observer.observe_system(test_wd)
        print(test_wd.observed)
        if test_wd.observed:
            # Then this is the new upper limit:
            upper_mdot = mdot
        else:
            lower_mdot = mdot
        converged = abs(upper_mdot - lower_mdot) < tolerance
    print(teff)
    print(upper_mdot)
    print(lower_mdot)
    # We want to convert this to log10(kg/Myr)
    toret = np.log10((lower_mdot/1000) * 1000000 * pc.seconds_per_year) # lower_mdot is the conservative choice
    return toret

def main():
    print('Remember to allow at least 0.4 dex leeway for errors')
    manager = mn.Manager()
    manager.publish_live_data(0, ti.TimescaleType.KoesterOvershoot) # The timescale choice shouldn't matter (to order-of-magnitude) as long as system has reached steady state
    min_mdot_dict = dict()
    for teff in np.linspace(3000, 20000, 10):
        min_mdot = find_mdot_cutoff_as_function_of_teff('DA', 'ELB_DT', teff)
        min_mdot_dict[teff] = min_mdot
    print(min_mdot_dict)
    for k, v in min_mdot_dict.items():
        print(str(k) + ',' + str(v))

if __name__ == '__main__':
    main()
