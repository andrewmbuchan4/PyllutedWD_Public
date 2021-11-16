#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import chemistry_info as ci
import geology_info as gi
import graph_factory as gf

def generate_mixes():
    geology_model = gi.GeologyModel(None, False)  # For now, not normalising the abundances because then N in mantle becomes negative!
    fcfs = [0.75, 0.99]
    mix_dict = dict()
    for fcf in fcfs:
        #Comment/Uncomment the next bits as appropriate:
        #        For Earth:
        #mix = geology_model.get_earth_mix(fcf)
        #        For Mars:
        #mix = geology_model.get_earth_mix(fcf, geology_model.mars_abundances)
        #        For what the model thinks is (roughly) Earth:
        #abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(geology_model.get_earth_differentiation_pressure(), geology_model.get_earth_oxygen_fugacity())
        #mix = geology_model.get_earth_mix(fcf, abundances)
        #        For what the model thinks is (roughly) Mars:
        #abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(geology_model.get_mars_differentiation_pressure(), geology_model.get_mars_oxygen_fugacity())
        #mix = geology_model.get_earth_mix(fcf, abundances)
        #        For Ideal HP Mantle (low fO2 to make Ni siderophilic as poss)
        #abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(60, -3)
        #mix = geology_model.get_earth_mix(fcf, abundances)
        #        For Ideal HP Core (high fO2 to make Cr/Si lithophilic as poss)
        #abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(60, -1)
        #mix = geology_model.get_earth_mix(fcf, abundances)
        #        For Ideal LP Mantle (low fO2 to make Ni siderophilic as poss)
        #abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(0, -3)
        #mix = geology_model.get_earth_mix(fcf, abundances)
        #        For Ideal LP Core (high fO2 to make Cr/Si lithophilic as poss)
        abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(0, -1)
        mix = geology_model.get_earth_mix(fcf, abundances)
        mix_dict[fcf] = mix
    return mix_dict

def generate_observations(mix_dict, high_pressure_idealisation=False, low_pressure_idealisation=False):
    idealisation_shifts = {
        'High Pressure': {
            'Mantle': {
                ci.Element.Ni: -0.1,
                ci.Element.Fe: -0.2,
                ci.Element.Cr: -0.25,
                ci.Element.Si: -0.2
            },
            'Core': {
                ci.Element.Ni: 0.2,
                ci.Element.Fe: 0.2,
                ci.Element.Cr: 0.3,
                ci.Element.Si: 0.3
            }
        },
        'Low Pressure': {
            'Mantle': {
                ci.Element.Ni: -1.2,
                ci.Element.Fe: -0.2,
                ci.Element.Cr: 0.1,
                ci.Element.Si: -0.2
            },
            'Core': {
                ci.Element.Ni: 0.2,
                ci.Element.Fe: 0.2,
                ci.Element.Cr: -0.2,
                ci.Element.Si: -1.2
            }
        }
    }
    # Generate synthetic observations, assuming an arbitrary amount of H/He
    HorHe = 10000
    error = 0.05  # Add an arbitrary error to each observation
    obs_dict = dict()
    for fcf, mix in mix_dict.items():
        obs_dict[fcf] = dict()
        for element in ci.usual_elements:
            obs_value = np.log10(mix[element]/HorHe)
            if high_pressure_idealisation:
                if fcf < 0.17:
                    obs_value += idealisation_shifts['High Pressure']['Mantle'].get(element, 0)
                else:
                    obs_value += idealisation_shifts['High Pressure']['Core'].get(element, 0)
            if low_pressure_idealisation:
                if fcf < 0.17:
                    obs_value += idealisation_shifts['Low Pressure']['Mantle'].get(element, 0)
                else:
                    obs_value += idealisation_shifts['Low Pressure']['Core'].get(element, 0)
            obs_dict[fcf][element] = (obs_value, error)
    return obs_dict

def print_synthetic_observations(obs_dict):
    for fcf, obs_set in obs_dict.items():
        strings_to_print = list()
        for element in ci.usual_elements:
            strings_to_print.append(str(obs_set[element][0]) + ',' + str(obs_set[element][1]))
        print(','.join(strings_to_print))

def plot_mixes(mix_dict):
    graph_fac = gf.GraphFactory()
    for fcf, mix in mix_dict.items():
        modified_mix_dict = dict()
        for el in [ci.Element.O, ci.Element.Si, ci.Element.Fe, ci.Element.Ni, ci.Element.Mg]:
            if mix[el] > 0.02:
                modified_mix_dict[el] = mix[el]
        modified_mix_dict[ci.Element.Placeholder] = 1 - sum(x[1] for x in modified_mix_dict.items())
        graph_fac.plot_mix_pie_chart(fcf, modified_mix_dict)

def main():
    mix_dict = generate_mixes()
    plot_mixes(mix_dict)
    obs_dict = generate_observations(mix_dict, False, False)
    print_synthetic_observations(obs_dict)

if __name__ == '__main__':
    main()
