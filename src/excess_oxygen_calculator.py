#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum
import copy
import numpy as np

import chemistry_info as ci
import geology_info as gi
import graph_factory as gf

class OxidationStrategy(Enum):
    default = 0
    conservative = 1

    def __str__(self):
        return self.name

class ExcessOxygenCalculator:

    def __init__(self):
        self.oxidation_dict = {
            OxidationStrategy.default: {
                ci.Element.Mg: 1, #MgO
                ci.Element.Si: 2, #SiO2
                ci.Element.Al: 1.5, #Al2O3
                ci.Element.Ca: 1, #CaO
                ci.Element.Fe: 1, #FeO
            },
            OxidationStrategy.conservative: {
                ci.Element.Mg: 1, #MgO
                ci.Element.Si: 2, #SiO2
                ci.Element.Al: 1.5, #Al2O3
                ci.Element.Ca: 1, #CaO
                #ci.Element.C: 2, #CO2
                ci.Element.Fe: 1.5, #Fe2O3
            }
        }
        self.species_names = {
            OxidationStrategy.default: {
                ci.Element.Mg: 'MgO',
                ci.Element.Si: 'SiO$_2$',
                ci.Element.Al: 'Al$_2$O$_3$',
                ci.Element.Ca: 'CaO',
                ci.Element.Fe: 'FeO',
            },
            OxidationStrategy.conservative: {
                ci.Element.Mg: 'MgO',
                ci.Element.Si: 'SiO$_2$',
                ci.Element.Al: 'Al$_2$O$_3$',
                ci.Element.Ca: 'CaO',
                ci.Element.Fe: 'Fe$_2$O$_3$',
                ci.Element.C: 'CO$_2$',
            }
        }
        dummy_gm = gi.GeologyModel()
        self.solar = dummy_gm.solar_abundances
        self.layer_priority = [gi.Layer.mantle, gi.Layer.crust, gi.Layer.bulk, gi.Layer.core]
        # Next line is for the excess O/Teff plot ONLY! - shouldnt be necessary now?
        #self.layer_priority = [gi.Layer.bulk, gi.Layer.mantle, gi.Layer.crust, gi.Layer.core]

    def get_favoured_layer(self, abundances):
        for layer in self.layer_priority:
            if abundances[ci.Element.O].get(layer, None) is not None:
                return layer
        return None

    def normalise_abundances(self, abundances):
        total_abundances = {
            gi.Layer.bulk: 0,
            gi.Layer.core: 0,
            gi.Layer.mantle: 0,
            gi.Layer.crust: 0
        }
        for layer in total_abundances.keys():
            for element in abundances.keys():
                try:
                    total_abundances[layer] += abundances[element][layer]
                except KeyError:
                    pass
        for layer in total_abundances.keys():
            for element in abundances.keys():
                try:
                    abundances[element][layer] /= total_abundances[layer]
                except KeyError:
                    pass
        return abundances

    def assign_oxygen(self, normalised_abundances, favoured_layer, ox_strat):
        toret = dict()
        for element, ratio in self.oxidation_dict[ox_strat].items():
            toret[element] = ratio*normalised_abundances[element][favoured_layer]
        return toret

    def calculate_excess_oxygen(self, normalised_abundances, favoured_layer, ox_strat):
        oxygen_assignations = self.assign_oxygen(normalised_abundances, favoured_layer, ox_strat)
        total_oxygen_assigned = sum(v for v in oxygen_assignations.values())
        excess_oxygen = normalised_abundances[ci.Element.O][favoured_layer] - total_oxygen_assigned
        fractional_excess_oxygen = excess_oxygen/normalised_abundances[ci.Element.O][favoured_layer]
        return excess_oxygen, fractional_excess_oxygen, oxygen_assignations

    def calculate_excess_oxygen_mass(self, number_abundance_dict, favoured_layer, fractional_excess_oxygen, total_mass=1):
        geo_model = gi.GeologyModel()
        mass_abundance_dict = geo_model.convert_number_abundances_to_mass(number_abundance_dict)
        o_mass_fraction = mass_abundance_dict[ci.Element.O][favoured_layer]
        excess_o_mass = fractional_excess_oxygen*o_mass_fraction*total_mass
        return excess_o_mass

    def calculate_water_abundance(self, number_abundance_dict, favoured_layer, fractional_excess_oxygen, total_mass=1):
        excess_o_mass = self.calculate_excess_oxygen_mass(number_abundance_dict, favoured_layer, fractional_excess_oxygen, total_mass)
        o_element_mass = ci.get_element_mass(ci.Element.O)
        h_element_mass = ci.get_element_mass(ci.Element.H)
        h2o_element_mass = (2*h_element_mass) + o_element_mass
        water_mass = (h2o_element_mass*excess_o_mass)/o_element_mass
        water_mass_fraction = water_mass/total_mass
        return water_mass, water_mass_fraction, excess_o_mass

    def calculate_sigma_significance(self, starting_abundances, favoured_layer, ox_strat, o_error):
        if o_error == 0:
            return None
        o_lower_bound = 0
        o_upper_bound = 10
        tolerance = 0.000001
        converged_o = None
        converged_abundances = None
        while converged_o is None:
            # Do a simple binary search
            test_o = 0.5*(o_lower_bound + o_upper_bound)
            test_abundances = copy.deepcopy(starting_abundances)
            test_abundances[ci.Element.O] = test_o
            reformatted_abundances = self.reformat_abundances(test_abundances, favoured_layer)
            normalised_abundances = self.normalise_abundances(reformatted_abundances)
            excess_oxygen, fractional_excess_oxygen, oxygen_assignations = self.calculate_excess_oxygen(normalised_abundances, favoured_layer, ox_strat)
            diff = abs(excess_oxygen)
            if diff < tolerance:
                converged_abundances = normalised_abundances
                converged_o = normalised_abundances[ci.Element.O][favoured_layer]
            else:
                if excess_oxygen > 0:
                    # Need to reduce the O abundance
                    o_upper_bound = test_o
                else:
                    o_lower_bound = test_o
        o_factor = converged_o/starting_abundances[ci.Element.O][favoured_layer]
        delta_log_O = np.log10(o_factor)
        sigma_significance = abs(delta_log_O/o_error)
        return sigma_significance

    def reformat_abundances(self, abundances, favoured_layer=gi.Layer.bulk):
        # Designed to take a {Element: Value} and return {Element: {Layer.bulk: Value}} if needed
        reformatted_dict = dict()
        for element, value in abundances.items():
            if isinstance(value, dict):
                reformatted_dict[element] = value
            else:
                reformatted_dict[element] = {favoured_layer: value}
        return reformatted_dict

    def check_for_nans(self, abundances, favoured_layer):
        any_nans = False
        for el, val in abundances.items():
            if val is None:
                any_nans = True
            else:
                if isinstance(val, dict):
                    test_val = val.get(favoured_layer)
                    if test_val is None or np.isnan(test_val):
                        any_nans = True
                else:
                    if np.isnan(val):
                        any_nans = True
        return any_nans

    def calculate_steady_state_composition(self, wd_abundances, predicted_composition, wd_timescales, reference_element, reference_layer):
        # Composition of data points, if they were in steady state (using predicted composition to fill in missing elements)
        # Firstly, assemble the composition we're going to use
        to_use = dict()
        for element, default in predicted_composition.items():
            if element in wd_abundances and wd_abundances[element] != 0:
                raw_abundance = 10**(wd_abundances[element] - wd_abundances[reference_element])
                scaling_factor = wd_timescales[reference_element]/wd_timescales.get(element)
                to_use[element] = {reference_layer: raw_abundance*scaling_factor}
            else:
                to_use[element] = {reference_layer: default[reference_layer]/predicted_composition[reference_element][reference_layer]}
        # Everything is now scaled to reference_element (unnormalised), and corrected for sinking where needed
        # Next scale everything by its sinking timescales
        normalised_abundances = self.normalise_abundances(to_use)
        return normalised_abundances

    def run_full_calculation(self, abundances_dict, total_mass=1, include_conservative_solarFe_composition=False, include_dataO_composition=False, observations=None, errors={ci.Element.O: 0}, wd_timescales=None):
        eoc_stat_dict = dict()
        additional_abundances = dict()
        o_error = errors.get(ci.Element.O, 0)
        for composition_name, abundances in abundances_dict.items():
            if abundances is not None:
                reformatted_abundances = self.reformat_abundances(abundances)
                normalised_abundances = self.normalise_abundances(reformatted_abundances)
                favoured_layer = self.get_favoured_layer(reformatted_abundances)
                any_nans = self.check_for_nans(abundances, favoured_layer)
                eoc_stat_dict[composition_name] = dict()
                for ox_strat in self.oxidation_dict.keys():
                    if any_nans:
                        eoc_stat_dict[composition_name][ox_strat] = None
                    else:
                        excess_oxygen, fractional_excess_oxygen, oxygen_assignations = self.calculate_excess_oxygen(normalised_abundances, favoured_layer, ox_strat)
                        water_mass, water_mass_fraction, excess_o_mass = self.calculate_water_abundance(normalised_abundances, favoured_layer, fractional_excess_oxygen, total_mass)
                        sigma_significance = self.calculate_sigma_significance(normalised_abundances, favoured_layer, ox_strat, o_error)
                        eoc_stat_dict[composition_name][ox_strat] = excess_oxygen, fractional_excess_oxygen, oxygen_assignations, water_mass, water_mass_fraction, excess_o_mass, normalised_abundances, self.species_names[ox_strat], sigma_significance, favoured_layer
                # NB: If we ever include conservative compositions that involve altering all elements AND not forcing a fit to data,
                # need to improve the logic in calculate_conservative_composition - see the TODO there. Also TODO: why are we manually setting everything to bulk here?
                if include_conservative_solarFe_composition:
                    solar_Fe_comp = self.calculate_conservative_composition(normalised_abundances, None, ci.Element.Mg, gi.Layer.bulk, [ci.Element.Fe])
                    #additional_abundances[composition_name + ', allow solar Fe'] = solar_Fe_comp
                if include_dataO_composition:
                    data_O_comp = self.calculate_conservative_composition(normalised_abundances, observations, ci.Element.Mg, gi.Layer.bulk, [ci.Element.O], True)
                    #additional_abundances[composition_name + ', data O'] = data_O_comp
                if include_dataO_composition:
                    data_comp = self.calculate_conservative_composition(normalised_abundances, observations, ci.Element.Mg, gi.Layer.bulk, None, True)
                    #additional_abundances[composition_name + ', data'] = data_comp
                if wd_timescales is not None:
                    ss_comp = self.calculate_steady_state_composition(observations, normalised_abundances, wd_timescales, ci.Element.Mg, gi.Layer.bulk)
                    additional_abundances[composition_name + ', data (SS)'] = ss_comp
        for composition_name, abundances in additional_abundances.items():
            favoured_layer = self.get_favoured_layer(abundances)
            eoc_stat_dict[composition_name] = dict()
            for ox_strat in self.oxidation_dict.keys():
                excess_oxygen, fractional_excess_oxygen, oxygen_assignations = self.calculate_excess_oxygen(abundances, favoured_layer, ox_strat)
                water_mass, water_mass_fraction, excess_o_mass = self.calculate_water_abundance(abundances, favoured_layer, fractional_excess_oxygen, total_mass)
                sigma_significance = self.calculate_sigma_significance(abundances, favoured_layer, ox_strat, o_error)
                eoc_stat_dict[composition_name][ox_strat] = excess_oxygen, fractional_excess_oxygen, oxygen_assignations, water_mass, water_mass_fraction, excess_o_mass, abundances, self.species_names[ox_strat], sigma_significance, favoured_layer
        return eoc_stat_dict

    def calculate_conservative_composition(self, nominal_model, observations, reference_element, reference_layer, elements_to_alter=None, force_data_fit=False):
        #if elements_to_alter is None, assume we need to alter all of them
        #nominal_model is linear, observations and solar are on a log scale
        # TODO: Assuming maximum abundance for all non-O elements is definitely not the most conservative composition in general
        # (because you have to normalise everything afterwards)
        # Need to prioritise maximising the more oxidising elements somehow
        toret = dict()
        assert reference_element in nominal_model
        reformatted_abundances = self.reformat_abundances(nominal_model, reference_layer)
        for element, nominal_entry in reformatted_abundances.items():
            possible_abundances = list()
            observed_abundance = None
            solar_abundance = None
            if observations is not None and element in observations and reference_element in observations and observations.get(element, 0) != 0:
                observed_abundance = 10**(observations[element] - observations[reference_element])
                possible_abundances.append(observed_abundance)
            if element in self.solar and reference_element in self.solar and not force_data_fit:
                solar_abundance = 10**(self.solar[element] - self.solar[reference_element])
                possible_abundances.append(solar_abundance)
            nominal_abundance = nominal_entry[reference_layer]/reformatted_abundances[reference_element][reference_layer]
            possible_abundances.append(nominal_abundance)
            if elements_to_alter is None or (elements_to_alter is not None and element in elements_to_alter):
                if force_data_fit:
                    if observed_abundance is None:
                        conservative_abundance = nominal_abundance
                    else:
                        conservative_abundance = observed_abundance
                else:
                    if element == ci.Element.O:
                        conservative_abundance = min(possible_abundances)
                    else:
                        conservative_abundance = max(possible_abundances)
            else:
                conservative_abundance = nominal_abundance
            toret[element] = {reference_layer: conservative_abundance}
        normalised_abundances = self.normalise_abundances(toret)
        return normalised_abundances

def main():
    eoc = ExcessOxygenCalculator()
    input_dict = {
        ci.Element.Al:{gi.Layer.bulk: 0.077449061001879},
        ci.Element.Ti:{gi.Layer.bulk: 0.002509694000061},
        ci.Element.Ca:{gi.Layer.bulk: 0.065982759001601},
        ci.Element.Ni:{gi.Layer.bulk: 0.054985632001334},
        ci.Element.Fe:{gi.Layer.bulk: 0.912114606022125},
        ci.Element.Cr:{gi.Layer.bulk: 0.013812800000335},
        ci.Element.Mg:{gi.Layer.bulk: 1.00000000002426},
        ci.Element.Si:{gi.Layer.bulk: 0.978014541023723},
        ci.Element.Na:{gi.Layer.bulk: 0.050125881001216},
        ci.Element.O: {gi.Layer.bulk: 15.1257355},
        ci.Element.C: {gi.Layer.bulk: 7.22713864317531},
        ci.Element.N: {gi.Layer.bulk: 1.95037272104731}
    }
    observations = {
        ci.Element.Ca: -6.314824132,
        ci.Element.Mg: -5.234735014,
        ci.Element.Si: -5.154189483,
        ci.Element.O: -4.234890122
    }
    mass = 1 # There is an inherent inconsistency in that we later assume that this mass includes H, although H is not present in the dict
    # Basically we assume M_hydrogen << M_object -> put in paper!
    dummy_gm = gi.GeologyModel()
    earth = dummy_gm.element_info
    print(earth)
    conservative_input_dict = eoc.calculate_conservative_composition(input_dict, observations, ci.Element.Si, gi.Layer.bulk)
    conservative_input_dict_fe_only = eoc.calculate_conservative_composition(input_dict, observations, ci.Element.Mg, gi.Layer.bulk, [ci.Element.Fe])
    print(conservative_input_dict)
    input_compositions = {
        'Model': input_dict,
        'Conservative': conservative_input_dict,
        'Conservative Fe': conservative_input_dict_fe_only
    }
    eoc_stats = eoc.run_full_calculation(input_compositions, mass, False, False, None, {ci.Element.O: 0.1})
    print(eoc_stats)
    print('Earth:')
    print(earth)
    eoc_stats = eoc.run_full_calculation({'Earth': {el: earth[el] for el in ci.usual_elements}}, mass)
    print(eoc_stats)
    graph_fac = gf.GraphFactory()
    graph_fac.make_excess_oxygen_plot(eoc_stats, 'test')

if __name__ == '__main__':
    main()
