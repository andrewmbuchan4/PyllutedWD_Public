#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections as cn
import csv
from enum import Enum
import itertools
import numpy as np
import scipy.interpolate as si
import scipy.stats as st

import abundance_model as am
import chemistry_info as ci
import geology_info as gi
import graph_factory as gf
import model_parameters as mp
import prior_functions as pf
import pwd_utils as pu
import synthetic_population as sp
import timescale_interpolator as ti
import white_dwarf_model as wdm

class ModellerType(Enum):
    SimpleFcfInterpolation = 0
    AnalyticApproximation = 1
    GridInterpolation = 2
    Null = 3

    def __str__(self):
        return self.name

# TODO add tests for this bit
class SyntheticGrid:
#Consists of a SyntheticPopulation and a dict of interpolators
    def __init__(self, synthetic_grid_filename, regular_grid=True):
        self.population = sp.SyntheticPopulation(None, None, None, synthetic_grid_filename)
        # Let's start by just considering 3 elements and expand later...
        self.input_elements = [ci.Element.Ca, ci.Element.Fe, ci.Element.Na]
        # Let's start by just considering 2 variables and expand later...
        self.variables_to_interpolate = [mp.ModelParameter.fragment_core_frac, mp.ModelParameter.formation_distance]
        self.boundaries = {
            mp.ModelParameter.fragment_core_frac: {pf.Limit.Lower: 0, pf.Limit.Upper: 1},
            mp.ModelParameter.formation_distance: {pf.Limit.Lower: -2, pf.Limit.Upper: 2}
        }
        if regular_grid:
            element_vals_dict = dict()
            for element in self.input_elements:
                element_vals_dict[element] = list(set(self.population.observed_abundances(element)))
                element_vals_dict[element].sort()

            grids = dict()
            for vti in self.variables_to_interpolate:
                shape = [len(element_vals_dict[element]) for element in self.input_elements]
                grids[vti] = np.full(shape, None)
            list_of_lists = [range(len(element_vals_dict[element])) for element in self.input_elements]
            for grid_indices in itertools.product(*list_of_lists):
                element_dictionary = dict()
                for i, element in enumerate(self.input_elements):
                    element_dictionary[element] = element_vals_dict[element][grid_indices[i]]
                matching_systems = self.population.find_systems_with_abundances(element_dictionary)
                if len(matching_systems) < 1:
                    print('Warning! Grid point had no corresponding reference system, skipping')
                elif len(matching_systems) == 1:
                    # We want to use the initial input values as the reference values for the interpolation grid
                    if matching_systems[0].pollution_properties is None:
                        print('Warning! Grid point was missing modelled properties, skipping')
                    else:
                        for vti in self.variables_to_interpolate:
                            grids[vti][grid_indices] = matching_systems[0].modelled_properties.get(vti, None)
                else:
                    raise NotImplementedError('Found >1 matching reference system. Which to use?')
            # 'linear', False, None means linear interpolation, no error if we go out of bounds, extrapolate in that case
            self.interpolators = dict()  # TODO: Can we just have one grid, and interpolate multiple variables on that grid?
            for vti in self.variables_to_interpolate:
                self.interpolators[vti] = si.RegularGridInterpolator(
                    (element_vals_dict[ci.Element.Ca], element_vals_dict[ci.Element.Fe], element_vals_dict[ci.Element.Na]),
                    grids[vti],
                    'linear',
                    False,
                    None
                )
        else:
            # This is not really needed/tested
            data_points = list()
            variable_vals_dict = dict()
            for vti in self.variables_to_interpolate:
                variable_vals_dict[vti] = list()
            for system in self.population:
                data_point = [system.observed_abundances[element] for element in self.input_elements]
                data_points.append(data_point)
                for vti in self.variables_to_interpolate:
                    variable_vals_dict[vti].append(system.pollution_properties[vti])
            self.interpolators = dict()
            for vti in self.variables_to_interpolate:
                self.interpolators[vti] = si.LinearNDInterpolator(
                    data_points,
                    variable_vals_dict[vti]
                )

    def interpolate_modelled_properties(self, system):
        list_to_sample = list()
        for element in self.input_elements:
            list_to_sample.append(system.observed_abundances[element])  # If we don't have this, error out
        point_to_sample = np.array([list_to_sample])
        interpolated_properties = dict()
        for vti in self.variables_to_interpolate:
            interpolated_value = self.interpolators[vti](point_to_sample)[0]

            if vti in self.boundaries:
                if self.boundaries[vti].get(pf.Limit.Lower, None) is not None:
                    interpolated_value = max(interpolated_value, self.boundaries[vti].get(pf.Limit.Lower, None))
                if self.boundaries[vti].get(pf.Limit.Upper, None) is not None:
                    interpolated_value = min(interpolated_value, self.boundaries[vti].get(pf.Limit.Upper, None))

            interpolated_properties[vti] = [interpolated_value]
        # Things to take into account somehow:
        # Effect of sinking
        # Some normalisation to control for pollution fraction. Maybe dimensions should be [X/Mg]?
        # Control for WD properties?
        return interpolated_properties

class Modeller:

    def __init__(self, modeller_type, modeller_args=None):
        self.modeller_type = modeller_type
        if self.modeller_type == ModellerType.SimpleFcfInterpolation:
            self.model_system = self.apply_simple_fcf_interpolation
        elif self.modeller_type == ModellerType.AnalyticApproximation:
            self.model_system = self.apply_analytic_approximation
            self.sample_across_stars = modeller_args[0]
            self.max_star = 957
            self.stellar_compositions = self.load_generic_float_data_csv('StellarCompositionsSortFE.csv')
            self.stars_to_sample = range(0, self.max_star + 1)
            self.default_star = list(self.stars_to_sample)
            #self.default_star = [478, 479]
            self.timescale_interpolator = ti.TimescaleInterpolator()
            self.geology_model = gi.GeologyModel()
        elif self.modeller_type == ModellerType.GridInterpolation:
            grid_filename = modeller_args[1]
            self.grid = SyntheticGrid(grid_filename)
            self.model_system = self.interpolate_on_synthetic_grid
        elif self.modeller_type == ModellerType.Null:
            self.model_system = self.null_model
        else:
            raise ValueError('Unrecognised modelling strategy ' + str(modeller_type))
        self.metallicity_elements = [ci.Element.Al, ci.Element.Ti, ci.Element.Ca, ci.Element.Mg]
        self.sinking_element_pairs_in_order_of_preference = [
            (ci.Element.Al, ci.Element.Ca),                    # We need a pair of elements that behave similarly in terms
            (ci.Element.Al, ci.Element.Ti),                    # of heating and partitioning as far as the model is concerned
            (ci.Element.Ca, ci.Element.Ti),
            (ci.Element.Ca, ci.Element.Mg),
            (ci.Element.Ti, ci.Element.Mg),
            (ci.Element.Al, ci.Element.Mg),
            (ci.Element.Ni, ci.Element.Fe),
            (ci.Element.Fe, ci.Element.Cr),
            (ci.Element.Ni, ci.Element.Cr)
        ]
        self.heating_element_pairs_in_order_of_preference = [  # The first element always needs to be the more volatile one
            (ci.Element.Na, ci.Element.Ca), #TODO: Should Ca actually be lower within the Na/O/Mg groups?
            (ci.Element.Na, ci.Element.Al),
            (ci.Element.Na, ci.Element.Ti),
            (ci.Element.Na, ci.Element.Mg),
            (ci.Element.O, ci.Element.Ca), # O is dubious here as it is non-monotonic - potentially move this bunch below the Mg bunch?
            (ci.Element.O, ci.Element.Al),
            (ci.Element.O, ci.Element.Ti),
            (ci.Element.O, ci.Element.Mg),
            (ci.Element.Mg, ci.Element.Ca),
            (ci.Element.Mg, ci.Element.Al),
            (ci.Element.Mg, ci.Element.Ti)
        ]
        self.partitioning_element_pairs_in_order_of_preference = [  # One of these should be lithophile, the other siderophile
            (ci.Element.Ca, ci.Element.Fe),
            (ci.Element.Mg, ci.Element.Fe),
            (ci.Element.Al, ci.Element.Fe),
            (ci.Element.Ti, ci.Element.Fe),
            (ci.Element.Ca, ci.Element.Ni), # Ni might actually be preferable to Fe as it is actually more siderophile
            (ci.Element.Mg, ci.Element.Ni),
            (ci.Element.Al, ci.Element.Ni),
            (ci.Element.Ti, ci.Element.Ni),
            (ci.Element.Ca, ci.Element.Cr),
            (ci.Element.Mg, ci.Element.Cr),
            (ci.Element.Al, ci.Element.Cr),
            (ci.Element.Ti, ci.Element.Cr)
        ]

    def model_populations(self, population_dict, overwrite=False):
        for pop_name, population in population_dict.items():
            print('Modelling population ' + pop_name)
            self.model_population(population, overwrite)
        #self.apply_kstest_across_pops(population_dict) This needs to be updated to handle list entries

    def model_population(self, population, overwrite=False):
        for system in population:
            if system.observed and system.observed_abundances is not None:
                if overwrite or system.modelled_properties is None:
                    print('Modelling system ' + str(system.id))
                    self.model_system(system)
                else:
                    print('System ' + str(system.id) + ' was already modelled, skipping')
        population.dump_to_csv()
        #self.apply_kstest_to_single_pop(population) This needs to be updated to handle list entries

    def apply_kstest_to_single_pop(self, population):
        for parameter in [mp.ModelParameter.fragment_core_frac, mp.ModelParameter.formation_distance, mp.ModelParameter.t_sinceaccretion]:
            input_distribution = population.input_values(parameter)
            output_distribution = population.modelled_values(parameter)
            bad_indices = list()
            other_bad_indices = list()
            for i, od in enumerate(output_distribution):
                if od is None or od in [np.nan]:
                    bad_indices.append(i)
            for bi in sorted(bad_indices, reverse=True):
                output_distribution.pop(bi) # Should we remove the corresponding numbers from input? Issue is that they won't necessarily be 1:1 at this point
            for i, od in enumerate(input_distribution):
                if od is None or od in [np.nan]:
                    other_bad_indices.append(i)
            for obi in sorted(other_bad_indices, reverse=True):
                input_distribution.pop(obi)
            if len(input_distribution) > 0 and len(output_distribution) > 0:
                result = st.ks_2samp(input_distribution, output_distribution)# if p < 0.05, samples are drawn from different distributions
            else:
                result = None
            population.io_ks_test_results[parameter] = result

    def apply_kstest_across_pops(self, population_dict):
        for pop_name, population in population_dict.items():
            for other_name, other_population in population_dict.items():
                for parameter in [mp.ModelParameter.fragment_core_frac, mp.ModelParameter.formation_distance, mp.ModelParameter.t_sinceaccretion]:
                    output_distribution = population.modelled_values(parameter)
                    other_output_distribution = other_population.modelled_values(parameter)
                    bad_indices = list()
                    other_bad_indices = list()
                    for i, od in enumerate(output_distribution):
                        if od is None or od in [np.nan]:
                            bad_indices.append(i)
                    for bi in sorted(bad_indices, reverse=True):
                        output_distribution.pop(bi)
                    for i, od in enumerate(other_output_distribution):
                        if od is None or od in [np.nan]:
                            other_bad_indices.append(i)
                    for obi in sorted(other_bad_indices, reverse=True):
                        other_output_distribution.pop(obi)
                    if len(output_distribution) > 0 and len(other_output_distribution) > 0:
                        result = st.ks_2samp(output_distribution, other_output_distribution) # if p < 0.05, samples are drawn from different distributions
                    else:
                        result = None
                    try:
                        population.pop_ks_test_results[other_name][parameter] = result
                    except KeyError:
                        population.pop_ks_test_results[other_name] = dict()
                        population.pop_ks_test_results[other_name][parameter] = result
                    try:
                        other_population.pop_ks_test_results[pop_name][parameter] = result
                    except KeyError:
                        other_population.pop_ks_test_results[pop_name] = dict()
                        other_population.pop_ks_test_results[pop_name][parameter] = result

    def load_generic_float_data_csv(self, input_filename):
        # TODO: Change this to a with clause to prevent ResourceWarnings
        with open(pu.get_path_to_data() + input_filename, encoding='utf-8') as generic_csv:
            generic_list = [row for row in csv.reader(generic_csv)]
            generic_array = np.asarray(generic_list)
        return generic_array.astype(np.float)

    def apply_simple_fcf_interpolation(self, system):
        MgHx = system.observed_abundances.get(ci.Element.Mg, None)
        FeHx = system.observed_abundances.get(ci.Element.Fe, None)
        if MgHx is None or FeHx is None:
            print('Warning! Could not calculate fcf for system ' + str(system))
        else:
            MgFe = MgHx - FeHx
            cnf = ((0 - MgFe) + 1) /3
            if cnf < 0:
                cnf = 0
            if cnf > 1:
                cnf = 1
            output_dict = {
                mp.ModelParameter.fragment_core_frac: [cnf]
            }
            system.set_modelled_properties(output_dict)

    def null_model(self, system):
        for model_parameter in mp.ModelParameter:
            if system.modelled_properties is None:
                system.modelled_properties = {model_parameter: [None]}
            else:
                system.modelled_properties[model_parameter] = system.modelled_properties.get(model_parameter, [None])

    def apply_analytic_approximation(self, system):
        HorHe = None
        if system.wd_properties[mp.WDParameter.spectral_type] == 'DB':
            HorHe = 'He'
        elif system.wd_properties[mp.WDParameter.spectral_type] == 'DA':
            HorHe = 'H'
        else:
            pass
        sinking_timescales = self.timescale_interpolator.get_wd_timescales(
            HorHe,
            system.wd_properties[mp.WDParameter.logg],
            system.wd_properties[mp.WDParameter.temperature],
            system.observed_abundances.get(ci.Element.Ca, -14)  # If Ca not present, we'll assume it's -14: this is the default case, dropoff starts around -12
        ) # These should be in years

        metallicity = self.estimate_metallicity(system, sinking_timescales)
        t_disc, t_sinceaccretion = self.estimate_sinking(system, sinking_timescales, metallicity)
        #if len(metallicity) > 1:
        #    t_disc = t_disc*len(metallicity)
        #    t_sinceaccretion = t_sinceaccretion*len(metallicity)
        d_formation = self.estimate_heating(system, sinking_timescales, t_disc, t_sinceaccretion, metallicity)
        fcf = self.estimate_fragment_core_fraction(system, sinking_timescales, t_disc, t_sinceaccretion, d_formation, metallicity)

        proposed_d_formation = [np.log10(d) if d is not None else None for d in d_formation]

        output_dict = {
            mp.ModelParameter.metallicity: metallicity,
            mp.ModelParameter.formation_distance: self.collapse_list_of_repeats(proposed_d_formation),
            mp.ModelParameter.fragment_core_frac: self.collapse_list_of_repeats(fcf),
            mp.ModelParameter.accretion_timescale: t_disc, # It's important not to collapse these
            mp.ModelParameter.t_sinceaccretion: t_sinceaccretion
        }

        system.set_modelled_properties(output_dict)

    def collapse_list_of_repeats(self, input_list):
        if len(input_list) <= 1:
            return input_list
        if all(val == input_list[0] for val in input_list):
            return [input_list[0]]
        else:
            return input_list

    def get_stellar_abundance(self, star, element):
        stellar_comp = self.stellar_compositions[star]
        element_index = ci.usual_elements.index(element)
        if element_index == 6:
            element_index = len(ci.usual_elements)
        elif element_index > 6:
            element_index -= 1
        try:
            initial_el = stellar_comp[element_index]
        except IndexError:
            initial_el = 1
        return initial_el

    def get_elements_for_sinking_estimate(self, system):
        for pair in self.sinking_element_pairs_in_order_of_preference:
            if system.observed_abundances.get(pair[0], None) is not None and system.observed_abundances.get(pair[1], None) is not None:
                return pair[0], pair[1]
        print('Warning! No suitable pair of elements was found for estimating sinking')
        return None, None

    def get_elements_for_heating_estimate(self, system):
        for pair in self.heating_element_pairs_in_order_of_preference:
            if system.observed_abundances.get(pair[0], None) is not None and system.observed_abundances.get(pair[1], None) is not None:
                return pair[0], pair[1]
        print('Warning! No suitable pair of elements was found for estimating heating')
        return None, None

    def get_elements_for_fcf_estimate(self, system):
        for pair in self.partitioning_element_pairs_in_order_of_preference:
            if system.observed_abundances.get(pair[0], None) is not None and system.observed_abundances.get(pair[1], None) is not None:
                return pair[0], pair[1]
        print('Warning! No suitable pair of elements was found for estimating fcf')
        return None, None

    def estimate_sinking(self, system, sinking_timescales, metallicity=[None]):
        if system.wd_properties[mp.WDParameter.spectral_type] == 'DA':
            # Then it really doesn't matter what the exact numbers are (and it's kind of meaningless anyway) - just put it in steady state
            t_Mg = sinking_timescales[ci.Element.Mg]
            t_disc = [20*t_Mg]*len(metallicity) # Multiply by len(metallicity) to guarantee that these 3 lists all have the same length
            t_sinceaccretion = [10*t_Mg/1000000]*len(metallicity)  # 20 and 10 are arbitrary, but will put this system in steady state. just as long as 5t_Mg (ish) < t_sinceaccretion < t_disc
        else:
            element1, element2 = self.get_elements_for_sinking_estimate(system)
            XHx1 = system.observed_abundances.get(element1, None)
            XHx2 = system.observed_abundances.get(element2, None)
            if XHx1 is None or XHx2 is None:
                return [None], [None]
            t_disc, t_sinceaccretion = self.calculate_sinking_from_element_pair(element1, element2, 10**(XHx1 - XHx2), sinking_timescales[element1], sinking_timescales[element2])
        return t_disc, t_sinceaccretion

    def calculate_sinking_from_element_pair(self, element1, element2, element_ratio, t_1, t_2):
        # element_ratio is in linear space
        # Assume initial composition is Earth like (according to geology_info)
        #initial_Ca = 0.011174350504526264
        #initial_Fe = 0.15006893982819816
        ratio_increases_with_time = t_1 > t_2

        test_t_disc = 3*max(t_1, t_2) # This is a bit experimental - ideally we want this to be another parameter we search for
        wd_timescales = [t_1, t_2]
        #non_zero_wd_timescales = [t_1, t_2]
        pollutionfraction = -6 # Exact value irrelevant: we only care about ratio not absolute quantity
        tolerance = 0.000001
        max_iterations = 100
        zero_time_threshold = 0.00001

        t_disc_toret = list()
        t_sinceaccretion_toret = list()
        stars_to_sample = self.default_star
        if self.sample_across_stars:
            stars_to_sample = self.stars_to_sample
        # Assume initial composition is Solar like (according to geology_info) (These are not normalised but doesn't matter)
        # NB: Using geology_info.solar_abundances seemed to underestimate Al/Ca relative to other stars?!? Why?
        # Also (V IMPORTANT) Al/Ca starts increasing (due to 'heating') at distances less than about 0.25AU (-0.6 in log units)
        # This represents a fundamental limit to this method: if distance < -0.6 or so, the model will think we're in
        # declining phase (to try to match the elevated Al/Ca) even if we're not, with various knock-on effects eg overestimating fcf
        #initial_XHx1 = self.geology_model.solar_ratiod_to_H[element1]
        #initial_XHx2 = self.geology_model.solar_ratiod_to_H[element2]
        #initial_CaHx = -1.22077587
        #initial_FeHx = -0.079526231

        # Update: Now assume the composition corresponds to one of the Brewer stars
        for star in stars_to_sample:
            initial_el1 = self.get_stellar_abundance(star, element1)
            initial_el2 = self.get_stellar_abundance(star, element2)
            initial_ratio = initial_el1/initial_el2
            #initial_ratio = 10**(initial_XHx1 - initial_XHx2)
            planetesimal_abundance = [initial_ratio, 1]  # This gets normalised later
            #non_zero_planetesimal_abundance = [initial_Ca, initial_Fe]
            lower_t_sinceaccretion_bound = 0 # Minimum possible value (in Myr)
            upper_t_sinceaccretion_bound = 200*(max(t_1, t_2)/1000000) # Maximum possible value (in Myr) Setting it equal to N times t_Ca
            iteration_count = 0
            while True:
                # Do a simple binary search
                test_t_sinceaccretion = 0.5*(lower_t_sinceaccretion_bound + upper_t_sinceaccretion_bound)
                #test_abundances = wdm.process_abundances(test_t_sinceaccretion, test_t_disc, planetesimal_abundance, non_zero_planetesimal_abundance, wd_timescales, non_zero_wd_timescales, pollutionfraction, True)
                test_abundances = wdm.process_abundances(test_t_sinceaccretion, test_t_disc, planetesimal_abundance, wd_timescales, pollutionfraction, True)
                test_ratio = 10**(test_abundances[0] - test_abundances[1])
                rel_diff = abs((test_ratio - element_ratio)/element_ratio)
                if rel_diff < tolerance:
                    print('Sinking calculation converged (' + str(iteration_count + 1) + ' iterations)')
                    t_disc_toret.append(test_t_disc)
                    t_sinceaccretion_toret.append(test_t_sinceaccretion)
                    break
                else:
                    if test_ratio > element_ratio:
                        # We guessed too high
                        if ratio_increases_with_time:
                            # Go earlier
                            upper_t_sinceaccretion_bound = test_t_sinceaccretion
                        else:
                            # Go later
                            lower_t_sinceaccretion_bound = test_t_sinceaccretion
                    else:
                        # We guessed too low
                        if ratio_increases_with_time:
                            # Go later
                            lower_t_sinceaccretion_bound = test_t_sinceaccretion
                        else:
                            # Go earlier
                            upper_t_sinceaccretion_bound = test_t_sinceaccretion
                iteration_count += 1
                if iteration_count > max_iterations:
                    # To prevent infinite loops
                    if test_t_sinceaccretion < zero_time_threshold:
                        print('Warning! Sinking calculation exceeded iteration limit (' + str(max_iterations) + '), but test t = ' + str(test_t_sinceaccretion) + ' so will approximate as 0')
                        t_disc_toret.append(test_t_disc)
                        t_sinceaccretion_toret.append(0)
                        break
                    else:
                        print('Warning! Sinking calculation exceeded iteration limit (' + str(max_iterations) + '), returning None, None')
                        t_disc_toret.append(None)
                        t_sinceaccretion_toret.append(None)
                        break
        return t_disc_toret, t_sinceaccretion_toret

    def estimate_metallicity(self, system, sinking_timescales):
        if system.wd_properties[mp.WDParameter.spectral_type] == 'DB':
            return [None]
        #Plan: put system into steady state, compare to all 958 stars for the best match (for elements less susceptible to heating, partitioning differences etc)
        elements_to_use = self.metallicity_elements
        steady_state_adjusted_abundances = cn.OrderedDict()
        for element in elements_to_use:
            uncorrected_abundance = system.observed_abundances.get(element, None)
            if uncorrected_abundance is not None:
                scaled_abundance = (10**uncorrected_abundance)/sinking_timescales[element]
                steady_state_adjusted_abundances[element] = scaled_abundance
        if len(steady_state_adjusted_abundances) < 2:
            return [None]
        comparison_results = dict()
        first_entry = next(iter(steady_state_adjusted_abundances.items()))
        base_el = first_entry[0]
        base_el_ss_abundance = first_entry[1]
        for star in self.stars_to_sample:
            stellar_abundances = [self.get_stellar_abundance(star, element) for element in steady_state_adjusted_abundances]  # TODO check this returns 1 for Mg
            # Now compare steady_state_adjusted_abundances to stellar_abundances
            # We'll call the first element el_base, then compare ratios of all other elements relative to el_base
            # Use L1 normalisation for now rather than L2 - don't want to be overly swayed by outliers
            el_index = 0
            total_penalty = 0
            for element, ss_abundance in steady_state_adjusted_abundances.items():
                if element == base_el:
                    el_index += 1
                    continue
                ss_ratio = ss_abundance/base_el_ss_abundance
                stellar_ratio = stellar_abundances[el_index]/stellar_abundances[0]
                penalty = ss_ratio/stellar_ratio if ss_ratio > stellar_ratio else stellar_ratio/ss_ratio
                penalty -= 1 # So that identical results give 0 penalty
                total_penalty += penalty
                el_index += 1
            comparison_results[star] = total_penalty
        potential_best_matches = list()
        current_best_score = np.inf
        for star, penalty in comparison_results.items():
            if penalty == current_best_score:
                potential_best_matches.append(star)
            if penalty < current_best_score:
                potential_best_matches = [star]
                current_best_score = penalty
        return potential_best_matches

    def estimate_fragment_core_fraction(self, system, sinking_timescales, t_disc, t_sinceaccretion, d_formation, metallicity=[None]):
        stars_to_sample = metallicity if not all(m is None for m in metallicity) else self.default_star
        if self.sample_across_stars:
            stars_to_sample = self.stars_to_sample
        element1, element2 = self.get_elements_for_fcf_estimate(system)
        XHx1 = system.observed_abundances.get(element1, None)
        XHx2 = system.observed_abundances.get(element2, None)
        if XHx1 is None or XHx2 is None:
            return [None]*len(stars_to_sample)
        observed_element_ratio = 10**(XHx1 - XHx2)
        fcf_toret = list()
        t_sinceaccretionyears = [t*1000000 if t is not None else None for t in t_sinceaccretion]
        if len(t_sinceaccretionyears) == 1:
            t_sinceaccretionyears = len(stars_to_sample)*t_sinceaccretionyears
        if len(t_disc) == 1:
            t_disc = len(stars_to_sample)*t_disc
        for star_count, star in enumerate(stars_to_sample):
            if t_sinceaccretionyears[star_count] is None or t_disc[star_count] is None:
                fcf_toret.append(None)
            else:
                scaling_el1 = wdm.calculate_buildup_scaling_factors(t_sinceaccretionyears[star_count], t_disc[star_count], sinking_timescales[element1])*wdm.calculate_sinkout_scaling_factors(t_sinceaccretionyears[star_count], t_disc[star_count], sinking_timescales[element1])
                scaling_el2 = wdm.calculate_buildup_scaling_factors(t_sinceaccretionyears[star_count], t_disc[star_count], sinking_timescales[element2])*wdm.calculate_sinkout_scaling_factors(t_sinceaccretionyears[star_count], t_disc[star_count], sinking_timescales[element2])

                uncompensated_fragment_element_ratio = observed_element_ratio * (scaling_el2/scaling_el1)

                # Now need to factor in heating: based on the d_formation we calculated, we skewed this ratio away from what it would be based
                # on the fcf distribution we put in at the start, need to compensate for this
                z_formation = 0.05
                t_formation = 1.5
                if d_formation[star_count] in [None, np.nan]:
                    # then we can't do anything
                    heating_ratio = 1
                else:
                    heated_abundances = am.get_all_abundances([element1, element2], d_formation[star_count], z_formation, t_formation)
                    heating_ratio = heated_abundances[element1]/heated_abundances[element2] # if this comes out as el1 > el2 then our ufer was high largely because of heating, and should be lowered (and vice versa)
                fragment_element_ratio = uncompensated_fragment_element_ratio/heating_ratio
                fcf = self.calculate_fcf_from_element_pair(element1, element2, fragment_element_ratio)
                fcf_toret.append(fcf)
        return fcf_toret

    def calculate_fcf_from_element_pair(self, element1, element2, fragment_element_ratio):
        # Input to this function should be in linear space, not log
        # Assume Earth-like composition, taken from geology_info
        el1_mantle = self.geology_model.element_info[element1][gi.Layer.mantle] #0.008885525909622828#
        el1_core = self.geology_model.element_info[element1][gi.Layer.core] #0.0#
        el2_mantle = self.geology_model.element_info[element2][gi.Layer.mantle] #0.030077124487241223#
        el2_core = self.geology_model.element_info[element2][gi.Layer.core] #0.9031235537406759#
        useful_term = el1_mantle - (fragment_element_ratio*el2_mantle)
        other_term = (el2_core*fragment_element_ratio) - el1_core
        fcf = useful_term / (useful_term + other_term)
        if fcf < 0 or fcf > 1:
            print('Warning! Unphysical fcf inferred')
            return np.nan
        return fcf

    def estimate_heating(self, system, sinking_timescales, t_disc, t_sinceaccretion, metallicity=[None]):
        element1, element2 = self.get_elements_for_heating_estimate(system)
        XHx1 = system.observed_abundances.get(element1, None)
        XHx2 = system.observed_abundances.get(element2, None)
        if XHx1 is None or XHx2 is None:
            target_ratio = None
        else:
            target_ratio = 10**(XHx1 - XHx2)
        d_formation = self.find_d_formation(element1, element2, target_ratio, t_disc, t_sinceaccretion, sinking_timescales.get(element1), sinking_timescales.get(element2), metallicity)
        return d_formation

    def find_d_formation(self, element1, element2, target_ratio, t_disc, t_sinceaccretion, t_1, t_2, metallicity=[None]):
        # Assume initial composition is Solar like (according to geology_info) (These are not normalised but doesn't matter)
        #initial_CaHx = -1.22077587
        #initial_NaHx = -1.3596521970000002
        #initial_XHx1 = self.geology_model.solar_abundances[element1]
        #initial_XHx2 = self.geology_model.solar_abundances[element2]
        #initial_XHx1 = self.geology_model.solar_ratiod_to_H[element1]
        #initial_XHx2 = self.geology_model.solar_ratiod_to_H[element2]
        #initial_ratio = 10**(initial_XHx1 - initial_XHx2)

        #initial_ratio = 0.048984876/0.058807196
        d_formation_toret = list()
        stars_to_sample = metallicity if not all(m is None for m in metallicity) else self.default_star
        if self.sample_across_stars:
            stars_to_sample = self.stars_to_sample
        z_formation = 0.05
        t_formation = 1.5
        tolerance = 0.000001
        max_iterations = 50
        pollutionfraction = -6 # This doesn't matter since we only care about ratios
        if len(t_sinceaccretion) == 1:
            t_sinceaccretion = len(stars_to_sample)*t_sinceaccretion
        if len(t_disc) == 1:
            t_disc = len(stars_to_sample)*t_disc
        wd_timescales = np.array([t_1, t_2])
        for star_count, star in enumerate(stars_to_sample):
            if target_ratio is None or t_sinceaccretion[star_count] is None or t_disc[star_count] is None or element1 is None or element2 is None:
                d_formation_toret.append(None)
            else:
                initial_el1 = self.get_stellar_abundance(star, element1)
                initial_el2 = self.get_stellar_abundance(star, element2)
                initial_ratio = initial_el1/initial_el2
                lower_d_formation_bound = 0 # Minimum possible value (in AU)
                upper_d_formation_bound = 3 # Maximum possible value (in AU)
                max_distance_threshold = upper_d_formation_bound - 0.0001
                iteration_count = 0
                while True:
                    # Do a simple binary search
                    test_d_formation = 0.5*(lower_d_formation_bound + upper_d_formation_bound)
                    test_abundances = am.get_all_abundances(ci.usual_elements, test_d_formation, z_formation, t_formation, star)
                    presunk_ratio = initial_ratio*(test_abundances[element1]/test_abundances[element2])
                    postsinking_abundances = wdm.process_abundances(t_sinceaccretion[star_count], t_disc[star_count], [presunk_ratio, 1], wd_timescales, pollutionfraction, True)
                    test_ratio = 10**(postsinking_abundances[0] - postsinking_abundances[1])
                    rel_diff = abs((test_ratio - target_ratio)/target_ratio)
                    if rel_diff < tolerance:
                        print('D formation calculation converged (' + str(iteration_count + 1) + ' iterations)')
                        d_formation_toret.append(test_d_formation)
                        break
                    else:
                        # This logic relies on the assumption that element1/element2 increases with d - so element 1 needs to be the more volatile
                        if test_ratio > target_ratio:
                            # We guessed too high
                            upper_d_formation_bound = test_d_formation
                        else:
                            # We guessed too low
                            lower_d_formation_bound = test_d_formation
                    iteration_count += 1
                    if iteration_count > max_iterations:
                        # To prevent infinite loops
                        if test_d_formation > max_distance_threshold:
                            print('Warning! D formation calculation exceeded iteration limit (' + str(max_iterations) + '), but d_formation was close to max, returning ' + str(upper_d_formation_bound))
                            d_formation_toret.append(upper_d_formation_bound)
                            break
                        else:
                            print('Warning! D formation calculation exceeded iteration limit (' + str(max_iterations) + '), returning None')
                            d_formation_toret.append(None)
                            break
        return d_formation_toret

    def interpolate_on_synthetic_grid(self, system):
        modelled_properties = self.grid.interpolate_modelled_properties(system)
        system.set_modelled_properties(modelled_properties)

def plot_modelling_functions():
    modeller = Modeller(ModellerType.AnalyticApproximation, [False, 'NA'])
    test_CaFe_list = np.linspace(0, 0.5, 1000)
    t_Ca = 2000000  # These are roughly representative (at least in terms of their ratio)
    t_Fe = 1500000
    output_vals_dict = {
        mp.ModelParameter.accretion_timescale: list(),
        mp.ModelParameter.t_sinceaccretion: list(),
        mp.ModelParameter.fragment_core_frac: list()
    }
    for CaFe in test_CaFe_list:
        t_disc, t_sinceaccretion = modeller.calculate_sinking_from_element_pair(ci.Element.Ca, ci.Element.Fe, CaFe, t_Ca, t_Fe)
        fcf = modeller.calculate_fcf_from_element_pair(ci.Element.Ca, ci.Element.Fe, CaFe)
        output_vals_dict[mp.ModelParameter.accretion_timescale].append(t_disc)
        output_vals_dict[mp.ModelParameter.t_sinceaccretion].append(t_sinceaccretion)
        output_vals_dict[mp.ModelParameter.fragment_core_frac].append(fcf)
    print(test_CaFe_list)
    print(output_vals_dict)
    graph_fac = gf.GraphFactory()
    graph_fac.plot_modeller_functions(test_CaFe_list, ci.Element.Ca, ci.Element.Fe, output_vals_dict, mp.ModelParameter.fragment_core_frac, mp.ModelParameter.t_sinceaccretion)

def main():
    plot_modelling_functions()

if __name__ == '__main__':
    main()
