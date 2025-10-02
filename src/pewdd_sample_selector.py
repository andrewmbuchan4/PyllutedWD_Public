#!/usr/bin/env python

from argparse import Namespace
from itertools import combinations

import csv
import numpy as np
import scipy

import chemistry_info as ci
import model_parameters as mp
import pewdd_interface as pi
import pwd_utils as pu
import solar_abundances as sa
import timescale_interpolator as ti

max_num_observable_sinking_timescales = 5
assumed_accretion_event_timescale = 100000 #yr
dummy_timescale = 0.01 #yr - doesn't matter much, just needs to be short relative to assumed_accretion_event_timescale
min_thermohaline_temp = 10000 # The minimum temperature that we will consider for thermohaline DAs
max_thermohaline_temp = 20502
synthetic_default_error = 0.2 # The error to be assumed for synthetic wds
min_elements_for_correlation = 4

simplified_D_estimation = False

num_els_string = 'num_modellable_elements'
timescale_pairs_string = 'timescale_pairs'
score_threshold_string = 'score_threshold'
thermohaline_dummy_override_string = 'THERMOHALINE_DUMMY'

sample_selection_settings = {
    'BVK_DA': {
        num_els_string: 2,
        score_threshold_string: 1,
        mp.WDParameter.atmospheric_type: ci.Element.H,
        timescale_pairs_string: [
            (ti.TimescaleType.KoesterNoOvershoot, ti.TimescaleType.BedardNoOvershoot)
        ]
    },
    'BVK_DB': {
        num_els_string: 4,  # Just to keep it manageable
        score_threshold_string: 2,
        mp.WDParameter.atmospheric_type: ci.Element.He,
        timescale_pairs_string: [
            (ti.TimescaleType.KoesterNoOvershoot, ti.TimescaleType.BedardNoOvershoot)
        ]
    },
    'OVERSHOOT_DA': {
        num_els_string: 2,
        score_threshold_string: 1,
        mp.WDParameter.atmospheric_type: ci.Element.H,
        timescale_pairs_string: [
            (ti.TimescaleType.Bedard3DOvershoot, ti.TimescaleType.BedardNoOvershoot)
        ]
    },
    'OVERSHOOT_DA_PATCHED': {
        num_els_string: 2,
        score_threshold_string: 1,
        mp.WDParameter.atmospheric_type: ci.Element.H,
        timescale_pairs_string: [
            (ti.TimescaleType.Bedard3DOvershootPatched, ti.TimescaleType.BedardNoOvershoot)
        ]
    },
    'OVERSHOOT_DA_CONTROL': {
        num_els_string: 2,
        score_threshold_string: -0.5,
        mp.WDParameter.atmospheric_type: ci.Element.H,
        timescale_pairs_string: [
            (ti.TimescaleType.Bedard3DOvershoot, ti.TimescaleType.BedardNoOvershoot)
        ]
    },
    'OVERSHOOT_DB': {
        num_els_string: 2,
        score_threshold_string: 1,
        mp.WDParameter.atmospheric_type: ci.Element.He,
        timescale_pairs_string: [
            (ti.TimescaleType.BedardVariableOvershoot, ti.TimescaleType.BedardNoOvershoot)
        ]
    },
    'FIXED_OVERSHOOT_DB': {
        num_els_string: 2,
        score_threshold_string: 1,
        mp.WDParameter.atmospheric_type: ci.Element.He,
        timescale_pairs_string: [
            (ti.TimescaleType.BedardVariableOvershoot, ti.TimescaleType.BedardOvershoot)
        ]
    },
    'THERMOHALINE_DA': {
        num_els_string: 2,
        score_threshold_string: 1,
        mp.WDParameter.atmospheric_type: ci.Element.H,
        timescale_pairs_string: [
            (ti.TimescaleType.Bedard3DOvershootPatched, thermohaline_dummy_override_string)
        ]
    }
}

cached_timescale_interpolator = ti.TimescaleInterpolator()

def get_cached_timescale_interpolator():
    return cached_timescale_interpolator

def find_shared_grid_range(atmospheric_type, pairs_of_timescale_types):
    min_teff = -np.inf
    max_teff = np.inf
    min_logg = -np.inf
    max_logg = np.inf
    for pair in pairs_of_timescale_types:
        members = [pair[0], pair[1]]
        for member in members:
            if member == thermohaline_dummy_override_string:
                min_teff = max(min_teff, min_thermohaline_temp)
                max_teff = min(max_teff, max_thermohaline_temp)
            else:
                teff_range = ti.teff_ranges[atmospheric_type][member]
                logg_range = ti.logg_ranges[atmospheric_type][member]
                if teff_range[0] > min_teff:
                    min_teff = teff_range[0]
                if teff_range[1] < max_teff:
                    max_teff = teff_range[1]
                if logg_range[0] > min_logg:
                    min_logg = logg_range[0]
                if logg_range[1] < max_logg:
                    max_logg = logg_range[1]
    return min_teff, max_teff, min_logg, max_logg

def flatten_timescale_pairs(timescale_pairs):
    toret = list()
    for combo in timescale_pairs:
        for index in [0, 1]:
            if combo[index] not in toret:
                toret.append(combo[index])
    return toret

def calculate_discrepancy_metric(wd, timescale_pairs, timescale_override_dict=None, D_override=None):
    # This estimates the number of sigma discrepancy due to difference in timescale values relative to the abundance uncertainties
    # So, broadly speaking, if the max_discrepency is above 1, that's a sign that we should maybe care
    starting_maximum_discrepancy = 0
    max_discrepency = starting_maximum_discrepancy
    try:
        detected_elements = wd.get_elements_present()
    except AttributeError:
        # For Synthetic white dwarfs
        # Really the solution to this problem is to make real and synthetic white dwarfs inherit from the same class and share basic functions like this
        if wd.observed_abundances is not None:
            detected_elements = [el for el, val in wd.observed_abundances.items() if val not in [None, np.nan, -np.inf]]
        else:
            detected_elements = list()
    elements_to_compare = [el for el in ci.usual_elements if el in detected_elements]
    most_discrepent_el1 = None
    most_discrepent_el2 = None
    most_discrepent_timescale_pair = None
    flattened_timescale_pairs = flatten_timescale_pairs(timescale_pairs)
    if D_override is None:
        D = estimate_declining_phase_depth(wd, flattened_timescale_pairs)
    else:
        D = D_override
    #print('D = ' + str(D))
    default_error = 0.2 # If there's no error on an element, we'd probably replace it with this when running, so do the same here
    for el1, el2 in combinations(elements_to_compare, 2):
        for combo in timescale_pairs:
            if timescale_override_dict is None:
                try:
                    #timescales1 = wd.get_timescale_values_dict(combo[0])
                    timescales1 = get_timescale_values_from_wd_or_dummy(wd, combo[0])
                except KeyError:
                    timescales1 = None
                try:
                    #timescales2 = wd.get_timescale_values_dict(combo[1])
                    timescales2 = get_timescale_values_from_wd_or_dummy(wd, combo[1])
                except KeyError:
                    timescales2 = None
            else:
                if combo[0] == thermohaline_dummy_override_string:
                    timescales1 = get_timescale_values_from_wd_or_dummy(wd, combo[0])
                else:
                    timescales1 = timescale_override_dict.get(combo[0])
                if combo[1] == thermohaline_dummy_override_string:
                    timescales2 = get_timescale_values_from_wd_or_dummy(wd, combo[1])
                else:
                    timescales2 = timescale_override_dict.get(combo[1])
            if timescales1 is not None and timescales2 is not None:
                try:
                    ratio1 = timescales1[el1]/timescales1[el2]
                    ratio2 = timescales2[el1]/timescales2[el2]
                except KeyError:
                    ratio1 = None
                    ratio2 = None
                if ratio1 is not None and ratio2 is not None:
                    try:
                        sigma_X = wd.get_abundance(el1).upper_error
                    except AttributeError:
                        sigma_X = synthetic_default_error
                    try:
                        sigma_Y = wd.get_abundance(el2).upper_error
                    except AttributeError:
                        sigma_Y = synthetic_default_error
                    if sigma_X <= 0.0:
                        sigma_X = default_error
                    if sigma_Y <= 0.0:
                        sigma_Y = default_error
                    metric = metric_function_for_given_element_pair_X_Y(ratio1, ratio2, D, sigma_X, sigma_Y)
                    if metric > max_discrepency:
                        max_discrepency = metric
                        most_discrepent_el1 = el1
                        most_discrepent_el2 = el2
                        most_discrepent_timescale_pair = combo
    return max_discrepency

def metric_function_for_given_element_pair_X_Y(ratio1, ratio2, D, sigma_X, sigma_Y):
    # This can be derived by considering the difference in log abundance (relative to errors) caused by swapping to a new set of timescales
    # And doing a neat substitution (t = D*tau) to remove the awkward difference-in-reciprocals term
    # And also making a few assumptions eg that the system reaches SS, then hits declining phase (such that the discrepencies just add)
    #ratio1 = tau_X_1/tau_Y_1
    #ratio2 = tau_X_2/tau_Y_2
    SS_term = np.log10(ratio1/ratio2) # Discrepancy due to SS phase
    declining_term = (D/np.log(10))*(ratio1 - ratio2)
    error_term = np.sqrt(sigma_X**2 + sigma_Y**2) # This has some implicit assumptions which I don't think we need to worry about...
    numerator = abs(SS_term + declining_term)
    toret = numerator/error_term
    return toret

#def rank_list(list_to_rank):
#    ranks = scipy.stats.rankdata(list_to_rank, method='max')
#    max_rank = max(ranks)
#    toret = [max_rank - (r - 1) for r in ranks]
#    return toret

#def rank_candidates_by_timescale_discrepancy(wd_candidates, timescale_pairs):
#    effects = list()
#    for wd in wd_candidates:
#        effect = estimate_effect_of_timescale_discrepancy(wd, timescale_pairs)
#        effects.append(effect)
#    ranks_toret = rank_list(effects)
#    return effects, ranks_toret

def estimate_declining_phase_depth(white_dwarf, timescale_types_to_try):
    # This is super hand wavy, but has some logic behind it! We need to take two things into account:
    corr = estimate_correlation_between_elements_and_timescales(white_dwarf, timescale_types_to_try)
    chance = estimate_a_priori_chance_of_declining_phase(white_dwarf, timescale_types_to_try)
    # How can we estimate D? Firstly, it should scale with a_priori_probability in some way, in order to recover D ~ 0 for systems with short sinking timescales (and low a_priori_probabilities)
    # It should also increase with the strength of correlation
    # This suggests something like D ~ a*c where a is the a_priori_probability and c is the correlation
    # Caveat: c can be negative, and in this case it seems (intuitively) that we shouldn't decrease the metric (arbitrary choice!)
    # Also, the maximum value of this product is 1, when actually D can go up to max_num_observable_sinking_timescales
    # So we arrive at this:
    if simplified_D_estimation:
        estimate = chance
    else:
        estimate = chance*max(0, corr)*max_num_observable_sinking_timescales
    print(white_dwarf.full_name())
    print('corr = ' + str(corr))
    print('chance = ' + str(chance))
    print('estimate = ' + str(estimate))
    print()
    return estimate

def estimate_a_priori_chance_of_declining_phase(white_dwarf, timescale_types_to_try):
    # The idea here is that the longer WD sinking timescales are, the greater the chance that they'll be in the declining phase
    reference_element = ci.Element.Mg # Arbitrary, but also unimportant
    all_chances = list()
    for timescale_type in timescale_types_to_try:
        timescales = get_timescale_values_from_wd_or_dummy(white_dwarf, timescale_type, [reference_element])
        reference_sinking_timescale = timescales[reference_element]
        observable_declining_phase_time = max_num_observable_sinking_timescales*reference_sinking_timescale
        max_observable_time = assumed_accretion_event_timescale + observable_declining_phase_time
        chance = observable_declining_phase_time/max_observable_time
        all_chances.append(chance)
    toret = np.mean(all_chances)
    return toret

def estimate_correlation_between_elements_and_timescales(white_dwarf, timescale_types_to_try):
    elements = ci.usual_elements
    all_corrs = list()
    for timescale_type in timescale_types_to_try:
        timescales = get_timescale_values_from_wd_or_dummy(white_dwarf, timescale_type, elements)
        # Makes sense to work in terms of log abundances since the log abundances change linearly in time (in declining phase)
        try:
            abundances = white_dwarf.get_abundance_values_dict(elements) # This includes upper bounds which must be removed later
        except AttributeError:
            abundances = dict() if white_dwarf.observed_abundances is None else white_dwarf.observed_abundances
        try:
            reference_element = list(abundances)[0] # It doesn't matter what the element is, as long as it's actually in this dict: correlation should be invariant to an offset (which is the effect of changing the reference element)
        except IndexError:
            return 0 # This means we don't have any abundances to work with!
        scaled_abundances = sa.scale_abundances_to_solar(abundances, reference_element)
        try:
            detected_elements = white_dwarf.get_elements_present()
        except AttributeError:
            # For Synthetic white dwarfs
            # Really the solution to this problem is to make real and synthetic white dwarfs inherit from the same class and share basic functions like this
            detected_elements = [el for el, val in white_dwarf.observed_abundances.items() if val not in [None, np.nan, -np.inf]]
        elements_to_compare = [el for el in elements if el in detected_elements]
        if len(elements_to_compare) < min_elements_for_correlation:
            return 0 # We don't care about these at all - any correlation here is pretty close to meaningless. We want to essentially use this as a 'bonus' to pick out a few extra system that really seem like they must be in the declining phase
        timescales_to_compare = [timescales[el] for el in elements_to_compare]
        abundances_to_compare = [scaled_abundances[el] for el in elements_to_compare]
        try:
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(timescales_to_compare, abundances_to_compare)
            all_corrs.append(r_value)
        except ValueError:
            # Occurs for the dummy timescales: they're all the same, so no correlation can be calculated
            pass
    toret = np.mean(all_corrs)
    return toret

def get_timescale_values_from_wd_or_dummy(white_dwarf, timescale_type, elements=ci.usual_elements):
    if timescale_type == thermohaline_dummy_override_string:
        toret = {el: dummy_timescale for el in elements}
    else:
        try:
            toret = white_dwarf.get_timescale_values_dict(timescale_type, elements)
        except AttributeError:
            timescale_interpolator = get_cached_timescale_interpolator()
            atm_type = ci.Element.He if white_dwarf.wd_properties[mp.WDParameter.spectral_type] == 'DB' else ci.Element.H
            try:
                ca_to_use = white_dwarf.observed_abundances.get(ci.Element.Ca, -15)
            except AttributeError:
                ca_to_use = -15
            toret = timescale_interpolator.get_wd_timescales(
                atm_type,
                white_dwarf.wd_properties[mp.WDParameter.logg],
                white_dwarf.wd_properties[mp.WDParameter.temperature],
                ca_to_use
            )[timescale_type]
    return toret

def compile_sample(wd_candidates, score_threshold, timescale_pairs):
    scores = list()
    indices_for_sample = list()
    for i, wd in enumerate(wd_candidates):
        #discrepancy_effect = discrepancy_effects[i]
        #declining_effect = declining_effects[i]
        #a_priori_probability = a_priori_probabilities[i]
        #final_score = overall_metric(discrepancy_effect, declining_effect, a_priori_probability)
        score = calculate_discrepancy_metric(wd, timescale_pairs)
        scores.append(score)
        wd.pewdd_metric = score

    #final_ranks = rank_list(scores)
    #for j, rank in enumerate(final_ranks):
    #    if rank <= sample_size:
    #        indices_for_sample.append(j)
    if score_threshold > 0:
        for j, s in enumerate(scores):
            if s >= score_threshold: # Broadly speaking, if this is above 1, then the discrepency due to timescales is greater than the error
                indices_for_sample.append(j)
    else:
        for j, s in enumerate(scores):
            if 0 < s <= -score_threshold: # Excluding 0 - there's literally no point looking at those
                indices_for_sample.append(j)
    sample = [wd_candidates[k] for k in indices_for_sample]
    print('Sample')
    for wd in sample:
        print(wd.full_name())
        print(wd.get_teff())
        print(wd.get_logg())
        print(wd.pewdd_metric)
        print()
    return sample

def output_sample_to_csv(name_of_sample, final_sample, manager):
    indices_to_print = [0] # Always print header
    for i, wd in enumerate(manager.white_dwarfs):
        if wd in final_sample:
            indices_to_print.append(i+1) #Add 1 for header!
    infile = pu.get_path_to_data() + 'PEWDD.csv'
    outfile = pu.get_path_to_data() + 'PEWDD_' + name_of_sample + '.csv'
    with open(outfile, 'w', newline='', encoding='utf-8') as of:
        to_write = csv.writer(of)
        with open(infile, encoding='utf-8') as in_f:
            for j, row in enumerate(csv.reader(in_f, delimiter=',')):
                if j in indices_to_print:
                    to_write.writerow(row)

def pick_out_sample(name_of_sample):
    min_teff, max_teff, min_logg, max_logg = find_shared_grid_range(
        sample_selection_settings[name_of_sample][mp.WDParameter.atmospheric_type],
        sample_selection_settings[name_of_sample][timescale_pairs_string]
    )
    manager = pi.load_pewdd_manager()
    filter1 = pi.filter_wds_by_num_modellable_detected_elements(manager.white_dwarfs, sample_selection_settings[name_of_sample][num_els_string])
    filter2 = pi.filter_wds_by_atmospheric_type(filter1, sample_selection_settings[name_of_sample][mp.WDParameter.atmospheric_type])
    wd_candidates = pi.filter_wds_by_teff_logg(filter2, min_teff, max_teff, min_logg, max_logg)
    for wd in wd_candidates:
        print(wd.full_name())

    final_sample = compile_sample(
        wd_candidates,
        sample_selection_settings[name_of_sample][score_threshold_string],
        sample_selection_settings[name_of_sample][timescale_pairs_string]
    )
    output_sample_to_csv(name_of_sample, final_sample, manager)
    return final_sample

def pick_out_all_dbs():
    manager = pi.load_pewdd_manager()
    dbs = pi.filter_wds_by_atmospheric_type(manager.white_dwarfs, ci.Element.He)
    toret = pi.filter_wds_by_teff_logg(dbs, -np.inf, np.inf, -np.inf, np.inf)
    return toret

def pick_out_all_das():
    manager = pi.load_pewdd_manager()
    das = pi.filter_wds_by_atmospheric_type(manager.white_dwarfs, ci.Element.H)
    toret = pi.filter_wds_by_teff_logg(das, -np.inf, np.inf, -np.inf, np.inf)
    return toret

def extract_D_values_from_test_pop(test_pop):
    toret = list()
    for wd in test_pop:
        D = get_declining_phase_depth_from_synthetic_wd(wd, timescale_override_dict, timescale_pairs)
        if D < max_num_observable_sinking_timescales:
            toret.append(D)
        print(toret)
        raise
    return toret

def get_declining_phase_depth_from_synthetic_wd(wd, timescale_override_dict, timescale_pairs):
    event_duration = wd.pollution_properties[mp.ModelParameter.accretion_timescale] # years
    time = 1000000*wd.pollution_properties[mp.ModelParameter.t_sinceaccretion] # in which we convert from Myr to yr
    delta_time = max(0, time - event_duration) # yr
    flattened_timescale_pairs = flatten_timescale_pairs(timescale_pairs)
    reference_element = ci.Element.Mg # arbitrary
    reference_sinking_time = np.mean([timescale_override_dict[timescale_type_to_try][reference_element] for timescale_type_to_try in flattened_timescale_pairs])
    return delta_time/reference_sinking_time

def get_proxy_pat(dm):
    # The idea here is that the percentage-above-threshold is a bit crude
    # It doesn't reflect situations where values are tightly clustered just below the threshold
    # We'll try to fix this by calculating a proxy p-a-t
    # Where every individual sample contributes a probability of being affected
    # And by assumption we'll just adopt an error function to describe that
    sigma = 0.5 # Again, an assumption
    scaled_discrepency_metric = (dm - 1)/sigma
    raw_erf = scipy.special.erf(scaled_discrepency_metric)
    toret = 0.5*(1+raw_erf)
    return toret

def plot_metric_for_synthetic_pop(name_of_pewdd_sample_to_plot_against):
    # Idea here is to, for each Teff and logg, calculate the expected metric value by averaging across a synthetic population of pollutants
    # We can use the population synthesis code to do this - but we will need to override the WD teff and log(g)

    grid_steps = 20
    weak_metric_threshold = 0.5
    metric_threshold = 1
    strong_metric_threshold = 2
    include_wd_markers = True

    atm_type = sample_selection_settings[name_of_pewdd_sample_to_plot_against][mp.WDParameter.atmospheric_type]
    timescale_pairs = sample_selection_settings[name_of_pewdd_sample_to_plot_against][timescale_pairs_string]

    all_timescale_types = flatten_timescale_pairs(timescale_pairs)
    if thermohaline_dummy_override_string in all_timescale_types:
        min_teff = 5000
        max_teff = 21000
        min_logg = 7.5
        max_logg = 8.5
    else:
        timescale_types_for_grid_range = [att for att in all_timescale_types if att != thermohaline_dummy_override_string]
        min_teff = max([ti.teff_ranges[atm_type][timescale_type][0] for timescale_type in timescale_types_for_grid_range])
        max_teff = min([ti.teff_ranges[atm_type][timescale_type][1] for timescale_type in timescale_types_for_grid_range])
        min_logg = max([ti.logg_ranges[atm_type][timescale_type][0] for timescale_type in timescale_types_for_grid_range])
        max_logg = min([ti.logg_ranges[atm_type][timescale_type][1] for timescale_type in timescale_types_for_grid_range])

    import synthetic_population as sp
    ca_to_use = -15

    Teff_values = np.linspace(min_teff, max_teff, grid_steps+1)
    logg_values = np.linspace(min_logg, max_logg, grid_steps+1)
    DM_values = np.zeros((len(logg_values), len(Teff_values)))
    PAT_values = np.zeros((len(logg_values), len(Teff_values)))
    PAWT_values = np.zeros((len(logg_values), len(Teff_values)))
    PAST_values = np.zeros((len(logg_values), len(Teff_values)))
    PROXYPAT_values = np.zeros((len(logg_values), len(Teff_values)))
    local_modellable_wd_density = np.zeros((len(logg_values), len(Teff_values)))

    #if atm_type == ci.Element.He:
        #D_values = [0, 0.1, 0.3, 0.6, 1, 2, 5] # Sampling appropriately here might be the best approach
        #D_values = [1] #...on second thoughts, mathematically it should be equivalent to just take the mean value alone! Which is about 1 (see fig 4.7 of my thesis)
        # On third thoughts it will make a difference for the number of times above a certain threshold, but we can just calculate it for each synthetic wd no?
    #else:
    #    D_values = [0]

    population_size = None
    wd_config_to_use = None
    pollution_config_to_use = None
    if atm_type == ci.Element.He:
        #synth_pop_file_name = 'popdump_TestPop1_TestObs_TestMod1.csv'
        synth_pop_file_name = 'popdump_ReferenceDB_RealisticObservererr0p2_NullModeller.csv' #<-- this is not fully self consistent because we used VO to generate the population
        #synth_pop_file_name = 'popdump_SyntheticHollandsTidalKO_HollandsObservererr0p2_NullModeller.csv'
        #synth_pop_file_name = 'popdump_SyntheticHollandsCollisionalKO_HollandsObservererr0p2_NullModeller.csv'
    else:
        synth_pop_file_name = 'popdump_ReferenceDA_RealisticObservererr0p2_NullModeller.csv' #<-- this is not fully self consistent because we used 3P to generate the population
        #synth_pop_file_name = 'popdump_DADeltaFcfPop_RealisticObservererr0p2_StandardModeller.csv'
    synth_pop_file = pu.get_path_to_pipeline_base_dir() + 'popdumps/' + synth_pop_file_name
    test_pop = sp.SyntheticPopulation(population_size, wd_config_to_use, pollution_config_to_use, synth_pop_file)

    #The distribution of D should be independent of Teff and logg. I should use fig 4.7 in my thesis.
    # Doesn't matter if D is 'too high' for warm DAs - the 'chance' parameter will filter them out later.
    # In principle it should also be exponentially decaying... can we sample, say 100 values from an exponentially decaying function?
    # The mean of an exponential distribution is 1/lambda (the first parameter that sets the decay)
    # We want the mean to be about 1 (see fig 4.7 in my thesis)
    # hence lambda = 1 (the second argument is how many values we want to sample)
    #D_values = np.random.exponential(1, 50)

    num_D_values = 50

    cached_timescale_interpolator = get_cached_timescale_interpolator()
    for i, logg in enumerate(logg_values):
        for j, Teff in enumerate(Teff_values):
            print('Calculating grid step ' + str((i*len(logg_values)) + j + 1) + '/' + str((grid_steps+1)*(grid_steps+1)))
            dm_values = list()
            proxy_pat_values = list()
            timescale_override_dict = cached_timescale_interpolator.get_wd_timescales(
                atm_type,
                logg,
                Teff,
                ca_to_use
            )
            max_Teff_dist = 2000
            max_logg_dist = 0.1
            first_iteration = True
            local_modellable_wd_density_val = 0
            while len(dm_values) < (5*num_D_values):
                for wd in test_pop:
                    # We should also filter the wd by Teff and logg, because these might correlate with detected elements!
                    Teff_dist = abs(wd.wd_properties[mp.WDParameter.temperature] - Teff)
                    logg_dist = abs(wd.wd_properties[mp.WDParameter.logg] - logg)
                    if Teff_dist < max_Teff_dist and logg_dist < max_logg_dist:
                        if first_iteration:
                            local_modellable_wd_density_val += 1
                        if wd.observed_abundances is not None and len(wd.observed_abundances) > 1:
                            lambda_val = estimate_a_priori_chance_of_declining_phase(wd, [tt for tt in all_timescale_types if tt not in [None, thermohaline_dummy_override_string]])
                            D_values = np.random.exponential(lambda_val, num_D_values)
                            for D_override in D_values:
                                dm = calculate_discrepancy_metric(wd, timescale_pairs, timescale_override_dict, D_override)
                                dm_values.append(dm)
                                proxy_pat = get_proxy_pat(dm)
                                proxy_pat_values.append(proxy_pat)
                max_Teff_dist += 1000
                max_logg_dist += 0.05
                first_iteration = False
            percentage_above_threshold = (np.asarray(dm_values) >= metric_threshold).sum()/len(dm_values)
            percentage_above_weak_threshold = (np.asarray(dm_values) >= weak_metric_threshold).sum()/len(dm_values)
            percentage_above_strong_threshold = (np.asarray(dm_values) >= strong_metric_threshold).sum()/len(dm_values)
            mean_dm_across_all_wds_and_Ds = np.mean(dm_values)
            mean_proxy_pat_across_all_wds_and_Ds = np.mean(proxy_pat_values)
            DM_values[i,j] = mean_dm_across_all_wds_and_Ds
            PAT_values[i,j] = percentage_above_threshold
            PAWT_values[i,j] = percentage_above_weak_threshold
            PAST_values[i,j] = percentage_above_strong_threshold
            PROXYPAT_values[i,j] = mean_proxy_pat_across_all_wds_and_Ds
            local_modellable_wd_density[i,j] = local_modellable_wd_density_val

    reference_systems = dict()
    full_sample = dict()
    if include_wd_markers:
        reference_systems = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pick_out_sample(name_of_pewdd_sample_to_plot_against)}
        if atm_type == ci.Element.H:
            full_sample = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pick_out_all_das()}
        elif atm_type == ci.Element.He:
            full_sample = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pick_out_all_dbs()}
        else:
            print(atm_type)
    local_modellable_wd_density /= local_modellable_wd_density.sum()
    import graph_factory as gf
    graph_fac = gf.GraphFactory()
    graph_fac.plot_discrepancy_metric(Teff_values, logg_values, DM_values, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample)
    graph_fac.plot_discrepancy_metric_as_percentage_above_threshold(Teff_values, logg_values, PAT_values, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample, metric_threshold)
    graph_fac.plot_discrepancy_metric_as_percentage_above_threshold(Teff_values, logg_values, PAST_values, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample, strong_metric_threshold)
    graph_fac.plot_discrepancy_metric_as_percentage_above_threshold(Teff_values, logg_values, PAWT_values, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample, weak_metric_threshold)
    graph_fac.plot_discrepancy_metric_as_proxy_percentage_above_threshold(Teff_values, logg_values, PROXYPAT_values, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample)
    graph_fac.plot_modellable_wd_density(Teff_values, logg_values, local_modellable_wd_density, atm_type)
    graph_fac.plot_discrepancy_metric(Teff_values, logg_values, DM_values*local_modellable_wd_density, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample, True)
    graph_fac.plot_discrepancy_metric_as_percentage_above_threshold(Teff_values, logg_values, PAT_values*local_modellable_wd_density, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample, metric_threshold, True)
    graph_fac.plot_discrepancy_metric_as_percentage_above_threshold(Teff_values, logg_values, PAST_values*local_modellable_wd_density, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample, strong_metric_threshold, True)
    graph_fac.plot_discrepancy_metric_as_percentage_above_threshold(Teff_values, logg_values, PAWT_values*local_modellable_wd_density, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample, weak_metric_threshold, True)
    graph_fac.plot_discrepancy_metric_as_proxy_percentage_above_threshold(Teff_values, logg_values, PROXYPAT_values*local_modellable_wd_density, timescale_pairs[0][0], timescale_pairs[0][1], atm_type, reference_systems, full_sample, True)

def main():
    #pick_out_sample('OVERSHOOT_DA_CONTROL')
    #pick_out_sample('BVK_DB')
    #pick_out_sample('BVK_DA')
    #pick_out_sample('OVERSHOOT_DB')
    #pick_out_sample('FIXED_OVERSHOOT_DB')
    #pick_out_sample('OVERSHOOT_DA')
    #pick_out_sample('THERMOHALINE_DA')
    #plot_metric_for_synthetic_pop('OVERSHOOT_DA')
    #plot_metric_for_synthetic_pop('FIXED_OVERSHOOT_DB')

    #plot_metric_for_synthetic_pop('BVK_DA')
    #plot_metric_for_synthetic_pop('OVERSHOOT_DA_PATCHED')
    #plot_metric_for_synthetic_pop('THERMOHALINE_DA')
    plot_metric_for_synthetic_pop('BVK_DB')
    plot_metric_for_synthetic_pop('OVERSHOOT_DB')

    pass

if __name__ == '__main__':
    main()
