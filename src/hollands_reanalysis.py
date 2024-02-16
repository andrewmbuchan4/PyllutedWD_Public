#!/usr/bin/env python
# -*- coding: utf-8 -*-

# A script to produce some new plots/results for my thesis, based
# on having a second go at the analysis of the Hollands sample,
# the first go being Harrison+ 2021
# This will use data from/results produced by the original_codebase
# but I'm putting this script in src because I want to leave original_codebase
# untouched as far as possible

from scipy import stats as stat
import csv
import numpy as np

import chemistry_info as ci
import graph_factory as gf
import hollands_abundances as ha
import model_analyser as ma
import white_dwarf_model_new as wdm
print('Using experimental wdm')

def extract_relevant_ratios(wd_names, relevant_wds, cafe, mgfe):
    wd_index = 0
    toret = list()
    for wd_name in wd_names:
        if wd_name in relevant_wds:
            toret.append((cafe[wd_index], mgfe[wd_index]))
        wd_index += 1
    return toret

def get_og_indices(wd_names):
    row_count = 0
    toret = dict()
    with open('../original_codebase/wd_data_1112.csv', encoding='utf-8') as timescales_csv:
        for row in csv.reader(timescales_csv):
            if row[0] in wd_names:
                toret[row_count]= row[0]
            row_count += 1
    return toret

def extract_hollands_timescales(wd_names, wd_indices):
    row_count = 0
    toret = dict()
    with open('../original_codebase/wd_timescales_1112.csv', encoding='utf-8') as timescales_csv:
        for row in csv.reader(timescales_csv):
            if row_count in wd_indices:
                toret[wd_indices[row_count]] = {
                    ci.Element.Al: float(row[0]),
                    ci.Element.Ti: float(row[1]),
                    ci.Element.Ca: float(row[2]),
                    ci.Element.Ni: float(row[3]),
                    ci.Element.Fe: float(row[4]),
                    ci.Element.Cr: float(row[5]),
                    ci.Element.Mg: float(row[6]),
                    ci.Element.Si: float(row[7]),
                    ci.Element.Na: float(row[8]),
                    ci.Element.O: float(row[9]),
                    ci.Element.C: float(row[10]),
                    ci.Element.N: float(row[11])
                }
            row_count += 1
    print(toret)
    return toret

def calculate_original_ratio(observed_ratio, t_1, t_2, t_sinceaccretion_myr, t_disc_myr):
    if observed_ratio is None:
        return None
    # Going to find what the ratio of element1 and element2 must have been in the pollutant in order to match what we observe
    # Without loss of generality, assume the abundance of element 1 is 1, and find the abundance of element 2 to match the observed ratio
    min_element2 = 0.0000001
    max_element2 = 100000000
    wd_timescales = [t_1, t_2]
    tolerance = 0.0001
    while True:
        print()
        print(observed_ratio)
        print(wd_timescales)
        print(t_sinceaccretion_myr)
        print(t_disc_myr)
        trial_element2 = (min_element2 + max_element2) / 2
        print(trial_element2)
        trial_planetesimal_abundances = [1, trial_element2]
        output_abundances = wdm.process_snapshot_abundances(t_sinceaccretion_myr, 1000000*t_disc_myr, trial_planetesimal_abundances, wd_timescales, -6) # On a log scale
        output_ratio = 10**(output_abundances[0] - output_abundances[1])
        print(output_ratio)
        error = abs((output_ratio - observed_ratio)/observed_ratio)
        if error < tolerance:
            print(84)
            print(1/trial_element2)
            return 1/trial_element2
        else:
            if output_ratio > observed_ratio:
                min_element2 = trial_element2
            else:
                max_element2 = trial_element2

def safe_exponentiate(base, power):
    try:
        return base**power
    except TypeError:
        return None

def plot_hollands_elel_ratios(graph_factory, correct_for_sinking=False, use_namg=False):
    manager = ha.load_manager()
    x_el_1 = ci.Element.Ca
    x_el_2 = ci.Element.Fe
    y_el_1 = ci.Element.Mg
    y_el_2 = ci.Element.Fe
    if use_namg:
        # Then all the 'mgfe' stuff is actually namg
        y_el_1 = ci.Element.Na
        y_el_2 = ci.Element.Mg
    cafe = [safe_exponentiate(10, elel) for elel in ha.get_hollands_el_el_values(x_el_1, x_el_2, manager)]
    mgfe = [safe_exponentiate(10, elel) for elel in ha.get_hollands_el_el_values(y_el_1, y_el_2, manager)]
    wd_names = ha.get_hollands_names(manager)
    print(wd_names)
    wd_og_indices = get_og_indices(wd_names)
    timescales = extract_hollands_timescales(wd_names, wd_og_indices)
    if correct_for_sinking:
        info_dict = ha.get_hollands_abundances_and_timescales_dict(manager)
        print(info_dict)
        corrected_cafe = list()
        corrected_mgfe = list()
        for i, wd_name in enumerate(wd_names):
            c_cafe = calculate_original_ratio(
                cafe[i],
                timescales[wd_name][x_el_1],
                timescales[wd_name][x_el_2],
                info_dict[wd_name]['Median time since accretion /Myr'],
                info_dict[wd_name]['Median accretion timescale /Myr']
            )
            c_mgfe = calculate_original_ratio(
                mgfe[i],
                timescales[wd_name][y_el_1],
                timescales[wd_name][y_el_2],
                info_dict[wd_name]['Median time since accretion /Myr'],
                info_dict[wd_name]['Median accretion timescale /Myr']
            )
            corrected_cafe.append(c_cafe)
            corrected_mgfe.append(c_mgfe)
        assert len(cafe) == len(corrected_cafe)
        assert len(mgfe) == len(corrected_mgfe)
    #print(len(wd_names))
    #   0           1                  2            3               4                        5                          6            7           8             9                10              11                    12                     13                          14              15                       16                 17
    #['System', 'Best fit model', 'Evidence', 'Chi Squared', 'Number of Elements', 'Chi Squared per data point', 'Good Fit?', 'Primitive?', 'Core Rich?', 'Mantle Rich?', 'Crust Rich?', 'Volatile Rich?', 'Volatile Depleted?', 'Moderate Volatile Depleted?', 'Build Up Phase?', 'Steady State Phase?', 'Declining Phase?', 'Model parameters']
    primitive_wds = list()
    core_rich_wds = list()
    mantle_rich_wds = list()
    crust_rich_wds = list()
    volatile_rich_wds = list()
    volatile_depleted_wds = list()
    moderate_volatile_depleted_wds = list()
    build_up_wds = list()
    steady_state_wds = list()
    declining_wds = list()
    row_no = 0
    with open('classifications.csv', encoding='utf-8') as best_fits_csv:
        for row in csv.reader(best_fits_csv):
            if row_no == 0:
                pass
            else:
                if row[7] == 'True':
                    primitive_wds.append(row[0])
                if row[8] == 'True':
                    core_rich_wds.append(row[0])
                if row[9] == 'True':
                    mantle_rich_wds.append(row[0])
                if row[10] == 'True':
                    crust_rich_wds.append(row[0])
                if row[11] == 'True':
                    volatile_rich_wds.append(row[0])
                if row[12] == 'True':
                    volatile_depleted_wds.append(row[0])
                if row[13] == 'True':
                    moderate_volatile_depleted_wds.append(row[0])
                if row[14] == 'True':
                    build_up_wds.append(row[0])
                if row[15] == 'True':
                    steady_state_wds.append(row[0])
                if row[16] == 'True':
                    declining_wds.append(row[0])
            row_no += 1
    print(primitive_wds)
    print(core_rich_wds)
    print(mantle_rich_wds)
    print(crust_rich_wds)
    print(volatile_rich_wds)
    print(volatile_depleted_wds)
    print(moderate_volatile_depleted_wds)
    print(build_up_wds)
    print(steady_state_wds)
    print(declining_wds)
    if correct_for_sinking:
        wd_data = dict()
        wd_data['Primitive'] = extract_relevant_ratios(wd_names, primitive_wds, corrected_cafe, corrected_mgfe)
        wd_data['Core rich'] = extract_relevant_ratios(wd_names, core_rich_wds, corrected_cafe, corrected_mgfe)
        #wd_data['Mantle rich'] = extract_relevant_ratios(wd_names, mantle_rich_wds, corrected_cafe, corrected_mgfe)
        wd_data['Crust rich'] = extract_relevant_ratios(wd_names, crust_rich_wds, corrected_cafe, corrected_mgfe)
        #wd_data['Volatile rich'] = extract_relevant_ratios(wd_names, volatile_rich_wds, corrected_cafe, corrected_mgfe)
        wd_data['VD'] = extract_relevant_ratios(wd_names, volatile_depleted_wds, corrected_cafe, corrected_mgfe)
        wd_data['MVD'] = extract_relevant_ratios(wd_names, moderate_volatile_depleted_wds, corrected_cafe, corrected_mgfe)
        #wd_data['Build-up'] = extract_relevant_ratios(wd_names, build_up_wds, corrected_cafe, corrected_mgfe)
        #wd_data['Steady State'] = extract_relevant_ratios(wd_names, steady_state_wds, corrected_cafe, corrected_mgfe)
        #wd_data['Declining'] = extract_relevant_ratios(wd_names, declining_wds, corrected_cafe, corrected_mgfe)
    else:
        wd_data = dict()
        wd_data['Primitive'] = extract_relevant_ratios(wd_names, primitive_wds, cafe, mgfe)
        wd_data['Core rich'] = extract_relevant_ratios(wd_names, core_rich_wds, cafe, mgfe)
        #wd_data['Mantle rich'] = extract_relevant_ratios(wd_names, mantle_rich_wds, cafe, mgfe)
        wd_data['Crust rich'] = extract_relevant_ratios(wd_names, crust_rich_wds, cafe, mgfe)
        #wd_data['Volatile rich'] = extract_relevant_ratios(wd_names, volatile_rich_wds, cafe, mgfe)
        wd_data['VD'] = extract_relevant_ratios(wd_names, volatile_depleted_wds, cafe, mgfe)
        wd_data['MVD'] = extract_relevant_ratios(wd_names, moderate_volatile_depleted_wds, cafe, mgfe)
        #wd_data['Build-up'] = extract_relevant_ratios(wd_names, build_up_wds, cafe, mgfe)
        #wd_data['Steady State'] = extract_relevant_ratios(wd_names, steady_state_wds, cafe, mgfe)
        #wd_data['Declining'] = extract_relevant_ratios(wd_names, declining_wds, cafe, mgfe)
    stellar_cafe = list()
    stellar_mgfe = list()
    for stellar_composition in manager.stellar_compositions:
        if use_namg:
            stellar_cafe.append(stellar_composition[2]/stellar_composition[4])
            stellar_mgfe.append(stellar_composition[7])
        else:
            stellar_cafe.append(stellar_composition[2]/stellar_composition[4])
            stellar_mgfe.append(1/stellar_composition[4])
    stellar_data = list(zip(stellar_cafe, stellar_mgfe))
    graph_name = 'hollands_'+ str(x_el_1) + str(x_el_2) + str(y_el_1) + str(y_el_2) +'_scatter'
    if correct_for_sinking:
        graph_name += '_corrected'
    plot_dict = graph_factory.make_hollands_elel_scatter_plot(wd_data, stellar_data, x_el_1, x_el_2, y_el_1, y_el_2, graph_name)
    return plot_dict

def analyse_population_material():
    graph_fac = gf.GraphFactory()
    plot_dict1 = plot_hollands_elel_ratios(graph_fac, False, False)
    plot_dict2 = plot_hollands_elel_ratios(graph_fac, True, False)
    plot_dict3 = plot_hollands_elel_ratios(graph_fac, False, True)
    plot_dict4 = plot_hollands_elel_ratios(graph_fac, True, True)
    graph_fac.multipanelise(
        [plot_dict1, plot_dict2],
        1,
        2,
        ['hollands_material_cafemgfe_multipanel.pdf', 'hollands_material_cafemgfe_multipanel.png'],
        5,
        10,
        None,
        None,
        False,
        False
    )
    graph_fac.multipanelise(
        [plot_dict3, plot_dict4],
        1,
        2,
        ['hollands_material_cafenamg_multipanel.pdf', 'hollands_material_cafenamg_multipanel.png'],
        5,
        10,
        None,
        None,
        False,
        False
    )

def estimate_accretion_phase(t_sinceaccretion_samples, accretion_timescale_samples, time_to_reach_steady_state_myr):
    # t_sinceaccretion will be in Myr, accretion_timescales AKA t_event AKA t_disc in log(yr)
    assert len(t_sinceaccretion_samples) == len(accretion_timescale_samples)
    i = 0
    bu_count = 0
    ss_count = 0
    bu_dec_count = 0
    ss_dec_count = 0
    delta_times = list()
    while i < len(accretion_timescale_samples):
        t_sinceaccretion = t_sinceaccretion_samples[i]
        t_disc = (10**accretion_timescale_samples[i])/1000000
        delta_times.append(t_sinceaccretion - t_disc)
        if t_disc > t_sinceaccretion:
            #Then it's either build-up or steady state
            if t_sinceaccretion > time_to_reach_steady_state_myr:
                # Steady state
                ss_count += 1
            else:
                # Build - up
                bu_count += 1
        else:
            # Then it's some flavour of declining phase
            if t_disc > time_to_reach_steady_state_myr:
                # Then it reached steady state before declining
                ss_dec_count += 1
            else:
                # It finished accreting while still in build-up
                bu_dec_count += 1
        i += 1
    median_delta_time = np.percentile(delta_times, 50)
    return {'bu': bu_count, 'ss': ss_count, 'bu_dec': bu_dec_count, 'ss_dec': ss_dec_count, 'total': i, 't_ss': time_to_reach_steady_state_myr, 'median_delta_time': median_delta_time}

def check_if_primitive(model_parameters):
    primitive_parameters = ['Stellar metallicity indices', 'Time since Accretion/Myrs', 'log(Pollution Fraction)', 'log(Accretion Event Timescale/Yrs)']
    for model_parameter in model_parameters:
        if model_parameter not in primitive_parameters:
            return False
    return True

def reclassify_systems():
    manager = ha.load_manager()
    wd_names = ha.get_hollands_names(manager)
    accretion_timescale_dict = dict()
    accretion_timescale_median_dict = dict()
    accretion_timescale_onesigma_dict = dict()

    wd_names = ha.get_hollands_names(manager)
    wd_og_indices = get_og_indices(wd_names)
    timescales = extract_hollands_timescales(wd_names, wd_og_indices)

    #["Stellar metallicity indices","Time since Accretion/Myrs", "log(Formation Distance/AU)","Feeding Zone Size/AU", "Parent Core Fraction", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],

    classification = dict()

    for wd_name in sorted(wd_names):
        stats, weightpost, log_likelihood, best_fit_param_values, model_parameters, model_number = ha.extract_posteriors(wd_name)
        print()
        print(model_parameters)
        print(model_number)
        print(best_fit_param_values)

        accretion_timescale_index = model_parameters.index('log(Accretion Event Timescale/Yrs)')
        accretion_timescale_samples = weightpost[:,accretion_timescale_index]
        t_sinceaccretion_index = model_parameters.index('Time since Accretion/Myrs')
        t_sinceaccretion_samples = weightpost[:,t_sinceaccretion_index]

        accretion_phase_estimate_dict = estimate_accretion_phase(t_sinceaccretion_samples, accretion_timescale_samples, 0.000005*timescales[wd_name][ci.Element.Mg])
        # I worry that this is horribly biased towards declining phase solutions, since we sample the accretion timescale in log space, so the alternative
        # below just looks at the single best fit solution. But I feel like really the bias isn't a problem - it just reflects the Bayesian nature of what we're doing
        # (ie the priors favour declining phase)
        #
        #accretion_phase_estimate_dict = estimate_accretion_phase(
        #    [best_fit_param_values[t_sinceaccretion_index]],
        #    [best_fit_param_values[accretion_timescale_index]],
        #    0.000005*timescales[wd_name][ci.Element.Mg]
        #)


        bu_count = accretion_phase_estimate_dict['bu']
        ss_count = accretion_phase_estimate_dict['ss']
        dec_count = accretion_phase_estimate_dict['bu_dec'] + accretion_phase_estimate_dict['ss_dec']
        if bu_count > ss_count and bu_count > dec_count:
            buildup_phase = True
            steady_state_phase = False
            declining_phase = False
        elif ss_count > bu_count and ss_count > dec_count:
            buildup_phase = False
            steady_state_phase = True
            declining_phase = False
        elif dec_count > bu_count and dec_count > ss_count:
            buildup_phase = False
            steady_state_phase = False
            declining_phase = True
        else:
            raise # Maybe we have a tie?!

        total_count = bu_count + ss_count + dec_count

        #bu_sigma = abs(stat.norm.ppf(bu_count/total_count))
        #ss_sigma = abs(stat.norm.ppf(ss_count/total_count))
        #dec_sigma = abs(stat.norm.ppf(dec_count/total_count))
        #
        #max_sigma = abs(stat.norm.ppf(1/total_count))
        #if np.isinf(bu_sigma):
        #    bu_sigma = max_sigma
        #if np.isinf(ss_sigma):
        #    ss_sigma = max_sigma
        #if np.isinf(dec_sigma):
        #    dec_sigma = max_sigma

        best_model_evidence = stats['global evidence']

        #best_chi_square_M1 = -2.0 * (log_likelihood - norm_log) # Probably easier to reuse old calculation rather than calculate norm_log

        primitive = check_if_primitive(model_parameters) # In the end, I just manually adjusted the original table in the cases where it was apparently primitive and (moderate) volatile depleted

        volatile_poor = "log(Formation Distance/AU)" in model_parameters # New logic: if the model invoked heating, its volatile poor

        classification[wd_name] = {
            'Model': model_parameters,
            'Evidence': best_model_evidence, # This parameter is basically just here to make me comfortable that we have the right model
            #'ChiSq': chi_sq, # not reanalysing this
            #'NumElements': num_elements, # not reanalysing this
            #'ChiSqPerDataPoint': chi_sq_per_data_point, # not reanalysing this
            #'GoodFit': good_fit, # not reanalysing this
            'Primitive': primitive,
            #'CoreRich': core_rich,  # not reanalysing this
            #'MantleRich': mantle_rich, # not reanalysing this
            #'CrustRich': crust_rich, # not reanalysing this
            'VolatileRich': not volatile_poor,
            'VolatilePoor': volatile_poor,
            # 'ModerateVolatileDepleted': moderate_volatile_depleted, Not separating VD and MVD in this analysis - can reuse from previous
            'BuildUpPhase': buildup_phase,
            'SteadyStatePhase': steady_state_phase,
            'DecliningPhase': declining_phase,
            #'BuildUpSigma': bu_sigma,
            #'SteadyStateSigma': ss_sigma,
            #'DecliningSigma': dec_sigma,
            'BuildUpPC': 100*(bu_count/total_count),
            'SteadyStatePC': 100*(ss_count/total_count),
            'DecliningPC': 100*(dec_count/total_count)
        }
    print(classification)
    with open('best_fits_hb20_revisited.csv', 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            'System',
            'Model parameters',
            'Evidence',
            'Primitive?',
            'Volatile Rich?',
            'Volatile Poor?',
            'Build Up Phase?',
            'Steady State Phase?',
            'Declining Phase?',
            #'Build Up Sigma',
            #'Steady State Sigma',
            #'Declining Sigma',
            'Build Up %',
            'Steady State %',
            'Declining %'
        ])
        for sys, vals in classification.items():
            to_write.writerow([
                sys,
                vals['Model'],
                vals['Evidence'],
                vals['Primitive'],
                vals['VolatileRich'],
                vals['VolatilePoor'],
                vals['BuildUpPhase'],
                vals['SteadyStatePhase'],
                vals['DecliningPhase'],
                #vals['BuildUpSigma'],
                #vals['SteadyStateSigma'],
                #vals['DecliningSigma'],
                vals['BuildUpPC'],
                vals['SteadyStatePC'],
                vals['DecliningPC']
            ])
    return classification

def extract_accretion_lifetime_info(accretion_phase_estimate_dict, threshold=0.67):
    lower_bounds = list()
    upper_bounds = list()
    for wd_name, phase_info in accretion_phase_estimate_dict.items():
        print()
        print(wd_name)
        if phase_info['bu']/phase_info['total'] > threshold:
            print('P(bu) = ' + str(phase_info['bu']/phase_info['total']))
            print('t_event should be more than t_sinceaccretion')
            raise # This seems to be irrelevant!
        elif phase_info['ss']/phase_info['total'] > threshold:
            print('P(ss) = ' + str(phase_info['ss']/phase_info['total']))
            print('t_event should be more than ' + str(phase_info['t_ss']) + ' Myr')
            lower_bounds.append(phase_info['t_ss'])
        elif phase_info['bu_dec']/phase_info['total'] > threshold:
            print('P(bu_dec) = ' + str(phase_info['bu_dec']/phase_info['total']))
            print('t_event should be less than ' + str(phase_info['t_ss']) + ' Myr')
            upper_bounds.append(phase_info['t_ss'])
        elif phase_info['ss_dec']/phase_info['total'] > threshold:
            print('P(ss_dec) = ' + str(phase_info['bu_dec']/phase_info['total']))
            print('t_event should be more than ' + str(phase_info['t_ss']) + ' Myr')
            lower_bounds.append(phase_info['t_ss'])
        else:
            print('Accretion phase unclear')
    print(lower_bounds)
    print(upper_bounds)
    # These two numbers define a range of possible values of t_event that would be consistent with all the (meaningful) constraints
    # There's no guarantee that this would even give a sensible range, but it happens to come out to 0.2-0.7 Myr
    # The constraining systems are SDSSJ1211+2326 and SDSSJ0046+2717
    print(max(lower_bounds))
    print(min(upper_bounds))

def analyse_accretion_lifetime():
    manager = ha.load_manager()
    wd_names = ha.get_hollands_names(manager)
    accretion_timescale_dict = dict()
    accretion_timescale_median_dict = dict()
    accretion_timescale_onesigma_dict = dict()

    wd_names = ha.get_hollands_names(manager)
    wd_og_indices = get_og_indices(wd_names)
    timescales = extract_hollands_timescales(wd_names, wd_og_indices)

    accretion_phase_estimate_dict = dict()

    for wd_name in wd_names:
        stats, weightpost, log_likelihood, best_fit_param_values, model_parameters, model_number = ha.extract_posteriors(wd_name)
        accretion_timescale_index = model_parameters.index('log(Accretion Event Timescale/Yrs)')
        accretion_timescale_samples = weightpost[:,accretion_timescale_index]
        #accretion_timescale_index = model_parameters.index('Time since Accretion/Myrs')

        accretion_timescale_median_dict[wd_name] = np.percentile(accretion_timescale_samples, 50)
        percentile_16 = np.percentile(accretion_timescale_samples, 16)
        percentile_84 = np.percentile(accretion_timescale_samples, 84)
        filter_1 = np.extract(accretion_timescale_samples <= percentile_84, accretion_timescale_samples)  # Remove anything above the 84th percentile
        filter_2 = np.extract(filter_1 >= percentile_16, filter_1) # Remove anything below the 16th percentile. What's left is the one sigma distribution
        accretion_timescale_onesigma_dict[wd_name] = filter_2
        accretion_timescale_dict[wd_name] = accretion_timescale_samples

        t_sinceaccretion_index = model_parameters.index('Time since Accretion/Myrs')
        t_sinceaccretion_samples = weightpost[:,t_sinceaccretion_index]
        accretion_phase_estimate_dict[wd_name] = estimate_accretion_phase(t_sinceaccretion_samples, accretion_timescale_samples, 0.000005*timescales[wd_name][ci.Element.Mg])
    extract_accretion_lifetime_info(accretion_phase_estimate_dict)

    print(accretion_timescale_dict)
    print(accretion_timescale_median_dict)
    print(accretion_timescale_onesigma_dict)
    print(accretion_phase_estimate_dict)
    scaled_median_delta_times = [accretion_phase_estimate_dict[wd_name]['median_delta_time']/(0.000001*timescales[wd_name][ci.Element.Mg]) for wd_name in accretion_phase_estimate_dict]
    print()
    print()
    print()
    print(accretion_phase_estimate_dict['SDSSJ0744+4649']['median_delta_time'])
    print(0.000001*timescales['SDSSJ0744+4649'][ci.Element.Mg])
    print(sorted(scaled_median_delta_times))
    model_analyser = ma.ModelAnalyser('graphs')

    half_delta_bin_size = 0.25
    delta_time_bins, delta_time_bar_centres = model_analyser.generate_bins_and_bar_centres(0, 5, half_delta_bin_size)

    median_delta_time_heights, dtbins2 = np.histogram(
        scaled_median_delta_times,
        delta_time_bins,
        density=True
    )
    model_analyser.graph_fac.make_histogram(
        delta_time_bar_centres,
        [median_delta_time_heights],
        ['Cool DZs'],
        'declining_phase',
        2*half_delta_bin_size,
        1.2,
        'Mg sinking timescales since accretion',
        'delta_times',
        None,
        None
    )

    raise
    all_medians = [accretion_timescale_median_dict[wd_name] for wd_name in accretion_timescale_median_dict]
    print(all_medians)

    individual_systems_to_plot_by_index = [20, 50, 71, 154, 180]
    individual_systems_to_plot_by_name = [wd_names[i] for i in individual_systems_to_plot_by_index]

    all_samples = list()
    for wd_name, samples in accretion_timescale_dict.items():
        all_samples.extend(samples)

    all_onesigma_samples = list()
    for wd_name, samples in accretion_timescale_onesigma_dict.items():
        all_onesigma_samples.extend(samples)

    small_half_bin_size = 0.1
    bins_small, bar_centres_small = model_analyser.generate_bins_and_bar_centres(0, 8, small_half_bin_size)
    big_half_bin_size = 0.2
    bins_big, bar_centres_big = model_analyser.generate_bins_and_bar_centres(0, 8, big_half_bin_size)

    medians_heights, bins2 = np.histogram(
        all_medians,
        bins_big,
        density=True
    )

    model_analyser.graph_fac.make_histogram(
        bar_centres_big,
        [medians_heights],
        ['Medians'],
        't_disc',
        big_half_bin_size*2,
        1.2,
        'log(Accretion Event Timescale/Yrs)',
        'medians',
        None,
        None
    )

    all_onesigma_heights, bins2 = np.histogram(
        all_onesigma_samples,
        bins_small,
        density=True
    )

    model_analyser.graph_fac.make_histogram(
        bar_centres_small,
        [all_onesigma_heights],
        ['Summed Posteriors, 1 Sigma only'],
        't_disc',
        small_half_bin_size*2,
        1.2,
        'log(Accretion Event Timescale/Yrs)',
        'summed_onesigma',
        None,
        None
    )

    all_heights, bins2 = np.histogram(
        all_samples,
        bins_small,
        density=True
    )

    model_analyser.graph_fac.make_histogram(
        bar_centres_small,
        [all_heights],
        ['Summed Posteriors'],
        't_disc',
        small_half_bin_size*2,
        1.2,
        'log(Accretion Event Timescale/Yrs)',
        'summed',
        None,
        None
    )

    individual_binned_names = list()
    individual_binned_samples = list()
    for name, samples in accretion_timescale_dict.items():
        if name in individual_systems_to_plot_by_name:
            system_heights, bins3 = np.histogram(
                samples,
                bins_small,
                density=True
            )
            individual_binned_names.append(name)
            individual_binned_samples.append(system_heights)

    model_analyser.graph_fac.make_histogram(
        bar_centres_small,
        individual_binned_samples,
        individual_binned_names,
        't_disc',
        small_half_bin_size*2,
        1.2,
        'log(Accretion Event Timescale/Yrs)',
        'individuals',
        None,
        None
    )

def main():
    reclassify_systems()
    #analyse_population_material()
    #analyse_accretion_lifetime()

if __name__ == '__main__':
    main()
