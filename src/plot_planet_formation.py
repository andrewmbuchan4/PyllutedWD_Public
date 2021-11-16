#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import chemistry_info as ci
import geology_info as gi
import graph_factory as gf

def generate_pressure_vals():
    return list(range(0, 61, 1))
    #return list(range(0, 10, 1)) + list(range(10, 105, 5))
    
def generate_fO2_vals():
    return[-3, -2.8, -2.6, -2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1]  # Reducing range to avoid non-convergence!

def generate_fcf_vals():
    return [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019] + [f / 100.0 for f in range(2, 101, 1)]

earth_Ds = { #Values from Rudge 2010 supplementary table 5
    ci.Element.W: 10**1.513,
    ci.Element.Ni: 10**1.418,
    ci.Element.P: 10**1.398,
    ci.Element.Co: 10**1.381,
    ci.Element.Pb: 10**1.159,
    ci.Element.Fe: 10**1.136,
    ci.Element.Cu: 10**0.801,
    ci.Element.V: 10**0.262,
    ci.Element.Cr: 10**0.195,
    ci.Element.Mn: 10**-0.155,
    ci.Element.Nb: 10**-0.276,
    ci.Element.Ta: 10**-0.611,
    ci.Element.Si: 10**-0.728,
    ci.Element.Zn: 10**-0.824,
    ci.Element.Ga: 10**-1.0,
    ci.Element.Ti: 0
}

colour_dict = {
    ci.Element.Ca: '#191970',
    ci.Element.Al: '#66CDAA',
    ci.Element.Fe: '#C0C0C0',
    ci.Element.Si: '#DF2E20',
    ci.Element.Mg: '#87CEFA',
    ci.Element.O: '#FFD700',
    ci.Element.C: '#7FFFD4',
    ci.Element.Cr: '#808000',
    ci.Element.Ni: '#22CD22',
    ci.Element.S: '#FFA500',
    ci.Element.Placeholder: '#000000'
}

def form_planets(geo_model, start_from_prev_result=False):
    # Each planet is defined by a pressure/fO2 combo
    planet_dict = dict()
    planet_cnf_dict = dict()
    planet_Ds_dict = dict()
    critical_fragment_core_fraction_dict = dict()
    dlogX_dP_dict = dict()
    for fO2 in generate_fO2_vals():
        Ds = None
        cnf = None
        prev_Ds = None
        prev_cnf = None
        prev_P = None
        for pressure in generate_pressure_vals():
            print()
            print(pressure)
            print(fO2)
            abundances, cnf, Ds, all_Ds = geo_model.form_a_planet_iteratively(pressure, fO2, None, None, cnf if start_from_prev_result else None, Ds if start_from_prev_result else None, True)
            #abundances, cnf, Ds, all_Ds = geo_model.form_a_planet_iteratively(pressure, fO2, None, None, cnf if start_from_prev_result else None, earth_Ds, True)
            if abundances is not None:
                planet_dict[(pressure, fO2)] = abundances
            if cnf is not None:
                planet_cnf_dict[(pressure, fO2)] = cnf
            if all_Ds is not None:
                planet_Ds_dict[(pressure, fO2)] = all_Ds
            critical_fragment_core_fraction_dict[(pressure, fO2)] = dict()
            dlogX_dP_dict[(pressure, fO2)] = dict()
            if abundances is None:
                print('Warning! Abundances = None, could distort dlogX_dP results')
            else:
                for element in abundances.keys():
                    if element == ci.Element.O:
                        continue  # This calculation is not applicable to O at present TODO: It will be soon after implementing James Badro's O parametrisation
                    if prev_Ds is None or prev_cnf is None or prev_P is None:
                        critical_fragment_core_fraction_dict[(pressure, fO2)][element] = None
                        for fcf in generate_fcf_vals():
                            if dlogX_dP_dict[(pressure, fO2)].get(fcf) is None:
                                dlogX_dP_dict[(pressure, fO2)][fcf] = dict()
                            dlogX_dP_dict[(pressure, fO2)][fcf][element] = None
                    else:
                        # Very crude estimates coming up...
                        dD_dP = (Ds[element] - prev_Ds[element])/(pressure - prev_P)
                        dcnf_dP = (cnf - prev_cnf)/(pressure - prev_P)
                        f_crit = geo_model.get_critical_fragment_core_fraction(Ds[element], cnf, dD_dP, dcnf_dP)
                        critical_fragment_core_fraction_dict[(pressure, fO2)][element] = f_crit
                        for fcf in generate_fcf_vals():
                            if dlogX_dP_dict[(pressure, fO2)].get(fcf) is None:
                                dlogX_dP_dict[(pressure, fO2)][fcf] = dict()
                            dlogX_dP = geo_model.get_dlogX_dV(Ds[element], cnf, fcf, dD_dP, dcnf_dP)
                            dlogX_dP_dict[(pressure, fO2)][fcf][element] = dlogX_dP
                            if pressure == 54 and element == ci.Element.Ni:
                                print(dlogX_dP)
            prev_Ds = Ds
            prev_cnf = cnf
            prev_P = pressure
    return planet_dict, planet_cnf_dict, planet_Ds_dict, critical_fragment_core_fraction_dict, dlogX_dP_dict
    
def form_planets_with_variable_sulfur(geo_model, pressure=None, fO2=None):
    # Each planet is defined by a bulk sulfur abundance
    planet_dict = dict()
    planet_cnf_dict = dict()
    planet_Ds_dict = dict()
    if pressure is None:
        pressure = geo_model.get_earth_differentiation_pressure()
    if fO2 is None:
        fO2 = geo_model.get_earth_oxygen_fugacity()
    geo_model.reinit()  # This ensures that the element_info is returned from an unknown state to the default
    default_abundances = dict()
    for el, abundances in geo_model.element_info.items(): 
        default_abundances[el] = abundances[gi.Layer.bulk]
    sulfur_vals = [s / 10000 for s in range(0, 501, 1)]
    #s_val = 0
    #while s_val < 0.05:
    #    sulfur_vals.append(s_val)
    #    s_val += 0.0005
    for s in sulfur_vals:
        new_abundances = dict()
        for el, abundance in default_abundances.items(): 
            new_abundances[el] = abundance
        new_abundances[ci.Element.S] = s
        geo_model.reinit(new_abundances)
        abundances, cnf, Ds, all_Ds = geo_model.form_a_planet_iteratively(pressure, fO2)
        if abundances is not None:
            planet_dict[s] = abundances
        if cnf is not None:
            planet_cnf_dict[s] = cnf
        if all_Ds is not None:
            planet_Ds_dict[s] = Ds
    return planet_dict, planet_cnf_dict, planet_Ds_dict, sulfur_vals

def plot_wd_comparisons(planet_dict, tag):
    wd_dict = {
        'TestCore65': {
            'cnf': 0.65,
            'observed_log_elfe_ratios': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elfe_ratio_errors': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elmg_ratios': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elmg_ratio_errors': {ci.Element.Cr: None, ci.Element.Ni: None}
        },
        'TestCore99': {
            'cnf': 0.99,
            'observed_log_elfe_ratios': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elfe_ratio_errors': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elmg_ratios': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elmg_ratio_errors': {ci.Element.Cr: None, ci.Element.Ni: None}
        },
        'TestCore20': {
            'cnf': 0.2,
            'observed_log_elfe_ratios': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elfe_ratio_errors': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elmg_ratios': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elmg_ratio_errors': {ci.Element.Cr: None, ci.Element.Ni: None}
        },
        'TestMantle99': {
            'cnf': 0.01,
            'observed_log_elfe_ratios': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elfe_ratio_errors': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elmg_ratios': {ci.Element.Cr: None, ci.Element.Ni: None},
            'observed_log_elmg_ratio_errors': {ci.Element.Cr: None, ci.Element.Ni: None}
        }
    }
    geo_model = gi.GeologyModel()
    earth_pressure = geo_model.get_earth_differentiation_pressure()
    earth_fO2 = geo_model.get_earth_oxygen_fugacity()
    earth_pressure_index = None
    earth_fO2_index = None
    for wd, wd_props in wd_dict.items():
        log_crfe_array = np.zeros((len(generate_fO2_vals()), len(generate_pressure_vals())))
        log_nife_array = np.zeros((len(generate_fO2_vals()), len(generate_pressure_vals())))
        log_crmg_array = np.zeros((len(generate_fO2_vals()), len(generate_pressure_vals())))
        log_nimg_array = np.zeros((len(generate_fO2_vals()), len(generate_pressure_vals())))
        for f_index, f in enumerate(generate_fO2_vals()):
            for p_index, p in enumerate(generate_pressure_vals()):
                if p == earth_pressure:
                    earth_pressure_index = p_index
                if f == earth_fO2:
                    earth_fO2_index = f_index
                try:
                    wd_abundances = geo_model.find_system_specific_abundances(planet_dict[(p, f)], wd_props['cnf'])
                    log_crfe_array[f_index][p_index] = np.log10(wd_abundances[ci.Element.Cr]/wd_abundances[ci.Element.Fe])
                    log_nife_array[f_index][p_index] = np.log10(wd_abundances[ci.Element.Ni]/wd_abundances[ci.Element.Fe])
                    log_crmg_array[f_index][p_index] = np.log10(wd_abundances[ci.Element.Cr]/wd_abundances[ci.Element.Mg])
                    log_nimg_array[f_index][p_index] = np.log10(wd_abundances[ci.Element.Ni]/wd_abundances[ci.Element.Mg])
                except KeyError:
                    log_crfe_array[f_index][p_index] = None
                    log_nife_array[f_index][p_index] = None
                    log_crmg_array[f_index][p_index] = None
                    log_nimg_array[f_index][p_index] = None
        wd_props['log_crfe_array'] = log_crfe_array
        wd_props['log_nife_array'] = log_nife_array
        wd_props['log_crmg_array'] = log_crmg_array
        wd_props['log_nimg_array'] = log_nimg_array  # TODO generate the 4 pairs dynamically...
    for wd, wd_props in wd_dict.items():
        variable_p_log_crfe = list()
        variable_p_log_nife = list()
        variable_p_log_crmg = list()
        variable_p_log_nimg = list()
        for p_index, p in enumerate(generate_pressure_vals()):
            variable_p_log_crfe.append(wd_props['log_crfe_array'][earth_fO2_index][p_index])
            variable_p_log_nife.append(wd_props['log_nife_array'][earth_fO2_index][p_index])
            variable_p_log_crmg.append(wd_props['log_crmg_array'][earth_fO2_index][p_index])
            variable_p_log_nimg.append(wd_props['log_nimg_array'][earth_fO2_index][p_index])
        wd_props['variable_p_log_crfe'] = variable_p_log_crfe
        wd_props['variable_p_log_nife'] = variable_p_log_nife
        wd_props['variable_p_log_crmg'] = variable_p_log_crmg
        wd_props['variable_p_log_nimg'] = variable_p_log_nimg
        variable_f_log_crfe = list()
        variable_f_log_nife = list()
        variable_f_log_crmg = list()
        variable_f_log_nimg = list()
        for f_index, f in enumerate(generate_fO2_vals()):
            variable_f_log_crfe.append(wd_props['log_crfe_array'][f_index][earth_pressure_index])
            variable_f_log_nife.append(wd_props['log_nife_array'][f_index][earth_pressure_index])
            variable_f_log_crmg.append(wd_props['log_crmg_array'][f_index][earth_pressure_index])
            variable_f_log_nimg.append(wd_props['log_nimg_array'][f_index][earth_pressure_index])
        wd_props['variable_f_log_crfe'] = variable_f_log_crfe
        wd_props['variable_f_log_nife'] = variable_f_log_nife
        wd_props['variable_f_log_crmg'] = variable_f_log_crmg
        wd_props['variable_f_log_nimg'] = variable_f_log_nimg
    graph_fac = gf.GraphFactory()
    for wd, wd_props in wd_dict.items():
        ignore_lack_of_observations = wd.startswith('Test')
        graph_fac.plot_log_elel_ratio(wd_props['variable_p_log_crfe'], generate_pressure_vals(), 'Pressure', wd, ci.Element.Cr, ci.Element.Fe, wd_props['cnf'], wd_props['observed_log_elfe_ratios'][ci.Element.Cr], wd_props['observed_log_elfe_ratio_errors'][ci.Element.Cr], tag, earth_fO2, colour_dict, ignore_lack_of_observations)
        graph_fac.plot_log_elel_ratio(wd_props['variable_p_log_nife'], generate_pressure_vals(), 'Pressure', wd, ci.Element.Ni, ci.Element.Fe, wd_props['cnf'], wd_props['observed_log_elfe_ratios'][ci.Element.Ni], wd_props['observed_log_elfe_ratio_errors'][ci.Element.Ni], tag, earth_fO2, colour_dict, ignore_lack_of_observations)
        graph_fac.plot_log_elel_ratio(wd_props['variable_f_log_crfe'], generate_fO2_vals(), 'fO2', wd, ci.Element.Cr, ci.Element.Fe, wd_props['cnf'], wd_props['observed_log_elfe_ratios'][ci.Element.Cr], wd_props['observed_log_elfe_ratio_errors'][ci.Element.Cr], tag, earth_pressure, colour_dict, ignore_lack_of_observations)
        graph_fac.plot_log_elel_ratio(wd_props['variable_f_log_nife'], generate_fO2_vals(), 'fO2', wd, ci.Element.Ni, ci.Element.Fe, wd_props['cnf'], wd_props['observed_log_elfe_ratios'][ci.Element.Ni], wd_props['observed_log_elfe_ratio_errors'][ci.Element.Ni], tag, earth_pressure, colour_dict, ignore_lack_of_observations)
        graph_fac.plot_3d_log_elel_ratio(wd_props['log_crfe_array'], generate_pressure_vals(), generate_fO2_vals(), wd, ci.Element.Cr, ci.Element.Fe, wd_props['cnf'], wd_props['observed_log_elfe_ratios'][ci.Element.Cr], wd_props['observed_log_elfe_ratio_errors'][ci.Element.Cr], tag)
        graph_fac.plot_3d_log_elel_ratio(wd_props['log_nife_array'], generate_pressure_vals(), generate_fO2_vals(), wd, ci.Element.Ni, ci.Element.Fe, wd_props['cnf'], wd_props['observed_log_elfe_ratios'][ci.Element.Ni], wd_props['observed_log_elfe_ratio_errors'][ci.Element.Ni], tag)
        graph_fac.plot_log_elel_ratio(wd_props['variable_p_log_crmg'], generate_pressure_vals(), 'Pressure', wd, ci.Element.Cr, ci.Element.Mg, wd_props['cnf'], wd_props['observed_log_elmg_ratios'][ci.Element.Cr], wd_props['observed_log_elmg_ratio_errors'][ci.Element.Cr], tag, earth_fO2, colour_dict, ignore_lack_of_observations)
        graph_fac.plot_log_elel_ratio(wd_props['variable_p_log_nimg'], generate_pressure_vals(), 'Pressure', wd, ci.Element.Ni, ci.Element.Mg, wd_props['cnf'], wd_props['observed_log_elmg_ratios'][ci.Element.Ni], wd_props['observed_log_elmg_ratio_errors'][ci.Element.Ni], tag, earth_fO2, colour_dict, ignore_lack_of_observations)
        graph_fac.plot_log_elel_ratio(wd_props['variable_f_log_crmg'], generate_fO2_vals(), 'fO2', wd, ci.Element.Cr, ci.Element.Mg, wd_props['cnf'], wd_props['observed_log_elmg_ratios'][ci.Element.Cr], wd_props['observed_log_elmg_ratio_errors'][ci.Element.Cr], tag, earth_pressure, colour_dict, ignore_lack_of_observations)
        graph_fac.plot_log_elel_ratio(wd_props['variable_f_log_nimg'], generate_fO2_vals(), 'fO2', wd, ci.Element.Ni, ci.Element.Mg, wd_props['cnf'], wd_props['observed_log_elmg_ratios'][ci.Element.Ni], wd_props['observed_log_elmg_ratio_errors'][ci.Element.Ni], tag, earth_pressure, colour_dict, ignore_lack_of_observations)
        graph_fac.plot_3d_log_elel_ratio(wd_props['log_crmg_array'], generate_pressure_vals(), generate_fO2_vals(), wd, ci.Element.Cr, ci.Element.Mg, wd_props['cnf'], wd_props['observed_log_elmg_ratios'][ci.Element.Cr], wd_props['observed_log_elmg_ratio_errors'][ci.Element.Cr], tag)
        graph_fac.plot_3d_log_elel_ratio(wd_props['log_nimg_array'], generate_pressure_vals(), generate_fO2_vals(), wd, ci.Element.Ni, ci.Element.Mg, wd_props['cnf'], wd_props['observed_log_elmg_ratios'][ci.Element.Ni], wd_props['observed_log_elmg_ratio_errors'][ci.Element.Ni], tag)
    return 0

def extract_stacked_abundances(stacked_named_elements, variable_p_vals_core, variable_p_vals_mantle, variable_fO2_vals_core, variable_fO2_vals_mantle, planet_cnf_dict, earth_fO2, earth_pressure):
    stacked_variable_p_vals = {gi.Layer.core: dict(), gi.Layer.mantle: dict()}
    bulk_stacked_variable_p_vals = {gi.Layer.core: dict(), gi.Layer.mantle: dict()}  # These will be weighted by core/mantle number fraction at a particular pressure/fO2
    weighted_stacked_variable_offset = [0]*len(variable_p_vals_core[ci.Element.Fe]) # Assume they all have the same number of data points...we should bail if this isn't true
    for layer in [gi.Layer.core, gi.Layer.mantle]:
        ref_vals = variable_p_vals_core if layer == gi.Layer.core else variable_p_vals_mantle
        stacked_variable_offset = [0]*len(ref_vals[ci.Element.Fe])  # Assume they all have the same number of data points... we should bail if this isn't true
        for sne in stacked_named_elements[layer]:
            data = [sum(a) if a[1] is not None else None for a in zip(stacked_variable_offset, ref_vals[sne])]  # If the abundance is None for one element, it should be None for all of them
            stacked_variable_p_vals[layer][sne] = data
            stacked_variable_offset = data
            if layer == gi.Layer.core:
                weighted_ref_vals = [r*planet_cnf_dict[(generate_pressure_vals()[p_index], earth_fO2)] if r is not None else None for p_index, r in enumerate(ref_vals[sne])]  # Multiply by cnf at that pressure. The catch is that each ref val corresponds to a different pressure
            else:
                weighted_ref_vals = [r*(1 - planet_cnf_dict[(generate_pressure_vals()[p_index], earth_fO2)]) if r is not None else None for p_index, r in enumerate(ref_vals[sne])]  # Multiply by 1 - cnf at that pressure
            weighted_data = [sum(a) if a[1] is not None else None for a in zip(weighted_stacked_variable_offset, weighted_ref_vals)]
            bulk_stacked_variable_p_vals[layer][sne] = weighted_data
            weighted_stacked_variable_offset = weighted_data
        for ele, abundance in ref_vals.items():
            if ele not in stacked_named_elements[layer]:
                stacked_variable_offset = [sum(a) if a[1] is not None else None for a in zip(stacked_variable_offset, abundance)]
                if layer == gi.Layer.core:
                    weighted_ref_vals = [r*planet_cnf_dict[(generate_pressure_vals()[p_index], earth_fO2)] if r is not None else None for p_index, r in enumerate(abundance)]  # Multiply by cnf at that pressure. The catch is that each ref val corresponds to a different pressure
                else:
                    weighted_ref_vals = [r*(1 - planet_cnf_dict[(generate_pressure_vals()[p_index], earth_fO2)]) if r is not None else None for p_index, r in enumerate(abundance)]  # Multiply by 1 - cnf at that pressure
                weighted_stacked_variable_offset = [sum(a) if a[1] is not None else None for a in zip(weighted_stacked_variable_offset, weighted_ref_vals)]
        stacked_variable_p_vals[layer][ci.Element.Placeholder] = stacked_variable_offset
        bulk_stacked_variable_p_vals[layer][ci.Element.Placeholder] = weighted_stacked_variable_offset
        
    stacked_variable_fO2_vals = {gi.Layer.core: dict(), gi.Layer.mantle: dict()}
    bulk_stacked_variable_fO2_vals = {gi.Layer.core: dict(), gi.Layer.mantle: dict()}  # These will be weighted by core/mantle number fraction at a particular pressure/fO2
    weighted_stacked_variable_offset = [0]*len(variable_fO2_vals_core[ci.Element.Fe]) # Assume they all have the same number of data points...we should bail if this isn't true
    for layer in [gi.Layer.core, gi.Layer.mantle]:
        ref_vals = variable_fO2_vals_core if layer == gi.Layer.core else variable_fO2_vals_mantle
        stacked_variable_offset = [0]*len(ref_vals[ci.Element.Fe])  # Assume they all have the same number of data points...we should bail if this isn't true
        for sne in stacked_named_elements[layer]:
            data = [sum(a) if a[1] is not None else None for a in zip(stacked_variable_offset, ref_vals[sne])]  # If the abundance is None for one element, it should be None for all of them
            stacked_variable_fO2_vals[layer][sne] = data
            stacked_variable_offset = data
            if layer == gi.Layer.core:
                weighted_ref_vals = [r*planet_cnf_dict[(earth_pressure, generate_fO2_vals()[fO2_index])] if r is not None else None for fO2_index, r in enumerate(ref_vals[sne])]  # Multiply by cnf at that pressure. The catch is that each ref val corresponds to a different pressure
            else:
                weighted_ref_vals = [r*(1 - planet_cnf_dict[(earth_pressure, generate_fO2_vals()[fO2_index])]) if r is not None else None for fO2_index, r in enumerate(ref_vals[sne])]  # Multiply by 1 - cnf at that pressure
            weighted_data = [sum(a) if a[1] is not None else None for a in zip(weighted_stacked_variable_offset, weighted_ref_vals)]
            bulk_stacked_variable_fO2_vals[layer][sne] = weighted_data
            weighted_stacked_variable_offset = weighted_data
        for ele, abundance in ref_vals.items():
            if ele not in stacked_named_elements[layer]:
                stacked_variable_offset = [sum(a) if a[1] is not None else None for a in zip(stacked_variable_offset, abundance)]
                if layer == gi.Layer.core:
                    weighted_ref_vals = [r*planet_cnf_dict[(earth_pressure, generate_fO2_vals()[fO2_index])] if r is not None else None for fO2_index, r in enumerate(abundance)]  # Multiply by cnf at that pressure. The catch is that each ref val corresponds to a different pressure
                else:
                    weighted_ref_vals = [r*(1 - planet_cnf_dict[(earth_pressure, generate_fO2_vals()[fO2_index])]) if r is not None else None for fO2_index, r in enumerate(abundance)]  # Multiply by 1 - cnf at that pressure
                weighted_stacked_variable_offset = [sum(a) if a[1] is not None else None for a in zip(weighted_stacked_variable_offset, weighted_ref_vals)]
        stacked_variable_fO2_vals[layer][ci.Element.Placeholder] = stacked_variable_offset
        bulk_stacked_variable_fO2_vals[layer][ci.Element.Placeholder] = weighted_stacked_variable_offset
    
    return stacked_variable_p_vals, bulk_stacked_variable_p_vals, stacked_variable_fO2_vals, bulk_stacked_variable_fO2_vals

def main():
    geo_model = gi.GeologyModel()
    planet_dict, planet_cnf_dict, planet_Ds_dict, critical_fragment_core_fraction_dict, dlogX_dP_dict = form_planets(geo_model)
    
    earth_pressure = geo_model.get_earth_differentiation_pressure()
    earth_fO2 = geo_model.get_earth_oxygen_fugacity()
    # planet_dict = {Planet_key: {Element1: {Layer1: number, ... } ...} ...}
    elements = planet_dict[(earth_pressure, earth_fO2)].keys()
    
    
    colour_dict = {
        ci.Element.Ca: '#006400',
        ci.Element.C: '#00ff00',
        ci.Element.Fe: '#C0C0C0',
        ci.Element.Si: '#ff0000',
        ci.Element.Ni: '#0000cd',
        ci.Element.O: '#FFD700',
        ci.Element.Al: '#e9967a',
        ci.Element.Mg: '#00ffff',
        ci.Element.Cr: '#ff24C3',
        ci.Element.Placeholder: '#000000'
    }

    variable_p_vals_core = dict()
    variable_p_vals_mantle = dict()
    variable_fO2_vals_core = dict()
    variable_fO2_vals_mantle = dict()
    variable_p_vals_core_rel_iron = dict()
    variable_p_vals_mantle_rel_iron = dict()
    variable_fO2_vals_core_rel_iron = dict()
    variable_fO2_vals_mantle_rel_iron = dict()
    observed_earth_core_abundances = dict()
    observed_earth_mantle_abundances = dict()
    observed_earth_core_abundances_rel_iron = dict()
    observed_earth_mantle_abundances_rel_iron = dict()
    variable_p_cnf_core = list()
    variable_fO2_cnf_core = list()
    variable_Ds_dict = dict()
    variable_fcrit_dict = dict()
    variable_dlogX_dP_dict = dict()
    observed_cnfs = {'Earth': geo_model.earth_layer_number_fractions[gi.Layer.core]}
    
    for element in elements:
        variable_p_vals_core[element] = list()
        variable_p_vals_mantle[element] = list()
        variable_fO2_vals_core[element] = list()
        variable_fO2_vals_mantle[element] = list()
        variable_p_vals_core_rel_iron[element] = list()
        variable_p_vals_mantle_rel_iron[element] = list()
        variable_fO2_vals_core_rel_iron[element] = list()
        variable_fO2_vals_mantle_rel_iron[element] = list()
        observed_earth_core_abundances[element] = list()  # The next four lists don't really need to be lists!
        observed_earth_mantle_abundances[element] = list()
        observed_earth_core_abundances_rel_iron[element] = list()
        observed_earth_mantle_abundances_rel_iron[element] = list()
        variable_fcrit_dict[element] = list()

    dlogX_dP_pressures = [1, 54]
    for element in elements:
        if element is ci.Element.O:
            continue
        variable_dlogX_dP_dict[element] = dict()
        for p in dlogX_dP_pressures:
            variable_dlogX_dP_dict[element][p] = list()
            for fcf in generate_fcf_vals():
                variable_dlogX_dP_dict[element][p].append(dlogX_dP_dict[(p, earth_fO2)][fcf][element])

    for element in elements:
        try:
            ca = geo_model.get_core_abundance(element)
            ma = geo_model.get_mantle_abundance(element)
            if ca is not None:
                observed_earth_core_abundances[element].append(geo_model.get_core_abundance(element))
                observed_earth_core_abundances_rel_iron[element].append(geo_model.get_core_abundance(element)/geo_model.get_core_abundance(ci.Element.Fe))
            if ma is not None:
                observed_earth_mantle_abundances[element].append(geo_model.get_mantle_abundance(element))
                observed_earth_mantle_abundances_rel_iron[element].append(geo_model.get_mantle_abundance(element)/geo_model.get_mantle_abundance(ci.Element.Fe))
        except (KeyError, TypeError) as e:
            # Then we don't have information for this element
            pass

    major_elements_for_conv = [ci.Element.Hf, ci.Element.U, ci.Element.Ta, ci.Element.Pb, ci.Element.Nb, ci.Element.Si, ci.Element.Mn, ci.Element.Zn, ci.Element.Ga, ci.Element.V, ci.Element.Cr,
    ci.Element.Cu, ci.Element.Ti, ci.Element.Fe, ci.Element.W, ci.Element.P, ci.Element.Co, ci.Element.Ni, ci.Element.O, ci.Element.C, ci.Element.S]
    for key, Ds_list in planet_Ds_dict.items():
        variable_Ds_dict[key] = dict()
        for element in major_elements_for_conv:
            variable_Ds_dict[key][element] = list()
        for Ds in Ds_list:
            for element in major_elements_for_conv:
                variable_Ds_dict[key][element].append(Ds[element])

    for p in generate_pressure_vals():
        for element in elements:
            if planet_dict.get((p, earth_fO2), None) is not None:
                variable_p_vals_core[element].append(planet_dict[(p, earth_fO2)][element][gi.Layer.core])
                variable_p_vals_mantle[element].append(planet_dict[(p, earth_fO2)][element][gi.Layer.mantle])
                variable_p_vals_core_rel_iron[element].append(planet_dict[(p, earth_fO2)][element][gi.Layer.core]/planet_dict[(p, earth_fO2)][ci.Element.Fe][gi.Layer.core])
                variable_p_vals_mantle_rel_iron[element].append(planet_dict[(p, earth_fO2)][element][gi.Layer.mantle]/planet_dict[(p, earth_fO2)][ci.Element.Fe][gi.Layer.mantle])
                try:
                    variable_fcrit_dict[element].append(critical_fragment_core_fraction_dict[(p, earth_fO2)][element])
                except KeyError:
                    # ignore Oxygen
                    pass
            else:
                variable_p_vals_core[element].append(None)
                variable_p_vals_mantle[element].append(None)
                variable_p_vals_core_rel_iron[element].append(None)
                variable_p_vals_mantle_rel_iron[element].append(None)
                try:
                    variable_fcrit_dict[element].append(None)
                except KeyError:
                    # ignore Oxygen
                    pass
        variable_p_cnf_core.append(planet_cnf_dict.get((p, earth_fO2), None))
    for f in generate_fO2_vals():
        for element in elements:
            try:
                variable_fO2_vals_core[element].append(planet_dict[(earth_pressure, f)][element][gi.Layer.core])
                variable_fO2_vals_mantle[element].append(planet_dict[(earth_pressure, f)][element][gi.Layer.mantle])
                variable_fO2_vals_core_rel_iron[element].append(planet_dict[(earth_pressure, f)][element][gi.Layer.core]/planet_dict[(earth_pressure, f)][ci.Element.Fe][gi.Layer.core])
                variable_fO2_vals_mantle_rel_iron[element].append(planet_dict[(earth_pressure, f)][element][gi.Layer.mantle]/planet_dict[(earth_pressure, f)][ci.Element.Fe][gi.Layer.mantle])
            except KeyError as e:
                variable_fO2_vals_core[element].append(None)
                variable_fO2_vals_mantle[element].append(None)
        variable_fO2_cnf_core.append(planet_cnf_dict.get((earth_pressure, f), None))
    
    stacked_named_elements_combined = { # We want to bundle up minor elements for readability for the bulk plot
        gi.Layer.core: [ci.Element.Ni, ci.Element.Si, ci.Element.O, ci.Element.Fe],
        gi.Layer.mantle: [ci.Element.Ca, ci.Element.Al, ci.Element.Fe, ci.Element.Mg, ci.Element.Si, ci.Element.O]
    }
    stacked_named_elements = {
        gi.Layer.core: [ci.Element.Cr, ci.Element.C, ci.Element.Ni, ci.Element.Si, ci.Element.O, ci.Element.Fe],
        gi.Layer.mantle: [ci.Element.C, ci.Element.Cr, ci.Element.Ni, ci.Element.Ca, ci.Element.Al, ci.Element.Fe, ci.Element.Mg, ci.Element.Si, ci.Element.O]
    }
    
    stacked_variable_p_vals_for_single_plots, bulk_stacked_variable_p_vals_for_single_plots, stacked_variable_fO2_vals_for_single_plots, bulk_stacked_variable_fO2_vals_for_single_plots = extract_stacked_abundances(
        stacked_named_elements,
        variable_p_vals_core,
        variable_p_vals_mantle,
        variable_fO2_vals_core,
        variable_fO2_vals_mantle,
        planet_cnf_dict,
        earth_fO2,
        earth_pressure
    )
    stacked_variable_p_vals_for_bulk_plots, bulk_stacked_variable_p_vals_for_bulk_plots, stacked_variable_fO2_vals_for_bulk_plots, bulk_stacked_variable_fO2_vals_for_bulk_plots = extract_stacked_abundances(
        stacked_named_elements_combined,
        variable_p_vals_core,
        variable_p_vals_mantle,
        variable_fO2_vals_core,
        variable_fO2_vals_mantle,
        planet_cnf_dict,
        earth_fO2,
        earth_pressure
    )
    #planet_dict_sulf, planet_cnf_dict_sulf, planet_Ds_dict_sulf, sulf_vals = form_planets_with_variable_sulfur(geo_model)
    #
    #variable_s_vals_core = dict()
    #variable_s_vals_mantle = dict()
    ##all_elements = planet_dict_sulf[sulf_vals[0]].keys()
    #important_s_elements = [ci.Element.Si, ci.Element.Cr, ci.Element.Ni, ci.Element.Ti, ci.Element.O, ci.Element.C, ci.Element.Fe]
    #for element in important_s_elements:
    #    variable_s_vals_core[element] = list()
    #    variable_s_vals_mantle[element] = list()
    #for element in important_s_elements:
    #    for s in sulf_vals:
    #        try:
    #            variable_s_vals_core[element].append(planet_dict_sulf[s][element][gi.Layer.core])
    #        except KeyError:
    #            variable_s_vals_core[element].append(None)
    #        try:
    #            variable_s_vals_mantle[element].append(planet_dict_sulf[s][element][gi.Layer.mantle])
    #        except KeyError:
    #            variable_s_vals_mantle[element].append(None)
    
    run_tag = geo_model.config_name
    
    graph_fac = gf.GraphFactory()
    #for ktp in [(54, -2), (66, -2), (67, -2), (68, -2), (69, -2), (70, -2), (71, -2), (72, -2), (73, -2), (74, -2), (75, -2), (90, -2), (154, -2)]:
    #    tag = str(ktp[0]) + '_' + str(ktp[1])
    #    try:
    #        graph_fac.make_convergence_plot(variable_Ds_dict[ktp], run_tag + '_' + tag, colour_dict)
    #    except KeyError as e:
    #        print(e)
    #        print(tag)
    #        print('Key not found, moving on to next plot')
    s_tag = 'no_sulf'
    #graph_fac.make_planet_formation_plot(variable_p_vals_core, 'Pressure', generate_pressure_vals(), 'core', observed_earth_core_abundances, run_tag + s_tag)
    #graph_fac.make_planet_formation_plot(variable_p_vals_mantle, 'Pressure', generate_pressure_vals(), 'mantle', observed_earth_mantle_abundances, run_tag + s_tag)
    #graph_fac.make_planet_formation_plot(variable_fO2_vals_core, 'fO2', generate_fO2_vals(), 'core', observed_earth_core_abundances, run_tag + s_tag)
    #graph_fac.make_planet_formation_plot(variable_fO2_vals_mantle, 'fO2', generate_fO2_vals(), 'mantle', observed_earth_mantle_abundances, run_tag + s_tag)
    #graph_fac.make_planet_cnf_plot(variable_p_cnf_core, 'Pressure', generate_pressure_vals(), observed_cnfs)
    #graph_fac.make_planet_cnf_plot(variable_fO2_cnf_core, 'fO2', generate_fO2_vals(), observed_cnfs)
    
    #graph_fac.make_planet_formation_plot(variable_p_vals_core_rel_iron, 'Pressure', generate_pressure_vals(), 'core', observed_earth_core_abundances_rel_iron, run_tag + s_tag, 'Iron')
    #graph_fac.make_planet_formation_plot(variable_p_vals_mantle_rel_iron, 'Pressure', generate_pressure_vals(), 'mantle', observed_earth_mantle_abundances_rel_iron, run_tag + s_tag, 'Iron')
    
    
    #graph_fac.make_earth_composition_plot(observed_earth_core_abundances, {'Model (' + str(earth_pressure) + ' GPa, IW' + str(earth_fO2) + ')': planet_dict[(earth_pressure, earth_fO2)]}, gi.Layer.core, 'core' + run_tag + s_tag)
    #graph_fac.make_earth_composition_plot(observed_earth_mantle_abundances, {'Model (' + str(earth_pressure) + 'GPa, IW' + str(earth_fO2) + ')': planet_dict[(earth_pressure, earth_fO2)]}, gi.Layer.mantle, 'mantle' + run_tag + s_tag)
    
    #graph_fac.make_planet_formation_stacked_plot(stacked_variable_p_vals_for_single_plots[gi.Layer.core], 'Pressure', generate_pressure_vals(), 'core', observed_earth_core_abundances, colour_dict, run_tag + s_tag)
    #graph_fac.make_planet_formation_stacked_plot(stacked_variable_p_vals_for_single_plots[gi.Layer.mantle], 'Pressure', generate_pressure_vals(), 'mantle', observed_earth_mantle_abundances, colour_dict, run_tag + s_tag)
    #graph_fac.make_planet_formation_combined_stacked_plot(bulk_stacked_variable_p_vals_for_bulk_plots, 'Pressure', generate_pressure_vals(), 'bulk', None, colour_dict, run_tag + s_tag)
    
    #graph_fac.make_planet_formation_stacked_plot(stacked_variable_fO2_vals_for_single_plots[gi.Layer.core], 'fO2', generate_fO2_vals(), 'core', observed_earth_core_abundances, colour_dict, run_tag + s_tag)
    #graph_fac.make_planet_formation_stacked_plot(stacked_variable_fO2_vals_for_single_plots[gi.Layer.mantle], 'fO2', generate_fO2_vals(), 'mantle', observed_earth_mantle_abundances, colour_dict, run_tag + s_tag)
    #graph_fac.make_planet_formation_combined_stacked_plot(bulk_stacked_variable_fO2_vals_for_bulk_plots, 'fO2', generate_fO2_vals(), 'bulk', None, colour_dict, run_tag + s_tag)
    
    #graph_fac.make_planet_formation_multipanel_stacked_plot(
    #    stacked_variable_p_vals_for_single_plots,
    #    bulk_stacked_variable_p_vals_for_single_plots,
    #    stacked_variable_fO2_vals_for_single_plots,
    #    bulk_stacked_variable_fO2_vals_for_single_plots,
    #    generate_pressure_vals(),
    #    generate_fO2_vals(),
    #    'multipanel',
    #    None,
    #    colour_dict,
    #    run_tag + s_tag
    #)
    
    #plot_wd_comparisons(planet_dict, run_tag) # Deprecated: Use plot_isolated_pressure_fO2_effect.py instead
    #print(planet_dict[(earth_pressure, earth_fO2)])
    #print(planet_cnf_dict[(earth_pressure, earth_fO2)])
    #print(planet_Ds_dict[(earth_pressure, earth_fO2)][-1:])
    #geo_model.tabulate_output(planet_dict[(earth_pressure, earth_fO2)], planet_cnf_dict[(earth_pressure, earth_fO2)], planet_Ds_dict[(earth_pressure, earth_fO2)][-1:][0], 'test_out')
    
    #print('Ni core')
    #print(variable_s_vals_core[ci.Element.Ni])
    #print('Ni mantle')
    #print(variable_s_vals_mantle[ci.Element.Ni])
    #print('Cr core')
    #print(variable_s_vals_core[ci.Element.Cr])
    #print('Cr mantle')
    #print(variable_s_vals_mantle[ci.Element.Cr])
    #graph_fac.make_planet_formation_plot(variable_s_vals_core, 'Sulfur Content', sulf_vals, str(gi.Layer.core), None, '54_-2')
    #graph_fac.make_planet_formation_plot(variable_s_vals_mantle, 'Sulfur Content', sulf_vals, str(gi.Layer.mantle), None, '54_-2')
    
    graph_fac.make_fcrit_plot(variable_fcrit_dict, generate_pressure_vals(), earth_fO2)
    graph_fac.make_dlogX_dP_plot(variable_dlogX_dP_dict, generate_fcf_vals(), earth_fO2)
    
if __name__ == '__main__':
    main()
