#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import chemistry_info as ci
import graph_factory as gf

# A newer, neater version of white_dwarf_model.py
# I use this version in the synthetic pipeline, and plan to use it in the main Bayesian code as well
# Problem being that this version is much slower for some reason that I haven't investigated yet (prohibitively slow for the Bayesian model) (TODO!)

# Function to take a set of planetesimal abundances and calculate their abundance in a wd atmosphere after a certain time
# Ideally we'd like planetesimal_abundance and wd_timescales to be numpy arrays but will cast from list if necessary
def process_abundances(t_sinceaccretion, t_disc, planetesimal_abundance, wd_timescales, pollution_level, snapshot=True):
    if isinstance(planetesimal_abundance, list):
        planetesimal_abundance = np.array(planetesimal_abundance)
    if isinstance(wd_timescales, list):
        wd_timescales = np.array(wd_timescales)
    if snapshot:
        return process_snapshot_abundances(
            t_sinceaccretion,
            t_disc,
            planetesimal_abundance,
            wd_timescales,
            pollution_level
        )
    else:
        return process_lifetime_abundances(
            t_sinceaccretion,
            t_disc,
            planetesimal_abundance,
            wd_timescales,
            pollution_level
        )

# Function to take a set of planetesimal abundances and calculate their abundance in a wd atmosphere after a certain time t_sinceaccretion is in Myr, t_disc is in yr
def process_snapshot_abundances(t_sinceaccretion, t_disc, planetesimal_abundance, wd_timescales, pollution_level):
    t_sinceaccretion_years = t_sinceaccretion*1000000
    buildup_scaling_factors = calculate_buildup_scaling_factors(t_sinceaccretion_years, t_disc, wd_timescales)
    sinkout_scaling_factors = calculate_sinkout_scaling_factors(t_sinceaccretion_years, t_disc, wd_timescales)
    proposed_abundances = planetesimal_abundance*buildup_scaling_factors*sinkout_scaling_factors
    # The summation is the bottleneck: bizarrely, manually summing like this seems to be faster than using the built in sum or np.sum functions
    summation = 0
    for pa in proposed_abundances:
        summation += pa
    normalised_abundances = proposed_abundances/summation
    toret = np.log10(normalised_abundances) + pollution_level
    return toret

def calculate_buildup_scaling_factors(t_sinceaccretion_years, t_disc, sinking_timescales): # Accounts for build-up while accretion ongoing
    # Sending in inverse timescales to avoid (slow) division operations - can just multiply by the inverse timescales instead. Means we only have to do the division once.
    if t_sinceaccretion_years <= 0:
        return np.ones_like(sinking_timescales)
    else:
        if isinstance(sinking_timescales, list):
            sinking_timescales = np.array(sinking_timescales)
        return sinking_timescales*(1 - np.exp(-(min(t_sinceaccretion_years, t_disc)/sinking_timescales)))

def calculate_sinkout_scaling_factors(t_sinceaccretion_years, t_disc, sinking_timescales): # Accounts for sinking out after accretion ends
    if t_sinceaccretion_years <= 0:
        return np.ones_like(sinking_timescales)
    else:
        if isinstance(sinking_timescales, list):
            sinking_timescales = np.array(sinking_timescales)
        return np.exp(min(t_disc - t_sinceaccretion_years, 0)/sinking_timescales)


#def process_snapshot_abundances(t_sinceaccretion, t_disc, planetesimal_abundance, wd_timescales, pollution_level):
#    t_sinceaccretion_years = t_sinceaccretion*1000000
#    proposed_abundances = list()
#    sum_of_abundances = 0
#    for i, pa in enumerate(planetesimal_abundance):
#        buildup_scaling_factor = calculate_buildup_scaling_factor(t_sinceaccretion_years, t_disc, wd_timescales[i])
#        sinkout_scaling_factor = calculate_sinkout_scaling_factor(t_sinceaccretion_years, t_disc, wd_timescales[i])
#        proposed_abundance = pa*buildup_scaling_factor*sinkout_scaling_factor
#        sum_of_abundances += proposed_abundance
#        proposed_abundances.append(proposed_abundance)
#    for i, pa in enumerate(proposed_abundances):
#        proposed_abundances[i] = pa/sum_of_abundances
#    proposed_toret = list()
#    for pa in proposed_abundances:
#        proposed_toret.append(np.log10(pa) + pollution_level)
#    return proposed_toret
#
#def calculate_buildup_scaling_factor(t_sinceaccretion_years, t_disc, sinking_timescale): # Accounts for build-up while accretion ongoing
#    if t_sinceaccretion_years <= 0:
#        return 1
#    else:
#        return sinking_timescale*(1 - np.exp(-(min(t_sinceaccretion_years, t_disc)/sinking_timescale)))
#
#def calculate_sinkout_scaling_factor(t_sinceaccretion_years, t_disc, sinking_timescale): # Accounts for sinking out after accretion ends
#    if t_sinceaccretion_years <= 0:
#        return 1
#    else:
#        return np.exp(min(t_disc - t_sinceaccretion_years, 0)/sinking_timescale)


# Function to take a set of planetesimal abundances and calculate the total amount of each element that entered a wd atmosphere prior to the current time
# This should return equivalent X:HorHe ratios which will be larger than the actual X:HorHe values (to capture all the material we no longer observe)
def process_lifetime_abundances(t_sinceaccretion, t_disc, planetesimal_abundance, wd_timescales, pollution_level):
    t_sinceaccretion_years = t_sinceaccretion*1000000
    if t_sinceaccretion_years == 0:
        compensation_factor = 1
    elif 0 < t_sinceaccretion_years < t_disc:
        compensation_factor = (t_sinceaccretion_years/wd_timescales) + np.exp(-t_sinceaccretion_years/wd_timescales)
    else:
        compensation_factor = (t_disc/wd_timescales) + np.exp((t_sinceaccretion_years-t_disc)/wd_timescales) + np.exp(-t_disc/wd_timescales)
    return np.log10(planetesimal_abundance*compensation_factor) + pollution_level  # Assuming this is going to tell us the proportion of the pollution fraction corresponding to each element

def process_absolute_abundance(element, atmosphere_element, sinking_timescale, accretion_rate, mass_cvz, t_obs, t_event): # All times in yr, except t_obs (Myr)
    # TODO: Should probably rewire the code to use this function, it's more versatile and exchanges pollution fraction for masses
    # The unit of accretion_rate is effectively M yr^-1 where M is the same units as mass_cvz
    # mass_cvz could be calculated using timescale_interpolator
    # accretion_rate would be a new parameter in place of pollution fraction
    # This should return the number abundance relative to Hx
    t_sinceaccretion_years = t_obs*1000000
    t_limit = min(t_sinceaccretion_years, t_event)
    prefactor = (ci.get_element_mass(atmosphere_element)*accretion_rate*sinking_timescale)/(ci.get_element_mass(element)*mass_cvz)
    factor1 = np.exp(-t_sinceaccretion_years/sinking_timescale)
    factor2 = np.exp(t_limit/sinking_timescale) - 1
    toret = prefactor*factor1*factor2
    return toret

def plot_accretion_example():
    el_X = ci.Element.Fe
    el_Y = ci.Element.Ca
    t_X = 0.5
    t_Y = 1
    planetesimal_abundances = [1, 1]
    wd_timescales = [t_X, t_Y]
    t_event = 15
    t_ss = 5*max(t_X, t_Y)
    min_t = 0
    max_t = t_event + t_ss
    number_of_steps = 1000
    t_range = np.linspace(min_t, max_t, number_of_steps + 1)
    accretion_rate = 25
    mass_cvz = 1

    N_X = [process_absolute_abundance(ci.Element.Fe, ci.Element.He, t_X, accretion_rate, mass_cvz, t/1000000, t_event) for t in t_range]
    N_Y = [process_absolute_abundance(ci.Element.Fe, ci.Element.He, t_Y, accretion_rate, mass_cvz, t/1000000, t_event) for t in t_range]

    pol_frac = -6 # irrelevant
    test_abundances = [process_abundances(t/1000000, t_event, planetesimal_abundances, wd_timescales, pol_frac, True) for t in t_range]
    test_ratios = [10**(test_abundance[0] - test_abundance[1]) for test_abundance in test_abundances]

    graph_fac = gf.GraphFactory()
    graph_fac.plot_accretion_phase_demo(t_range, test_ratios, N_X, N_Y)

def main():
    plot_accretion_example()

if __name__ == '__main__':
    main()
