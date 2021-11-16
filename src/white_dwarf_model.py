#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

# Function to take a set of planetesimal abundances and calculate their abundance in a wd atmosphere after a certain time
def process_abundances(t_sinceaccretion, t_disc, planetesimal_abundance, non_zero_planetesimal_abundance, wd_timescales, non_zero_wd_timescales, pollutionfraction, snapshot=True):
    if snapshot:
        return process_snapshot_abundances(
            t_sinceaccretion,
            t_disc,
            planetesimal_abundance,
            non_zero_planetesimal_abundance,
            wd_timescales,
            non_zero_wd_timescales,
            pollutionfraction
        )
    else:
         return process_lifetime_abundances(
            t_sinceaccretion,
            t_disc,
            planetesimal_abundance,
            wd_timescales,
            pollutionfraction
        )
    
# Function to take a set of planetesimal abundances and calculate their abundance in a wd atmosphere after a certain time    
def process_snapshot_abundances(t_sinceaccretion, t_disc, planetesimal_abundance, non_zero_planetesimal_abundance, wd_timescales, non_zero_wd_timescales, pollutionfraction):
    t_sinceaccretion_years = t_sinceaccretion*1000000
    if t_sinceaccretion_years == 0:
        numerators = planetesimal_abundance
        denominators = non_zero_planetesimal_abundance
    elif 0 < t_sinceaccretion_years < t_disc:
        exp_term = (1-np.exp(-(t_sinceaccretion_years)/wd_timescales))
        exp_termd = (1-np.exp(-(t_sinceaccretion_years)/non_zero_wd_timescales))
        numerators = wd_timescales*planetesimal_abundance*exp_term
        denominators = non_zero_wd_timescales*non_zero_planetesimal_abundance*exp_termd
    else:
        exp_term1 = np.exp(-((t_sinceaccretion_years)-(t_disc))/wd_timescales)
        exp_term2 = (1-np.exp(-(t_disc)/wd_timescales))
        exp_term1d = np.exp(-((t_sinceaccretion_years)-(t_disc))/non_zero_wd_timescales)
        exp_term2d = (1-np.exp(-(t_disc)/non_zero_wd_timescales))
        numerators = wd_timescales*planetesimal_abundance*exp_term1*exp_term2
        denominators = non_zero_wd_timescales*non_zero_planetesimal_abundance*exp_term1d*exp_term2d
    premodel = numerators/np.sum(denominators)
    model = np.log10(premodel)
    return model + pollutionfraction

# Function to take a set of planetesimal abundances and calculate the total amount of each element that entered a wd atmosphere prior to the current time
# This should return equivalent X:HorHe ratios which will be larger than the actual X:HorHe values (to capture all the material we no longer observe)
def process_lifetime_abundances(t_sinceaccretion, t_disc, planetesimal_abundance, wd_timescales, pollutionfraction):
    t_sinceaccretion_years = t_sinceaccretion*1000000
    if t_sinceaccretion_years == 0:
        compensation_factor = 1
    elif 0 < t_sinceaccretion_years < t_disc:
        compensation_factor = (t_sinceaccretion_years/wd_timescales) + np.exp(-t_sinceaccretion_years/wd_timescales)
    else:
        compensation_factor = (t_disc/wd_timescales) + np.exp((t_sinceaccretion_years-t_disc)/wd_timescales) + np.exp(-t_disc/wd_timescales)
    return np.log10(planetesimal_abundance*compensation_factor) + pollutionfraction  # Assuming this is going to tell us the proportion of the pollution fraction corresponding to each element
