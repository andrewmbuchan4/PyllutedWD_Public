#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum
import numpy as np

import disc_model as dm
import live_data as ld
import model_parameters as mp

class Limit(Enum):
    Lower = 0
    Upper  = 1

# The order of parameters matters: some have priors which depend on the value of previously assigned parameters
parameters_to_cycle_through = [
    mp.ModelParameter.metallicity,
    mp.ModelParameter.formation_distance,
    mp.ModelParameter.feeding_zone_size,
    mp.ModelParameter.parent_core_frac,
    mp.ModelParameter.parent_crust_frac,
    mp.ModelParameter.fragment_core_frac,
    mp.ModelParameter.fragment_crust_frac,
    mp.ModelParameter.pollution_frac,
    mp.ModelParameter.accretion_timescale,
    mp.ModelParameter.t_sinceaccretion,
    mp.ModelParameter.pressure,
    mp.ModelParameter.oxygen_fugacity
]

def get_fragment_crust_frac_upper_limit(cube, parameter_indices):
    return 1 - cube[parameter_indices[mp.ModelParameter.fragment_core_frac]]

def get_pollution_frac_lower_limit(cube, parameter_indices):
    return np.log10(np.sum(10**(ld._live_non_zero_wd_abundances))) - 1.5
    
def get_pollution_frac_upper_limit(cube, parameter_indices):
    return np.log10(np.sum(10**(ld._live_non_zero_wd_abundances))) + 1.5
    
def get_t_sinceaccretion_upper_limit(cube, parameter_indices):
    return ((12*ld._live_t_mg)+(10**(cube[parameter_indices[mp.ModelParameter.accretion_timescale]])))/1000000

prior_limits_dict = {
    'Default': {
        mp.ModelParameter.metallicity: {
            Limit.Lower: 0,
            Limit.Upper: 958
        },
        mp.ModelParameter.formation_distance: {
            Limit.Lower: -2,
            Limit.Upper: np.log10(dm.S_disc(1.5))
        },
        mp.ModelParameter.feeding_zone_size: {
            Limit.Lower: 0,
            Limit.Upper: 0.15
        },
        mp.ModelParameter.parent_core_frac: {
            Limit.Lower: 0,
            Limit.Upper: 0.19
        },
        mp.ModelParameter.parent_crust_frac: {
            Limit.Lower: 0,
            Limit.Upper: 0.25
        },
        mp.ModelParameter.fragment_core_frac: {
            Limit.Lower: 0,
            Limit.Upper: 1
        },
        mp.ModelParameter.fragment_crust_frac: {
            Limit.Lower: 0,
            Limit.Upper: get_fragment_crust_frac_upper_limit
        },
        mp.ModelParameter.pollution_frac: {
            Limit.Lower: get_pollution_frac_lower_limit,
            Limit.Upper: get_pollution_frac_upper_limit
        },
        mp.ModelParameter.accretion_timescale: {
            Limit.Lower: 0,
            Limit.Upper: 8
        },
        mp.ModelParameter.t_sinceaccretion: {
            Limit.Lower: 0,
            Limit.Upper: get_t_sinceaccretion_upper_limit
        },
        mp.ModelParameter.pressure: {
            Limit.Lower: 0,
            Limit.Upper: 60
        },
        mp.ModelParameter.oxygen_fugacity: {
            Limit.Lower: -3,
            Limit.Upper: -1
        }
    },
    'TestPrior': {
        mp.ModelParameter.formation_distance: {
            Limit.Lower: -1.5
        }
    },
    'HighPressure': {
        mp.ModelParameter.pressure: {
            Limit.Lower: 45,
            Limit.Upper: 60
        }
    },
    'RaisedPressure': {
        mp.ModelParameter.pressure: {
            Limit.Lower: 15,
            Limit.Upper: 60
        }
    },
    'LowPressure': {
        mp.ModelParameter.pressure: {
            Limit.Lower: 0,
            Limit.Upper: 15
        }
    }
}


def universal_prior(cube):
    
    parameter_indices = mp.parameter_indices(ld._live_model)
    
    for parameter in parameters_to_cycle_through:
        if mp.model_uses_parameter(ld._live_model, parameter):
            try:
                upper_limit_raw = prior_limits_dict[ld._live_prior][parameter][Limit.Upper]
            except KeyError:
                upper_limit_raw = prior_limits_dict['Default'][parameter][Limit.Upper]
            try:
                lower_limit_raw = prior_limits_dict[ld._live_prior][parameter][Limit.Lower]
            except KeyError:
                lower_limit_raw = prior_limits_dict['Default'][parameter][Limit.Lower]
            try:
                # Some limits are given by executing a function - if the dictionary entry is callable, we should call it
                upper_limit = upper_limit_raw(cube, parameter_indices)
            except TypeError:
                upper_limit = upper_limit_raw
            try:
                lower_limit = lower_limit_raw(cube, parameter_indices)
            except TypeError:
                lower_limit = lower_limit_raw
            scaling_term = upper_limit - lower_limit
            cube[parameter_indices[parameter]] = (scaling_term*cube[parameter_indices[parameter]]) + lower_limit
    
    return cube
