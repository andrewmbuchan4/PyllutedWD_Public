#!/usr/bin/env python

from enum import Enum
import numpy as np

import chemistry_info as ci
import disc_model as dm
import live_data as ld
import model_parameters as mp
import physical_constants as pc

class Limit(Enum):
    Lower = 0
    Upper = 1

# The order of parameters matters: some have priors which depend on the value of previously assigned parameters
parameters_to_cycle_through = [
    mp.ModelParameter.metallicity,
    mp.ModelParameter.formation_distance,
    mp.ModelParameter.feeding_zone_size,
    mp.ModelParameter.parent_core_frac,
    mp.ModelParameter.parent_crust_frac,
    mp.ModelParameter.fragment_core_frac,
    mp.ModelParameter.fragment_crust_frac,
    mp.ModelParameter.fragment_mass,
    mp.ModelParameter.accretion_timescale,
    mp.ModelParameter.t_sinceaccretion,
    mp.ModelParameter.pressure,
    mp.ModelParameter.oxygen_fugacity
]

def get_fragment_crust_frac_upper_limit(cube, parameter_indices):
    return 1 - cube[parameter_indices[mp.ModelParameter.fragment_core_frac]]

# The pollution fraction functions are deprecated but here for posterity I guess
def get_pollution_frac_lower_limit(cube, parameter_indices):
    return ld._live_white_dwarf.estimate_minimum_pollution_fraction(ci.usual_elements) - 0.5 # you can't be too far below this limit!

def get_pollution_frac_upper_limit(cube, parameter_indices):
    return ld._live_white_dwarf.estimate_maximum_pollution_fraction(ci.usual_elements) + 1.5 # The 1.5 is a safety factor: your reference element could happen to be underabundant. need a bigger margin at the top end

def get_t_sinceaccretion_upper_limit(cube, parameter_indices):
    return ((12*ld._live_t_mg)+(10**(cube[parameter_indices[mp.ModelParameter.accretion_timescale]])))/1000000  # NB the variable we're indexing into here is the accretion timescale, not the time since accretion!

#def get_feeding_zone_size_upper_limit(cube, parameter_indices):
    #
    #return cube[parameter_indices[mp.ModelParameter.accretion_timescale]]

def get_minimum_distance():
    min_distance_AU = 0.0955 # Roughly the distance where all the abundances start dropping to zero rapidly. Weird behaviour happens further in than this!
    return np.log10(min_distance_AU)

prior_limits_dict = {
    'Default': {
        mp.ModelParameter.metallicity: {
            Limit.Lower: 0,
            Limit.Upper: 958
        },
        mp.ModelParameter.formation_distance: {
            #Limit.Lower: -2,
            Limit.Lower: get_minimum_distance(),
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
        #mp.ModelParameter.pollution_frac: { # This should perhaps not be uniform but skewed towards low values (as in the synthetic population code)
        #    Limit.Lower: get_pollution_frac_lower_limit,
        #    Limit.Upper: get_pollution_frac_upper_limit
        #},
        #mp.ModelParameter.fragment_mass: { # Putting this on a log scale, so this is log_10 kg
        #    Limit.Lower: 10,
        #    Limit.Upper: 25  # TODO: This is probably not a very realistic prior! Maybe look at collisional cascade distributions?
        #},
        #mp.ModelParameter.fragment_mass: { # this is in kg
        #    Limit.Lower: 1E10, # Loosely based on the lower mass asteroids John considered. Might want to make this higher though really!
        #    Limit.Upper: pc.M_Earth,
        #    'Power': 11/6 # Collisional cascade (Dohnanyi 1969)
        #},
        mp.ModelParameter.fragment_mass: { # Putting this on a log scale, so this is log_10 kg
            Limit.Lower: 13.5, # Bringing this in line with the synthetic pipeline - no particular significance to this though. Just a necessary truncation
            Limit.Upper: np.log10(pc.M_Earth),
            'Power': 11/6,
            'Log': True
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
    },
    'NarrowTdisc': {
        mp.ModelParameter.accretion_timescale: {
            Limit.Lower: 4,
            Limit.Upper: 7
        }
    },
    'LongTdisc': {
        mp.ModelParameter.accretion_timescale: {
            Limit.Lower: 5,
            Limit.Upper: 8
        }
    },
    'RestrictedFormation': {
        mp.ModelParameter.formation_distance: {
            Limit.Lower: get_minimum_distance(),
            Limit.Upper: np.log10(dm.S_disc(1.5))
        }
        #mp.ModelParameter.feeding_zone_size: { # May also want this to be a function of the formation distance?
        #    Limit.Lower: 0, # Arguably, this shouldn't be allowed? Or, if it is allowed, it should be the default value in model_parameters.
        #    Limit.Upper: get_feeding_zone_size_upper_limit
        #}
    }
}

def universal_prior(cube):
    # Including a more generalised scaling for non-uniform distributions
    # In general, this transformation is the inverse of the cumulative distribution
    # For a variable m, distributed according to dn = m^-alpha dm between m1 and m2
    # According to my maths
    # the inverse cumulative function at a particular value of m=m_test is
    # m_test = ((ncum,m_test*(m2^beta - m1^beta)) + m1^beta)^(1/beta)
    # where ncum,m_test is the cumulative count at m=m_test and beta = 1 - alpha
    # Sanity check: when alpha = 0 (i.e. a uniform distribution) this reduces to
    # m_test = (ncum,m_test*(m2 - m1)) + m1
    # which is the standard uniform distribution transformation
    # with ncum,m_test being the pre-transform cube value and m_test being the post-transform value
    # So I think the transformation we need is
    # cube[index] = ((cube[index]*(m2^beta - m1^beta)) + m1^beta)^(1/beta)
    # which is a generalisation of what we already have, but with some extra scalings

    # Then we can take it a step further and suppose that we are dealing with the log of a variable
    # which is distributed according to dn = m^-alpha dm
    # So we have a variable l, where l = log10(m) and
    # cube[index] = log10((cube[index]*((10**l2)^beta - (10**l1)^beta) + (10**l1)^beta)^(1/beta))

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
                power_raw = prior_limits_dict[ld._live_prior][parameter]['Power'] # This will play the role of alpha
            except KeyError:
                power_raw = 0
            try:
                # Some limits are given by executing a function - if the dictionary entry is callable, we should call it
                upper_limit = upper_limit_raw(cube, parameter_indices)
            except TypeError:
                upper_limit = upper_limit_raw
            try:
                lower_limit = lower_limit_raw(cube, parameter_indices)
            except TypeError:
                lower_limit = lower_limit_raw
            try:
                power = power_raw(cube, parameter_indices)
            except TypeError:
                power = power_raw
            try:
                log_scale = prior_limits_dict[ld._live_prior][parameter]['Log']
            except KeyError:
                log_scale = False
            beta = 1 - power
            if log_scale:
                linear_upper_limit = 10**upper_limit
                scaled_upper_limit = linear_upper_limit**beta
                linear_lower_limit = 10**lower_limit
                scaled_lower_limit = linear_lower_limit**beta
                scaling_term = scaled_upper_limit - scaled_lower_limit
                cube[parameter_indices[parameter]] = np.log10(((scaling_term*cube[parameter_indices[parameter]]) + scaled_lower_limit)**(1/beta))
            else:
                scaled_upper_limit = upper_limit**beta
                scaled_lower_limit = lower_limit**beta
                scaling_term = scaled_upper_limit - scaled_lower_limit
                cube[parameter_indices[parameter]] = ((scaling_term*cube[parameter_indices[parameter]]) + scaled_lower_limit)**(1/beta)
    return cube
