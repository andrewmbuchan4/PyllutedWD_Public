#!/usr/bin/env python

# A newer, neater version of white_dwarf_model.py which uses mass

import numpy as np
import sys

import chemistry_info as ci
import physical_constants as pc
import pwd_utils as pu
import thermohaline_interpolator as thi

therm_factor_interpolator = thi.ThermohalineInterpolator()

def calculate_abundance_by_number(t, t_event, Hx, M_cvz, mu_X, M_X, tau_X, consider_thermohaline=False, Teff=None, logg=None):
    critical_threshold = np.log10(0.5) # This must be negative for the logic of the function to work
    mu_Hx = ci.get_element_mass(Hx)
    therm_factor = 0
    thermohaline_regime = False
    declining_phase = t > t_event
    if Hx != ci.Element.H:
        consider_thermohaline = False
    if consider_thermohaline:
        # This doesn't necessarily mean will we will be in a thermohaline regime, it just means we will consider this possibility
        unit_conversion_factor = 1000*pc.M_Sun/pc.seconds_per_year
        log_Mdot = np.log10(unit_conversion_factor*np.sum(M_X)/t_event)
        therm_factor = therm_factor_interpolator((log_Mdot, Teff, logg)) # A thermohaline mixing correction factor
        if np.isnan(therm_factor) or therm_factor is None:
            therm_factor = 0 # In this case, assume no thermohaline mixing - seems like the most conservative thing to do, since the most likely cause of this happening is low temperature.
        if therm_factor < critical_threshold: # Note the therm_factor gets more negative the more relevant thermohaline mixing is
            thermohaline_regime = True
    # Basically this is equation 1.8 in my thesis, with M_dot expressed as M/t_event, and reformulated to avoid enormous exponentials as detailed below
    # mu_X - expecting a numpy array of ci.get_element_mass(X)
    # M_X - expecting a numpy array of the mass of each element in the accreting body (what units? doesn't matter as long as it's the same as M_cvz)
    # tau_X - expecting a numpy array of sinking timescales for each element (years)
    # M_cvz - mass of convective zone (what units? doesn't matter as long as it's the same as M_X)
    # t_event - accretion timescale event (years)
    # What units is this returning? It's number ratio relative to Hx (on a log10 scale)

    #The linear space version of this is of the form:
    #N = P * e^A * (e^B-1)

    #so N = P * e^A * e^B - P * e^A
    #let Z = P * e^A * e^B, Y = P * e^A s.t. N = Z - Y
    #then

    #ln(Z) = ln(P) + A + B
    #ln(Y) = ln(P) + A

    #so we calculate these terms, then do

    #N = e^ln(Z) - e^ln(Y)

    #And hopefully we avoid doing any wild exponentials. Main point is to dodge any e^A or e^B, because they can be extreme

    if thermohaline_regime and declining_phase:
        t_to_use = t_event
    else:
        t_to_use = t

    prefactor = (mu_Hx*M_X*tau_X)/(mu_X*M_cvz*t_event)

    t_l = min(t_to_use, t_event) # The 'limiting' time
    lnP = np.log(prefactor)

    A = -t_to_use/tau_X
    B = t_l/tau_X

    lnZ = lnP + A + B
    lnY = lnP + A

    N_X = np.exp(lnZ) - np.exp(lnY)
    uncorrected_N_X = np.log10(N_X)
    corrected_N_X = uncorrected_N_X + therm_factor # therm_factor is usually negative, so the corrected_N_X is smaller.

    if thermohaline_regime:
        total_N_X = np.sum(10**corrected_N_X)
        no_differential_sinking_N_X = M_X/mu_X # Without differential sinking, the number abundances are in proportion to M_X, but scaled by mu_X to switch from mass to number abundances
        scaling_factor = total_N_X/np.sum(no_differential_sinking_N_X)
        rescaled_N_X = scaling_factor*no_differential_sinking_N_X
        logscale_rescaled_N_X = np.log10(rescaled_N_X)
        if not declining_phase:
            return logscale_rescaled_N_X
        else:
            declining_phase_correction_factor = (t-t_event)/(np.log(10)*tau_X) # correction factor to account for declining phase
            logscale_rescaled_N_X -= declining_phase_correction_factor
            return logscale_rescaled_N_X
    else:
        return corrected_N_X

# The next couple of functions are basically deprecated and only used in the synthetic modeller
def calculate_buildup_scaling_factors(t_sinceaccretion_years, t_disc, sinking_timescales):
    if t_sinceaccretion_years <= 0:
        return np.ones_like(sinking_timescales)
    else:
        if isinstance(sinking_timescales, list):
            sinking_timescales = np.array(sinking_timescales)
        return sinking_timescales*(1 - np.exp(-(min(t_sinceaccretion_years, t_disc)/sinking_timescales)))

def calculate_sinkout_scaling_factors(t_sinceaccretion_years, t_disc, sinking_timescales):
    if t_sinceaccretion_years <= 0:
        return np.ones_like(sinking_timescales)
    else:
        if isinstance(sinking_timescales, list):
            sinking_timescales = np.array(sinking_timescales)
        return np.exp(min(t_disc - t_sinceaccretion_years, 0)/sinking_timescales)

def example():
    Hx = ci.Element.H
    mu_X = np.array([6, 8, 10])
    M_X = np.array([0.000000000000000001, 0.00000000000000001, 0.00000000000000002]) # This is effectively now in Solar masses
    tau_X = np.array([31, 11, 17])
    M_cvz = 0.000000000000003 # This needs to also be in Solar masses now
    t = 11000
    t_event = 10000
    Teff = 12101
    logg = 7.95
    N_X = calculate_abundance_by_number(t, t_event, Hx, M_cvz, mu_X, M_X, tau_X, True, Teff, logg)
    print(N_X)

def main():
    example()

if __name__ == '__main__':
    main()
