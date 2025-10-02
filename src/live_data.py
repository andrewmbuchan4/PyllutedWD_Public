#!/usr/bin/env python

import csv
import numpy as np

# This file is a bit of a hack. It exists because we need to pass external information into the prior function,
# but the prior function can only contain one argument (cube) otherwise pymultinest complains.
# So the idea here is that the manager will update the following variables, which can then be accessed by the prior

_live_model = None
_live_prior = None
_live_all_wd_errors = None
_live_all_wd_abundances = None
_live_t_mg = None
_live_stellar_compositions = None
#_live_non_zero_wd_errors = None
#_live_non_zero_wd_abundances = None
#_live_non_zero_wd_timescales = None
_live_all_wd_timescales = None
_geo_model = None
_live_elements_present = None
_live_q = None
_live_mass = None
_live_type = None
_live_white_dwarf = None
_live_timescale_type = None

# These get set by publish_live_model in manager.py
_live_model = None
_live_prior = None
_live_enhancement_model = None
_live_consider_thermohaline = None

# These last few variables are particularly egregious - they're ultimately here because the complete_model function needs to work for both Bayesian and Synthetic code
# From the point of view of the Bayesian code, we would ideally just send the white dwarf object in as an argument
# But in the Synthetic code, there is no white dwarf. Setting these global variables to the relevant values is a compromise

_live_Hx = None
_live_M_cvz = None
_live_teff = None
_live_logg = None
