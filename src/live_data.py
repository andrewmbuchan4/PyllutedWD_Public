#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np

# This file is a bit of a hack. It exists because we need to pass external information into the prior function,
# but the prior function can only contain one argument (cube) otherwise pymultinest complains.
# So the idea here is that the manager will update the following variables, which can then be accessed by the prior
# TODO: Make a more elegant design that doesn't have these pseudo-global variables!
# Maybe some sort of static function in the Manager?
# Or maybe make the pollution model class into a subclass of a pymultinest solver, then use self.wd_abundances etc
# Then you still have to 'publish' the data but at least it's contained in an object
# Or some hybrid solution where the pollution model subclass has a self.live_data which is given a LiveData object containing the live data

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
