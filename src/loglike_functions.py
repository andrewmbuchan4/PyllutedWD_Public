#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import complete_model as cm
import live_data as ld
import model_parameters as mp

def evaluate_log_likelihood(model_result, min_likelihood=mp.minimum_likelihood):
    if model_result is None:
        #like = min_likelihood  # For testing purposes only
        like = 1.1*min_likelihood  # ...So these points should be ignored (no information gained)
    else:
        like = ld._live_white_dwarf.log_likelihood(model_result, min_likelihood)
        # First check if we violate any bounds:
        #for el, el_result in model_result.items():
        #    upper_bound = ld._live_upper_bounds[el]
        #    if upper_bound is not None and el_result >= upper_bound:
        #        #return min_likelihood   # For testing purposes only
        #        return 0.9*min_likelihood   # ...So these points should not be ignored
        #    lower_bound = ld._live_lower_bounds[el]
        #    if lower_bound is not None and el_result <= lower_bound:
        #        #return min_likelihood   # For testing purposes only
        #        return 0.9*min_likelihood   # ...So these points should not be ignored
        #non_zero_model_values = list()
        #for el in ld._live_elements_present:
        #    non_zero_model_values.append(model_result[el])
        #model_values = np.asarray(non_zero_model_values)
        #neg = 1/((2*np.pi)*((ld._live_non_zero_wd_errors)**2))
        #chi2 = ((ld._live_non_zero_wd_abundances - model_values)**2)*(neg)*2*np.pi
        #like = -0.5*np.sum(chi2 - np.log(neg))
        #if like == -np.inf:
        #    #like = min_likelihood # For testing purposes only
        #    like = 0.9*min_likelihood

        #norm = -0.5*len(ld._live_elements_present)*np.log(2.*np.pi) - len(ld._live_elements_present)*np.log(ld._live_non_zero_wd_errors)

        # chi-squared
        #chisq = np.sum(((ld._live_non_zero_wd_abundances - model_values)/(ld._live_non_zero_wd_errors))**2)

        #alt_like = bonus_norm - 0.5*bonus_chisq

    return like

def universal_loglike(cube):

    parameter_indices = mp.parameter_indices(ld._live_model)

    input_values = list()
    for param in mp.get_model_params_in_order():
        to_append = cube[parameter_indices[param]] if mp.model_uses_parameter(ld._live_model, param) else mp.default_values[ld._enhancement_model][param]
        input_values.append(to_append)

    # Could rethink the structure of this part. (call model.execute() and go into this function: ?)

    # Another (possibly major) improvement to make is that the relative abundances are usually better known
    # than the absolute - in other words, there's a correlation between the errors on the abundances
    # that is not being taken into account, and we're making fitting the data unnecessarily hard on the model
    # The chi squared calculation below should really take the error correlation into account somehow
    # Or the whole calculation could be restructured so that we fit the absolute abundance of e.g. Ca/Hx,
    # and then the abundances of everything else relative to Ca (error on X/Ca should be smaller than that of X/Hx)

    model_result, diagnostics_for_post_processing = cm.complete_model_calculation(
        input_values[0], #fe_star
        input_values[1], #t_sinceaccretion
        input_values[2], #d_formation
        input_values[3], #z_formation
        input_values[4], #N_c
        input_values[5], #N_o
        input_values[6], #f_c
        input_values[7], #f_o
        input_values[8], #pollutionfraction
        10**input_values[9], #t_disc
        input_values[10], #pressure
        input_values[11], #fO2 (oxygen fugacity relative to Iron Wuestite buffer, in log units)
        ld._enhancement_model
    )

    min_likelihood = mp.minimum_likelihood              # Any points with likelihood less than this are ignored
                                                        # --> use for errors, but not for bounds violations
                                                        #( = -1e90 by default in the model, -1e100 by default in pymultinest)

    like = evaluate_log_likelihood(model_result, min_likelihood)
    return like

    #if model_result is None:
    #    #like = min_likelihood  # For testing purposes only
    #    like = 1.1*min_likelihood  # ...So these points should be ignored (no information gained)
    #else:
    #    # First check if we violate any bounds:
    #    for el, el_result in model_result.items():
    #        upper_bound = ld._live_upper_bounds[el]
    #        if upper_bound is not None and el_result >= upper_bound:
    #            #return min_likelihood   # For testing purposes only
    #            return 0.9*min_likelihood   # ...So these points should not be ignored
    #        lower_bound = ld._live_lower_bounds[el]
    #        if lower_bound is not None and el_result <= lower_bound:
    #            #return min_likelihood   # For testing purposes only
    #            return 0.9*min_likelihood   # ...So these points should not be ignored
    #    non_zero_model_values = list()
    #    for el in ld._live_elements_present:
    #        non_zero_model_values.append(model_result[el])
    #    model_values = np.asarray(non_zero_model_values)
    #    neg = 1/((2*np.pi)*((ld._live_non_zero_wd_errors)**2))
    #    chi2 = ((ld._live_non_zero_wd_abundances - model_values)**2)*(neg)*2*np.pi
    #    like = -0.5*np.sum(chi2 - np.log(neg))
    #    if like == -np.inf:
    #        #like = min_likelihood # For testing purposes only
    #        like = 0.9*min_likelihood
    #return like
