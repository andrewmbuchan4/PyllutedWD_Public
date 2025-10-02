#!/usr/bin/env python

import numpy as np

import complete_model as cm
import live_data as ld
import model_parameters as mp


def evaluate_log_likelihood(model_result, min_likelihood=mp.minimum_likelihood):
    if model_result is None:
        # like = min_likelihood  # For testing purposes only
        like = 1.1 * min_likelihood
        # ...so these points should be ignored (no information gained)
    else:
        like = ld._live_white_dwarf.log_likelihood(model_result, min_likelihood)
    return like


def universal_loglike(cube):

    parameter_indices = mp.parameter_indices(ld._live_model)
    input_values = list()
    for param in mp.get_model_params_in_order():
        to_append = (
            cube[parameter_indices[param]]
            if mp.model_uses_parameter(ld._live_model, param)
            else mp.default_values[ld._live_enhancement_model][param]
        )
        input_values.append(to_append)

    model_result, diagnostics_for_post_processing = cm.complete_model_calculation(
        input_values[0],  # fe_star
        input_values[1],  # t_sinceaccretion
        input_values[2],  # d_formation
        input_values[3],  # z_formation
        input_values[4],  # N_c
        input_values[5],  # N_o
        input_values[6],  # f_c
        input_values[7],  # f_o
        input_values[8],  # log_fragment mass
        10 ** input_values[9],  # t_disc
        input_values[10],  # pressure
        input_values[11],
        # ^ fO2 (oxygen fugacity relative to Iron Wuestite buffer, in log units)
        ld._live_enhancement_model,
        ld._live_consider_thermohaline,
    )

    min_likelihood = mp.minimum_likelihood
    # Any points with likelihood less than this are ignored
    # --> use for errors, but not for bounds violations
    # ( = -1e90 by default in the model, -1e100 by default in pymultinest)

    like = evaluate_log_likelihood(model_result, min_likelihood)
    return like
