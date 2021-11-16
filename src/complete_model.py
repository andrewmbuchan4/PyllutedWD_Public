#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import math
import numpy as np

import abundance_model as am
import chemistry_info as ci
import disc_model as dm
import enhancement_model as em
import geology_info as gi
import live_data as ld
import white_dwarf_model as wdm

# TODOs:
# Move enhancement_model argument somewhere else. Maybe make it a member of the PollutionModel class?
# It's a bit awkward to force the user to manually ensure that t_disc is linear (i.e. not log(t_disc)) given that's how we return it at the end
def complete_model_calculation(fe_star, t_sinceaccretion, d_formation, z_formation, N_c, N_o, f_c, f_o, pollutionfraction, t_disc, pressure, fO2, enhancement_model='NonEarthlike', t_formation=1.5, normalise_abundances=True, snapshot_wd_atm=True):

    if enhancement_model == 'Earthlike':
        normalise_abundances = False # Slight hack for testing purposes
    
    diagnostics = dict()
    elements = ci.usual_elements
    # This limit on fe_star exists because outside of this range,
    # ld._live_stellar_compositions[int(round(fe_star))] will give a KeyError (there are 958 compositions)
    floored_fe_star = math.floor(fe_star)
    if 0 <= floored_fe_star <= 957:
        # TODO: Remove all magic numbers! Do foo.Mg or foo['Mg'] or something instead of foo[6]
        linear_d_formation = 10**(d_formation)
        
        abundances = am.get_all_abundances(elements, linear_d_formation, z_formation, t_formation, fe_star)

        disc_abundances = dict()
        for el_index, element in enumerate(elements):
            if el_index == 6:  # TODO: get this by indexing elements?
                # Mg is special
                disc_abundances[element] = abundances[element]
            else:
                disc_abundances[element] = abundances[element] * ld._live_stellar_compositions[floored_fe_star][el_index - 1 if el_index > 6 else el_index]
        
        diagnostics['DiscAbundances'] = disc_abundances
        
        # This is to speed up performance by making sure we only need to fully initialise the geo_model (and by extension the partitioning model) once
        if ld._geo_model is None:
            ld._geo_model = gi.GeologyModel(disc_abundances)
        else:
            ld._geo_model.reinit(disc_abundances)
        
        enhancement_model = em.EnhancementModel(enhancement_model)
        enhancements_dict, enhancements_diagnostics = enhancement_model.find_enhancements(
            ld._geo_model,
            disc_abundances,
            elements,
            N_c,
            N_o,
            f_c,
            f_o,
            pressure,
            fO2,
            normalise_abundances
        )
        
        if enhancements_dict is None:
            return None, None
        
        diagnostics['Enhancements'] = enhancements_diagnostics
        
        enhancements = list()
        for element in elements:
            enhancements.append(enhancements_dict[element])

        terms = [None]*len(elements)
        t = 0
        denominator = 0
        while t < len(terms):
            if t == 6:
                terms[t] = 1
                # Magnesium is special! (Because we ratio to Mg)
            else:
                stellar_index = t
                if t > 6:
                    stellar_index -= 1
                terms[t] = enhancements[t]/enhancements[6]
            if ld._live_all_wd_abundances[t] != 0:
                denominator += terms[t]
            t += 1
        elements_present = list()
        inv_denominator = 1/denominator
        planetesimal_abundance = np.zeros(len(elements))
        non_zero_planetesimal_abundance = np.zeros(len(ld._live_non_zero_wd_abundances))
        nzp_index = 0
        for i, wd_abundance in enumerate(ld._live_all_wd_abundances):
            abundance = terms[i]*inv_denominator
            if wd_abundance != 0:
                non_zero_planetesimal_abundance[nzp_index] = abundance
                nzp_index += 1
            planetesimal_abundance[i] = abundance
            elements_present.append(elements[i])

        result = wdm.process_abundances(
            t_sinceaccretion,
            t_disc,
            planetesimal_abundance,
            non_zero_planetesimal_abundance,
            ld._live_all_wd_timescales,
            ld._live_non_zero_wd_timescales,
            pollutionfraction,
            snapshot_wd_atm
        )
    else:
        raise # This should never happen!
    elements_present_dict = collections.OrderedDict(zip(elements_present, result))
    toret = collections.OrderedDict()
    for element in elements:
        toret[element] = elements_present_dict.get(element)
    return toret, diagnostics
