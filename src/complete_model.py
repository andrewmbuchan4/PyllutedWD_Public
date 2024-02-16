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


# The only difference between this file and complete_model.py is that this imports
# white_dwarf_model_new (not white_dwarf_model) as wdm


#OrderedDict([
#(<Element.Al: 13>, -8.402920519628442),
#(<Element.Ti: 22>, -9.951585801900478),
#(<Element.Ca: 20>, -8.583525128145789),
#(<Element.Ni: 28>, -8.742000215870458),
#(<Element.Fe: 26>, -7.516201997514827),
#(<Element.Cr: 24>, -9.314063107265891),
#(<Element.Mg: 12>, -7.259119623605825),
#(<Element.Si: 14>, -7.300817124742995),
#(<Element.Na: 11>, -8.758592778738409),
#(<Element.O: 8>, -6.027977076736203),
#(<Element.C: 6>, -6.37123220460242),
#(<Element.N: 7>, -7.199452361943866)])

# TODOs:
# Move enhancement_model argument somewhere else. Maybe make it a member of the PollutionModel class?
# It's a bit awkward to force the user to manually ensure that t_disc is linear (i.e. not log(t_disc)) given that's how we return it at the end
# And we should really model relative abundances as far as possible! The uncertainty on those should be smaller than the absolute values
def complete_model_calculation(fe_star, t_sinceaccretion, d_formation, z_formation, N_c, N_o, f_c, f_o, pollutionfraction, t_disc, pressure, fO2, enhancement_model='NonEarthlike', t_formation=1.5, normalise_abundances=True, snapshot_wd_atm=True):

    #print()
    #print('cm42')
    #print(fe_star)
    #print(t_sinceaccretion)
    #print(d_formation)
    #print(z_formation)
    #print(N_c)
    #print(N_o)
    #print(f_c)
    #print(f_o)
    #print(pollutionfraction)
    #print(t_disc)
    #print(pressure)
    #print(fO2)
    #print(enhancement_model)


    #if enhancement_model == 'Earthlike':
    #    normalise_abundances = False # Slight hack for testing purposes

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

        planetesimal_abundance = list()
        for element in elements:
            planetesimal_abundance.append(enhancements_dict[element])
        result = wdm.process_abundances(
            t_sinceaccretion,
            t_disc,
            planetesimal_abundance,
            ld._live_all_wd_timescales,
            pollutionfraction,
            snapshot_wd_atm
        )
    else:
        raise ValueError("Metallicity must be between 0 and 958")
    elements_present_dict = collections.OrderedDict(zip(elements, result))
    toret = collections.OrderedDict()
    for element in elements:
        toret[element] = elements_present_dict.get(element)
    return toret, diagnostics
