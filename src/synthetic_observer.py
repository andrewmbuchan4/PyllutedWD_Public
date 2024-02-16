#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum
import numpy as np

import chemistry_info as ci
import detection_thresholds as dt
import model_parameters as mp
import synthetic_bandpass as sb

class ObservationType(Enum):
    NoCut = 0
    CaMgFeCutoff = 1
    HollandsColourCut = 2
    IndividualElementCutoff = 3
    TeffIndividualElementCutoff = 4

    def __str__(self):
        return self.name

class Observer:

    def __init__(self, observation_type, error_dict=dict(), threshold_type='Default'):
        self.error_dict = error_dict
        self.observation_type = observation_type
        self.threshold_offset = 0
        try:
            self.threshold_offset += threshold_type # This means you can set threshold_type to -1 and it'll use the default thresholds but all reduced by 1 dex
            self.threshold_type = 'Default'
        except TypeError:
            self.threshold_type = threshold_type
        self.threshold_def_dict = dt.threshold_bank[self.threshold_type]
        if self.observation_type == ObservationType.NoCut:
            self.observation_function = self.apply_NoCut
        elif self.observation_type == ObservationType.CaMgFeCutoff:
            self.observation_function = self.apply_CaMgFeCutoff
        elif self.observation_type == ObservationType.HollandsColourCut:
            self.observation_function = self.apply_HollandsColourCut
        elif self.observation_type == ObservationType.IndividualElementCutoff:
            self.observation_function = self.apply_IndividualElementCutoff
        elif self.observation_type == ObservationType.TeffIndividualElementCutoff:
            self.observation_function = self.apply_TeffDependentIndividualElementCutoff
        else:
            raise ValueError('Unrecognised observation type ' + str(observation_type))

    def observe_populations(self, population_dict, overwrite=False):
        observed_populations = dict()
        for pop_name, population in population_dict.items():
            self.observe_population(population, overwrite)

    def observe_population(self, population, overwrite=False):
        for system in population:
            if overwrite or system.observed is None:
                print('Trying to observe system ' + str(system.id))
                self.observe_system(system)
            else:
                print('System ' + str(system.id) + ' already assessed for observability, skipping')

    def observe_system(self, system):
        abundances_as_observed = self.observation_function(system)
        if len(abundances_as_observed) > 0:
            system.observed = True
            system.set_observed_abundances(abundances_as_observed)
        else:
            system.observed = False

    def apply_errors(self, abundances):
        new_pollution_abundances = dict()
        for element, abundance in abundances.items():
            error = self.error_dict.get(element, 0)
            new_abundance = abundance
            if error > 0:
                new_abundance = np.random.normal(abundance, error)
            new_pollution_abundances[element] = new_abundance
        return new_pollution_abundances

    def apply_CaMgFeCutoff(self, system):
        CaHx = system.pollution_abundances.get(ci.Element.Ca, -np.inf)
        MgHx = system.pollution_abundances.get(ci.Element.Mg, -np.inf)
        FeHx = system.pollution_abundances.get(ci.Element.Fe, -np.inf)
        cutoff = -9
        if (CaHx < cutoff) or (MgHx < cutoff) or (FeHx < cutoff):
            return dict()
        return self.apply_errors(system.pollution_abundances)

    def apply_HollandsColourCut(self, system):
        u = system.bandpass_magnitudes.get(sb.Bandpass.u)
        g = system.bandpass_magnitudes.get(sb.Bandpass.g)
        r = system.bandpass_magnitudes.get(sb.Bandpass.r)
        if u is None or g is None or r is None:
            return dict()
        y = u - g
        x = g - r
        y_in_bounds = y > 0.5 and y < 3.8  # All this is eyeballed from Fig 4 of Hollands 2017
        #        left edge defined by y = 3x + 1.5, right edge defined by y=3x-0.5 roughly
        # rearranging, we require x > y/3 - 0.5, x < y/3 + 1/6
        x_in_bounds = x > (y - 1.5)/3 and x < (y + 0.5)/3
        if y_in_bounds and x_in_bounds:
            return self.apply_errors(system.pollution_abundances)
        return dict()

    def apply_NoCut(self, system):
        return self.apply_errors(system.pollution_abundances)

    def apply_IndividualElementCutoff(self, system):
        # Numbers here are arbitrary - this function is unused!
        cutoff_thresholds = {  # Abundances have to be above this in order to detect element
            ci.Element.Al: -7 + self.threshold_offset,
            ci.Element.Ti: -7 + self.threshold_offset,
            ci.Element.Ca: -9 + self.threshold_offset,
            ci.Element.Ni: -7 + self.threshold_offset,
            ci.Element.Fe: -8 + self.threshold_offset,
            ci.Element.Cr: -7 + self.threshold_offset,
            ci.Element.Mg: -8 + self.threshold_offset,
            ci.Element.Si: -7 + self.threshold_offset,
            ci.Element.Na: -7 + self.threshold_offset,
            ci.Element.O: -7 + self.threshold_offset,
            ci.Element.C: -7 + self.threshold_offset,
            ci.Element.N: -7 + self.threshold_offset
        }
        # Firstly, apply random noise:
        noisy_abundances = self.apply_errors(system.pollution_abundances)
        # Now remove any elements which fall below the cutoff threshold
        toret = dict()
        for element, abundance in noisy_abundances.items():
            if abundance > cutoff_thresholds.get(element, np.inf):
                toret[element] = abundance
        return toret

    def calculate_element_cutoff_threshold(self, element, spectral_type, teff):
        # Assume there is a straight line defined by threshold = m*teff + c which sets the threshold
        m_c_tuple = self.threshold_def_dict[spectral_type][element]
        threshold = (m_c_tuple[0]*teff) + m_c_tuple[1] + self.threshold_offset # this is just y = mx + c, plus the threshold offset
        return threshold

    def apply_TeffDependentIndividualElementCutoff(self, system):
        cutoff_thresholds = {  # Abundances have to be above this in order to detect element
            ci.Element.Al: self.calculate_element_cutoff_threshold(ci.Element.Al, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.Ti: self.calculate_element_cutoff_threshold(ci.Element.Ti, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.Ca: self.calculate_element_cutoff_threshold(ci.Element.Ca, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.Ni: self.calculate_element_cutoff_threshold(ci.Element.Ni, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.Fe: self.calculate_element_cutoff_threshold(ci.Element.Fe, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.Cr: self.calculate_element_cutoff_threshold(ci.Element.Cr, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.Mg: self.calculate_element_cutoff_threshold(ci.Element.Mg, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.Si: self.calculate_element_cutoff_threshold(ci.Element.Si, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.Na: self.calculate_element_cutoff_threshold(ci.Element.Na, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.O: self.calculate_element_cutoff_threshold(ci.Element.O, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.C: self.calculate_element_cutoff_threshold(ci.Element.C, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature]),
            ci.Element.N: self.calculate_element_cutoff_threshold(ci.Element.N, system.wd_properties[mp.WDParameter.spectral_type], system.wd_properties[mp.WDParameter.temperature])
        }
        # Firstly, apply random noise:
        noisy_abundances = self.apply_errors(system.pollution_abundances)
        # Now remove any elements which fall below the cutoff threshold
        toret = dict()
        for element, abundance in noisy_abundances.items():
            if abundance > cutoff_thresholds.get(element, np.inf):
                toret[element] = abundance
        return toret
