#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is more of a fragment abundance calculator now, rather than an enhancement model!

from enum import Enum
import chemistry_info as ci
import geology_info as gi

# For now, ignore the crust. Treat parent_crust_number_fraction and fragment_crust_number_fraction as if they're always 0

class EnhancementModel():
    
    def __init__(self, model_type):
        self.model_type = model_type
        self.known_models = {
            'Earthlike': self.find_enhancements_earthlike,
            'NonEarthlike': self.find_enhancements_nonearthlike,
            'MantleOnly': self.find_enhancements_nonearthlike
        }
        self.model = self.known_models.get(self.model_type)
        if self.model is None:
            try:
                model_str = str(model_type)
            except:
                model_str = '[Could not cast model to string]'
            raise ValueError('Unknown EnhancementModel. Input was: ' + model_str + '. Known models: ' + ', '.join(self.known_models.keys()))
    
        self.earthlike_info = {
            ci.Element.C: {'default_e': 0, 'fill_core_first': False},
            ci.Element.N: {'default_e': 0, 'fill_core_first': False},
            ci.Element.O: {'default_e': 0, 'fill_core_first': True},
            ci.Element.Na: {'default_e': 0, 'fill_core_first': True},
            ci.Element.Mg: {'default_e': 0.0001, 'fill_core_first': True},
            ci.Element.Al: {'default_e': 0, 'fill_core_first': True},
            ci.Element.Si: {'default_e': 0, 'fill_core_first': True},
            ci.Element.Ca: {'default_e': 0, 'fill_core_first': True},
            ci.Element.Ti: {'default_e': 0, 'fill_core_first': True},
            ci.Element.Cr: {'default_e': 0, 'fill_core_first': True},
            ci.Element.Fe: {'default_e': 0, 'fill_core_first': True},
            ci.Element.Ni: {'default_e': 0, 'fill_core_first': False}
        }
    
    def find_enhancements(
        self,
        geo_model,
        disc_abundances,
        elements,
        parent_core_number_fraction,
        parent_crust_number_fraction,
        fragment_core_number_fraction,
        fragment_crust_number_fraction,
        pressure,
        fO2,
        normalise_abundances=True,
        extra_output=False
    ):  # TODO: Figure out something a bit more elegant than sending every parameter then ignoring half of them
        return self.model(
            geo_model,
            disc_abundances,
            elements,
            parent_core_number_fraction,
            parent_crust_number_fraction,
            fragment_core_number_fraction,
            fragment_crust_number_fraction,
            pressure,
            fO2,
            normalise_abundances,
            extra_output
        )

    def find_enhancements_nonearthlike(self, geo_model, disc_abundances, elements, parent_core_number_fraction, parent_crust_number_fraction, fragment_core_number_fraction, fragment_crust_number_fraction, pressure=54, fO2=-2, normalise_abundances=True, extra_output=False):  # 54/-2 from Fischer+ 2015
        if fragment_core_number_fraction is None:
            # This means no differentiation took place
            return disc_abundances, None
        # The geo_model should have been initialised using disc_abundances
        number_abundances, parent_core_number_fraction, Ds, all_Ds = geo_model.form_a_planet_iteratively(pressure, fO2)
        if number_abundances is None:
            return None, None
        enhancements = dict()
        fragment_mantle_number_fraction = 1 - fragment_core_number_fraction
        for element in elements:
            # The fragment abundances in each Layer are the same as the parent:
            fragment_core_abundance = number_abundances[element][gi.Layer.core]
            fragment_mantle_abundance = number_abundances[element][gi.Layer.mantle]
            # But the bulk is different if it has different mantle/core fractions:
            enhancements[element] = (fragment_mantle_number_fraction*fragment_mantle_abundance) + (fragment_core_number_fraction*fragment_core_abundance)
        diagnostics_dict = {
            'Abundances': number_abundances,
            'ParentCoreNumberFraction': parent_core_number_fraction,
            'Ds': Ds
        }
        return enhancements, diagnostics_dict
    
    def find_enhancements_earthlike(self, geo_model_ignore, disc_abundances, elements, parent_core_number_fraction, parent_crust_number_fraction, fragment_core_number_fraction, fragment_crust_number_fraction, pressure, fO2, normalise_abundances=True, extra_output=False):
        if (parent_core_number_fraction is None) or (parent_crust_number_fraction is None) or (fragment_core_number_fraction is None) or (fragment_crust_number_fraction is None):
            # This means no differentiation took place
            return disc_abundances
        parent_is_physical = 0.19 >= parent_core_number_fraction >= 0 and parent_crust_number_fraction >= 0 and parent_core_number_fraction + parent_crust_number_fraction <= 1
        fragment_is_physical = fragment_core_number_fraction >= 0 and fragment_crust_number_fraction >= 0 and fragment_core_number_fraction + fragment_crust_number_fraction <= 1
        fragment_is_physical = fragment_is_physical and not (fragment_core_number_fraction > parent_core_number_fraction and fragment_crust_number_fraction > parent_crust_number_fraction)
        fragment_is_physical = fragment_is_physical and not (fragment_core_number_fraction > (parent_core_number_fraction/(1-parent_crust_number_fraction)) and fragment_crust_number_fraction > 0.01)
        fragment_is_physical = fragment_is_physical and not ((parent_core_number_fraction/(1-parent_crust_number_fraction)) >= fragment_core_number_fraction >= parent_core_number_fraction and fragment_core_number_fraction - parent_core_number_fraction > (parent_crust_number_fraction-fragment_crust_number_fraction)*(parent_core_number_fraction/(1-parent_crust_number_fraction)))
        
        enhancements = dict()
        if not (fragment_is_physical and parent_is_physical):
            return None, {'ParentCoreNumberFraction': parent_core_number_fraction}
        else:
            # Ignore the input geo_model because it could be non-Earth-like, instead make a new one:
            earthlike_geo_model = gi.GeologyModel(None, normalise_abundances)
            fragment_mantle_number_fraction = 1 - (fragment_core_number_fraction + fragment_crust_number_fraction)
            parent_mantle_number_fraction = 1 - (parent_core_number_fraction + parent_crust_number_fraction)
            for element in elements:
                parent_bulk_abundance = earthlike_geo_model.get_bulk_abundance(element)
                parent_crust_abundance = earthlike_geo_model.get_crust_abundance(element)
                parent_core_abundance = earthlike_geo_model.get_core_abundance(element)
                parent_mantle_abundance = (parent_bulk_abundance - ((parent_crust_abundance*parent_crust_number_fraction) + (parent_core_abundance*parent_core_number_fraction)))/parent_mantle_number_fraction
                
                if parent_mantle_abundance >= 0:
                    fragment_mantle_abundance = parent_mantle_abundance
                    fragment_crust_abundance = parent_crust_abundance
                    fragment_core_abundance = parent_core_abundance
                else:
                    fragment_mantle_abundance = 0
                    if self.earthlike_info[element]['fill_core_first']:
                        fragment_crust_abundance = (parent_bulk_abundance - (parent_core_abundance*parent_core_number_fraction))/parent_crust_number_fraction
                        fragment_core_abundance = parent_core_abundance
                    else:
                        fragment_crust_abundance = parent_crust_abundance
                        fragment_core_abundance = (parent_bulk_abundance - (parent_crust_abundance*parent_crust_number_fraction))/parent_core_number_fraction
                
                fragment_bulk_abundance_term1 = fragment_crust_number_fraction*fragment_crust_abundance
                fragment_bulk_abundance_term2 = fragment_mantle_number_fraction*fragment_mantle_abundance
                fragment_bulk_abundance_term3 = fragment_core_number_fraction*fragment_core_abundance
                fragment_bulk_abundance =  fragment_bulk_abundance_term1 + fragment_bulk_abundance_term2 + fragment_bulk_abundance_term3
                
                enhancements[element] = disc_abundances[element]*(fragment_bulk_abundance/parent_bulk_abundance)
        diagnostics_dict = {'ParentCoreNumberFraction': parent_core_number_fraction}
        return enhancements, diagnostics_dict
