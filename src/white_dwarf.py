#!/usr/bin/env python
# -*- coding: utf-8 -*-

from copy import deepcopy
from enum import Enum

import collections as cn
import numpy as np

import chemistry_info as ci
import model_parameters as mp
import pwd_utils as pu
import solar_abundances as sa

class WhiteDwarfDataPointType(Enum):
    measurement = 0
    upper_bound  = 1
    lower_bound = 2
    label = 3 # i.e., something non-numerical that can be represented as a string

class WhiteDwarfDataPoint():

    def __init__(self, data_point_type, value, upper_error=None, included=True, lower_error=None):
        #data_point_type should be a WhiteDwarfDataPointType
        # By default, lower_error is assumed to be the same as upper_error such that upper_error is the only one you need to care about
        self.data_point_type = data_point_type
        if self.data_point_type != WhiteDwarfDataPointType.label:
            self.value = float(value)
        else:
            self.value = value
        if self.data_point_type == WhiteDwarfDataPointType.measurement:
            if upper_error is not None:
                self.upper_error = float(upper_error) # only applicable if data_point_type is measurement, in which case this is the 1 sigma error
                if lower_error is None:
                    self.lower_error = self.upper_error # Can set this separately if errors are asymmetrical
                else:
                    self.lower_error = float(lower_error)
            else:
                raise ValueError("WhiteDwarfDataPointType is measurement, but no error was supplied")
        else:
            self.upper_error = None
            self.lower_error = None
        self.included = included

    def __str__(self):
        if self.data_point_type == WhiteDwarfDataPointType.measurement:
            if self.upper_error == self.lower_error:
                return str(self.value) + '±' + str(self.upper_error)
            else:
                return str(self.value) + '^{+' + str(self.upper_error) + '}_{-' + str(self.lower_error) + '}'
        elif self.data_point_type == WhiteDwarfDataPointType.upper_bound:
            return '<' + str(self.value)
        elif self.data_point_type == WhiteDwarfDataPointType.lower_bound:
            return '>' + str(self.value)
        elif self.data_point_type == WhiteDwarfDataPointType.label:
            return str(self.value)
        else:
            raise ValueError("Unknown WhiteDwarfDataPointType")

    def __lt__(self, other):
        try:
            return self.value < other.value
        except AttributeError: # e.g. if we compared to a float
            return self.value < other

    def __gt__(self, other):
        try:
            return self.value > other.value
        except AttributeError: # e.g. if we compared to a float
            return self.value > other

    def __eq__(self, other):
        try:
            return self.value == other.value
        except AttributeError: # e.g. if we compared to a float
            return self.value == other

    def __neq__(self, other):
        return not self == other

    def __repr__(self):
        return str(self)

class WhiteDwarfAbundanceData():

    def __init__(self, abundance_data_dict):
        self.abundance_data_dict = abundance_data_dict # expecting {ci.Element: WhiteDwarfDataPoint, ... }
        for element, abundance_data in self.abundance_data_dict.items():
            if abundance_data.value > 0:
                abundance_data.value = -abundance_data.value # These values are always negative - now you don't need to keep typing minus!

    def get_plottable_points(self, elements_to_plot, reference_element, included=True, scale_to_solar=False): # Switch included to False to get the excluded points
        if reference_element not in self.abundance_data_dict:
            raise ValueError('Requested abundances relative to ' + str(reference_element) + ', which is absent from the WhiteDwarfAbundanceData')
        values_toret = list()
        upper_errors_toret = list()
        lower_errors_toret = list()
        upper_bounds_bools_toret = list()
        lower_bounds_bools_toret = list()
        elements_toret = list()
        error_arrow_size = 0.1
        for el in elements_to_plot:
            try:
                el_data = self.abundance_data_dict[el]
            except KeyError:
                values_toret.append(np.nan)
                upper_errors_toret.append(np.nan)
                lower_errors_toret.append(np.nan)
                upper_bounds_bools_toret.append(False)
                lower_bounds_bools_toret.append(False)
                elements_toret.append(el)
                continue
            if el != reference_element:
                if el_data.included == included:
                    if el_data.data_point_type == WhiteDwarfDataPointType.measurement:
                        values_toret.append(self.get_relative_value(el, reference_element, scale_to_solar))
                        upper_errors_toret.append(el_data.upper_error)
                        lower_errors_toret.append(el_data.lower_error)
                        upper_bounds_bools_toret.append(False)
                        lower_bounds_bools_toret.append(False)
                        elements_toret.append(el)
                    elif el_data.data_point_type == WhiteDwarfDataPointType.upper_bound:
                        values_toret.append(self.get_relative_value(el, reference_element, scale_to_solar))
                        upper_errors_toret.append(0)
                        lower_errors_toret.append(error_arrow_size)
                        upper_bounds_bools_toret.append(True)
                        lower_bounds_bools_toret.append(False)
                        elements_toret.append(el)
                    elif el_data.data_point_type == WhiteDwarfDataPointType.lower_bound:
                        values_toret.append(self.get_relative_value(el, reference_element, scale_to_solar))
                        upper_errors_toret.append(error_arrow_size)
                        lower_errors_toret.append(0)
                        upper_bounds_bools_toret.append(False)
                        lower_bounds_bools_toret.append(True)
                        elements_toret.append(el)
                    else:
                        raise ValueError("Unknown WhiteDwarfDataPointType")
                else:
                    values_toret.append(np.nan)
                    upper_errors_toret.append(np.nan)
                    lower_errors_toret.append(np.nan)
                    upper_bounds_bools_toret.append(False)
                    lower_bounds_bools_toret.append(False)
                    elements_toret.append(el)
        return elements_toret, values_toret, upper_errors_toret, lower_errors_toret, upper_bounds_bools_toret, lower_bounds_bools_toret

    def get_relative_value(self, element, reference_element, scale_to_solar):
        wd_rel_abundance = self.abundance_data_dict[element].value - self.abundance_data_dict[reference_element].value
        if scale_to_solar:
            solar_rel_abundance = sa.get_solar_relative_abundance(element, reference_element)
            return wd_rel_abundance - solar_rel_abundance
        else:
            return wd_rel_abundance

    def get_log_ratio(self, element1, element2):
        el1 = self.abundance_data_dict.get(element1)
        el2 = self.abundance_data_dict.get(element2)
        if el1 is None or el2 is None:
            ratio = None
        elif (not el1.included) or (not el2.included):
            ratio = None
        else:
            ratio = el1.value - el2.value
        return ratio

    def get_measurements_as_lists(self, return_included=True, return_excluded=False, elements=ci.all_elements, preserve_length=False):
        els_toret = list()
        values_toret = list()
        upper_errors_toret = list()
        lower_errors_toret = list()
        for el in elements:
            try:
                el_data = self.abundance_data_dict[el]
            except KeyError:
                if preserve_length:
                    els_toret.append(None)
                    values_toret.append(None)
                    upper_errors_toret.append(None)
                    lower_errors_toret.append(None)
                    continue
                else:
                    continue
            if el_data.data_point_type == WhiteDwarfDataPointType.measurement and ((el_data.included and return_included) or (not el_data.included and return_excluded)):
                els_toret.append(el)
                values_toret.append(el_data.value)
                upper_errors_toret.append(el_data.upper_error)
                lower_errors_toret.append(el_data.lower_error)
            else:
                if preserve_length:
                    els_toret.append(None)
                    values_toret.append(None)
                    upper_errors_toret.append(None)
                    lower_errors_toret.append(None)
                else:
                    pass
        return els_toret, values_toret, upper_errors_toret, lower_errors_toret

    def get_upper_bounds_as_lists(self, return_included=True, return_excluded=False, elements=ci.all_elements, preserve_length=False):
        els_toret = list()
        values_toret = list()
        for el in elements:
            try:
                el_data = self.abundance_data_dict[el]
            except KeyError:
                if preserve_length:
                    els_toret.append(None)
                    values_toret.append(None)
                    continue
                else:
                    continue
            if el_data.data_point_type == WhiteDwarfDataPointType.upper_bound and ((el_data.included and return_included) or (not el_data.included and return_excluded)):
                els_toret.append(el)
                values_toret.append(el_data.value)
            else:
                if preserve_length:
                    els_toret.append(None)
                    values_toret.append(None)
                else:
                    pass
        return els_toret, values_toret

    def get_lower_bounds_as_lists(self, return_included=True, return_excluded=False, elements=ci.all_elements, preserve_length=False):
        els_toret = list()
        values_toret = list()
        for el in elements:
            try:
                el_data = self.abundance_data_dict[el]
            except KeyError:
                if preserve_length:
                    els_toret.append(None)
                    values_toret.append(None)
                    continue
                else:
                    continue
            if el_data.data_point_type == WhiteDwarfDataPointType.lower_bound and ((el_data.included and return_included) or (not el_data.included and return_excluded)):
                els_toret.append(el)
                values_toret.append(el_data.value)
            else:
                if preserve_length:
                    els_toret.append(None)
                    values_toret.append(None)
                else:
                    pass
        return els_toret, values_toret

    def get_abundance_arrays(self, elements=ci.usual_elements):
        abundances = list()
        errors = list()
        upper_bounds = list()
        lower_bounds = list()
        for element in elements:
            data_point = self.abundance_data_dict.get(element)
            if data_point is None:
                abundances.append(None)
                errors.append(None)
                upper_bounds.append(None)
                lower_bounds.append(None)
            else:
                abundances.append(data_point.value if data_point.included and data_point.data_point_type == WhiteDwarfDataPointType.measurement else None)
                errors.append(data_point.upper_error if data_point.included and data_point.data_point_type == WhiteDwarfDataPointType.measurement else None)
                upper_bounds.append(data_point.value if data_point.included and data_point.data_point_type == WhiteDwarfDataPointType.upper_bound else None)
                lower_bounds.append(data_point.value if data_point.included and data_point.data_point_type == WhiteDwarfDataPointType.lower_bound else None)
        return np.array(abundances), np.array(errors), np.array(upper_bounds), np.array(lower_bounds)

    def get_elements_present(self):
        toret = list()
        for element, el_data in self.abundance_data_dict.items():
            if el_data.data_point_type == WhiteDwarfDataPointType.measurement and el_data.included and el_data.value not in [None, np.nan]:
                toret.append(element)
        return toret

    def get_modellable_detected_elements(self, detectable_elements=ci.usual_elements):
        toret = list()
        for el in detectable_elements:
            data_point = self.abundance_data_dict.get(el)
            if data_point is not None and data_point.data_point_type == WhiteDwarfDataPointType.measurement and data_point.upper_error > 0 and data_point.lower_error > 0:
                toret.append(el)
        return toret

    def __str__(self):
        toret = ''
        for el in ci.all_elements:
            try:
                el_data = self.abundance_data_dict[el]
            except KeyError:
                continue
            individual_str = str(el) + ': ' + str(el_data)
            if not el_data.included:
                individual_str += ' (not included)'
            toret += individual_str + '\n'
        return toret

    def __repr__(self):
        return str(self)

    def get_all_elements(self):
        return sorted(list(self.abundance_data_dict.keys()))

    def get_abundance(self, element):
        return self.abundance_data_dict.get(element)

    def __eq__(self, other):
        return self.abundance_data_dict == other.abundance_data_dict

    def __neq__(self, other):
        return not self == other

class WhiteDwarfPropertyData():

    def __init__(self, property_dict):
        self.property_dict = property_dict # expecting {mp.WDParameter: WhiteDwarfDataPoint, ... }
        # I'm here making the decision that a White Dwarf has a 'spectral_type' ('DA' or 'DB' for now)
        # And, separately, an 'atmospheric_type' (either ci.Element.H or ci.Element.He, accordingly)
        # I will further assume that DA <==> ci.Element.H and DB <==> ci.Element.He (not quite true in reality - this is a simplification which I think makes sense here)
        self.fill_in_spectral_and_atmospheric_types()
        self.check_spectral_and_atmospheric_types()

    def fill_in_spectral_and_atmospheric_types(self):
        if self.get_property_data(mp.WDParameter.spectral_type) is None and self.get_property_data(mp.WDParameter.atmospheric_type) is not None:
            if self.get_property_data(mp.WDParameter.atmospheric_type).value == ci.Element.H:
                self.property_dict[mp.WDParameter.spectral_type] = WhiteDwarfDataPoint(WhiteDwarfDataPointType.label, 'DA')
            elif self.get_property_data(mp.WDParameter.atmospheric_type).value == ci.Element.He:
                self.property_dict[mp.WDParameter.spectral_type] = WhiteDwarfDataPoint(WhiteDwarfDataPointType.label, 'DB')
            else:
                pass
                #raise ValueError('Must have a valid spectral type or atmospheric type! Received ' + str(self.get_property_data(mp.WDParameter.spectral_type)) + ' and ' + str(self.get_property_data(mp.WDParameter.atmospheric_type)))
        elif self.get_property_data(mp.WDParameter.spectral_type) is not None and self.get_property_data(mp.WDParameter.atmospheric_type) is None:
            if self.get_property_data(mp.WDParameter.spectral_type).value == 'DA':
                self.property_dict[mp.WDParameter.atmospheric_type] = WhiteDwarfDataPoint(WhiteDwarfDataPointType.label, ci.Element.H)
            elif self.get_property_data(mp.WDParameter.spectral_type).value == 'DB':
                self.property_dict[mp.WDParameter.atmospheric_type] = WhiteDwarfDataPoint(WhiteDwarfDataPointType.label, ci.Element.He)
            else:
                pass
                #raise ValueError('Must have a valid spectral type or atmospheric type! Received ' + self.get_property_data(mp.WDParameter.spectral_type) + ' and ' + self.get_property_data(mp.WDParameter.atmospheric_type))
        else:
            pass
            # Check the values are consistent

    def check_spectral_and_atmospheric_types(self):
        if self.get_property_data(mp.WDParameter.spectral_type) is None and self.get_property_data(mp.WDParameter.atmospheric_type) is None:
            # Allow this as a special case - we don't demand that the user must supply these properties
            return
        if self.get_property_data(mp.WDParameter.spectral_type) is None:
            # Something strange has happened
            raise ValueError('Spectral/atmospheric types were not valid')
        if self.get_property_data(mp.WDParameter.atmospheric_type) is None:
            raise ValueError('Spectral/atmospheric types were not valid')
        valid_da_combo = self.get_property_data(mp.WDParameter.spectral_type).value == 'DA' and self.get_property_data(mp.WDParameter.atmospheric_type).value == ci.Element.H
        valid_db_combo = self.get_property_data(mp.WDParameter.spectral_type).value == 'DB' and self.get_property_data(mp.WDParameter.atmospheric_type).value == ci.Element.He
        if not (valid_da_combo or valid_db_combo):
            raise ValueError('Spectral type and atmospheric types are not consistent! Received ' + str(self.get_property_data(mp.WDParameter.spectral_type).value) + ' and ' + str(self.get_property_data(mp.WDParameter.atmospheric_type).value))

    def __str__(self):
        toret = ''
        for wd_property in mp.WDParameter:
            if wd_property not in [mp.WDParameter.logq, mp.WDParameter.timescale_type, mp.WDParameter.consider_thermohaline]: # logq is associated with atmosphere data, it's not treated as a property of the white dwarf itself. Similarly, the timescale_type  (and consider_thermohaline) is only a WDParameter for the sake of the synthetic WDs - normally it would be associated with the atmosphere data
                property_value = self.property_dict.get(wd_property, '---')#for wd_property, property_value in self.property_dict.items():
                individual_str = str(wd_property) + ': ' + str(property_value)
                if self.property_dict.get(wd_property) is not None:
                    unit = mp.wd_parameter_units.get(wd_property, '')
                    if unit != '':
                        individual_str += ' ' + unit
                toret += individual_str + '\n'
        return toret

    def __repr__(self):
        return str(self)

    def get_property_data(self, property_to_get):
        # Could put default values here? eg 0.6 for mass, 8 for log(g)
        return self.property_dict.get(property_to_get)

    def __eq__(self, other):
        return self.property_dict == other.property_dict

    def __neq__(self, other):
        return not self == other

class WhiteDwarfAtmosphereData():
    # To contain logq and sinking timescales.
    def __init__(self, timescale_and_logq_dict):
        # expecting input to be a dict containing both logq and all the timescales, as this is what is calculated by TimescaleInterpolator
        self.logq_string = 'logq'
        try:
            self.logq = WhiteDwarfDataPoint(WhiteDwarfDataPointType.measurement, timescale_and_logq_dict[self.logq_string], 0) # For now, assume zero error on this
        except KeyError:
            self.logq = WhiteDwarfDataPoint(WhiteDwarfDataPointType.measurement, timescale_and_logq_dict[mp.WDParameter.logq], 0) # For now, assume zero error on this
        self.timescale_dict = {element: timescale_and_logq_dict[element] for element in timescale_and_logq_dict if element not in [self.logq_string, mp.WDParameter.logq]}

    def __str__(self):
        toret = 'log(q): ' + str(self.logq) + '\n'
        for element, timescale in self.timescale_dict.items():
            toret += str(element) + ': ' + str(timescale) + ' yr\n'
        toret += '\n'
        return toret

    def __repr__(self):
        return str(self)

    def get_timescales_as_array(self, elements=ci.usual_elements):
        return np.array([self.timescale_dict[el] for el in elements])

    def get_timescale_values_dict(self, elements=ci.usual_elements):
        toret = cn.OrderedDict()
        for el in elements:
            toret[el] = self.timescale_dict.get(el)
        return toret

    def __eq__(self, other):
        return self.logq == other.logq and self.timescale_dict == other.timescale_dict

    def __neq__(self, other):
        return not self == other

class WhiteDwarfAtmosphereDataset():
    # A dict of WhiteDwarfAtmosphereData
    def __init__(self, atmosphere_data_dict):
        self.atmosphere_data_dict = atmosphere_data_dict # expecting {ti.TimescaleType.something: WhiteDwarfAtmosphereData}

    def __str__(self):
        toret = ''
        for key, add in self.atmosphere_data_dict.items():
            toret += str(key) + ' Sinking Timescales:\n' + str(add)
        toret += '\n'
        return toret

    def __repr__(self):
        return str(self)

    def get_available_timescale_types(self):
        return list(self.atmosphere_data_dict.keys())

    def get_timescales_as_array(self, timescale_type, elements=ci.usual_elements):
        try:
            return self.atmosphere_data_dict[timescale_type].get_timescales_as_array(elements)
        except KeyError:
            return None

    def get_timescale_values_dict(self, timescale_type, elements=ci.usual_elements):
        return self.atmosphere_data_dict[timescale_type].get_timescale_values_dict(elements)

    def get_logq(self, timescale_type):
        return self.atmosphere_data_dict[timescale_type].logq

    def __eq__(self, other):
        if self.atmosphere_data_dict.keys() != other.atmosphere_data_dict.keys():
            return False
        for key in self.atmosphere_data_dict.keys():
            if self.atmosphere_data_dict[key] != other.atmosphere_data_dict[key]:
                return False
        return True

    def __neq__(self, other):
        return not self == other

class WhiteDwarf():

    def __init__(self, name_tuple_or_str, properties, abundances, timescales=dict()):
        if isinstance(name_tuple_or_str, tuple):
            self.system_name = name_tuple_or_str[0] # A string
            self.variant = name_tuple_or_str[1] # A string
            if self.variant == '':
                self.variant = None
        else:
            #Assume name_tuple_or_str is a string
            self.system_name = name_tuple_or_str
            self.variant = None
        self.abbreviated_name_as_used = None
        self.properties = deepcopy(properties) # Should be a WhiteDwarfPropertyData. Create a deepcopy because we want this to be unique
        self.abundances = deepcopy(abundances) # Should be a WhiteDwarfAbundanceData. Create a deepcopy because we want this to be unique
        # TODO : add sinking timescales to this
        if len(self.abundances.abundance_data_dict) > 0:
            if self.properties.get_property_data(mp.WDParameter.atmospheric_type) is None:
                raise ValueError('Must supply an atmospheric type if supplying abundances') # otherwise the abundances mean nothing
        self.abundances.abundance_data_dict[self.properties.get_property_data(mp.WDParameter.atmospheric_type).value] = WhiteDwarfDataPoint(WhiteDwarfDataPointType.measurement, 0, 0, False)
        self.atmosphere_data = self.unpack_timescales(timescales)

    def full_name(self):
        if self.variant is None:
            return self.system_name
        return self.system_name + '_' + self.variant

    def abbreviated_name(self, chars_to_remove):
        toret = self.full_name()
        max_length = len(toret) - chars_to_remove
        abbreviated_system_names = self.get_abbreviated_system_names()
        abbreviated_variants = self.get_abbreviated_variants()
        system_name_index = 0
        variant_index = 0
        while len(toret) > max_length:
            if len(abbreviated_variants[variant_index]) < 3 or variant_index == len(abbreviated_variants) - 1:
                abbreviate_system_name = True
            else:
                abbreviate_system_name = False
            if abbreviate_system_name:
                system_name_index += 1
            else:
                variant_index += 1
            if abbreviated_variants[variant_index] is None:
                toret = abbreviated_system_names[system_name_index]
            else:
                toret = abbreviated_system_names[system_name_index] + '_' + abbreviated_variants[variant_index]
        self.abbreviated_name_as_used = toret
        return toret

    def get_abbreviated_system_names(self):
        special_prefixes = ['USNO-A2.0', 'USNO-B1.0'] # <-- These mess around with the rest of the logic, so we'll handle them separately
        # What we're going to try here is look for coordinates that we can truncate. Identify these by presence of a + or -
        # Assuming we have only one of these, and it only occurs once
        coordinate_separators = ['+', '-']
        valid_chars = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '.', '?']
        min_coord_length = 2
        special_prefix = ''
        sys_name_to_use = self.system_name
        for sp in special_prefixes:
            if self.system_name.startswith(sp):
                special_prefix = sp
                sys_name_to_use = self.system_name.split(sp)[1]
        cs_present = None
        for cs in coordinate_separators:
            if cs in sys_name_to_use:
                cs_present = cs
        if cs_present is None:
            return [self.system_name]
        # Now we need to crawl through the string to identify the coordinate values
        cs_index = sys_name_to_use.index(cs_present)
        test_index = cs_index - 1
        coord1 = '' # Before the cs
        while test_index >= 0:
            test_char = sys_name_to_use[test_index]
            if test_char in valid_chars:
                coord1 += test_char
            else:
                break
            test_index -= 1
        coord1 = ''.join(reversed(coord1))
        prefix = sys_name_to_use.split(coord1)[0]
        test_index = cs_index + 1
        coord2 = '' # After the cs
        while test_index < len(sys_name_to_use):
            test_char = sys_name_to_use[test_index]
            if test_char in valid_chars:
                coord2 += test_char
            else:
                break
            test_index += 1
        coords = [coord1, coord2]
        if len(coord1) < 1 or len(coord2) < 1:
            # Then this won't work
            return [self.system_name]
        prefix = sys_name_to_use.split(coord1)[0]
        suffix = sys_name_to_use.split(coord2)[1]
        toret = [special_prefix + sys_name_to_use]
        # Sometimes the suffix may include special abbreviable words
        abbreviable_words = ['Phot', 'Opt', 'Spec', 'UV', 'uv', 'phot', 'spec', 'opt', 'NoOvershoot', 'Overshoot'] #case sensitive
        # There's a slight problem with this logic because we end up with no option to abbreviated the above strings and not also abbreviate the coordinate - but we should probably default to this if possible!
        for aw in abbreviable_words:
            suffix = suffix.replace(aw, aw[0])
        while len(coords[0]) > min_coord_length and len(coords[1]) > min_coord_length:
            #Now we need to know how many decimal places they each have, again assuming a dot will only appear (max) once...
            dp_counts = [0, 0]
            for i in [0, 1]:
                if '.' in coords[i]:
                    dp_counts[i] = len(coords[i]) - (coords[i].index('.') + 1)
            if dp_counts[0] > dp_counts[1]:
                to_remove = [1, 0]
            elif dp_counts[1] > dp_counts[0]:
                to_remove = [0, 1]
            else:
                to_remove = [1, 1]
            for i, coord in enumerate(coords):
                if to_remove[i] > 0:
                    coord = coord[:-to_remove[i]]
                if coord.endswith('.'):
                    coord = coord[:-1]
                coords[i] = coord
            #Now we actually perform the abbreviation:
            abbreviation = special_prefix + prefix + coords[0] + cs_present + coords[1] + suffix
            toret.append(abbreviation)
        return toret

    def get_abbreviated_variants(self):
        numbers = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
        allowable_letters_at_end = ['A', 'a', 'b']
        allowable_characters = numbers + allowable_letters_at_end
        #special_variants_dict = {
        #    'Gaiarerun': 'GR'
        #}
        special_suffixes = ['GR']
        suffix = ''
        toret = [self.variant]
        if self.variant is None:
            return toret
        #if self.variant in special_variants_dict:
        #    toret.append(special_variants_dict[self.variant])
        #    return toret
        test_index = len(self.variant) - 1
        while test_index >= 0:
            test_char = self.variant[test_index]
            if test_char in allowable_characters:
                suffix += test_char
                allowable_characters = numbers
            else:
                break
            test_index -= 1
        suffix = ''.join(reversed(suffix))
        if len(suffix) < 1:
            for ss in special_suffixes:
                if self.variant.endswith(ss):
                    suffix = ss
                    break
            else:
                # TODO in future: Rather than return, at this point we could sequentially remove characters (similar to the last bit of this function, but with no suffix)
                return toret
        prefix = self.variant.split(suffix)[0]
        if suffix.startswith('19') or suffix.startswith('20'):
            suffix = suffix[2:]
            toret.append(prefix + suffix)
        while len(prefix) > 1:
            prefix = prefix[:-1]
            toret.append(prefix + suffix)
        return toret

    def unpack_timescales(self, timescales):
        if isinstance(timescales, WhiteDwarfAtmosphereDataset):
            return timescales # They are already in the right format (assuming that everything inside is also correct...maybe we should check? Probably overkill)
        toret = dict()
        for key, value in timescales.items():
            if value is None:
                pass
            elif isinstance(value, WhiteDwarfAtmosphereData):
                toret[key] = value # Looks like it's already in the right format
            else:
                toret[key] = WhiteDwarfAtmosphereData(value)
        return WhiteDwarfAtmosphereDataset(toret)

    def get_plottable_points(self, elements_to_plot, reference_element=None, included=True, scale_to_solar=False): # Switch included to False to get the excluded points
        if reference_element is None:
            reference_element = self.get_atmospheric_type().value # This should be set to zero in __init__
        elements_plotted, values_toret, upper_errors_toret, lower_errors_toret, upper_bounds_bools_toret, lower_bounds_bools_toret = self.abundances.get_plottable_points(elements_to_plot, reference_element, included, scale_to_solar)
        return elements_plotted, reference_element, values_toret, upper_errors_toret, lower_errors_toret, upper_bounds_bools_toret, lower_bounds_bools_toret

    def get_all_elements(self):
        return self.abundances.get_all_elements()

    def get_all_elements_overlapping_with_fit_dict(self, fit_dict):
        elements_toret = list()
        for el in self.get_all_elements():
            add_to_elements_toret = True
            for fit_key, fit_data in fit_dict.items():
                if el not in fit_data:
                    add_to_elements_toret = False
            if add_to_elements_toret:
                elements_toret.append(el)
        return elements_toret

    def get_spectral_type(self):
        return self.properties.get_property_data(mp.WDParameter.spectral_type)

    def get_atmospheric_type(self):
        return self.properties.get_property_data(mp.WDParameter.atmospheric_type)

    def get_teff(self):
        return self.properties.get_property_data(mp.WDParameter.temperature)

    def get_logg(self):
        return self.properties.get_property_data(mp.WDParameter.logg)

    def get_logq(self, timescale_type):
        return self.atmosphere_data.get_logq(timescale_type)

    def get_logq_in_solar_masses(self, timescale_type):
        logq = self.get_logq(timescale_type).value
        mass = self.get_mass().value
        return mass * (10**logq)

    def get_mass(self):
        return self.properties.get_property_data(mp.WDParameter.mass)

    def get_distance(self):
        return self.properties.get_property_data(mp.WDParameter.distance)

    def get_abundance(self, element):
        return self.abundances.get_abundance(element)

    def get_abundance_arrays(self, elements=ci.usual_elements):
        return self.abundances.get_abundance_arrays(elements)

    def get_timescales_as_array(self, timescale_type, elements=ci.usual_elements):
        return self.atmosphere_data.get_timescales_as_array(timescale_type, elements)

    def get_errors_for_present_elements_as_array(self, elements=ci.usual_elements):
        errors_to_use = list()
        for el in elements:
            relevant_data = self.get_abundance(el)
            if relevant_data is not None and relevant_data.included:
                if relevant_data.data_point_type == WhiteDwarfDataPointType.measurement and relevant_data.value not in [None, np.nan]:
                    errors_to_use.append(relevant_data.upper_error)
        return np.array(errors_to_use)

    def log_likelihood(self, model_result, min_likelihood):
        data_abundances_to_use = list()
        model_abundances_to_use = list()
        errors_to_use = list()
        # First check if we violate any bounds:
        for el, el_result in model_result.items():
            relevant_data = self.get_abundance(el)
            if relevant_data is not None and relevant_data.included:
                if relevant_data.data_point_type == WhiteDwarfDataPointType.upper_bound and el_result > relevant_data.value:
                    return 0.9*min_likelihood
                elif relevant_data.data_point_type == WhiteDwarfDataPointType.lower_bound and el_result < relevant_data.value:
                    return 0.9*min_likelihood
                elif relevant_data.data_point_type == WhiteDwarfDataPointType.measurement and relevant_data.value not in [None, np.nan]:
                    data_abundances_to_use.append(relevant_data.value)
                    model_abundances_to_use.append(el_result)
                    errors_to_use.append(relevant_data.upper_error)
                else:
                    pass
        data_abundances_to_use = np.array(data_abundances_to_use)
        model_abundances_to_use = np.array(model_abundances_to_use)
        errors_to_use = np.array(errors_to_use)
        neg = 1/((2*np.pi)*(errors_to_use**2))
        chi2 = ((data_abundances_to_use - model_abundances_to_use)**2)*(neg)*2*np.pi
        like = -0.5*np.sum(chi2 - np.log(neg))
        if np.isnan(like) or like == -np.inf:
            like = 0.9*min_likelihood
        return like

    def estimate_minimum_pollution_fraction(self, elements_to_consider=ci.all_elements):
        list_of_abundances_for_pol_frac = list()
        for el in elements_to_consider:
            try:
                el_data = self.abundances.abundance_data_dict[el]
            except KeyError:
                continue
            if el_data.data_point_type == WhiteDwarfDataPointType.measurement and el_data.included and el not in ci.atmospheric_types:
                list_of_abundances_for_pol_frac.append(el_data.value)
            elif el_data.data_point_type == WhiteDwarfDataPointType.lower_bound and el_data.included and el not in ci.atmospheric_types:
                list_of_abundances_for_pol_frac.append(el_data.value)
            else:
                pass
        return np.log10(sum([10**(a) for a in list_of_abundances_for_pol_frac]))

    def estimate_maximum_pollution_fraction(self, elements_to_consider=ci.all_elements):
        # Assume solar composition for missing elements, what's approximately the highest pollution fraction we could get for this system?
        safety_margin = 1  # guarantees that, even if the below calculation goes wrong somehow, we always return at least the minimum possible pollution fraction + this margin
        list_of_abundances_for_estimate = list()
        reference_element_priority = [
            ci.Element.Mg,
            ci.Element.Si,
            ci.Element.Ca,
            ci.Element.Fe, #Somewhat arbitrary below here
            ci.Element.Ni,
            ci.Element.Cr,
            ci.Element.Al,
            ci.Element.Ti,
            ci.Element.Na,
            ci.Element.O,
            ci.Element.C,
            ci.Element.N
        ]
        for element in ci.all_elements:
            if element not in reference_element_priority and element not in ci.atmospheric_types:
                reference_element_priority.append(element)
        ref_el = None
        ref_el_value = None
        for test_ref_el in reference_element_priority:
            try:
                el_data = self.abundances.abundance_data_dict[test_ref_el]
            except KeyError:
                continue
            if el_data.data_point_type == WhiteDwarfDataPointType.measurement and el_data.included and test_ref_el not in ci.atmospheric_types and el_data.value is not None:
                ref_el = test_ref_el
                ref_el_value = el_data.value
                break
        if ref_el is None:
            raise ValueError('White dwarf ' + self.name + ' has no measurements to normalise solar abundances against')
        list_of_abundances_for_pol_frac = list()
        for el in elements_to_consider:
            solar_rel_abundance = sa.get_solar_relative_abundance(el, ref_el)
            if solar_rel_abundance is None:
                solar_abundance = -np.inf # Effectively we ignore this element if we have to add a default abundance (we don't know what it is)
            else:
                solar_abundance = ref_el_value + sa.get_solar_relative_abundance(el, ref_el)
            try:
                el_data = self.abundances.abundance_data_dict[el]
                if el_data.data_point_type == WhiteDwarfDataPointType.measurement and el_data.included and el not in ci.atmospheric_types:
                    list_of_abundances_for_pol_frac.append(el_data.value)
                elif el_data.data_point_type == WhiteDwarfDataPointType.upper_bound and el_data.included and el not in ci.atmospheric_types:
                    # Add the minimum of either the upper bound or the solar abundance. Sometimes the upper bounds are really high
                    list_of_abundances_for_pol_frac.append(min(el_data.value, solar_abundance))
                elif el_data.data_point_type == WhiteDwarfDataPointType.lower_bound and el_data.included and el not in ci.atmospheric_types:
                    list_of_abundances_for_pol_frac.append(max(el_data.value, solar_abundance))
                else:
                    if el not in ci.atmospheric_types:
                        list_of_abundances_for_pol_frac.append(solar_abundance) # Won't this add the other of H and He?!
            except KeyError:
                if el not in ci.atmospheric_types:
                    list_of_abundances_for_pol_frac.append(solar_abundance)
        estimate = np.log10(sum([10**(a) for a in list_of_abundances_for_pol_frac]))
        return max(estimate, self.estimate_minimum_pollution_fraction(elements_to_consider) + safety_margin)

    def get_elements_present(self):
        return self.abundances.get_elements_present()

    def get_measurements_as_lists(self, return_included=True, return_excluded=False, elements=ci.all_elements, preserve_length=False):
        return self.abundances.get_measurements_as_lists(return_included, return_excluded, elements, preserve_length)

    def get_upper_bounds_as_lists(self, return_included=True, return_excluded=False, elements=ci.all_elements, preserve_length=False):
        return self.abundances.get_upper_bounds_as_lists(return_included, return_excluded, elements, preserve_length)

    def get_lower_bounds_as_lists(self, return_included=True, return_excluded=False, elements=ci.all_elements, preserve_length=False):
        return self.abundances.get_lower_bounds_as_lists(return_included, return_excluded, elements, preserve_length)

    def get_modellable_detected_elements(self, detectable_elements=ci.usual_elements):
        return self.abundances.get_modellable_detected_elements(detectable_elements)
        #toret = 0
        #for el in detectable_elements:
        #    abundance = self.get_abundance(el)

    def get_log_ratio(self, element1, element2):
        return self.abundances.get_log_ratio(element1, element2)

    def get_abundance_values_dict(self, elements=ci.usual_elements):
        toret = dict()
        for element in elements:
            abundance = self.get_abundance(element)
            if abundance is None:
                pass
            else:
                try:
                    if abundance.included:
                        toret[element] = abundance.value
                    else:
                        pass
                except:
                    pass
        return toret

    def get_error_values_dict(self, elements=ci.usual_elements):
        toret = dict()
        for element in elements:
            abundance = self.get_abundance(element)
            if abundance is None:
                pass
            else:
                try:
                    if abundance.included:
                        toret[element] = abundance.upper_error
                    else:
                        pass
                except:
                    pass
        return toret

    def get_timescale_values_dict(self, timescale_type, elements=ci.usual_elements):
        return self.atmosphere_data.get_timescale_values_dict(timescale_type, elements)

    def get_available_timescale_types(self):
        return self.atmosphere_data.get_available_timescale_types()

    def __str__(self):
        toret = '\n'
        toret += '--- ' + self.full_name() + ' ---\n\n'
        toret += str(self.properties)
        toret += '\n'
        toret += str(self.abundances)
        toret += '\n'
        toret += str(self.atmosphere_data)
        return toret

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return self.full_name() == other.full_name() and self.properties == other.properties and self.abundances == other.abundances and self.atmosphere_data == other.atmosphere_data

    def __neq__(self, other):
        return not self == other
