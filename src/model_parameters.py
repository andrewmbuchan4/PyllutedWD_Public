#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum
import numpy as np

minimum_likelihood = -1.0e90

class ModelParameter(Enum):
    metallicity = 0
    t_sinceaccretion = 1
    formation_distance = 2
    feeding_zone_size = 3
    parent_core_frac = 4
    parent_crust_frac = 5
    fragment_core_frac = 6
    fragment_crust_frac = 7
    pollution_frac = 8  # TODO: Rename this to pollution level
    accretion_timescale = 9
    pressure = 10
    oxygen_fugacity = 11

    def __str__(self):
        word_list = [word.capitalize() for word in self.name.split('_')]
        subbed_word_list = list()
        for word in word_list:
            if word == 'Frac':
                subbed_word_list.append('Fraction')
            elif word == 'T':
                subbed_word_list.append('Time')
            elif word == 'Sinceaccretion':
                subbed_word_list.append('Since')
                subbed_word_list.append('Accretion')
            else:
                subbed_word_list.append(word)
                if word in ['Core', 'Crust']:
                    subbed_word_list.append('Number')
        return ' '.join(subbed_word_list)

    def units(self):
        unit_dict = {
            'pressure': ('', 'GPa'),
            'oxygen_fugacity': ('', 'ΔIW'),
            't_sinceaccretion': ('', 'Myr'),
            'formation_distance': ('log(', 'AU)'),
            'feeding_zone_size': ('', 'AU)'),
            'accretion_timescale': ('log(', 'yr)'),
        }
        return unit_dict.get(self.name, ('', ''))

    def full_unit_string(self):
        pre_unit, post_unit = self.units()
        if post_unit != '':
            post_unit = ' /' + post_unit
        return pre_unit + str(self) + post_unit

class WDParameter(Enum):
    spectral_type = 0
    atmospheric_type = 1
    temperature = 2
    logg = 3
    mass = 4
    distance = 5
    logq = 6

    def __str__(self):
        return ' '.join([word.capitalize() for word in self.name.split('_')])

wd_parameter_units = {
    #WDParameter.spectral_type: None,
    WDParameter.temperature: 'K',
    #WDParameter.logg: None,
    WDParameter.mass: 'M_{Solar}',
    WDParameter.distance: 'pc'
    #WDParameter.atmospheric_type: None
}

wd_parameters_for_synthesis = [  # this pre-dates atmospheric_type!
    WDParameter.spectral_type,
    WDParameter.temperature,
    WDParameter.logg,
    WDParameter.mass,
    WDParameter.distance
]

model_parameter_strings = {
    ModelParameter.metallicity: "Stellar metallicity indices",
    ModelParameter.t_sinceaccretion: "Time since Accretion/Myrs",
    ModelParameter.formation_distance: "log(Formation Distance/AU)",
    ModelParameter.feeding_zone_size: "Feeding Zone Size/AU",
    ModelParameter.parent_core_frac: "Parent Core Fraction",
    ModelParameter.parent_crust_frac: "Parent Crust Fraction",
    ModelParameter.fragment_core_frac: "Fragment Core Fraction",
    ModelParameter.fragment_crust_frac: "Fragment Crust Fraction",
    ModelParameter.pollution_frac: "log(Pollution Fraction)",
    ModelParameter.accretion_timescale: "log(Accretion Event Timescale/Yrs)",
    ModelParameter.pressure: "Pressure /GPa",  # NB: This used to be log!
    ModelParameter.oxygen_fugacity: "Oxygen Fugacity /ΔIW"
}

# These must begin with 'Hierarchy'
# Warning: some logic will probably break if you have more than 10 levels
# manager.py assumes each level is just 1 digit
hierarchy_definitions_dict = {
    'Hierarchy_Basic': {
        0: [
            ModelParameter.metallicity,
            ModelParameter.t_sinceaccretion,
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale
        ],
        1: [
            ModelParameter.formation_distance
        ]
    },
    'Hierarchy_Default': {
        0: [
            ModelParameter.metallicity,
            ModelParameter.t_sinceaccretion,
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale
        ],
        1: [
            ModelParameter.formation_distance
        ],
        2: [
            ModelParameter.feeding_zone_size
        ],
        3: [
            ModelParameter.fragment_core_frac,
            ModelParameter.pressure,
            ModelParameter.oxygen_fugacity
        ]
    },
    'Hierarchy_Test': {
        0: [
            ModelParameter.metallicity,
            ModelParameter.t_sinceaccretion,
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale,
            ModelParameter.fragment_core_frac,
            ModelParameter.pressure,
            ModelParameter.oxygen_fugacity
        ],
        1: [
            ModelParameter.formation_distance  #This is a dummy - it's irrelevent
        ]
    },
    'Hierarchy_Test2': {
        0: [
            ModelParameter.metallicity,
            ModelParameter.t_sinceaccretion,
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale,
            ModelParameter.formation_distance,
            ModelParameter.feeding_zone_size
        ],
        1: [
            ModelParameter.fragment_core_frac,
            ModelParameter.pressure,
            ModelParameter.oxygen_fugacity
        ]
    },
    'Hierarchy_Test3': {
        0: [
            ModelParameter.metallicity,
            ModelParameter.t_sinceaccretion,
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale
        ],
        1: [
            ModelParameter.formation_distance,
            ModelParameter.feeding_zone_size,
            ModelParameter.fragment_core_frac,
            ModelParameter.pressure,
            ModelParameter.oxygen_fugacity
        ]
    },
    'Hierarchy_OM': {
        0: [
            ModelParameter.metallicity,
            ModelParameter.t_sinceaccretion,
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale
        ],
        1: [
            ModelParameter.formation_distance
        ],
        2: [
            ModelParameter.feeding_zone_size,
        ],
        3: [
            ModelParameter.fragment_core_frac,
            ModelParameter.fragment_crust_frac  # This isn't as flexible as it should be. fcf should be tried on its own, then fof added if fcf is added...
        ],
        4: [
            ModelParameter.parent_core_frac
        ],
        5: [
            ModelParameter.parent_crust_frac
        ]
    },
    'Hierarchy_OM_reduced': {
        0: [
            ModelParameter.metallicity,
            ModelParameter.t_sinceaccretion,
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale
        ],
        1: [
            ModelParameter.formation_distance
        ],
        2: [
            ModelParameter.feeding_zone_size,
        ],
        3: [
            ModelParameter.fragment_core_frac,
            ModelParameter.fragment_crust_frac
        ]
    },
    'Hierarchy_Earth': {
        0: [
            #ModelParameter.metallicity,  # Not used - invoke SolarCompositions.csv to fix initial composition at solar
            ModelParameter.t_sinceaccretion, #In principle these timescale ones could also be removed
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale,
            ModelParameter.formation_distance,
            #ModelParameter.fragment_core_frac, # Not used - we already know this in principle, so use correct default values
            ModelParameter.pressure,
            ModelParameter.oxygen_fugacity
        ],
        1: [
            ModelParameter.feeding_zone_size
        ]
    },
    # This one (below) should be used for Solar System objects, when the initial composition can be assumed to be solar and therefore removed as a parameter
    # Can also be used for other systems where you can fix the stellar composition somehow (eg wide binaries) in which case this becomes a misnomer (sorry)
    'Hierarchy_Meteorite': {
        0: [
            #ModelParameter.metallicity,  # Not used - invoke SolarCompositions.csv to fix initial composition at solar
            ModelParameter.t_sinceaccretion, #In principle these timescale ones could also be removed
            ModelParameter.pollution_frac,
            ModelParameter.accretion_timescale
        ],
        1: [
            ModelParameter.formation_distance
        ],
        2: [
            ModelParameter.feeding_zone_size
        ],
        3: [
            ModelParameter.fragment_core_frac,
            ModelParameter.pressure,
            ModelParameter.oxygen_fugacity
        ]
    }
}

model_definitions_dict = {
    'Model_24': {
        ModelParameter.metallicity: True,
        ModelParameter.t_sinceaccretion: True,
        ModelParameter.formation_distance: True,
        ModelParameter.feeding_zone_size: True,
        ModelParameter.parent_core_frac: True,
        ModelParameter.parent_crust_frac: True,
        ModelParameter.fragment_core_frac: True,
        ModelParameter.fragment_crust_frac: True,
        ModelParameter.pollution_frac: True,
        ModelParameter.accretion_timescale: True,
        ModelParameter.pressure: False,
        ModelParameter.oxygen_fugacity: False
    },
    'Model_4_equiv': {
        ModelParameter.metallicity: True,
        ModelParameter.t_sinceaccretion: True,
        ModelParameter.formation_distance: False,
        ModelParameter.feeding_zone_size: False,
        ModelParameter.parent_core_frac: False,
        ModelParameter.parent_crust_frac: False,
        ModelParameter.fragment_core_frac: True,
        ModelParameter.fragment_crust_frac: False,
        ModelParameter.pollution_frac: True,
        ModelParameter.accretion_timescale: True,
        ModelParameter.pressure: True,
        ModelParameter.oxygen_fugacity: True
    },
    'Model_Andy': {
        ModelParameter.metallicity: True,
        ModelParameter.t_sinceaccretion: True,
        ModelParameter.formation_distance: True,
        ModelParameter.feeding_zone_size: True,
        ModelParameter.parent_core_frac: True,
        ModelParameter.parent_crust_frac: True,
        ModelParameter.fragment_core_frac: True,
        ModelParameter.fragment_crust_frac: True,
        ModelParameter.pollution_frac: True,
        ModelParameter.accretion_timescale: True,
        ModelParameter.pressure: True,
        ModelParameter.oxygen_fugacity: False
    },
    'Model_Andy_No_Crust': {
        ModelParameter.metallicity: True,
        ModelParameter.t_sinceaccretion: True,
        ModelParameter.formation_distance: True,
        ModelParameter.feeding_zone_size: True,
        ModelParameter.parent_core_frac: True,
        ModelParameter.parent_crust_frac: False,
        ModelParameter.fragment_core_frac: True,
        ModelParameter.fragment_crust_frac: False,
        ModelParameter.pollution_frac: True,
        ModelParameter.accretion_timescale: True,
        ModelParameter.pressure: True,
        ModelParameter.oxygen_fugacity: False
    },
    'Model_Full': {
        ModelParameter.metallicity: True,
        ModelParameter.t_sinceaccretion: True,
        ModelParameter.formation_distance: True,
        ModelParameter.feeding_zone_size: True,
        ModelParameter.parent_core_frac: True,
        ModelParameter.parent_crust_frac: True,
        ModelParameter.fragment_core_frac: True,
        ModelParameter.fragment_crust_frac: True,
        ModelParameter.pollution_frac: True,
        ModelParameter.accretion_timescale: True,
        ModelParameter.pressure: True,
        ModelParameter.oxygen_fugacity: True
    },
    'Model_Full_No_Crust': {
        ModelParameter.metallicity: True,
        ModelParameter.t_sinceaccretion: True,
        ModelParameter.formation_distance: True,
        ModelParameter.feeding_zone_size: True,
        ModelParameter.parent_core_frac: False,
        ModelParameter.parent_crust_frac: False,
        ModelParameter.fragment_core_frac: True,
        ModelParameter.fragment_crust_frac: False,
        ModelParameter.pollution_frac: True,
        ModelParameter.accretion_timescale: True,
        ModelParameter.pressure: True,
        ModelParameter.oxygen_fugacity: True
    },
    'Model_Mantle_Only': {
        ModelParameter.metallicity: True,
        ModelParameter.t_sinceaccretion: True,
        ModelParameter.formation_distance: True,
        ModelParameter.feeding_zone_size: True,
        ModelParameter.parent_core_frac: False,
        ModelParameter.parent_crust_frac: False,
        ModelParameter.fragment_core_frac: False,
        ModelParameter.fragment_crust_frac: False,
        ModelParameter.pollution_frac: True,
        ModelParameter.accretion_timescale: True,
        ModelParameter.pressure: True,
        ModelParameter.oxygen_fugacity: True
    }
}

def model_uses_parameter(model, parameter):
    return model_definitions_dict[model].get(parameter, False)

def get_model_params_in_order():
    return [m for m in ModelParameter]

def parameter_indices(model):
    toret = dict()
    index = 0
    for param in get_model_params_in_order():
        if model_uses_parameter(model, param):
            toret[param] = index
            index += 1
    return toret

default_values = {
    'NonEarthlike': {
        ModelParameter.metallicity: 479, #Assume average Fe/H, i.e. go halfway through the 958 indices
        ModelParameter.t_sinceaccretion: 0,
        ModelParameter.formation_distance: 2,
        ModelParameter.feeding_zone_size: 0.05,
        ModelParameter.parent_core_frac: 0.17,
        ModelParameter.parent_crust_frac: 0.01,
        ModelParameter.fragment_core_frac: None,
        ModelParameter.fragment_crust_frac: 0,
        ModelParameter.pollution_frac: -6.5, # This is a made-up but not unreasonable value
        ModelParameter.accretion_timescale: 5.6,
        ModelParameter.pressure: 54,
        ModelParameter.oxygen_fugacity: -2
    },
    'Earthlike': {
        ModelParameter.metallicity: 479, #Assume average Fe/H, i.e. go halfway through the 958 indices
        ModelParameter.t_sinceaccretion: 0,
        ModelParameter.formation_distance: 2,
        ModelParameter.feeding_zone_size: 0.05,
        ModelParameter.parent_core_frac: 0.17,
        ModelParameter.parent_crust_frac: 0.01,
        ModelParameter.fragment_core_frac: 0.17,
        ModelParameter.fragment_crust_frac: 0,
        ModelParameter.pollution_frac: -6.5, # This is a made-up but not unreasonable value
        ModelParameter.accretion_timescale: 5.6,
        ModelParameter.pressure: 54,
        ModelParameter.oxygen_fugacity: -2
    },
    'MantleOnly': {  # The same as NEL, but we assume no core (unless this is overridden by the prior)
        ModelParameter.metallicity: 479, #Assume average Fe/H, i.e. go halfway through the 958 indices
        ModelParameter.t_sinceaccretion: 0,
        ModelParameter.formation_distance: 2,
        ModelParameter.feeding_zone_size: 0.05,
        ModelParameter.parent_core_frac: 0.17,
        ModelParameter.parent_crust_frac: 0.01,
        ModelParameter.fragment_core_frac: 0,
        ModelParameter.fragment_crust_frac: 0,
        ModelParameter.pollution_frac: -6.5, # This is a made-up but not unreasonable value
        ModelParameter.accretion_timescale: 5.6,
        ModelParameter.pressure: 54,
        ModelParameter.oxygen_fugacity: -2
    },
    'EarthMantle': {
        ModelParameter.metallicity: 0, #We should use SolarComposition.csv, which has only 1 entry (the Sun) at index 0
        ModelParameter.t_sinceaccretion: 0,
        ModelParameter.formation_distance: 2,
        ModelParameter.feeding_zone_size: 0.05,
        ModelParameter.parent_core_frac: 0.17,
        ModelParameter.parent_crust_frac: 0.01,
        ModelParameter.fragment_core_frac: 0, # If we feed it mantle, we can just tell it it's mantle
        ModelParameter.fragment_crust_frac: 0,
        ModelParameter.pollution_frac: -6.5, # This is a made-up but not unreasonable value
        ModelParameter.accretion_timescale: 5.6,
        ModelParameter.pressure: 54,
        ModelParameter.oxygen_fugacity: -2
    },
    'Meteorite': {
        ModelParameter.metallicity: 0, #We should use SolarComposition.csv, which has only 1 entry (the Sun) at index 0
        ModelParameter.t_sinceaccretion: 0,
        ModelParameter.formation_distance: 2,
        ModelParameter.feeding_zone_size: 0.05,
        ModelParameter.parent_core_frac: 0.17,
        ModelParameter.parent_crust_frac: 0.01,
        ModelParameter.fragment_core_frac: 0.17,
        ModelParameter.fragment_crust_frac: 0,
        ModelParameter.pollution_frac: -6.5, # This is a made-up but not unreasonable value
        ModelParameter.accretion_timescale: 5.6,
        ModelParameter.pressure: 54,
        ModelParameter.oxygen_fugacity: -2
    }
}
