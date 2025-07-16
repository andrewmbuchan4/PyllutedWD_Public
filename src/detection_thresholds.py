#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import chemistry_info as ci

min_mdot_bank = {
    # Calculate this using the find_mdot_cutoff_as_function_of_teff function in find_synthetic_wd_mass_cutoff.py
    'Hollands': {
        'DB': {
            'gradient': 0.00102,
            'yintercept': 15.23103026,
            'leeway': 5 # Seems a bit high - but I've calibrated it and seems to work
        }
    },
    'Default': {
        'DA': { # These values are not wrong, but are a bit inefficient because it means we will likely forward model a lot of systems with no detectable pollution
            'gradient': 0,
            'yintercept': -np.inf,
            'leeway': 0
        }
    },
    'Realistic': {
         # These values are not wrong, but are a bit inefficient because it means we will likely forward model a lot of systems with no detectable pollution
        'DA': {
            'gradient': 0,
            'yintercept': -np.inf,
            'leeway': 0
        },
        'DB': {
            'gradient': 0,
            'yintercept': -np.inf,
            'leeway': 0
        }
    },
    'ELB_DT': {
        'DA': {
            'gradient': -0.00004963316712,
            'yintercept': 18.68687806,
            'leeway': 4 # This is calibrated VERY roughly!
        },
        'DB': {
            'gradient': 0.0001503674854,
            'yintercept': 14.18686497,
            'leeway': 0.5 # This is calibrated VERY roughly!
        }
    }
}

threshold_bank = { # Used for defining detection thresholds if observation_type == ObservationType.TeffIndividualElementCutoff
    'Default': {
        'DA': {
            ci.Element.Al: (0.000384039372, -(14.1531876-1.337685817)),  # 1st element is the gradient, 2nd element is the y-intercept
            ci.Element.Ti: (0.000384039372, -(14.1531876-0.579460997)),
            ci.Element.Ca: (0.000384039372, -14.1531876),
            ci.Element.Ni: (0.000384039372, -(14.1531876-1.907685817)),
            ci.Element.Fe: (0.000384039372, -(14.1531876-2.12545136)),
            ci.Element.Cr: (0.000384039372, -(14.1531876-1.17545136)),
            ci.Element.Mg: (0.000384039372, -(14.1531876-0.424720522)),
            ci.Element.Si: (0.000384039372, -(14.1531876-0.873705764)),
            ci.Element.Na: (0.000384039372, -(14.1531876-1.532413312)),
            ci.Element.O: (0.000384039372, -(14.1531876-2.405138313)),
            ci.Element.C: (0.000384039372, -(14.1531876-2.405138313)),  # Duplicating O
            ci.Element.N: (0.000384039372, -(14.1531876-2.405138313)),  # Duplicating O
            ci.Element.S: (0.000384039372, -(14.1531876-2.405138313)) #Duplicating O

        },
        'DB': {
            ci.Element.Al: (0.000416915147, -(16.6297136-0.955322979)),
            ci.Element.Ti: (0.000416915147, -(16.6297136-0.025126899)),
            ci.Element.Ca: (0.000416915147, -16.6297136),
            ci.Element.Ni: (0.000416915147, -(16.6297136-1.119121373)),
            ci.Element.Fe: (0.000416915147, -(16.6297136-1.705460667)),
            ci.Element.Cr: (0.000416915147, -(16.6297136-0.338201724)),
            ci.Element.Mg: (0.000416915147, -(16.6297136-1.312608751)),
            ci.Element.Si: (0.000416915147, -(16.6297136-0.808628868)),
            ci.Element.Na: (0.000416915147, -(16.6297136-2.142770512)),
            ci.Element.O: (0.000416915147, -(16.6297136-2.103844062)),
            ci.Element.C: (0.000416915147, -(16.6297136-2.103844062)),  # Duplicating O
            ci.Element.N: (0.000416915147, -(16.6297136-2.103844062)),  # Duplicating O
            ci.Element.S: (0.000416915147, -(16.6297136-2.103844062))  # Duplicating O
        }
    },
    '560mA': {
        'DA': {
            ci.Element.Al: (0.00038403937, -(12.5511276-1.337685817)),  # 1st element is the gradient, 2nd element is the y-intercept
            ci.Element.Ti: (0.00038403937, -(12.5511276-0.579460997)),
            ci.Element.Ca: (0.00038403937, -12.5511276),
            ci.Element.Ni: (0.00038403937, -(12.5511276-1.907685817)),
            ci.Element.Fe: (0.00038403937, -(12.5511276-2.12545136)),
            ci.Element.Cr: (0.00038403937, -(12.5511276-1.17545136)),
            ci.Element.Mg: (0.00038403937, -(12.5511276-0.424720522)),
            ci.Element.Si: (0.00038403937, -(12.5511276-0.873705764)),
            ci.Element.Na: (0.00038403937, -(12.5511276-1.532413312)),
            ci.Element.O: (0.00038403937, -(12.5511276-2.405138313)),
            ci.Element.C: (0.00038403937, -(12.5511276-2.405138313)),  # Duplicating O
            ci.Element.N: (0.00038403937, -(12.5511276-2.405138313))  # Duplicating O
        },
        'DB': {
            ci.Element.Al: (0.000416915147, -(15.0276536-0.955322979)),
            ci.Element.Ti: (0.000416915147, -(15.0276536-0.025126899)),
            ci.Element.Ca: (0.000416915147, -15.0276536),
            ci.Element.Ni: (0.000416915147, -(15.0276536-1.119121373)),
            ci.Element.Fe: (0.000416915147, -(15.0276536-1.705460667)),
            ci.Element.Cr: (0.000416915147, -(15.0276536-0.338201724)),
            ci.Element.Mg: (0.000416915147, -(15.0276536-1.312608751)),
            ci.Element.Si: (0.000416915147, -(15.0276536-0.808628868)),
            ci.Element.Na: (0.000416915147, -(15.0276536-2.142770512)),
            ci.Element.O: (0.000416915147, -(15.0276536-2.103844062)),
            ci.Element.C: (0.000416915147, -(15.0276536-2.103844062)),  # Duplicating O
            ci.Element.N: (0.000416915147, -(15.0276536-2.103844062))  # Duplicating O
        }
    },
    'v2': {
        'DA': {
            ci.Element.Al: (0.00017, -11),  # 1st element is the gradient, 2nd element is the y-intercept
            ci.Element.Ti: (0.00013, -11),
            ci.Element.Ca: (0.0002, -11.4),
            ci.Element.Ni: (0.000188, -11.4),
            ci.Element.Fe: (0.000215, -10.5),
            ci.Element.Cr: (0.00015, -10.1),
            ci.Element.Mg: (0.00012, -9.5),
            ci.Element.Si: (0, -8.3),
            ci.Element.Na: (0.0002, -10.5),
            ci.Element.O: (0, -6.2),
            ci.Element.C: (0, -5),
            ci.Element.N: (0, -5)
        },
        'DB': {
            ci.Element.Al: (0.00017, -11),
            ci.Element.Ti: (0.00013, -12.5),
            ci.Element.Ca: (0.00025, -13),
            ci.Element.Ni: (0.000188, -12.5),
            ci.Element.Fe: (0.000194, -11.1),
            ci.Element.Cr: (0.0002, -12.5),
            ci.Element.Mg: (0.00018, -11),
            ci.Element.Si: (0, -8.3),
            ci.Element.Na: (0.0003, -13),
            ci.Element.O: (0, -7),
            ci.Element.C: (0, -5),
            ci.Element.N: (0, -5)
        }
    },
    'Hollands': { # For present purposes I only care about the Ca, Fe and Mg thresholds. Everything else is set to arbitrary levels (makes no difference). Could expand later!
        'DB': {
            ci.Element.Al: (0, 0),
            ci.Element.Ti: (0, 0),
            ci.Element.Ca: (0.00101, -16.6),
            ci.Element.Ni: (0, 0),
            ci.Element.Fe: (0.0011, -16.2),
            ci.Element.Cr: (0, 0),
            ci.Element.Mg: (0.00102, -15.2),
            ci.Element.Si: (0, 0),
            ci.Element.Na: (0, 0),
            ci.Element.O: (0, 0),
            ci.Element.C: (0, 0),
            ci.Element.N: (0, 0)
        }
    },
    'HollandsCr': { # The same as above but now with Cr detectable, to see if it matters
        'DB': {
            ci.Element.Al: (0, 0),
            ci.Element.Ti: (0, 0),
            ci.Element.Ca: (0.00101, -16.6),
            ci.Element.Ni: (0, 0),
            ci.Element.Fe: (0.0011, -16.2),
            ci.Element.Cr: (0.00098, -15.8),
            ci.Element.Mg: (0.00102, -15.2),
            ci.Element.Si: (0, 0),
            ci.Element.Na: (0, 0),
            ci.Element.O: (0, 0),
            ci.Element.C: (0, 0),
            ci.Element.N: (0, 0)
        }
    },
    'ELB_DT': { #These are the 'real', not the conservative, limits. Interpolated between 10,000 and 20,000 K
        'DA': {
            ci.Element.C: (-0.0001, -7),
            ci.Element.N: (-0.00025, -4),
            ci.Element.O: (0.0001, -9),
            ci.Element.Mg: (0.00005, -6.5),
            ci.Element.Al: (-0.00005, -5),
            ci.Element.Si: (-0.0002, -3),
            ci.Element.P: (-0.00005, -7.5),
            ci.Element.S: (-0.00005, -6.5),
            ci.Element.Ca: (0.00025, -11),
            ci.Element.Fe: (0, -6.5),
            ci.Element.Ni: (0, -7),
            ci.Element.Cu: (0.00005, -10),
            ci.Element.Cr: (0, 0),
            ci.Element.Na: (0, 0),
            ci.Element.Ti: (0, 0)
        },
        'DB': {
            ci.Element.C: (0.00015, -12.5),
            ci.Element.N: (0.00005, -10),
            ci.Element.O: (0.0003, -13.5),
            ci.Element.Mg: (0.0002, -10.5),
            ci.Element.Al: (0.00005, -8),
            ci.Element.Si: (0, -8),
            ci.Element.P: (0.00015, -12.5),
            ci.Element.S: (0.0001, -11),
            ci.Element.Ca: (0.0003, -13),
            ci.Element.Fe: (0.00015, -10),
            ci.Element.Ni: (0.0002, -11.5),
            ci.Element.Cu: (0.00025, -14.5),
            ci.Element.Cr: (0, 0),
            ci.Element.Na: (0, 0),
            ci.Element.Ti: (0, 0)
        }
    }
    #'v1': {
    #    'DA': {
    #        ci.Element.Al: (0.0001714, -10.86),  # 1st element is the gradient, 2nd element is the y-intercept, to 4sf
    #        ci.Element.Ti: (0.0002824, -11.81),
    #        ci.Element.Ca: (0.0002, -11.4),
    #        ci.Element.Ni: (0.0001882, -11.04),
    #        ci.Element.Fe: (0.0002171, -10.49),
    #        ci.Element.Cr: (0.0002222, -11.11),
    #        ci.Element.Mg: (0.000176, -9.84),
    #        ci.Element.Si: (0.000072, -7.94),
    #        ci.Element.Na: (0, -8.1),
    #        ci.Element.O: (0, -5.6),
    #        ci.Element.C: (0, -7.6),
    #        ci.Element.N: (0, -6)
    #    },
    #    'DB': {
    #        ci.Element.Al: (0.0001684, -10.04),
    #        ci.Element.Ti: (0.00012, -12.4),
    #        ci.Element.Ca: (0.000208, -12.36),
    #        ci.Element.Ni: (0.00016, -10.8),
    #        ci.Element.Fe: (0.0001942, -10.97),
    #        ci.Element.Cr: (0.0001371, -11.29),
    #        ci.Element.Mg: (0.0001828, -10.31),
    #        ci.Element.Si: (0, -8.1),
    #        ci.Element.Na: (0.0003636, -12.18),
    #        ci.Element.O: (0, -6.7),
    #        ci.Element.C: (0, -8),
    #        ci.Element.N: (0, -8.2)
    #    }
    #}
}
