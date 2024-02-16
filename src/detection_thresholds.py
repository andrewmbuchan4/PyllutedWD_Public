#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chemistry_info as ci

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
            ci.Element.C: (0.000384039372, -(14.1531876-2.405138313)),  # Duplicating O, but since C is not used in modelling this choice makes no difference
            ci.Element.N: (0.000384039372, -(14.1531876-2.405138313))  # Duplicating O, but since C is not used in modelling this choice makes no difference
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
            ci.Element.C: (0.000416915147, -(16.6297136-2.103844062)),  # Duplicating O, but since C is not used in modelling this choice makes no difference
            ci.Element.N: (0.000416915147, -(16.6297136-2.103844062))  # Duplicating O, but since C is not used in modelling this choice makes no difference
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
            ci.Element.C: (0.00038403937, -(12.5511276-2.405138313)),  # Duplicating O, but since C is not used in modelling this choice makes no difference
            ci.Element.N: (0.00038403937, -(12.5511276-2.405138313))  # Duplicating O, but since C is not used in modelling this choice makes no difference
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
            ci.Element.C: (0.000416915147, -(15.0276536-2.103844062)),  # Duplicating O, but since C is not used in modelling this choice makes no difference
            ci.Element.N: (0.000416915147, -(15.0276536-2.103844062))  # Duplicating O, but since C is not used in modelling this choice makes no difference
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
# The Default set above are calculated from MWDD (03/08/22) looking at all available WDs with pollution recorded and specifying 2 coordinates that the threshold seems to pass through (from visual inspection, leaving a ~0.5 dex margin) These are:
# Actually from WDMS, seems to be more complete!?
# Al: DA: (5000, -10) and (22500, -7) DB: (5000, -9.2) and (24000, -6)
# Ti: DA: (5000, -10.4) and (13500, -8) DB: (10000, -11.2) and (15000, -10.6)
# Ca: DA: (5000, -10.4) and (24000, -6.6) DB: (7500, -10.8) and (20000, -8.2) where WDJ1922+0233 is an outlier
# Ni: DA: (5500, -10) and (22500, -6.8) DB: (5000, -10) and (17500, -8) where WD1425+540 is an outlier
# Fe: DA: (5000, -9.4) and (22500, -5.6) DB: (5000, -10) and (22500, -6.6)
# Cr: DA: (5000, -10) and (23000, -6) DB: (5000, -10.6) and (22500, -8.2) where GD378 is an outlier
# Mg: DA: (2500, -9.4) and (15000, -7.2) DB: (5000, -9.4) and (22500, -6.2) where WD1425+540 is an outlier
# Si: DA: (7500, -7.4) and (20000, -6.5) where WD2230-125 is an outlier DB: (5000, -8.1) and (24000, -8.1)
# Na: DA: (5000, -8.1) and (24000, -8.1) (only 1 data point for this!) DB: (6000, -10) and (17000, -6) where WDJ2147-4035 is an outlier
# O: DA: (5000, -5.6) and (22500, -5.6) DB: (5000, -6.7) and (24000, -6.7)
# C: DA: (5000, -7.6) and (22500, -7.6) DB: (5000, -8) and (24000, -8)
# N: DA: (5000, -6) and (22500, -6) no data points! basing this on lowest upper bound DB: (5000, -8.2) and (24000, -8.2) only one data point!
