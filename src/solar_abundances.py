#!/usr/bin/env python

import chemistry_info as ci

solar_ratiod_to_H = {
    # Asplund 2021
    ci.Element.H: 12.00,
    ci.Element.He: 10.914,
    ci.Element.Li: 0.96,
    ci.Element.Be: 1.38,
    ci.Element.B: 2.70,
    ci.Element.C: 8.46,
    ci.Element.N: 7.83,
    ci.Element.O: 8.69,
    ci.Element.F: 4.40,
    ci.Element.Ne: 8.06,
    ci.Element.Na: 6.22,
    ci.Element.Mg: 7.55,
    ci.Element.Al: 6.43,
    ci.Element.Si: 7.51,
    ci.Element.P: 5.41,
    ci.Element.S: 7.12,
    ci.Element.Cl: 5.31,
    ci.Element.Ar: 6.38,
    ci.Element.K: 5.07,
    ci.Element.Ca: 6.30,
    ci.Element.Sc: 3.14,
    ci.Element.Ti: 4.97,
    ci.Element.V: 3.90,
    ci.Element.Cr: 5.62,
    ci.Element.Mn: 5.42,
    ci.Element.Fe: 7.46,
    ci.Element.Co: 4.94,
    ci.Element.Ni: 6.20,
    ci.Element.Cu: 4.18,
    ci.Element.Zn: 4.56,
    ci.Element.Ga: 3.02,
    ci.Element.Ge: 3.62,
    ci.Element.As: None,
    ci.Element.Se: None,
    ci.Element.Br: None,
    ci.Element.Kr: 3.12,
    ci.Element.Rb: 2.32,
    ci.Element.Sr: 2.83,
    ci.Element.Y: 2.21,
    ci.Element.Zr: 2.59,
    ci.Element.Nb: 1.47,
    ci.Element.Mo: 1.88,
    ci.Element.Ru: 1.75,
    ci.Element.Rh: 0.78,
    ci.Element.Pd: 1.57,
    ci.Element.Ag: 0.96,
    ci.Element.Cd: None,
    ci.Element.In: 0.80,
    ci.Element.Sn: 2.02,
    ci.Element.Sb: None,
    ci.Element.Te: None,
    ci.Element.I: None,
    ci.Element.Xe: 2.22,
    ci.Element.Cs: None,
    ci.Element.Ba: 2.27,
    ci.Element.La: 1.11,
    ci.Element.Ce: 1.58,
    ci.Element.Pr: 0.75,
    ci.Element.Nd: 1.42,
    ci.Element.Sm: 0.95,
    ci.Element.Eu: 0.52,
    ci.Element.Gd: 1.08,
    ci.Element.Tb: 0.31,
    ci.Element.Dy: 1.10,
    ci.Element.Ho: 0.48,
    ci.Element.Er: 0.93,
    ci.Element.Tm: 0.11,
    ci.Element.Yb: 0.85,
    ci.Element.Lu: 0.10,
    ci.Element.Hf: 0.85,
    ci.Element.Ta: None,
    ci.Element.W: 0.79,
    ci.Element.Re: None,
    ci.Element.Os: 1.35,
    ci.Element.Ir: None,
    ci.Element.Pt: None,
    ci.Element.Au: 0.91,
    ci.Element.Hg: None,
    ci.Element.Tl: 0.92,
    ci.Element.Pb: 1.95,
    ci.Element.Bi: None,
    ci.Element.Th: 0.03,
    ci.Element.U: None,
}

for el in solar_ratiod_to_H:
    # These are given with an offset of 12, which we need to remove
    try:
        solar_ratiod_to_H[el] -= 12
    except TypeError:
        # For Nones
        pass

# The 1% and 99% percentiles from
# https://docs.google.com/spreadsheets/d/1B2WXfFtx3KhP4Ucl2kmJGmhq73ldJqqZiofG8m10n48

upper_X_ratiod_to_solar = {
    ci.Element.H: {
        # This is the reference element
        ci.Element.C: 0.3745,
        ci.Element.N: 0.5445,
        ci.Element.O: 0.45,
        ci.Element.Na: 0.57,
        ci.Element.Mg: 0.34,
        ci.Element.Al: 0.44,
        ci.Element.Si: 0.38,
        ci.Element.Ca: 0.4,
        ci.Element.Ti: 0.3741,
        ci.Element.V: 0.3545,
        ci.Element.Cr: 0.3841,
        ci.Element.Mn: 0.4841,
        ci.Element.Fe: 0.4,
        ci.Element.Ni: 0.45,
        ci.Element.Y: 0.5141,
    },
    # ci.Element.He:{
    #    ci.Element.C: 0.3745,
    #    ci.Element.N: 0.5445,
    #    ci.Element.O: 0.45,
    #    ci.Element.Na: 0.57,
    #    ci.Element.Mg: 0.34,
    #    ci.Element.Al: 0.44,
    #    ci.Element.Si: 0.38,
    #    ci.Element.Ca: 0.4,
    #    ci.Element.Ti: 0.3741,
    #    ci.Element.V: 0.3545,
    #    ci.Element.Cr: 0.3841,
    #    ci.Element.Mn: 0.4841,
    #    ci.Element.Fe: 0.4,
    #    ci.Element.Ni: 0.45,
    #    ci.Element.Y: 0.5141
    # },
    ci.Element.Mg: {
        ci.Element.C: 0.2,
        ci.Element.N: 0.3841,
        ci.Element.O: 0.37,
        ci.Element.Na: 0.25,
        ci.Element.Mg: 0,
        ci.Element.Al: 0.1341,
        ci.Element.Si: 0.11,
        ci.Element.Ca: 0.2282,
        ci.Element.Ti: 0.18,
        ci.Element.V: 0.2082,
        ci.Element.Cr: 0.16,
        ci.Element.Mn: 0.1941,
        ci.Element.Fe: 0.17,
        ci.Element.Ni: 0.15,
        ci.Element.Y: 0.42,
    },
}
lower_X_ratiod_to_solar = {
    ci.Element.H: {
        ci.Element.C: -0.42,
        ci.Element.N: -0.54,
        ci.Element.O: -0.2145,
        ci.Element.Na: -0.5141,
        ci.Element.Mg: -0.39,
        ci.Element.Al: -0.45,
        ci.Element.Si: -0.4041,
        ci.Element.Ca: -0.42,
        ci.Element.Ti: -0.3341,
        ci.Element.V: -0.4045,
        ci.Element.Cr: -0.5541,
        ci.Element.Mn: -0.7582,
        ci.Element.Fe: -0.5182,
        ci.Element.Ni: -0.4982,
        ci.Element.Y: -0.5482,
    },
    # ci.Element.He:{
    #    ci.Element.C: -0.42,
    #    ci.Element.N: -0.54,
    #    ci.Element.O: -0.2145,
    #    ci.Element.Na: -0.5141,
    #    ci.Element.Mg: -0.39,
    #    ci.Element.Al: -0.45,
    #    ci.Element.Si: -0.4041,
    #    ci.Element.Ca: -0.42,
    #    ci.Element.Ti: -0.3341,
    #    ci.Element.V: -0.4045,
    #    ci.Element.Cr: -0.5541,
    #    ci.Element.Mn: -0.7582,
    #    ci.Element.Fe: -0.5182,
    #    ci.Element.Ni: -0.4982,
    #    ci.Element.Y: -0.5482
    # },
    ci.Element.Mg: {
        ci.Element.C: -0.2141,
        ci.Element.N: -0.26,
        ci.Element.O: -0.0841,
        ci.Element.Na: -0.1841,
        ci.Element.Mg: 0,
        ci.Element.Al: -0.26,
        ci.Element.Si: -0.08,
        ci.Element.Ca: -0.08,
        ci.Element.Ti: -0.02,
        ci.Element.V: -0.15,
        ci.Element.Cr: -0.23,
        ci.Element.Mn: -0.4541,
        ci.Element.Fe: -0.22,
        ci.Element.Ni: -0.16,
        ci.Element.Y: -0.2982,
    },
}


def get_solar_relative_abundance(element1, element2):
    try:
        return solar_ratiod_to_H[element1] - solar_ratiod_to_H[element2]
    except (TypeError, KeyError):
        return None


def scale_abundances_to_solar(
    abundance_dict, reference_element, normalisation_element=None, descale=False
):
    # The difference between the reference_element and normalisation_element is subtle
    # but important.
    # The reference_element is ultimately going to be the element that appears in the
    # denominator on the y axis of the various composition plots.
    # The problem I ran into is that, when the reference_element is the atmospheric type
    # of the WD (ie Hx).
    # (ie we're essentially dealing with absolute abundances)
    # The vertical position can (not always) be meaningless. It's meaningless if, for
    # example, we are plotting the composition of accreted material.
    # Only the relative abundances have any meaning - the vertical offset is arbitrary.
    # So I just fixed it to match the Mg abundance. BUT
    # This actually introduced further confusion because the uncertainty remains
    # non-zero, which makes no sense - if the Mg abundance is always set to a certain
    # value, how could there be any uncertainty?
    # What I'd done was show the uncertainty on the absolute abundances in a context
    # where we should only care about the relative abundances
    # The solution is that the scaling to Mg needs to happen at sampling time.
    # So we end up with an (optional) normalisation_element (Mg in this case) to go
    # alongside the reference_element (Hx).
    # descale is a separate option for if you have a set of abundances that has already
    # been scaled, and you want to undo that scaling!
    if abundance_dict is None:
        return None
    scaled_abundances = dict()
    sign = 1 if descale else -1
    for element, abundance in abundance_dict.items():
        try:
            scaled_abundances[element] = (
                abundance - abundance_dict[reference_element]
            ) + (sign * get_solar_relative_abundance(element, reference_element))
        except KeyError:
            # Assume that this happened because the reference_element wasn't one of the
            # ones we were modelling, so is the one that they were being modelled
            # relative to
            # i.e. abundance doesn't need any scaling apart from the solar scaling
            if normalisation_element is None:
                scaled_abundances[element] = abundance + (
                    sign * get_solar_relative_abundance(element, reference_element)
                )
            else:
                scaled_abundances[element] = (
                    abundance - abundance_dict[normalisation_element]
                ) + (
                    sign * get_solar_relative_abundance(element, normalisation_element)
                )
    return scaled_abundances


def main():
    print(get_solar_relative_abundance(ci.Element.O, ci.Element.Ca))


if __name__ == "__main__":
    main()
