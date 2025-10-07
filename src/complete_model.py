#!/usr/bin/env python

import collections
import math
import numpy as np

import abundance_model as am
import atmosphere_model as atm
import chemistry_info as ci
import enhancement_model as em
import geology_info as gi
import live_data as ld
import physical_constants as pc


# t_sinceaccretion in Myr, t_disc in yr
def complete_model_calculation(
    fe_star,
    t_sinceaccretion,
    d_formation,
    z_formation,
    N_c,
    N_o,
    f_c,
    f_o,
    log_fragment_mass,
    t_disc,
    pressure,
    fO2,
    enhancement_model="NonEarthlike",
    consider_thermohaline=False,
    t_formation=1.5,
    normalise_abundances=True,
):

    diagnostics = dict()
    elements = ci.writeable_elements
    # This limit on fe_star exists because outside of this range,
    # ld._live_stellar_compositions[int(round(fe_star))] will give a KeyError (there are
    # 958 compositions)
    floored_fe_star = math.floor(fe_star)
    if 0 <= floored_fe_star <= 957:
        linear_d_formation = 10 ** (d_formation)
        abundances = am.get_all_abundances(
            elements, linear_d_formation, z_formation, t_formation, fe_star
        )
        disc_abundances = dict()
        for el_index, element in enumerate(elements):
            try:
                if el_index == 6:
                    # Mg is special
                    disc_abundances[element] = abundances[element]
                else:
                    disc_abundances[element] = (
                        abundances[element]
                        * ld._live_stellar_compositions[floored_fe_star][
                            el_index - 1 if el_index > 6 else el_index
                        ]
                    )
            except IndexError:
                # If the element is not one where we have stellar data, assume absent
                disc_abundances[element] = 0.0

        diagnostics["DiscAbundances"] = disc_abundances
        # This is to speed up performance by making sure we only need to fully
        # initialise the geo_model (and by extension the partitioning model) once
        if ld._geo_model is None:
            ld._geo_model = gi.GeologyModel(disc_abundances)
        else:
            ld._geo_model.reinit(disc_abundances)

        enhancement_model = em.EnhancementModel(enhancement_model)
        enhancements_dict, enhancements_diagnostics = (
            enhancement_model.find_enhancements(
                ld._geo_model,
                disc_abundances,
                elements,
                N_c,
                N_o,
                f_c,
                f_o,
                pressure,
                fO2,
                normalise_abundances,
            )
        )

        diagnostics["Enhancements"] = enhancements_diagnostics

        if enhancements_dict is None:
            return None, diagnostics

        planetesimal_abundance = list()
        for element in elements:
            toappend = enhancements_dict[element]
            if np.isnan(toappend):
                # This is important: later normalisations will fail if nans are present
                planetesimal_abundance.append(0)
            else:
                planetesimal_abundance.append(toappend)

        Hx = ld._live_Hx
        # Hx = ld._live_white_dwarf.get_atmospheric_type().value
        mu_X = np.array([ci.get_element_mass(el) for el in elements])
        planetesimal_abundance_arr = np.array(planetesimal_abundance)

        if np.isnan(planetesimal_abundance_arr).all():
            return None, diagnostics

        relative_mass_fractions = mu_X * planetesimal_abundance_arr
        total_rmf = np.linalg.norm(relative_mass_fractions)
        if total_rmf == 0.0:
            return None, diagnostics
        # Can get odd output if fragment_mass is an int, so cast to float
        fragment_mass = 10 ** float(log_fragment_mass)
        M_X = (relative_mass_fractions / total_rmf) * (fragment_mass / pc.M_Sun)
        # ^ The mass of the cvz is in units of solar mass (not kg), so this converts to
        # consistent units
        tau_X = ld._live_all_wd_timescales
        M_cvz = ld._live_M_cvz
        # M_cvz = ld._live_white_dwarf.get_logq_in_solar_masses(ld._live_timescale_type)

        # Teff = ld._live_white_dwarf.get_teff().value
        # logg = ld._live_white_dwarf.get_logg().value
        Teff = ld._live_teff
        logg = ld._live_logg
        result = atm.calculate_abundance_by_number(
            1000000 * t_sinceaccretion,
            t_disc,
            Hx,
            M_cvz,
            mu_X,
            M_X,
            tau_X,
            consider_thermohaline,
            Teff,
            logg,
        )

    else:
        raise ValueError("Metallicity must be between 0 and 958")
    elements_present_dict = collections.OrderedDict(zip(elements, result))
    toret = collections.OrderedDict()
    for element in elements:
        toret[element] = elements_present_dict.get(element)
    return toret, diagnostics


def example():
    import manager as mn  # Just in the example, we don't need it for the actual usage
    import timescale_interpolator as ti
    from argparse import Namespace

    manager = mn.Manager(
        Namespace(
            wd_data_filename="WDInputData.csv",
            stellar_compositions_filename="StellarCompositionsSortFE.csv",
            n_live_points=0,  # This argument shouldn't matter
            pollution_model_names=["Model_24"],
            enhancement_model="NonEarthlike",
        )
    )
    manager.publish_live_data(0, ti.TimescaleType.KoesterOvershoot)
    # ^ This is necessary because the complete_model ends up checking the live stellar
    # abundances

    fe_star_cl = 0
    t_sinceaccretion_cl = 13.7762587485591173
    feeding_zone_size_cl = 0.000237425605399639312
    d_formation_cl = -1.43213572254115262
    t_event_cl = 10 ** (2.28103474755766999)
    pollution_frac_cl = -5.85996906635276815

    fe_star = 478
    t_sinceaccretion = 1
    d_formation = -0.3
    z_formation = 0.05
    N_c = 0.2
    N_o = 0.01
    f_c = 0.8
    f_o = 0.0001
    pollutionfraction = -5
    t_disc = 3000000
    pressure = 21
    fO2 = -2
    enhancement_model = "NonEarthlike"
    t_formation = 1.5
    normalise_abundances = True
    snapshot_wd_atm = True
    fragment_mass = 1.95e21
    # ^ This turns out to be roughly equivalent to pollutionfraction = -5 in this case
    # (I calibrated it so that the old_result and the new_result are basically the same)

    old_result = complete_model_calculation_old(
        fe_star_cl,
        t_sinceaccretion_cl,
        d_formation_cl,
        feeding_zone_size_cl,
        0.17,
        0.01,
        0.17,
        0,
        pollution_frac_cl,
        t_event_cl,
        54,
        -2,
        enhancement_model,
        t_formation,
        normalise_abundances,
        snapshot_wd_atm,
    )
    new_result = complete_model_calculation(
        fe_star_cl,
        t_sinceaccretion_cl,
        d_formation_cl,
        feeding_zone_size_cl,
        0.17,
        0.01,
        0.17,
        0,
        1e21,
        t_event_cl,
        54,
        -2,
        enhancement_model,
        t_formation,
        normalise_abundances,
    )
    print(old_result)
    print(new_result)

    old_result = complete_model_calculation_old(
        fe_star,
        t_sinceaccretion,
        d_formation,
        z_formation,
        N_c,
        N_o,
        f_c,
        f_o,
        pollutionfraction,
        t_disc,
        pressure,
        fO2,
        enhancement_model,
        t_formation,
        normalise_abundances,
        snapshot_wd_atm,
    )

    new_result = complete_model_calculation(
        fe_star,
        t_sinceaccretion,
        d_formation,
        z_formation,
        N_c,
        N_o,
        f_c,
        f_o,
        fragment_mass,
        t_disc,
        pressure,
        fO2,
        enhancement_model,
        t_formation,
        normalise_abundances,
    )

    print(old_result)
    print(new_result)


def main():
    example()


if __name__ == "__main__":
    main()
