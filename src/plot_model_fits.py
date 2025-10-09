#!/usr/bin/env python

import chemistry_info as ci
import complete_model as cm
import graph_factory as gf
import manager as mn
import timescale_interpolator as ti


def plot_fits():
    graph_fac = gf.GraphFactory()

    elements_to_plot = ci.usual_elements
    timescale_type = ti.TimescaleType.KoesterOvershoot

    manager = mn.Manager()

    # Special logic bc this particular system is being run with a dedicated nebular
    # composition
    manager.wd_data_filename = "WDJ0916InputData.csv"
    manager.stellar_compositions_filename = "WDJ0916StellarCompositionV2.csv"
    manager.load_global_data()

    system_fits = {
        # - Stellar metallicity indices
        # - Time since Accretion/Myrs
        # - log(Formation Distance/AU)
        # - Feeding Zone Size/AU
        # - PCF PcrustF Fragment Core Fraction FCrustF
        # - log(Fragment Mass /kg)
        # - log(Accretion Event Timescale/Yrs)
        # - Pressure /GPa
        # - Oxygen Fugacity /ΔIW
        "SDSSJ0916+2540": {
            "MedianVals": [
                0,
                0.467526553916363,
                -0.848113622562445,
                0.088902947433765,
                None,
                None,
                None,
                None,
                19.3866588618516,
                2.37329810536358,
                None,
                None,
            ],
            "SS": [
                0,
                5,
                -0.848113622562445,
                0.088902947433765,
                None,
                None,
                None,
                None,
                19.3866588618516,
                7,
                None,
                None,
            ],
            "NoHeating": [
                0,
                0.467526553916363,
                2,
                0.088902947433765,
                None,
                None,
                None,
                None,
                19.3866588618516,
                2.37329810536358,
                None,
                None,
            ],
        }
    }
    print(system_fits)

    for system_name, fits in system_fits.items():
        print(system_name)
        white_dwarf_index, white_dwarf = manager.get_white_dwarf_by_name(system_name)
        print(white_dwarf_index)
        print(white_dwarf)
        manager.publish_live_data(white_dwarf_index, timescale_type)
        outputs = dict()
        for fit_name, parameter_values in fits.items():
            print(fit_name)
            print(parameter_values)

            fe_star = parameter_values[0]
            t_sinceaccretion = parameter_values[1]
            d_formation = parameter_values[2]
            z_formation = parameter_values[3]
            N_c = parameter_values[4]
            N_o = parameter_values[5]
            f_c = parameter_values[6]
            f_o = parameter_values[7]
            fragment_mass = parameter_values[8]
            t_disc = 10 ** parameter_values[9]
            pressure = parameter_values[10]
            fO2 = parameter_values[11]
            enhancement_model = "NonEarthlike"
            t_formation = 1.5
            normalise_abundances = True
            result = cm.complete_model_calculation(
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
            print(result)
            outputs[fit_name] = result[0]
        model_name = None
        extra_text_dict = None
        video = False
        fit_dict_prescaled = False
        hack_legend_to_only_show_pressure = False
        graph_fac.make_composition_plot_mk3(
            white_dwarf,
            elements_to_plot,
            outputs,
            None,
            model_name,
            extra_text_dict,
            video,
            fit_dict_prescaled,
            hack_legend_to_only_show_pressure,
        )
        graph_fac.make_composition_plot_mk3(
            white_dwarf,
            elements_to_plot,
            outputs,
            ci.Element.Mg,
            model_name,
            extra_text_dict,
            video,
            fit_dict_prescaled,
            hack_legend_to_only_show_pressure,
        )


def main():
    plot_fits()


if __name__ == "__main__":
    main()
