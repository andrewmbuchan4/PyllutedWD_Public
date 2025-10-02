import csv
import numpy as np

import abundance_model as am
import disc_model as dm
import chemistry_info as ci
import graph_factory as gf
import live_data as ld
import solar_abundances as sa


def load_generic_float_data_csv(input_filename):
    with open("../data/" + input_filename) as generic_csv:
        generic_list = [row for row in csv.reader(generic_csv)]
    generic_array = np.asarray(generic_list)
    return generic_array.astype(float)


def get_abundances(d_formation, z_formation, t_formation, stellar_metallicity_index):
    ld._live_stellar_compositions = load_generic_float_data_csv(
        "StellarCompositionsSortFE.csv"
    )
    abundances = am.get_all_abundances(
        ci.usual_elements,
        d_formation,
        z_formation,
        t_formation,
        stellar_metallicity_index,
    )
    return abundances


def output_el_el_values(
    element1, element2, z_formation, t_formation, stellar_metallicity_index
):
    solar_ratio = sa.get_solar_relative_abundance(element1, element2)
    o_ca_values = list()
    log_elel_values = list()
    normalised_elel_values = list()
    d_formation_vals = np.linspace(10**-1, 10**1, 101)
    T_vals = list()
    for d_formation in d_formation_vals:
        abundances = get_abundances(
            d_formation, z_formation, t_formation, stellar_metallicity_index
        )
        disc_abundances = dict()
        for el_index, element in enumerate(ci.usual_elements):
            if el_index == 6:
                # Mg is special
                disc_abundances[element] = abundances[element]
            else:
                disc_abundances[element] = (
                    abundances[element]
                    * ld._live_stellar_compositions[stellar_metallicity_index][
                        el_index - 1 if el_index > 6 else el_index
                    ]
                )
        try:
            o_ca = disc_abundances[element1] / disc_abundances[element2]
        except ZeroDivisionError:
            o_ca = np.nan
        o_ca_values.append(o_ca)
        log_elel_values.append(np.log10(o_ca))
        normalised_elel_values.append(np.log10(o_ca) - solar_ratio)
        T_vals.append(dm.T_disc(d_formation, t_formation))
    d_formation_title = "d_formation /AU "
    T_title = "          T /K           "
    OCa_title = "       " + str(element1) + "/" + str(element2) + "            "
    log_elel_title = "   log(" + str(element1) + "/" + str(element2) + ")           "
    normalised_elel_title = (
        "       [" + str(element1) + "/" + str(element2) + "]           "
    )
    print(
        d_formation_title
        + " | "
        + T_title
        + " | "
        + OCa_title
        + " | "
        + log_elel_title
        + " | "
        + normalised_elel_title
    )
    for d_index, d in enumerate(d_formation_vals):
        print(
            str(d).ljust(len(d_formation_title))
            + " | "
            + str(T_vals[d_index]).ljust(len(T_title))
            + " | "
            + str(o_ca_values[d_index]).ljust(len(OCa_title))
            + " | "
            + str(log_elel_values[d_index]).ljust(len(log_elel_title))
            + " | "
            + str(normalised_elel_values[d_index]).ljust(len(normalised_elel_title))
        )


def plot_element_dropout(
    stellar_metallicity_index,
    z_formation,
    t_formation=1.5,
    elements_to_plot=ci.usual_elements,
):
    element_abundance_vals = {element: list() for element in elements_to_plot}
    d_formation_vals = np.linspace(0.01, 0.3, 101)
    T_vals = list()
    for d_formation in d_formation_vals:
        abundances = get_abundances(
            d_formation, z_formation, t_formation, stellar_metallicity_index
        )
        disc_abundances = dict()
        for el_index, element in enumerate(ci.usual_elements):
            if el_index == 6:
                # Mg is special
                disc_abundances[element] = abundances[element]
            else:
                disc_abundances[element] = (
                    abundances[element]
                    * ld._live_stellar_compositions[stellar_metallicity_index][
                        el_index - 1 if el_index > 6 else el_index
                    ]
                )
        for element in elements_to_plot:
            element_abundance_vals[element].append(disc_abundances[element])
        T_vals.append(dm.T_disc(d_formation, t_formation))
    graph_fac = gf.GraphFactory()
    graph_fac.plot_element_dropout(
        d_formation_vals, T_vals, element_abundance_vals, z_formation
    )
    graph_fac.plot_element_dropout_ratio(
        d_formation_vals,
        T_vals,
        element_abundance_vals,
        ci.Element.Ca,
        ci.Element.Mg,
        z_formation,
    )
    graph_fac.plot_element_dropout_ratio(
        d_formation_vals,
        T_vals,
        element_abundance_vals,
        ci.Element.Ti,
        ci.Element.Mg,
        z_formation,
    )
    graph_fac.plot_element_dropout_ratio(
        d_formation_vals,
        T_vals,
        element_abundance_vals,
        ci.Element.Ca,
        ci.Element.Ti,
        z_formation,
    )


def main():
    # plot_element_dropout(0, 0.088, 1.5, [ci.Element.Al, ci.Element.Ti, ci.Element.Ca, ci.Element.Fe, ci.Element.Mg, ci.Element.Si])
    z_formation = 0  # AU
    t_formation = 1.5  # Myr
    stellar_metallicity_index = 0
    # stellar_metallicity_index = 478
    output_el_el_values(
        ci.Element.N, ci.Element.Si, z_formation, t_formation, stellar_metallicity_index
    )


if __name__ == "__main__":
    main()
