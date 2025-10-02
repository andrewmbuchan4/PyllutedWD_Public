#!/usr/bin/env python

from argparse import Namespace
import numpy as np

import chemistry_info as ci
import complete_model as cm
import graph_factory as gf
import manager as mn
import model_parameters as mp
import pwd_utils as pu
import solar_abundances as sa
import timescale_interpolator as ti
import white_dwarf as wd

video_configurations = {
    'volatile_depletion': {
        'png_dir': 'volatile_depletion_video/',
        'interpolation_steps': [10, 200],
        'leg_descriptors': [
            ('log(distance)', 2),
            ('log(distance)', 2)
        ],
        'list_of_coordinates': [
            #[fe_star, t_sinceaccretion (Myr), d_formation, z_formation, N_c, N_o, f_c, f_o, log_fragment_mass, t_disc (yr), pressure, fO2]
            #[478, 5, 2, 0.05, None, None, 0.17, None, 17, 10000000, 45, -2], # Nominal
            [478, 5, 2, 0.1, None, None, 0.17, None, 17, 10000000, 45, -2],
            [478, 5, 0, 0.1, None, None, 0.17, None, 17, 10000000, 45, -2],
            [478, 5, -2, 0.1, None, None, 0.17, None, 17, 10000000, 45, -2]
        ]
    },
    'differentiation': {
        'png_dir': 'differentiation_video/',
        'interpolation_steps': [200],
        'leg_descriptors': [
            ('core fraction', 6)
        ],
        'list_of_coordinates': [
            #[fe_star, t_sinceaccretion (Myr), d_formation, z_formation, N_c, N_o, f_c, f_o, log_fragment_mass, t_disc (yr), pressure, fO2]
            [478, 5, 0, 0.1, None, None, 0, None, 17, 10000000, 54, -1.3],
            [478, 5, 0, 0.1, None, None, 1, None, 17, 10000000, 54, -1.3]
        ]
    },
    'GD61': {
        'png_dir': 'GD61_param_variation_video/',
        'interpolation_steps': [20, 20, 20, 20, 20, 20],
        'leg_descriptors': [
            ('metallicity', 0),
            ('log(distance)', 2),
            ('core fraction', 6),
            ('pressure', 10),
            ('fO2', 11),
            ('time', 1)
        ],
        'list_of_coordinates': [
            #[fe_star, t_sinceaccretion (Myr), d_formation, z_formation, N_c, N_o, f_c, f_o, log_fragment_mass, t_disc (yr), pressure, fO2]
            [478, 0.001, -2, 0.11, None, None, 0.17, None, 18.62112, 10**3.17, 1, -2],
            [382.06, 0.001, -2, 0.11, None, None, 0.17, None, 18.62112, 10**3.17, 1, -2],
            [382.06, 0.001, 0.38, 0.11, None, None, 0.17, None, 18.62112, 10**3.17, 1, -2],
            [382.06, 0.001, 0.38, 0.11, None, None, 0.03, None, 18.62112, 10**3.17, 1, -2],
            [382.06, 0.001, 0.38, 0.11, None, None, 0.03, None, 18.62112, 10**3.17, 44.93, -2],
            [382.06, 0.001, 0.38, 0.11, None, None, 0.03, None, 18.62112, 10**3.17, 44.93, -2.65],
            [382.06, 0.2, 0.38, 0.11, None, None, 0.03, None, 18.62112, 10**3.17, 44.93, -2.65]
        ]
    },
    'GD362': {
        'png_dir': 'GD362_param_variation_video/',
        'interpolation_steps': [3, 5, 5, 3],
        'wd_type': ci.Element.He,
        'Teff': 10057,
        'logg': 7.95,
        'wd_abundances': {
            ci.Element.Al: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.4, 0.2),
            ci.Element.Ti: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.95, 0.1),
            ci.Element.Ca: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.24, 0.1),
            ci.Element.Ni: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.07, 0.15),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.65, 0.1),
            ci.Element.Cr: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.41, 0.1),
            ci.Element.Mg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.98, 0.25),
            ci.Element.Si: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.84, 0.3),
            ci.Element.Na: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.79, 0.2),
            ci.Element.O: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -5.14, None, False),
            ci.Element.N: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -4.14, None, False)
        },
        'leg_descriptors': [
            ('log(distance)', 2),
            ('log(distance)', 2),
            ('log(distance)', 2),
            ('log(distance)', 2)
        ],
        'list_of_coordinates': [
            #[fe_star, t_sinceaccretion (Myr), d_formation, z_formation, N_c, N_o, f_c, f_o, log_fragment_mass, t_disc (yr), pressure, fO2]
#           Stellar metallicity indices	Time since Accretion/Myrs	log(Formation Distance/AU)	Feeding Zone Size/AU	log(Pollution Fraction)	log(Accretion Event Timescale/Yrs)
#Median:
#568.697671995007	0.226503762288808	-0.811497911365654	0.120301945952567	-5.18123294496173	4.46747153477622
#Mode:
#399.52	0.1	-0.7	0.14	-5.17	6.59

            [568.697671995007, 0.226503762288808, 2, 0.120301945952567, None, None, None, None, 21, 10**4.46747153477622, 54, -2],
            [568.697671995007, 0.226503762288808, 0, 0.120301945952567, None, None, None, None, 21, 10**4.46747153477622, 54, -2],
            [568.697671995007, 0.226503762288808, -0.811497911365654, 0.120301945952567, None, None, None, None, 17, 10**4.46747153477622, 54, -2],
            [568.697671995007, 0.226503762288808, 0, 0.120301945952567, None, None, None, None, 21, 10**4.46747153477622, 54, -2],
            [568.697671995007, 0.226503762288808, 2, 0.120301945952567, None, None, None, None, 21, 10**4.46747153477622, 54, -2]
        ]
    },
    'WD0059+257': {
        'png_dir': 'WD0059+257_param_variation_video/',
        'interpolation_steps': [70, 30, 70, 70, 30, 70],
        'wd_type': ci.Element.H,
        'Teff': 20113,
        'logg': 7.89,
        'wd_abundances': {
            ci.Element.Al: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.97),
            ci.Element.Ca: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.71, 0.19),
            ci.Element.Ni: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.8, 0.24),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.54, 0.16),
            ci.Element.Cr: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.59),
            ci.Element.Mg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.84, 0.14),
            ci.Element.Si: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.26, 0.24),
            ci.Element.Na: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.79, 0.2),
            ci.Element.O: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.72, 0.14)
        },
        'leg_descriptors': [
            ('core fraction', 6),
            ('core fraction', 6),
            ('core fraction', 6),
            ('core fraction', 6),
            ('core fraction', 6),
            ('core fraction', 6)
        ],
        'list_of_coordinates': [

#            #[fe_star, t_sinceaccretion (Myr), d_formation, z_formation, N_c, N_o, f_c, f_o, log_fragment_mass, t_disc (yr), pressure, fO2]
#Stellar metallicity indices	Time since Accretion/Myrs	log(Formation Distance/AU)	Fragment Core Fraction	log(Fragment Mass /kg)	log(Accretion Event Timescale/Yrs)	Pressure /GPa	Oxygen Fugacity /ΔIW

#Median:
#459.152406077556	9.51060362842148E-07	-0.225163316389305	0.452100129532203	3295380772673.62	0.305741302429945	16.2438353981353	-1.75587458913119

            [459.152406077556, 9.51060362842148E-07, -0.225163316389305, 0.05, None, None, 0, None, 20.5, 10**0.305741302429945, 16.2438353981353, -1.75587458913119],
            [459.152406077556, 9.51060362842148E-07, -0.225163316389305, 0.05, None, None, 0.452100129532203, None, 20.5, 10**0.305741302429945, 16.2438353981353, -1.75587458913119],
            [459.152406077556, 9.51060362842148E-07, -0.225163316389305, 0.05, None, None, 0.452100129532203, None, 20.5, 10**0.305741302429945, 16.2438353981353, -1.75587458913119],
            [459.152406077556, 9.51060362842148E-07, -0.225163316389305, 0.05, None, None, 0.8, None, 20.5, 10**0.305741302429945, 16.2438353981353, -1.75587458913119],
            [459.152406077556, 9.51060362842148E-07, -0.225163316389305, 0.05, None, None, 0.452100129532203, None, 20.5, 10**0.305741302429945, 16.2438353981353, -1.75587458913119],
            [459.152406077556, 9.51060362842148E-07, -0.225163316389305, 0.05, None, None, 0.452100129532203, None, 20.5, 10**0.305741302429945, 16.2438353981353, -1.75587458913119],
            [459.152406077556, 9.51060362842148E-07, -0.225163316389305, 0.05, None, None, 0, None, 20.5, 10**0.305741302429945, 16.2438353981353, -1.75587458913119]
        ]
    },
    'J1227GR': {
        'png_dir': 'J1227GR_heating_variation/',
        'interpolation_steps': [50, 50],
        'wd_type': ci.Element.He,
        'Teff': 7946.72,
        'logg': 8.069893,
        'wd_abundances': {
            ci.Element.Al: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -7),
            ci.Element.Ti: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -9.7, 0.3),
            ci.Element.Ca: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -8.7, 0.2),
            ci.Element.Ni: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -8.9, 0.2),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.5, 0.2),
            ci.Element.Cr: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -9),
            ci.Element.Mg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.3, 0.2),
            ci.Element.Si: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -7.1),
            ci.Element.Na: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -8.3, 0.2),
            ci.Element.O: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -4)
        },
        'leg_descriptors': [
            ('formdist', 2),
            ('formdist', 2)
        ],
        'list_of_coordinates': [

#            #[fe_star, t_sinceaccretion (Myr), d_formation, z_formation, N_c, N_o, f_c, f_o, log_fragment_mass, t_disc (yr), pressure, fO2]
#Stellar metallicity indices	Time since Accretion/Myrs	log(Formation Distance/AU)	Fragment Core Fraction	log(Fragment Mass /kg)	log(Accretion Event Timescale/Yrs)	Pressure /GPa	Oxygen Fugacity /ΔIW

#Median:

            [695.511722949721, 0.935056105420901, -1, 0.080417918508727, None, None, 0.17, None, 18.9744928319028, 3.02727336546314, 54, -2],
            [695.511722949721, 0.935056105420901, 0.233852560787135, 0.080417918508727, None, None, 0.17, None, 18.9744928319028, 3.02727336546314, 54, -2],
            [695.511722949721, 0.935056105420901, 2, 0.080417918508727, None, None, 0.17, None, 18.9744928319028, 3.02727336546314, 54, -2]

        ]
    }
}

def interpolate_between_coords(coord1, coord2, interpolation_fraction):
    assert len(coord1) == len(coord2)
    toret = list()
    for i, entry1 in enumerate(coord1):
        entry2 = coord2[i]
        try:
            interpolated_entry = (entry1*(1-interpolation_fraction)) + (entry2*interpolation_fraction)
            print(interpolated_entry)
            toret.append(interpolated_entry)
        except TypeError:
            toret.append(None)
    return toret

def produce_pngs(wd_name, png_dir, list_of_coordinates, leg_descriptors, interpolation_steps, wd_type, Teff, logg, CaHx, timescale_type, wd_abundance_data_raw=dict(), elements_to_plot=ci.usual_elements):
    graph_fac = gf.GraphFactory(pu.get_path_to_default_graphs() + png_dir)
    wd_property_data_raw1 = {
        mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, wd_type),
        mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, Teff, 0),
        mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, logg, 0)
    }
    property_data1 = wd.WhiteDwarfPropertyData(wd_property_data_raw1)

    #wd_abundance_data_raw = dict()
    #wd_abundance_data_raw = {
    #    ci.Element.Ca: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, CaHx, 0.2)
    #}

    abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)

    tim = ti.TimescaleInterpolator()
    timescales = tim.get_wd_timescales(wd_type, logg, Teff, CaHx)

    white_dwarf = wd.WhiteDwarf(wd_name, property_data1, abundance_data, timescales)

    mu_X = np.array([ci.get_element_mass(el) for el in elements_to_plot])
    M_X = np.array([(10**sa.solar_ratiod_to_H[el]) * (ci.get_element_mass(el)) for el in elements_to_plot])
    tau_X = white_dwarf.get_timescales_as_array(timescale_type, elements_to_plot)
    #tau_X = np.array([5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000]) # for testing. should remove differential sinking

    reference_element = ci.Element.Mg
    extra_text_dict = None

    max_frames = sum(interpolation_steps)
    print(max_frames)
    digits_needed = int(np.ceil(np.log10(max_frames)))

    frame = 0
    leg = 0
    for coord_index, coord in enumerate(list_of_coordinates):
        try:
            next_coord = list_of_coordinates[coord_index + 1]
        except IndexError:
            #Done!
            break
        if frame == 0:
            step = 0
        else:
            step = 1 # Otherwise we duplicate the previous frame
        interpolation_steps_to_use = interpolation_steps[leg] - 1 if leg == 0 else interpolation_steps[leg]
        while step <= interpolation_steps_to_use:
            frame += 1
            print('Outputting frame ' + str(frame) + '/' + str(max_frames))
            interpolation_fraction = step/interpolation_steps_to_use
            coords_to_use = interpolate_between_coords(coord, next_coord, interpolation_fraction)
            N_X, ignore = cm.complete_model_calculation(
                coords_to_use[0],
                coords_to_use[1],
                coords_to_use[2],
                coords_to_use[3],
                coords_to_use[4],
                coords_to_use[5],
                coords_to_use[6],
                coords_to_use[7],
                coords_to_use[8],
                coords_to_use[9],
                coords_to_use[10],
                coords_to_use[11],
                'NonEarthlike',
                False
            )
            fit_desc = leg_descriptors[leg]
            fit_dict = {fit_desc[0] + ' = {:.3f}'.format(coords_to_use[fit_desc[1]]): {el: N_X[el] for el in elements_to_plot}}
            model_name = str(frame).rjust(digits_needed, '0') + '_'
            graph_fac.make_composition_plot_mk3(white_dwarf, elements_to_plot, fit_dict, reference_element, model_name, extra_text_dict, True, False)
            step += 1
        leg += 1
    print('Produced ' + str(frame) + ' frames')

def main():
    configuration = 'J1227GR'
    png_dir = video_configurations[configuration]['png_dir']
    interpolation_steps = video_configurations[configuration]['interpolation_steps']
    leg_descriptors = video_configurations[configuration]['leg_descriptors']
    list_of_coordinates = video_configurations[configuration]['list_of_coordinates']
    wd_type = video_configurations[configuration]['wd_type']
    Teff = video_configurations[configuration]['Teff']
    logg = video_configurations[configuration]['logg']
    wd_abundance_data_raw = video_configurations[configuration]['wd_abundances']
    #Teff = 17280 # For GD61 (move this to config dict!)
    #logg = 8.2 # For GD61
    CaHx = -8.7 # Use -15 by default here!
    timescale_type = ti.TimescaleType.BedardNoOvershoot
    manager = mn.Manager()
    manager.publish_live_data(0, timescale_type)
    elements_to_plot = [
        ci.Element.Al,
        ci.Element.Ti,
        ci.Element.Ca,
        ci.Element.Ni,
        ci.Element.Fe,
        ci.Element.Cr,
        ci.Element.Mg,
        ci.Element.Si,
        ci.Element.Na,
        ci.Element.O
    ]
    produce_pngs(configuration, png_dir, list_of_coordinates, leg_descriptors, interpolation_steps, wd_type, Teff, logg, CaHx, timescale_type, wd_abundance_data_raw, elements_to_plot)
    print('Now run this in the output directory (' + pu.get_path_to_default_graphs() + png_dir + ')')
    print('(Might need minor modifications to match the filename format)')
    print('ffmpeg -framerate 25 -i  Example_%03d_composition_rel_He.png  -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" ' + configuration + '_video.mp4')

if __name__ == '__main__':
    main()
