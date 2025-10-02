#!/usr/bin/env python

from argparse import Namespace
import numpy as np

import atmosphere_model as atm
import chemistry_info as ci
import graph_factory as gf
import model_parameters as mp
import pwd_utils as pu
import solar_abundances as sa
import timescale_interpolator as ti
import white_dwarf as wd

def produce_pngs(png_dir, wd_type, Teff, logg, CaHx, t_kyr_min, t_kyr_max, timestep, t_event, timescale_type, M_cvz, elements_to_plot=ci.usual_elements):
    graph_fac = gf.GraphFactory(pu.get_path_to_default_graphs() + png_dir)
    wd_property_data_raw1 = {
        mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.He),
        mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, Teff, 0),
        mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, logg, 0)
    }
    property_data1 = wd.WhiteDwarfPropertyData(wd_property_data_raw1)

    wd_abundance_data_raw = dict()
    #wd_abundance_data_raw = {
    #    ci.Element.Ca: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, CaHx, 0.2)
    #}

    abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)

    tim = ti.TimescaleInterpolator()
    timescales = tim.get_wd_timescales(wd_type, logg, Teff, CaHx)

    white_dwarf = wd.WhiteDwarf('Example', property_data1, abundance_data, timescales)

    mu_X = np.array([ci.get_element_mass(el) for el in elements_to_plot])
    M_X = np.array([(10**sa.solar_ratiod_to_H[el]) * (ci.get_element_mass(el)) for el in elements_to_plot])
    tau_X = white_dwarf.get_timescales_as_array(timescale_type, elements_to_plot)
    #tau_X = np.array([5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000, 5000000]) # for testing. should remove differential sinking

    reference_element = None
    extra_text_dict = None

    t_kyr = t_kyr_min
    max_frames = ((t_kyr_max - t_kyr_min)/timestep) + 1
    digits_needed = int(np.ceil(np.log10(max_frames)))
    frame = 0
    while t_kyr <= t_kyr_max:
        print('Time: ' + str(t_kyr) + ' kyr')
        print('Frame ' + str(frame + 1) + '/' + str(int(max_frames)))
        t = t_kyr*1000
        N_X = atm.calculate_abundance_by_number(t, t_event, wd_type, M_cvz, mu_X, M_X, tau_X)
        t_Myr = t_kyr/1000
        fit_dict = {'t = ' + '{:.2f}'.format(t_Myr) + ' Myr': {el: N_X[i] for i, el in enumerate(elements_to_plot)}}
        model_name = str(frame).rjust(digits_needed, '0') + '_'
        graph_fac.make_composition_plot_mk3(white_dwarf, elements_to_plot, fit_dict, reference_element, model_name, extra_text_dict, True, False)
        t_kyr += timestep
        frame += 1
    print('Produced ' + str(frame) + ' frames')

def main():
    t_kyr_min = 0 #kyr
    t_kyr_max = 100000 #kyr
    timestep = 50 #kyr
    t_event = 25000000 #yr
    wd_type = ci.Element.He
    Teff = 10000
    logg = 8
    timescale_type = ti.TimescaleType.KoesterOvershoot
    CaHx = -15
    M_cvz = 10000 #0.01*np.linalg.norm(M_X) # An arbitrary normalisation factor for present purposes
    png_dir = 'ds_example_video_ko/'
    produce_pngs(png_dir, wd_type, Teff, logg, CaHx, t_kyr_min, t_kyr_max, timestep, t_event, timescale_type, M_cvz)
    print('Now run this in the output directory')
    print('ffmpeg -framerate 25 -i  Example_%04d_composition_rel_He.png  -c:v libx264 -r 30 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" out.mp4')

if __name__ == '__main__':
    main()
