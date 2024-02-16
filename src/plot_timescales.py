#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import chemistry_info as ci
import graph_factory as gf
import timescale_interpolator as ti

def get_logg_values():
    return [7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9]

def get_Teff_values():
    tlist = [
        2000, 2250, 2500, 2750,
        3000, 3250, 3500, 3750,
        4000, 4250, 4500, 4750,
        5000, 5250, 5500, 5750,
        6000, 6250, 6500, 6750,
        7000, 7250, 7500, 7750,
        8000, 8250, 8500, 8750,
        9000, 9250, 9500, 9750,
        10000, 10250, 10500, 10750,
        11000, 11250, 11500, 11750,
        12000, 12250, 12500, 12750,
        13000, 13250, 13500, 13750,
        14000, 14250, 14500, 14750,
        15000, 15250, 15500, 15750,
        16000,
        16250, 16500, 16750,
        17000, 17250, 17500, 17750,
        18000, 18250, 18500, 18750,
        19000, 19250, 19500, 19750,
        20000, 20250, 20500, 20750,
        21000, 21250, 21500, 21750,
        22000, 22250, 20500

    ]
    return tlist

def get_CaHe_values():
    return [-6.5, -7, -7.5, -8, -8.5, -9, -9.5, -10, -10.5, -11, -11.5, -12, -12.5, -13, -13.5, -14, -14.5, -15, -15.5, -16]

def get_variable_vals(variable):
    if variable == 'logg':
        return get_logg_values()
    if variable == 'Teff':
        return get_Teff_values()
    if variable == 'CaHe':
        return get_CaHe_values()
    return None

def get_elements_for_2D_plot():
    #return [ci.Element.Ca, ci.Element.Fe, ci.Element.Mg, ci.Element.O]
    #return [ci.Element.C, ci.Element.Si]
    return [ci.Element.Ca, ci.Element.Fe]

def get_elements_for_model_comparison(): # needs to be a subset of get_elements_for_2D_plot().
    return get_elements_for_2D_plot()
    #return [ci.Element.Ca, ci.Element.Mg]
    #return [ci.Element.C, ci.Element.Si]
    #return [ci.Element.Fe, ci.Element.Ca]

def generate_timescales():
    timescale_interpolator = ti.TimescaleInterpolator()
    logg_vals = get_logg_values()
    Teff_vals = get_Teff_values()
    CaHe_vals = get_CaHe_values()
    timescale_vals = {
        'H': dict(),
        'He': dict()
    }
    elements_we_care_about = get_elements_for_2D_plot()
    for element in elements_we_care_about:
        timescale_vals['H'][element] = np.zeros((len(logg_vals), len(Teff_vals)))
        timescale_vals['He'][element] = np.zeros((len(logg_vals), len(Teff_vals), len(CaHe_vals)))
    for i in range(len(logg_vals)):
        logg = logg_vals[i]
        for j in range(len(Teff_vals)):
            Teff = Teff_vals[j]
            H_timescales = timescale_interpolator.get_wd_timescales('H', logg, Teff, None)
            for element in elements_we_care_about:
                timescale_vals['H'][element][i,j] = H_timescales[element]
            for k in range(len(CaHe_vals)):
                CaHe = CaHe_vals[k]
                He_timescales = timescale_interpolator.get_wd_timescales('He', logg, Teff, CaHe)
                for element in elements_we_care_about:
                    timescale_vals['He'][element][i,j,k] = He_timescales[element]
    return timescale_vals

def extract_2D_vals_to_plot(timescale_vals, HorHe, logg, Teff, CaHe, elements_to_extract):
    # Whichever of these variables is None is the one we need to plot against
    num_variables = 0
    variable = None
    if logg is None:
        num_variables += 1
        variable = 'logg'
    if Teff is None:
        num_variables += 1
        variable = 'Teff'
    if CaHe is None:
        num_variables += 1
        variable = 'CaHe'
    if num_variables != 1:
        print('Error: Exactly 1 of logg, Teff and CaHe must be None. (This indicates the variable to plot against)')
    assert num_variables == 1
    if HorHe == 'H' and variable == 'CaHe':
        print('Error: H does not vary with CaHe')
        return None
    assert logg is None or logg in get_logg_values()
    assert Teff is None or Teff in get_Teff_values()
    assert CaHe is None or CaHe in get_CaHe_values()
    try:
        logg_index = get_logg_values().index(logg)
    except ValueError:
        logg_index = None
    try:
        Teff_index = get_Teff_values().index(Teff)
    except ValueError:
        Teff_index = None
    try:
        CaHe_index = get_CaHe_values().index(CaHe)
    except ValueError:
        CaHe_index = None
    to_plot = dict()
    for element in elements_to_extract:
        to_plot[element] = list()
        to_slice_into = timescale_vals[HorHe][element]
        if HorHe == 'H':
            if variable == 'logg':
                to_plot[element] = to_slice_into[:,Teff_index]
            if variable == 'Teff':
                to_plot[element] = to_slice_into[logg_index,:]
        elif HorHe == 'He':
            if variable == 'logg':
                to_plot[element] = to_slice_into[:,Teff_index,CaHe_index]
            if variable == 'Teff':
                to_plot[element] = to_slice_into[logg_index,:,CaHe_index]
            if variable == 'CaHe':
                to_plot[element] = to_slice_into[logg_index,Teff_index,:]
        else:
            print('Error: HorHe must be H or He')
            return None
    return to_plot, variable

def plot_2D_timescale(timescale_vals, HorHe, logg, Teff, CaHe):
    to_plot, variable = extract_2D_vals_to_plot(timescale_vals, HorHe, logg, Teff, CaHe, get_elements_for_2D_plot())
    graph_factory = gf.GraphFactory()
    graph_factory.plot_2D_timescales(to_plot, variable, get_variable_vals(variable), HorHe, logg, Teff, CaHe)

def plot_model_comparison(timescale_vals, HorHe, logg, Teff, CaHe):
    to_plot_new, variable = extract_2D_vals_to_plot(timescale_vals, HorHe, logg, Teff, CaHe, get_elements_for_model_comparison())
    if variable == 'logg':
        # These are at (or near) Teff = 6250 according to Simon (but we take the timescales themselves from Hollands) (also excluding those where logg was just assumed to be 8 by default)
        # System      logg     Ca       Fe
        #J0108-0537 7.256 1230268.771 1230268.771
        #J0201+2015 8.091 1348962.883 1348962.883
        #J0851+1543 8.12 1202264.435 1230268.771
        #J0939+4136 8.422 1122018.454 1148153.621
        #J0948+3008 7.21 1513561.248 1380384.265
        #J1024+4531 8.356 1258925.412 1258925.412
        #J1230+3143 7.676 1584893.192 1548816.619
        #J1314+3748 8.379 1071519.305 1071519.305
        #J1336+3547 8.057 1318256.739 1348962.883
        #J1356+2416 8.01 1584893.192 1412537.545
        #J1404+3620 8.216 1584893.192 1348962.883
        #J1430-0151 8.076 831763.7711 870963.59
        #J1448+1047 8.173 1445439.771 1445439.771
        #J1524+4049 8.203 1230268.771 1230268.771
        #J2340+0817 8.141 1288249.552 1174897.555
        old_vals_dict = {
            7.256:( 1230268.771, 1230268.771),
            8.091:( 1348962.883, 1348962.883),
            8.12: (1202264.435, 1230268.771) ,
            8.422:( 1122018.454, 1148153.621),
            7.21: (1513561.248, 1380384.265) ,
            8.356: (1258925.412, 1258925.412),
            7.676: (1584893.192, 1548816.619),
            8.379: (1071519.305, 1071519.305),
            8.057: (1318256.739, 1348962.883),
            8.01: (1584893.192 ,1412537.545 ),
            8.216: (1584893.192, 1348962.883),
            8.076: (831763.7711, 870963.59  ),
            8.173: (1445439.771, 1445439.771),
            8.203: (1230268.771, 1230268.771),
            8.141: (1288249.552, 1174897.555)
        }
        x_vals_old = list(old_vals_dict.keys())
        x_vals_old.sort()
        ca_vals = list()
        fe_vals = list()
        for x_val in x_vals_old:
            ca_vals.append(old_vals_dict[x_val][0])
            fe_vals.append(old_vals_dict[x_val][1])
        to_plot_old = {
            ci.Element.Ca: np.array(ca_vals),
            ci.Element.Fe: np.array(fe_vals)
        }
        print(x_vals_old)
        print(to_plot_old)
    elif variable == 'Teff':
        # These are at (or near) logg = 8 according to Simon (but we take the timescales themselves from Hollands) (also excluding those where logg was just assumed to be 8 by default)
        # System   Teff     Ca       Fe
        #J0158-0942 6115 2089296.131 1513561.248
        #J0252+0054 7478 1479108.388 1513561.248
        #J0807+4930 5172 851138.0382 851138.0382
        #J0842+1406 7075 1380384.265 1445439.771
        #J0901+0752 7263 954992.586 1000000
        #J0937+5228 6660 1288249.552 1318256.739
        #J1017+3447 6089 1819700.859 1513561.248
        #J1144+1218 5320 1737800.829 1288249.552
        #J1158+0454 5344 1023292.992 977237.221
        #J1158+1845 6696 1230268.771 1318256.739
        #J1158+5942 6046 1318256.739 1288249.552
        #J1224+2838 4991 2884031.503 1548816.619
        #J1254+3551 6417 1513561.248 1548816.619
        #J1259+3112 5664 2187761.624 1548816.619
        #J1356+2416 6173 1584893.192 1412537.545
        #J1356+0236 7662 1380384.265 1479108.388
        #J1401+3659 5931 2398832.919 1659586.907
        #J1421+1843 6517 1096478.196 1148153.621
        #J1616+3303 6491 1148153.621 1174897.555
        #J2109-0039 6132 1230268.771 1230268.771
        #J2231+0906 5679 2454708.916 1659586.907
        old_vals_dict = {
            6115: ( 2089296.131, 1513561.248),
            7478: ( 1479108.388, 1513561.248),
            5172: ( 851138.0382, 851138.0382),
            7075: ( 1380384.265, 1445439.771),
            7263: ( 954992.586, 1000000     ),
            6660: ( 1288249.552, 1318256.739),
            6089: ( 1819700.859, 1513561.248),
            5320: ( 1737800.829, 1288249.552),
            5344: ( 1023292.992, 977237.221 ),
            6696: ( 1230268.771, 1318256.739),
            6046: ( 1318256.739, 1288249.552),
            4991: ( 2884031.503, 1548816.619),
            6417: ( 1513561.248, 1548816.619),
            5664: ( 2187761.624, 1548816.619),
            6173: ( 1584893.192, 1412537.545),
            7662: ( 1380384.265, 1479108.388),
            5931: ( 2398832.919, 1659586.907),
            6517: ( 1096478.196, 1148153.621),
            6491: ( 1148153.621, 1174897.555),
            6132: ( 1230268.771, 1230268.771),
            5679: ( 2454708.916, 1659586.907)
        }
        x_vals_old = list(old_vals_dict.keys())
        x_vals_old.sort()
        ca_vals = list()
        fe_vals = list()
        for x_val in x_vals_old:
            ca_vals.append(old_vals_dict[x_val][0])
            fe_vals.append(old_vals_dict[x_val][1])
        to_plot_old = {
            ci.Element.Ca: np.array(ca_vals),
            ci.Element.Fe: np.array(fe_vals)
        }
        print(x_vals_old)
        print(to_plot_old)
    else:
        print('Old timescales only varied with Teff and logg')
        x_vals_old = None
        to_plot_old = None
    if HorHe != 'He':
        print('Old timescales only present for He')
        x_vals_old = None
        to_plot_old = None
    graph_factory = gf.GraphFactory()
    graph_factory.plot_timescale_model_comparison(to_plot_new, to_plot_old, variable, get_variable_vals(variable), x_vals_old, HorHe, logg, Teff, CaHe, get_elements_for_model_comparison())

def main():
    timescale_vals = generate_timescales()
    print(timescale_vals)
    #plot_2D_timescale(timescale_vals, 'He', None, 8000, -9.5)
    #plot_2D_timescale(timescale_vals, 'He', 8, None, -9.5)
    #plot_2D_timescale(timescale_vals, 'He', 8, 8000, None)
    #plot_2D_timescale(timescale_vals, 'H', None, 8000, -9.5)
    #plot_2D_timescale(timescale_vals, 'H', 8, None, -9.5)

    #plot_model_comparison(timescale_vals, 'H', None, 6250, -9.5) # At the moment, we're limited to Teff == 6250 only  (the only val. I got old model values for)
    #plot_model_comparison(timescale_vals, 'He', None, 6250, -9.5) # At the moment, we're limited to Teff == 6250 only  (the only val. I got old model values for)
    plot_model_comparison(timescale_vals, 'He', 8, 6500, None)
    plot_model_comparison(timescale_vals, 'He', 8, None, -9.5)
    plot_model_comparison(timescale_vals, 'He', None, 6500, -9.5)
    #plot_model_comparison(timescale_vals, 'H', 8, None, -9.5)  # At the moment, we're limited to logg == 8 only

if __name__ == '__main__':
    main()
