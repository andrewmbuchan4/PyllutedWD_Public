#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import chemistry_info as ci
import geology_info as gi
import graph_factory as gf

def main():
    geo_model = gi.GeologyModel()
    pressure_vals = [0.001, 0.01, 0.1, 0.5] + list(range(1, 61, 1))
    #pressure_vals = [20, 35, 50, 65, 80, 95]
    series_to_plot = dict()
    series_to_plot['Earth-like Fe'] = {'Pressure': list(), 'Radius': list(), 'Mass': list()}
    #series_to_plot['Conservative'] = {'Pressure': list(), 'Radius': list(), 'Mass': list()}
    #series_to_plot['High Fe'] = {'Pressure': list(), 'Radius': list(), 'Mass': list()}
    #print('Pressure /GPa,Radius /km,Mass /M_Earth')
    highFe = {ci.Element.Fe: {gi.Layer.bulk: 0.4005}, ci.Element.O: {gi.Layer.bulk: 0.5995}}
    for p in pressure_vals:
        r1, m1 = geo_model.calculate_radius_and_mass(p, geo_model.element_info)
        series_to_plot['Earth-like Fe']['Pressure'].append(p)
        series_to_plot['Earth-like Fe']['Radius'].append(r1)
        series_to_plot['Earth-like Fe']['Mass'].append(m1)
        #r2, m2 = geo_model.calculate_radius_and_mass(p, geo_model.element_info, True)
        #series_to_plot['Conservative']['Pressure'].append(p)
        #series_to_plot['Conservative']['Radius'].append(r2)
        #series_to_plot['Conservative']['Mass'].append(m2)
        #r3, m3 = geo_model.calculate_radius_and_mass(p, highFe)
        #series_to_plot['High Fe']['Pressure'].append(p)
        #series_to_plot['High Fe']['Radius'].append(r3)
        #series_to_plot['High Fe']['Mass'].append(m3)
    #    print(str(p) + ',' + str(r) + ',' + str(m))
    graph_fac = gf.GraphFactory()
    #graph_fac.plot_size_v_pressure(series_to_plot, 'Radius')
    #graph_fac.plot_size_v_pressure(series_to_plot, 'Mass')
    graph_fac.plot_size_v_pressure(series_to_plot)

if __name__ == '__main__':
    main()
