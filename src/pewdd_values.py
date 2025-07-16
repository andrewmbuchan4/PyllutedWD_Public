#!/usr/bin/env python
# -*- coding: utf-8 -*-

from argparse import Namespace

import csv
import numpy as np
import scipy

import chemistry_info as ci
import manager as mn
import solar_abundances as sa
import timescale_interpolator as ti

def get_common_thermohaline_el_values(list_of_elements):
    # Returns lists of element abundances for all wds that have ALL the specified elements, such that the length of each list is the same
    manager = load_manager('PEWDD_thermohaline.csv')
    toret = list()
    for element in list_of_elements:
        toret.append(list())
    i = 0
    for wd in manager.white_dwarfs:
        elements_present = wd.get_elements_present()
        include = True
        for element in list_of_elements:
            if element not in elements_present:
                include = False
        if include:
            for j, element in enumerate(list_of_elements):
                toret[j].append(wd.get_abundance(element).value)
        i += 1
    return toret

def filter_wds_by_num_modellable_detected_elements(white_dwarfs, num_elements, detectable_elements=ci.usual_elements):
    toret = list()
    for wd in white_dwarfs:
        m_d_els = wd.get_modellable_detected_elements(detectable_elements)
        if len(m_d_els) >= num_elements:
            toret.append(wd)
    return toret

def filter_wds_by_atmospheric_type(white_dwarfs, atmospheric_type):
    toret = list()
    for wd in white_dwarfs:
        atm_type = wd.get_atmospheric_type()
        if atm_type.value == atmospheric_type:
            toret.append(wd)
    return toret

def filter_wds_by_teff_logg(white_dwarfs, min_teff, max_teff, min_logg, max_logg):
    toret = list()
    for wd in white_dwarfs:
        teff = wd.get_teff()
        logg = wd.get_logg()
        if teff is not None and logg is not None and teff.value <= max_teff and teff.value >= min_teff and logg.value <= max_logg and logg.value >= min_logg:
            toret.append(wd)
    return toret

def load_manager(filename_override=None):
    manager = mn.Manager()
    if filename_override is not None:
        manager.wd_data_filename = filename_override
        manager.stellar_compositions_filename = 'StellarCompositionsSortFE.csv'
        manager.load_global_data()
    return manager

def load_pewdd_manager():
    manager = load_manager('PEWDD.csv')
    return manager

def load_pewdd_thermohaline_manager():
    manager = load_manager('PEWDD_thermohaline.csv')
    return manager

def main():
    print('Placeholder')

if __name__ == '__main__':
    main()
