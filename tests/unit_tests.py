#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: Find a less hacky way to import from /src/
# Maybe add an __init__.py to /src/

import csv
import numpy as np
import os
import unittest
import shutil
import sys

from argparse import Namespace
from pathlib import Path

def get_path_to_self():
    # Path(__file__) is the path of this file
    # resolve() returns the absolute path
    return Path(__file__).resolve()

def get_path_to_tests():
    # parents[0] is the path of the current directory
    return str(get_path_to_self().parents[0])
    
def get_path_to_parent():
    # parents[1] is the path of the parent directory
    return str(get_path_to_self().parents[1])

def get_path_to_base_dir():
    return '<your_filepath_here>'

def get_path_to_src():
    return get_path_to_parent() + '/src/'
    
def get_path_to_original_src():
    return get_path_to_parent() + '/original_codebase/'
    
def get_path_to_data():
    return get_path_to_parent() + '/data/'
    
def get_path_to_feni():
    return get_path_to_parent() + '/feni_src/feni/'

def get_path_to_test_output():
    return get_path_to_base_dir() + '/tests'
    
def get_path_to_chains():
    return get_path_to_tests() + '/chains/'

sys.path.append(get_path_to_src())
sys.path.append(get_path_to_original_src())

import abundance_model as am
import chemistry_info as ci
import complete_model as cm
import disc_model as dm
import enhancement_model as em
import geology_info as gi
import live_data as ld
import manager as mn
import model_parameters as mp
import original_complete_model as ocm
import original_enhancement_model as oe
import original_pressure_model as op
import partition_model as pam
import PWDCodeModel24Only as oc
import timescale_interpolator as ti
import white_dwarf_model as wdm

class TimescaleInterpolatorTests(unittest.TestCase):

    def test_init(self):
        test_ti = ti.TimescaleInterpolator()
        self.assertEqual(6.392, test_ti.timescale_data['He'][8][3000][-7][ci.Element.C])
        self.assertEqual(6.32, test_ti.timescale_data['He'][8][3000][-7][ci.Element.N])
        self.assertEqual(-5.533, test_ti.timescale_data['He'][8][3000][-7]['logq'])
        self.assertEqual(6.464, test_ti.timescale_data['He'][8][3000][-7.5][ci.Element.C])
        self.assertEqual(6.389, test_ti.timescale_data['He'][8][3250][-7.5][ci.Element.C])
        self.assertEqual(5.021, test_ti.timescale_data['He'][8.5][3250][-7.5][ci.Element.C])
        with self.assertRaises(KeyError):
            not_real = test_ti.timescale_data['He'][8][3001][-7][ci.Element.C]
        with self.assertRaises(KeyError):
            not_real = test_ti.timescale_data['H'][8][3000][-7][ci.Element.C]  # -7 is meaningless to a H interpolator: CaHe is not a variable
        self.assertEqual(None, test_ti.wd_data)
        self.assertEqual(['H', 'He'], list(test_ti.interpolators.keys()))
        self.assertEqual(2, len(test_ti.interpolators['H']['logq'].grid))  # H is 2D (Teff, logg)
        self.assertEqual(3, len(test_ti.interpolators['He']['logq'].grid)) # He is 3D (Teff, logg, CaHe)
        self.assertEqual([
             'logq',
             ci.Element.Li,
             ci.Element.Be,
             ci.Element.B, 
             ci.Element.C, 
             ci.Element.N, 
             ci.Element.O, 
             ci.Element.F, 
             ci.Element.Ne, 
             ci.Element.Na, 
             ci.Element.Mg, 
             ci.Element.Al, 
             ci.Element.Si, 
             ci.Element.P, 
             ci.Element.S, 
             ci.Element.Cl, 
             ci.Element.Ar, 
             ci.Element.K, 
             ci.Element.Ca, 
             ci.Element.Sc, 
             ci.Element.Ti, 
             ci.Element.V, 
             ci.Element.Cr, 
             ci.Element.Mn,
             ci.Element.Fe, 
             ci.Element.Co, 
             ci.Element.Ni, 
             ci.Element.Cu, 
             ci.Element.Zn
             ], list(test_ti.interpolators['He'].keys()))
    
    def test_interpolation(self):
        test_ti = ti.TimescaleInterpolator()
        wd_data = test_ti.load_wd_data()
        self.assertEqual(460, len(wd_data))
        self.assertEqual({'CaHe': -9.1, 'Teff': 6466, 'Type': 'He', 'logg': 8.257}, wd_data['SDSSJ0002+3209'])
        self.assertEqual({'CaHe': -6.11, 'Teff': 20409, 'Type': 'H', 'logg': 7.95}, wd_data['WD1929+011'])
        self.assertEqual(None, test_ti.extract_timescales('He', 8.1, 3123, 0))
        
        test_extraction = test_ti.extract_timescales('He', 8.1, 3123, -7.2)
        expected_test_extraction = {
            'logq': -5.827288959999999,
            ci.Element.Li: 6.174627040000002, 
            ci.Element.Be: 6.10975712, 
            ci.Element.B: 6.077452320000001, 
            ci.Element.C: 6.101025120000003, 
            ci.Element.N: 6.02990816, 
            ci.Element.O: 5.969210240000002, 
            ci.Element.F: 5.85992624, 
            ci.Element.Ne: 5.858646240000001, 
            ci.Element.Na: 5.783110240000001, 
            ci.Element.Mg: 5.773388960000001, 
            ci.Element.Al: 5.717176160000001, 
            ci.Element.Si: 5.712894240000001, 
            ci.Element.P: 5.658746240000001, 
            ci.Element.S: 5.654281440000002, 
            ci.Element.Cl: 5.596773440000001, 
            ci.Element.Ar: 5.5272793600000005, 
            ci.Element.K: 5.558265440000001, 
            ci.Element.Ca: 5.555288320000001, 
            ci.Element.Sc: 5.487475520000001, 
            ci.Element.Ti: 5.456463360000002, 
            ci.Element.V: 5.4264633600000005, 
            ci.Element.Cr: 5.422941440000002, 
            ci.Element.Mn: 5.396876640000001, 
            ci.Element.Fe: 5.394555360000001, 
            ci.Element.Co: 5.3685515200000005, 
            ci.Element.Ni: 5.37839152, 
            ci.Element.Cu: 5.336287360000001, 
            ci.Element.Zn: 5.325647360000001
        }
        self.assertEqual(expected_test_extraction, test_extraction)
        test_H_extraction = test_ti.extract_timescales('H', 8.1, 3123)
        expected_test_H_extraction = {
            'logq': -4.338032799999992, 
            ci.Element.He: 7.0982055999999965, # This is the only value I've manually checked
            ci.Element.Li: 6.913929599999997, 
            ci.Element.Be: 6.902627999999993, 
            ci.Element.B: 6.925250399999989, 
            ci.Element.C: 6.956859999999987, 
            ci.Element.N: 6.938161599999998, 
            ci.Element.O: 6.922964799999991, 
            ci.Element.F: 6.866155199999999, 
            ci.Element.Ne: 6.868450399999988, 
            ci.Element.Na: 6.826142399999992, 
            ci.Element.Mg: 6.819034399999995, 
            ci.Element.Al: 6.788231199999993, 
            ci.Element.Si: 6.790129600000001, 
            ci.Element.P: 6.753221599999994, 
            ci.Element.S: 6.752418399999996, 
            ci.Element.Cl: 6.706405599999994, 
            ci.Element.Ar: 6.655497599999994, 
            ci.Element.K: 6.67939599999999, 
            ci.Element.Ca: 6.679494399999996, 
            ci.Element.Sc: 6.627786399999998, 
            ci.Element.Ti: 6.608586399999994, 
            ci.Element.V: 6.5901863999999994, 
            ci.Element.Cr: 6.582478399999995, 
            ci.Element.Mn: 6.564576799999997, 
            ci.Element.Fe: 6.559570400000003, 
            ci.Element.Co: 6.542370399999989, 
            ci.Element.Ni: 6.557475199999995, 
            ci.Element.Cu: 6.521268799999995, 
            ci.Element.Zn: 6.51356719999999
        }
        self.assertEqual(expected_test_H_extraction, test_H_extraction)
        test_extrapolation = test_ti.extract_timescales('He', 9.1, 6123, -6.2)
        expected_extrapolation = {
            'logq': -8.82152, 
            ci.Element.Li: 3.2080931200000116, 
            ci.Element.Be: 3.128967200000009, 
            ci.Element.B: 3.077587040000001, 
            ci.Element.C: 3.0766246400000092, 
            ci.Element.N: 3.0176424000000033, 
            ci.Element.O: 2.964045120000005, 
            ci.Element.F: 2.8880315200000055, 
            ci.Element.Ne: 2.8726579200000053, 
            ci.Element.Na: 2.817027840000005, 
            ci.Element.Mg: 2.8029563200000034, 
            ci.Element.Al: 2.7603198400000046, 
            ci.Element.Si: 2.750029120000004, 
            ci.Element.P: 2.711339040000002, 
            ci.Element.S: 2.7043976000000036, 
            ci.Element.Cl: 2.6590310400000057, 
            ci.Element.Ar: 2.609198240000005, 
            ci.Element.K: 2.631071040000009, 
            ci.Element.Ca: 2.6231803199999995, 
            ci.Element.Sc: 2.5678310400000015, 
            ci.Element.Ti: 2.543594560000006, 
            ci.Element.V: 2.5164310400000076, 
            ci.Element.Cr: 2.5182696, 
            ci.Element.Mn: 2.4914296000000107, 
            ci.Element.Fe: 2.490613760000006, 
            ci.Element.Co: 2.466450240000004, 
            ci.Element.Ni: 2.479563040000002, 
            ci.Element.Cu: 2.447452320000007, 
            ci.Element.Zn: 2.4322123200000014
        }
        self.assertEqual(expected_extrapolation, test_extrapolation)
        wd_timescales_dict = test_ti.return_wd_timescales_as_dict(test_extraction)
        expected_wd_timescales_dict = {
            'logq': -5.827288959999999, 
            ci.Element.Al: 521406.1627526086, 
            ci.Element.Ti: 286064.1007553645, 
            ci.Element.Ca: 359160.2953306345, 
            ci.Element.Ni: 238996.48847621836, 
            ci.Element.Fe: 248059.2121373417, 
            ci.Element.Cr: 264814.30405580055, 
            ci.Element.Mg: 593456.5943745527, 
            ci.Element.Si: 516290.6261463279, 
            ci.Element.Na: 606890.3612551007, 
            ci.Element.O: 931558.7297104555, 
            ci.Element.C: 1261900.5219890976, 
            ci.Element.N: 1071292.7356341446
        }
        self.assertEqual(expected_wd_timescales_dict, wd_timescales_dict)
        wd_timescales_list = test_ti.return_wd_timescales_as_list(test_extraction)
        self.assertEqual([-5.827288959999999, 521406.1627526086, 286064.1007553645, 359160.2953306345, 238996.48847621836, 248059.2121373417, 264814.30405580055, 593456.5943745527, 516290.6261463279, 606890.3612551007, 931558.7297104555, 1261900.5219890976, 1071292.7356341446], wd_timescales_list)
        self.assertEqual(wd_timescales_dict, test_ti.get_wd_timescales('He', 8.1, 3123, -7.2))
        processed_wd_data = test_ti.process_wd_data()
        self.assertEqual(len(processed_wd_data), len(wd_data))
        test_real_extraction = test_ti.extract_timescales('He', 8.257, 6466, -9.1)
        expected_real_timescales = {
            'logq': -5.659591270399999, 
            ci.Element.Li: 6.145003929600002, 
            ci.Element.Be: 6.103735929600001, 
            ci.Element.B: 6.091647584000001, 
            ci.Element.C: 6.134778784000001, 
            ci.Element.N: 6.075668384, 
            ci.Element.O: 6.0251084224, 
            ci.Element.F: 5.921529337600001, 
            ci.Element.Ne: 5.9297075456, 
            ci.Element.Na: 5.859139584000001, 
            ci.Element.Mg: 5.857115545600001, 
            ci.Element.Al: 5.8047173376000005, 
            ci.Element.Si: 5.807701145600001, 
            ci.Element.P: 5.7563097920000015, 
            ci.Element.S: 5.756492499200001, 
            ci.Element.Cl: 5.701100384000001, 
            ci.Element.Ar: 5.6317805376, 
            ci.Element.K: 5.669069337600001, 
            ci.Element.Ca: 5.669908384, 
            ci.Element.Sc: 5.602536384000002, 
            ci.Element.Ti: 5.573314745600001, 
            ci.Element.V: 5.544337145600002, 
            ci.Element.Cr: 5.544308384000001, 
            ci.Element.Mn: 5.519113984000001, 
            ci.Element.Fe: 5.519977184000001, 
            ci.Element.Co: 5.494999584, 
            ci.Element.Ni: 5.508116384000001, 
            ci.Element.Cu: 5.465944384000001, 
            ci.Element.Zn: 5.4567685376
        }
        self.assertEqual(expected_real_timescales, test_real_extraction)
        self.assertEqual(expected_real_timescales, processed_wd_data['SDSSJ0002+3209'])

class ManagerTests(unittest.TestCase):
    
    def setUp(self):
        self.real_wd_file = 'BlouinConglomNewTimescales.csv'
        self.real_stellar_comps_file = 'StellarCompositionsSortFE.csv'
        self.test_args_none = Namespace(
            wd_data_filename=None,
            stellar_compositions_filename=None,
            n_live_points=None,
            enhancement_model=None,
            base_dir=get_path_to_test_output(),
            pollution_model_names=None
        )
        self.test_args_invalid_sc = Namespace(
            wd_data_filename=None,
            stellar_compositions_filename='foo',
            n_live_points=None,
            enhancement_model='NonEarthlike',
            base_dir=get_path_to_test_output(),
            pollution_model_names=['bar', 'baz']
        )
        self.test_args_invalid_pm = Namespace(
            wd_data_filename=None,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=None,
            enhancement_model='NonEarthlike',
            base_dir=get_path_to_test_output(),
            pollution_model_names=['bar', 'baz']
        )
        self.test_args_invalid_em = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=1500,
            enhancement_model='Dummy',
            base_dir=get_path_to_test_output(),
            pollution_model_names=['Model_24']
        )
        self.test_args_golden = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=1500,
            enhancement_model='NonEarthlike',
            base_dir=get_path_to_test_output(),
            pollution_model_names=['Model_24']
        )
        self.test_args_invalid_hierarchy = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=1500,
            enhancement_model='NonEarthlike',
            base_dir=get_path_to_test_output(),
            pollution_model_names=['Hierarchy_NOTREAL']
        )
        self.test_args_hierarchy = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=1500,
            enhancement_model='NonEarthlike',
            base_dir=get_path_to_test_output(),
            pollution_model_names=['Hierarchy_Basic']
        )
    
    def test_init(self):       
        with self.assertRaises(TypeError):
            test_mn = mn.Manager()
        
        test_mn = mn.Manager(self.test_args_none)
        self.assertIsNone(test_mn.wd_data_filename)
        self.assertIsNone(test_mn.stellar_compositions_filename)
        self.assertIsNone(test_mn.model_names)
        self.assertIsNone(test_mn.wd_abundances)
        self.assertIsNone(test_mn.wd_errors)
        self.assertIsNone(test_mn.wd_timescales)
        self.assertIsNone(test_mn.stellar_compositions)
        self.assertEqual(0, len(test_mn.models.keys()))
        
        with self.assertRaises(FileNotFoundError):
            test_mn = mn.Manager(self.test_args_invalid_sc)
        
        with self.assertRaises(KeyError):
            test_mn = mn.Manager(self.test_args_invalid_pm)
        with self.assertRaises(KeyError):
            test_mn = mn.Manager(self.test_args_invalid_hierarchy)

        test_mn = mn.Manager(self.test_args_golden)
        self.assertEqual(self.real_stellar_comps_file, test_mn.stellar_compositions_filename)
        self.assertEqual(['Model_24'], test_mn.model_names)
        self.assertEqual(False, test_mn.use_hierarchy)
        self.assertEqual(None, test_mn.parameter_hierarchy)
        self.assertEqual(1500, test_mn.n_live_points)
        self.assertEqual(958, len(test_mn.stellar_compositions))
        self.assertEqual(1, len(test_mn.models.keys()))
        # TODO: Should add another check going over the full contents of the file
        
        test_mn = mn.Manager(self.test_args_hierarchy)
        self.assertEqual(self.real_stellar_comps_file, test_mn.stellar_compositions_filename)
        self.assertEqual(None, test_mn.model_names)
        self.assertEqual(True, test_mn.use_hierarchy)
        self.assertEqual(
            {
                0: [
                    mp.ModelParameter.metallicity,
                    mp.ModelParameter.t_sinceaccretion,
                    mp.ModelParameter.pollution_frac,
                    mp.ModelParameter.accretion_timescale
                ],
                1: [mp.ModelParameter.formation_distance]
            },
            test_mn.parameter_hierarchy)
        self.assertEqual(1500, test_mn.n_live_points)
        self.assertEqual(958, len(test_mn.stellar_compositions))
        self.assertEqual(0, len(test_mn.models.keys()))
        
    def test_load_compositions(self):
        test_mn = mn.Manager(self.test_args_none)
        self.assertIsNone(test_mn.stellar_compositions_filename)
        self.assertIsNone(test_mn.stellar_compositions)
        
        test_mn.load_compositions()
        
        self.assertIsNone(test_mn.stellar_compositions_filename)
        self.assertIsNone(test_mn.stellar_compositions)
        
        with self.assertRaises(FileNotFoundError):
            test_mn.load_compositions('foo')
        
        test_mn.load_compositions(self.real_stellar_comps_file)
        
        self.assertEqual(self.real_stellar_comps_file, test_mn.stellar_compositions_filename)
        self.assertEqual(958, len(test_mn.stellar_compositions))
        
        with self.assertWarns(UserWarning):
            test_mn.load_compositions(self.real_stellar_comps_file)

    def test_load_methods(self):
        test_mn = mn.Manager(self.test_args_none)
        
        self.assertIsNone(test_mn.model_names)
        self.assertEqual(0, len(test_mn.models.keys()))
        
        test_mn.load_models()
        
        self.assertIsNone(test_mn.model_names)
        self.assertEqual(0, len(test_mn.models.keys()))
        
        with self.assertRaises(KeyError):
            test_mn.load_models(['bar', 'baz'])
        
        test_mn.load_models(['Model_24'])
        self.assertEqual(test_mn.model_names, ['Model_24'])
        self.assertEqual(1, len(test_mn.models.keys()))

        with self.assertWarns(UserWarning):
            test_mn.load_models(['Model_24'])

    def test_load_values(self):
        test_mn = mn.Manager(self.test_args_hierarchy)
        example_H_index = 222  # GALEX1931+0117
        self.assertEqual(test_mn.wd_names[example_H_index], 'GALEX1931+0117')
        self.assertEqual(test_mn.wd_types['GALEX1931+0117'], ci.Element.H)
        self.assertEqual(test_mn.wd_masses['GALEX1931+0117'], 0.596)
        # IMPORTANT: For purposes of updating John's code, I filled in all the timescales/logqs in the csv file.
        # Therefore, none of these are actually being overridden any more. Once I've added another system to the data that doesn't have timescales
        # these tests should be updated to check that those missing timescales do actually get replaced with freshly calculated ones
        self.assertAlmostEqual(test_mn.wd_qs['GALEX1931+0117'], -16.192102000000006)  # This one has been overridden from the actual input value of 0
        expected_H_timescales = {
            ci.Element.Al: 0.013820902396929203,
            ci.Element.Ti: 0.007857160053722625,
            ci.Element.Ca: 0.008551305064535511,
            ci.Element.Ni: 0.0058163160846718374,
            ci.Element.Fe: 0.006352677375621979,
            ci.Element.Cr: 0.006908510150094222,
            ci.Element.Mg: 0.01442371724114125,
            ci.Element.Si: 0.012453791771387156,
            ci.Element.Na: 0.004759581852912534,
            ci.Element.O: 0.00700801140594409,
            ci.Element.C: 0.01423467644313019,
            ci.Element.N: 0.008844139822659266
        }
        for key, t in expected_H_timescales.items():
            self.assertAlmostEqual(test_mn.wd_timescales['GALEX1931+0117'][key], expected_H_timescales[key]) # These have been overridden from the actual input value of 0
        expected_H_abundances = {
            ci.Element.Al: -6.2,
            ci.Element.Ti: 0,
            ci.Element.Ca: 0,
            ci.Element.Ni: -6.7,
            ci.Element.Fe: -4.5,
            ci.Element.Cr: -6.1,
            ci.Element.Mg: 0.0,
            ci.Element.Si: -4.75,
            ci.Element.Na: 0.0,
            ci.Element.O: -4.1,  # We ignore C, N for now so C gets overridden
            ci.Element.C: 0.0,
            ci.Element.N: 0.0
        }
        self.assertEqual(test_mn.wd_abundances['GALEX1931+0117'], expected_H_abundances)
        expected_H_errors = {
            ci.Element.Al: 0.2,
            ci.Element.Ti: 0.0,
            ci.Element.Ca: 0.0,
            ci.Element.Ni: 0.3,
            ci.Element.Fe: 0.3,
            ci.Element.Cr: 0.3,
            ci.Element.Mg: 0.0,
            ci.Element.Si: 0.2,
            ci.Element.Na: 0.0,
            ci.Element.O: 0.3,
            ci.Element.C: 0.0,
            ci.Element.N: 0.0
        }
        self.assertEqual(test_mn.wd_errors['GALEX1931+0117'], expected_H_errors)
        example_index = 226  # GD61
        self.assertEqual(test_mn.wd_names[example_index], 'GD61')
        self.assertEqual(test_mn.wd_types['GD61'], ci.Element.He)
        self.assertEqual(test_mn.wd_masses['GD61'], 0.71)
        self.assertEqual(test_mn.wd_qs['GD61'], -6.82)
        expected_timescales = {
            ci.Element.Al: 231415.26328517616, 
            ci.Element.Ti: 153586.13093388823, 
            ci.Element.Ca: 188014.78395681301, 
            ci.Element.Ni: 136332.68876863373, 
            ci.Element.Fe: 138099.46623484572, 
            ci.Element.Cr: 146336.276401041, 
            ci.Element.Mg: 254659.56916522002, 
            ci.Element.Si: 236936.51701234345, 
            ci.Element.Na: 253316.80396232716, 
            ci.Element.O: 333297.4569355848, 
            ci.Element.C: 401783.4096568484, 
            ci.Element.N: 366430.8246804239
        }
        for key, t in expected_H_timescales.items():
            self.assertAlmostEqual(test_mn.wd_timescales['GD61'][key], expected_timescales[key], 4) # These have been overridden from the actual input value of 0
        expected_abundances = {
            ci.Element.Al: 0.0,
            ci.Element.Ti: 0.0,
            ci.Element.Ca: -7.9,
            ci.Element.Ni: 0.0,
            ci.Element.Fe: -7.6,
            ci.Element.Cr: 0.0,
            ci.Element.Mg: -6.69,
            ci.Element.Si: -6.82,
            ci.Element.Na: 0.0,
            ci.Element.O: -5.95,
            ci.Element.C: 0.0,
            ci.Element.N: 0.0
        }
        self.assertEqual(test_mn.wd_abundances['GD61'], expected_abundances)
        expected_errors = {
            ci.Element.Al: 0.0,
            ci.Element.Ti: 0.0,
            ci.Element.Ca: 0.06,
            ci.Element.Ni: 0.0,
            ci.Element.Fe: 0.07,
            ci.Element.Cr: 0.0,
            ci.Element.Mg: 0.05,
            ci.Element.Si: 0.04,
            ci.Element.Na: 0.0,
            ci.Element.O: 0.04,
            ci.Element.C: 0.0,
            ci.Element.N: 0.0
        }
        self.assertEqual(test_mn.wd_errors['GD61'], expected_errors)
        expected_upper_bounds = {
            ci.Element.Al: -7.8,
            ci.Element.Ti: -8.6,
            ci.Element.Ca: None,
            ci.Element.Ni: -8.8,
            ci.Element.Fe: None,
            ci.Element.Cr: -8.0,
            ci.Element.Mg: None,
            ci.Element.Si: None,
            ci.Element.Na: -6.8,
            ci.Element.O: None,
            ci.Element.C: None,  # These ones at the end are still None because we're ignoring O C and N in the current version
            ci.Element.N: None
        }
        self.assertEqual(test_mn.wd_abundance_upper_bounds['GD61'], expected_upper_bounds)
        expected_lower_bounds = {  # TODO At some point it would be neat to have a system that actually has lower bounds
            ci.Element.Al: None,
            ci.Element.Ti: None,
            ci.Element.Ca: None,
            ci.Element.Ni: None,
            ci.Element.Fe: None,
            ci.Element.Cr: None,
            ci.Element.Mg: None,
            ci.Element.Si: None,
            ci.Element.Na: None,
            ci.Element.O: None,
            ci.Element.C: None,
            ci.Element.N: None
        }
        self.assertEqual(test_mn.wd_abundance_lower_bounds['GD61'], expected_lower_bounds)

class EnhancementTests(unittest.TestCase):
    
    def test_enhancements(self):
        elements = [
            ci.Element.Al,
            ci.Element.Ti,
            ci.Element.Ca,
            ci.Element.Ni,
            ci.Element.Fe,
            ci.Element.Cr,
            ci.Element.Mg,
            ci.Element.Si,
            ci.Element.Na,
            ci.Element.O,
            ci.Element.C,
            ci.Element.N
        ]
        test_disc_abundances = {
            ci.Element.Al: 0.1,
            ci.Element.Ti: 0.2,
            ci.Element.Ca: 0.3,
            ci.Element.Ni: 0.4,
            ci.Element.Fe: 0.5,
            ci.Element.Cr: 0.6,
            ci.Element.Mg: 0.7,
            ci.Element.Si: 0.8,
            ci.Element.Na: 0.9,
            ci.Element.O: 0.10,
            ci.Element.C: 0.11,
            ci.Element.N: 0.12
        }
        test_N_c = 0.17
        test_N_o = 0.1
        test_f_c = 0.6
        test_f_o = 0.005
        test_pressure = 30
        test_fO2 = -2
        geo_model = gi.GeologyModel()
        with self.assertRaises(TypeError):
            enhancement_model = em.EnhancementModel()
        with self.assertRaises(ValueError):
            enhancement_model = em.EnhancementModel('Dummy')
        enhancement_model1 = em.EnhancementModel('Earthlike')
        enhancement_model2 = em.EnhancementModel('NonEarthlike')
        enhancements1, ignore = enhancement_model1.find_enhancements(
            geo_model,
            test_disc_abundances,
            elements,
            test_N_c,
            test_N_o,
            test_f_c,
            test_f_o,
            test_pressure,
            test_fO2,
            False  # Not normalising (so that it will exactly match the original E_Al etc values)
        )
        enhancements2, ignore = enhancement_model2.find_enhancements(
            geo_model,
            test_disc_abundances,
            elements,
            test_N_c,
            test_N_o,
            test_f_c,
            test_f_o,
            test_pressure,
            test_fO2
        )
        self.assertEqual(12, len(enhancements1))
        self.assertEqual(12, len(enhancements2))
        expected_enhancements1 = {
            ci.Element.Al: 0.03353491807681976,
            ci.Element.Ti: 0.01446450809464508,
            ci.Element.Ca: 0.10323266195769829,
            ci.Element.Ni: 1.3349102755010887,
            ci.Element.Fe: 1.5734363502203488,
            ci.Element.Cr: 1.4482813550631326,
            ci.Element.Mg: 0.37007595957732214,
            ci.Element.Si: 0.676948841654683,
            ci.Element.Na: 0.1027164242338653,
            ci.Element.O: 0.04799630310971293,
            ci.Element.C: 0.3545932838588374,
            ci.Element.N: 0.4173428274782206
        }

        expected_enhancements2 = {
            ci.Element.Al: 0.007293931394696239,
            ci.Element.Ti: 0.000209760118540284,
            ci.Element.Ca: 0.005291198990178664,
            ci.Element.Ni: 0.03041874633682828,
            ci.Element.Fe: 0.5011740367392136,
            ci.Element.Cr: 0.004261300654280641,
            ci.Element.Mg: 0.07857423349502184,
            ci.Element.Si: 0.09169358252255623,
            ci.Element.Na: 0.0009710940033330876,
            ci.Element.O: 0.2506409926951525,
            ci.Element.C: 0.006146858954405018,
            ci.Element.N: 2.21339830538792e-05
        }
        self.assertEqual(expected_enhancements1, enhancements1)
        self.assertEqual(expected_enhancements2, enhancements2)
        

#TODO: See if it's possible to put this in the GeologyTests class, or otherwise not have this lying around
def get_mock_partition_coefficient(element, pressure=None, fO2=None):
    mock_pc_dict = {
        ci.Element.Ni: 50,
        ci.Element.Cr: 0.6
    }
    to_add = 0
    # This bit is just to ensure that we have a different D for Earth vs our mock value
    if pressure == 60:
        to_add = 3
    return mock_pc_dict[element] + to_add

class AbundanceTests(unittest.TestCase):
    
    def load_generic_float_data_csv(self, input_filename):
        # TODO: Change this to a with clause to prevent ResourceWarnings
        generic_csv = open(get_path_to_data() + input_filename)
        generic_list =  [row for row in csv.reader(generic_csv)]
        generic_array = np.asarray(generic_list)
        return generic_array.astype(np.float)
    
    def test_abundances(self):
        linear_d_formation = 0.5
        z_formation = 0.1
        t_formation = 1.5
        fe_star = 156
        ld._live_stellar_compositions = self.load_generic_float_data_csv('StellarCompositionsSortFE.csv')
        elements = [
            ci.Element.Al,
            ci.Element.Ti,
            ci.Element.Ca,
            ci.Element.Ni,
            ci.Element.Fe,
            ci.Element.Cr,
            ci.Element.Mg,
            ci.Element.Si,
            ci.Element.Na,
            ci.Element.O,
            ci.Element.C,
            ci.Element.N
        ]
        calculated_results = am.get_all_abundances(elements, linear_d_formation, z_formation, t_formation, fe_star)
        #calculated_results = dict()
        #calculated_results[ci.Element.Al] = am.Al(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.Ti] = am.Ti(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.Ca] = am.Ca(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.Ni] = am.Ni(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.Fe] = am.Fe(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.Cr] = am.Cr(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.Mg] = am.Mg(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.Si] = am.Si(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.Na] = am.Na(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.O] = am.O(
        #    dm.T_disc(linear_d_formation, t_formation),
        #    t_formation,
        #    linear_d_formation,
        #    fe_star,
        #    z_formation,
        #    calculated_results[ci.Element.Al],
        #    calculated_results[ci.Element.Ti],
        #    calculated_results[ci.Element.Ca],
        #    calculated_results[ci.Element.Ni],
        #    calculated_results[ci.Element.Fe],
        #    calculated_results[ci.Element.Cr],
        #    calculated_results[ci.Element.Mg],
        #    calculated_results[ci.Element.Si],
        #    calculated_results[ci.Element.Na]
        #)
        #calculated_results[ci.Element.C] = am.C(linear_d_formation, z_formation, t_formation)
        #calculated_results[ci.Element.N] = am.Nz(linear_d_formation, z_formation, t_formation)

        golden_results = {
            ci.Element.Al: 0.9999810892034007,
            ci.Element.Ti: 0.9999456663756952,
            ci.Element.Ca: 0.9993663361988419,
            ci.Element.Ni: 0.99509988719541,
            ci.Element.Fe: 0.9942026481672976,
            ci.Element.Cr: 0.9927085137755711,
            ci.Element.Mg: 0.994646712899971, 
            ci.Element.Si: 0.9941031003626686,
            ci.Element.Na: 0.91691422301259,
            ci.Element.O: 0.20384101469485028, 
            ci.Element.C: 0.0,
            ci.Element.N: 0.0
        }
        for element in golden_results.keys():
            required_dp = 7
            self.assertAlmostEqual(golden_results[element], calculated_results[element], required_dp)

class PamelaTests(unittest.TestCase):
    
    def test_init_and_gammas(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        self.assertEqual(2.7, pamela.nbot)
        self.assertEqual(3.0, 10**pamela.log10_gammaFe_sil)
        self.assertEqual([
            ci.Element.Hf,
            ci.Element.U,
            ci.Element.Ta,
            ci.Element.Pb,
            ci.Element.Nb,
            ci.Element.Si,
            ci.Element.Mn,
            ci.Element.Zn,
            ci.Element.Ga,
            ci.Element.V,
            ci.Element.Cr,
            ci.Element.Cu,
            #ci.Element.Ti,
            ci.Element.Fe,
            ci.Element.W,
            ci.Element.P,
            ci.Element.Co,
            ci.Element.Ni,
            ci.Element.O,
            ci.Element.C,
            ci.Element.S
        ], pamela.partitioners)
        self.assertEqual([
            ci.Element.Hf,
            ci.Element.U,
            ci.Element.Ta,
            ci.Element.Pb,
            ci.Element.Nb,
            ci.Element.Si,
            ci.Element.Mn,
            ci.Element.Zn,
            ci.Element.Ga,
            ci.Element.V,
            ci.Element.Cr,
            ci.Element.Cu,
            ci.Element.Fe,
            ci.Element.W,
            ci.Element.P,
            ci.Element.Co,
            ci.Element.Ni,
            ci.Element.O,
            ci.Element.C,
            ci.Element.S,
            ci.Element.N,
            ci.Element.Na,
            ci.Element.Mg,
            ci.Element.Al,
            ci.Element.Ti,
            ci.Element.Ca
        ], pamela.ele_set)
        self.assertEqual(pamela.v[ci.Element.O], -2)  # Just picking 1 example
        self.assertEqual(pamela.source[ci.Element.O], 'si')
        self.assertEqual(pamela.eps[1][2], -7.13737093741606)  # Should correspond to O/Si
        
        self.assertEqual([
        ci.Element.C, 
        ci.Element.O, 
        ci.Element.Si, 
        ci.Element.P,
        ci.Element.S,
        ci.Element.Ti,
        ci.Element.V,
        ci.Element.Cr,
        ci.Element.Mn, 
        ci.Element.Co,
        ci.Element.Ni,
        ci.Element.Ga, 
        ci.Element.Ge,
        ci.Element.Zr,
        ci.Element.Nb,
        ci.Element.Mo,
        ci.Element.Hf,
        ci.Element.Ta,
        ci.Element.W,
        ci.Element.Re,
        ci.Element.Cu,
        ci.Element.As,
        ci.Element.Zn
        ], pamela.logg0name)
        expected_logg0 = [0, 0, 0, 0, 0, -5.52, 0, 0, 0.36, 0, 0, 0, 0, -3.3, -1.61, 0, 0, -3.22, 0, 0, 2.14943391, 0, 0]
        self.assertEqual(len(expected_logg0), len(pamela.logg0))
        for exp, g0 in zip(expected_logg0, pamela.logg0):
            self.assertAlmostEqual(exp, g0)
        self.assertEqual([
        ci.Element.Fe,
        ci.Element.C,
        ci.Element.O,
        ci.Element.Si,
        ci.Element.P,
        ci.Element.S,
        ci.Element.Ti,
        ci.Element.V,
        ci.Element.Cr,
        ci.Element.Mn,
        ci.Element.Co,
        ci.Element.Ni,
        ci.Element.Ga,
        ci.Element.Ge,
        ci.Element.Zr,
        ci.Element.Nb,
        ci.Element.Mo,
        ci.Element.Hf,
        ci.Element.Ta,
        ci.Element.W,
        ci.Element.Re,
        ci.Element.Cu,
        ci.Element.As,
        ci.Element.Zn
        ], pamela.alloy_els)
        expected_C = [
        0.792963331154798,
        0.00862426542303857,
        0,
        0.106772634240816,
        0.00334429043701544,
        0.0306888135084857,
        0,
        0.00015250586745139,
        0.00896486718175481,
        0.000282826236838207,
        0.00219711477432245,
        0.0458868579200079,
        0,
        0.0000142601943827252,
        0,
        0,
        0.00000269924046268804,
        0,
        0,
        0.000000132412544707829,
        0.0000000639729121935093,
        0.000101880972047922,
        0.00000345646312151692,
        0
        ]
        self.assertEqual(len(expected_C), len(pamela.alloy))
        for exp, C in zip(expected_C, pamela.alloy):
            self.assertAlmostEqual(exp, C)
        
        old_gammas = pamela.calculate_gammas()
        self.assertEqual(0.19087291612663068, old_gammas[ci.Element.O])
        self.assertEqual(5.938117878057908, old_gammas[ci.Element.Si])
        self.assertEqual(0.8762221554491502, old_gammas[ci.Element.Fe])
        self.assertEqual(0.6472702432098675, old_gammas[ci.Element.Mn])
        new_gammas = pamela.calculate_gammas(pamela.T0, np.array(expected_C))
        self.assertEqual(new_gammas, old_gammas)
    
    def test_fischer_eps(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        pamela_es = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat',
            get_path_to_feni() + 'data/e_param_fischer_epsilon_update.dat'
        )
        self.assertIsNone(pamela.fischer_es.get((ci.Element.O, ci.Element.Si)))
        self.assertIsNone(pamela_es.fischer_es.get((ci.Element.O, ci.Element.Si)))
        self.assertEqual(-0.036, pamela_es.fischer_es[(ci.Element.V, ci.Element.C)])
        eps = pamela.construct_eps_table(3300, pamela.alloy[1:])
        eps_e = pamela_es.construct_eps_table(3300, pamela_es.alloy[1:])
        self.assertEqual(-7.13737093741606*(1873/3300), eps[1][2])  # O - Si default
        self.assertEqual(-5.94995266482236*(1873/3300), eps[6][0])  # V - C default
        self.assertEqual(-0.22918003981901158, eps_e[6][0])  # V - C gets overriden by es
        gammas = pamela.calculate_gammas(3300)
        gammas_e = pamela_es.calculate_gammas(3300)
        self.assertEqual(0.3906320799707601, gammas[ci.Element.O])
        self.assertEqual(1.186036978787397, gammas[ci.Element.V])
        self.assertEqual(1.2185224417219145, gammas_e[ci.Element.V]) # I'm assuming this is OK, I didn't check by hand...TODO

    def test_convert_abundance_dict_to_layer_composition(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        test_abundances = {ci.Element.Fe: {gi.Layer.core: 0.8}, ci.Element.Ni: {gi.Layer.core: 0.15}, ci.Element.C: {gi.Layer.core: 0.05}}
        test_result = pamela.convert_abundance_dict_to_layer_composition(gi.Layer.core, None)
        self.assertIsNone(test_result)
        test_result = pamela.convert_abundance_dict_to_layer_composition(gi.Layer.core, test_abundances)
        exp_result = {
            ci.Element.Hf: 0,
            ci.Element.U: 0,
            ci.Element.Ta: 0,
            ci.Element.Pb: 0,
            ci.Element.Nb: 0,
            ci.Element.Si: 0,
            ci.Element.Mn: 0,
            ci.Element.Zn: 0,
            ci.Element.Ga: 0,
            ci.Element.V: 0,
            ci.Element.Cr: 0,
            ci.Element.Cu: 0,
            ci.Element.Ti: 0,
            ci.Element.Fe: 0.8,
            ci.Element.W: 0, 
            ci.Element.P: 0, 
            ci.Element.Co: 0,
            ci.Element.Ni: 0.15,
            ci.Element.O: 0,
            ci.Element.C: 0.05,
            ci.Element.S: 0,
            ci.Element.N: 0,
            ci.Element.Na: 0,
            ci.Element.Mg: 0,
            ci.Element.Al: 0,
            ci.Element.Ca: 0
        } 
        self.assertEqual(len(exp_result), len(test_result))
        for el, el_val in exp_result.items():
            self.assertEqual(el_val, test_result[el])

    def test_convert_core_composition_to_alloy_composition(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        test_abundances = {ci.Element.Fe: 0.8, ci.Element.Ni: 0.15, ci.Element.C: 0.05}
        test_result = pamela.convert_core_composition_to_alloy_composition(test_abundances)
        exp_result = [0.8, 0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.assertEqual(len(exp_result), len(test_result))
        for er, tr in zip(exp_result, test_result):
            self.assertEqual(er, tr)
        
    def test_get_all_partition_coefficients(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        Ds = pamela.get_all_partition_coefficients(54, -2)
        exp_Ds = {
            ci.Element.Hf: 0.0,
            #ci.Element.U: 0,
            ci.Element.Ta: 0.017342227950395302,
            #ci.Element.Pb: 0,
            ci.Element.Nb: 0.5237944840651807,
            ci.Element.Si: 0.6839893208888564,
            ci.Element.Mn: 5.053769920186087,
            ci.Element.Zn: 0.45473960145421244,
            ci.Element.Ga: 1.08342627642091,
            ci.Element.V: 1.5196568184310144,
            ci.Element.Cr: 3.7650839227354473,
            ci.Element.Cu: 9.496099758919135,
            ci.Element.Ti: 0,
            ci.Element.Fe: 32.044056639589776,
            ci.Element.W: 149.64496964858944,
            ci.Element.P: 327.79381520839615,
            ci.Element.Co: 61.87014806303029,
            ci.Element.Ni: 59.003334996057504,
            ci.Element.O: 0.0989546319848992,
            ci.Element.C: 2633.570873664594,
            ci.Element.S: 131.37684416329742,
            ci.Element.N: 0,
            ci.Element.Na: 0,
            ci.Element.Mg: 0,
            ci.Element.Al: 0,
            ci.Element.Ca: 0
        }
        self.assertEqual(len(exp_Ds.keys()), len(Ds.keys()))
        for k, v in exp_Ds.items():
            self.assertAlmostEqual(v, Ds[k], 10)
        # Now run with some dummy initial abundances
        test_abundances = {
            ci.Element.C: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.01},
            ci.Element.N: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.01},
            ci.Element.O: {gi.Layer.mantle: 0.31, gi.Layer.core: 0.01},
            ci.Element.Na: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.01},
            ci.Element.Mg: {gi.Layer.mantle: 0.2, gi.Layer.core: 0.01},
            ci.Element.Al: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.01},
            ci.Element.Si: {gi.Layer.mantle: 0.2, gi.Layer.core: 0.01},
            ci.Element.Ca: {gi.Layer.mantle: 0.2, gi.Layer.core: 0.01},
            ci.Element.Ti: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.01},
            ci.Element.Cr: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.01},
            ci.Element.Fe: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.79},
            ci.Element.Ni: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.1},
            ci.Element.S: {gi.Layer.mantle: 0.01, gi.Layer.core: 0.01}
        }
        Ds = pamela.get_all_partition_coefficients(54, -2, test_abundances)
        exp_Ds = {
            ci.Element.Hf: 0,
            #ci.Element.U: 0,
            ci.Element.Ta: 0.09111458357273382,
            #ci.Element.Pb: 0,
            ci.Element.Nb: 0.5702315612044024,
            ci.Element.Si: 0.4062020447208237,
            ci.Element.Mn: 4.215817726905471,
            ci.Element.Zn: 0.4261041584435114,
            ci.Element.Ga: 1.015201755628057,
            ci.Element.V: 2.0028197287226903,
            ci.Element.Cr: 3.602352753002689,
            ci.Element.Cu: 10.454414015555049,
            ci.Element.Ti: 0,
            ci.Element.Fe: 30.026207842607352,
            ci.Element.W: 142.68467170875175,
            ci.Element.P: 587.1232514592169,
            ci.Element.Co: 58.00656313371489,
            ci.Element.Ni: 58.182664983343535,
            ci.Element.O: 0.17223566410074867,
            ci.Element.C: 4260.104869037981,
            ci.Element.S: 78.10203716903459,
            ci.Element.N: 0,
            ci.Element.Na: 0,
            ci.Element.Mg: 0,
            ci.Element.Al: 0,
            ci.Element.Ca: 0
        } 
        self.assertEqual(len(exp_Ds.keys()), len(Ds.keys()))
        for k, v in exp_Ds.items():
            self.assertAlmostEqual(v, Ds[k], 10)
        # Now replicate a sulfur calculation from here https://www.researchgate.net/publication/308404683_Sulfur_partitioning_calculator
        s_core_composition = {
            #See sheet 2 of https://docs.google.com/spreadsheets/d/1Jl5jEcamzo5bmx7APwfZvPIy4dj4cikFf54h_dDudZA
            ci.Element.Si: 0.001450960757,
            ci.Element.C: 0.1997948494,
            ci.Element.Fe: 0.6937848662,
            ci.Element.Ni: 0.029545891,
            ci.Element.O: 0.003679070409,
            ci.Element.S: 0.0717443623
        }
        s_mantle_composition = {
            ci.Element.S: 0.05,
            ci.Element.Mg: 0.336,
            ci.Element.Ca: 0.013,
            ci.Element.Na: 0.022,
            ci.Element.K: 0.002,
            ci.Element.Ti: 0.002,
            ci.Element.Fe: 0.029,
            ci.Element.O: 0.0000001  #This is a bit of a hack to get the inputs we need
        }
        Ds = pamela.calculate_D(5, 2073, -2, 2.7, 3, s_core_composition, s_mantle_composition)
        self.assertTrue(abs(10.57 - Ds[ci.Element.S])/10.57 < 0.01)
        s_core_composition = {
            #See sheet 2 of https://docs.google.com/spreadsheets/d/1Jl5jEcamzo5bmx7APwfZvPIy4dj4cikFf54h_dDudZA/edit#gid=1680331727
            ci.Element.Si: 0.01819771067,
            ci.Element.C: 0.2120095876,
            ci.Element.Fe: 0.6798027024,
            ci.Element.Ni: 0.03475429324,
            ci.Element.O: 0,
            ci.Element.S: 0.05523570611
        }
        s_mantle_composition = {
            ci.Element.S: 0.05,
            ci.Element.Mg: 0.354,
            ci.Element.Ca: 0.014,
            ci.Element.Na: 0.024,
            ci.Element.K: 0.002,
            ci.Element.Ti: 0.002,
            ci.Element.Fe: 0.013,
            ci.Element.O: 0.0000001  #This is a bit of a hack to get the inputs we need
        }
        Ds = pamela.calculate_D(6, 2073, -2, 2.7, 3, s_core_composition, s_mantle_composition)
        self.assertTrue(abs(4.19 - Ds[ci.Element.S])/4.19 < 0.05)
        
    def test_peridotite_liquidus(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        self.assertEqual(2020, pamela.peridotite_liquidus(0))
        self.assertEqual(-8422, pamela.peridotite_liquidus(-100))
        self.assertEqual(2228.84, pamela.peridotite_liquidus(2))
        self.assertEqual(3648.6016088871866, pamela.peridotite_liquidus(50))
        self.assertEqual(15587.101608887186, pamela.peridotite_liquidus(500))
        
    def test_calculate_logDFe(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        self.assertEqual(-0.1989700043360188, pamela.calculate_logDFe(np.log10(2), -7, 4))
        
    def test_calculate_logKD_app(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        self.assertEqual((-4.987777777777778, 0.4946204444444445), pamela.calculate_logKD_app(200, 4500, 3.2, ci.Element.W))
    
    def test_convert_composition_to_mass(self):
        pamela = pam.PartitionModel(
            get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            get_path_to_feni() + 'data/composition.dat'
        )
        test_comp = {
            ci.Element.At: 0.6,  #210
            ci.Element.Mc: 0.4   #290
        }
        exp_comp = pamela.convert_composition_to_mass(test_comp)
        self.assertEqual(126/242, exp_comp[ci.Element.At])
        self.assertEqual(116/242, exp_comp[ci.Element.Mc])

class GeologyTests(unittest.TestCase):
    
    def test_getters_and_setters(self):
        geology_model = gi.GeologyModel()
        self.assertEqual(geology_model.get_earth_differentiation_pressure(), 54)
        self.assertEqual(geology_model.get_earth_oxygen_fugacity(), -2)
        self.assertEqual(geology_model.get_earth_mean_molecular_weight(), 26)
        self.assertEqual(geology_model.get_earth_layer_number_fraction(gi.Layer.core), 0.17)
        self.assertEqual(geology_model.get_bulk_abundance(ci.Element.C), 0.0015917333226106879)
        self.assertEqual(geology_model.get_bulk_abundance(ci.Element.V), 0.00005395470924151831)
        self.assertEqual(geology_model.get_core_abundance(ci.Element.Si), 0.11082406539923671)
        self.assertEqual(geology_model.get_mantle_abundance(ci.Element.Si), 0.15814987756236518)
        self.assertEqual(geology_model.get_crust_abundance(ci.Element.Si), 0.18401433679401685)
        self.assertEqual(geology_model.get_mu(ci.Element.Pb), 1)
        total_bulk = 0
        total_core = 0
        total_mantle = 0
        total_crust = 0
        for el in ci.Element:
            total_bulk += geology_model.get_bulk_abundance(el) if geology_model.get_bulk_abundance(el) is not None else 0
            total_core += geology_model.get_core_abundance(el) if geology_model.get_core_abundance(el) is not None else 0
            total_mantle += geology_model.get_mantle_abundance(el) if geology_model.get_mantle_abundance(el) is not None else 0
            total_crust += geology_model.get_crust_abundance(el) if geology_model.get_crust_abundance(el) is not None else 0
        self.assertAlmostEqual(total_bulk, 1, 10)
        self.assertAlmostEqual(total_core, 1, 10)
        self.assertAlmostEqual(total_mantle, 1, 2)  # Can be more lenient here since mantle is calculated from the others
        self.assertAlmostEqual(total_crust, 1, 10)
        
    def test_form_a_planet_iteratively(self):
        # This is more of an integration test right now - should revise this!
        geology_model = gi.GeologyModel()
        abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(54, -2)
        self.assertEqual(len(abundances), 24)
        self.assertTrue(w_met > 0)
        self.assertTrue(w_met < 1)
        self.assertEqual(len(Ds), 24)
        self.assertEqual(all_Ds, None)  # We didn't tell it to store these
        
    def test_find_system_specific_abundances(self):
        geology_model = gi.GeologyModel()
        abundances = {'A': {gi.Layer.core: 0.5, gi.Layer.mantle: 0.25}, 'B': {gi.Layer.core: 0.5, gi.Layer.mantle: 0.75}}
        test_output = geology_model.find_system_specific_abundances(abundances, 0.3)
        required_dp = 10
        self.assertAlmostEqual(test_output['A'], 0.325, required_dp)
        self.assertAlmostEqual(test_output['B'], 0.675, required_dp)
        
    def test_find_system_specific_abundances_with_scalings(self):
        geology_model = gi.GeologyModel({'A': 5, 'B': 25})
        abundances = {'A': {gi.Layer.core: 0.5, gi.Layer.mantle: 0.25}, 'B': {gi.Layer.core: 0.5, gi.Layer.mantle: 0.75}}
        test_output = geology_model.find_system_specific_abundances(abundances, 0.3)
        required_dp = 10
        self.assertAlmostEqual(test_output['A'], 0.325, required_dp)
        self.assertAlmostEqual(test_output['B'], 0.675, required_dp)
    
    def test_pressure_radius_maths(self):
        geology_model = gi.GeologyModel()
        test_input_dict = {
            ci.Element.Fe: {
                gi.Layer.bulk: 0.3,
                gi.Layer.mantle: 0.1,
                gi.Layer.core: 0.9
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.4,
                gi.Layer.mantle: 0.2,
                gi.Layer.core: 0.1
            },
            ci.Element.Mg: {
                gi.Layer.bulk: 0.3,
                gi.Layer.mantle: 0.7
            }
        }
        expected_mass_dict = {
            ci.Element.Fe: {
                gi.Layer.bulk: 0.5502946335310696,
                gi.Layer.mantle: 0.21647194722030566,
                gi.Layer.core: 0.9691498715783141,
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.21020476537711125,
                gi.Layer.mantle: 0.12403383234229279,
                gi.Layer.core: 0.030850128421685913
            },
            ci.Element.Mg: {
                gi.Layer.bulk: 0.23950060109181923,
                gi.Layer.mantle: 0.6594942204374017,
                gi.Layer.core: 0.0
            }
        }
        
        mass_dict = geology_model.convert_number_abundances_to_mass(test_input_dict)
        self.assertEqual(expected_mass_dict, mass_dict)
        
        self.assertEqual(geology_model.calculate_radius_and_mass(geology_model.get_earth_differentiation_pressure(), dict()), (None, None))
        self.assertEqual(geology_model.calculate_radius_and_mass(geology_model.get_earth_differentiation_pressure(), test_input_dict), (6611.463949084282, 1.3620921110074933))  # Earth-like values, but with some dummy functions
        self.assertEqual(geology_model.calculate_radius_and_mass(135, test_input_dict, True), (6662.893525790423, 1.3962334556231406))  # Earth-like values, but finding a lower bound that corresponds to interpreting the Pressure as CMB
    
class ModelParameterTests(unittest.TestCase):
    
    def test_get_model_params_in_order(self):
        model_params_in_order = [
            mp.ModelParameter.metallicity,
            mp.ModelParameter.t_sinceaccretion,
            mp.ModelParameter.formation_distance,
            mp.ModelParameter.feeding_zone_size,
            mp.ModelParameter.parent_core_frac,
            mp.ModelParameter.parent_crust_frac,
            mp.ModelParameter.fragment_core_frac,
            mp.ModelParameter.fragment_crust_frac,
            mp.ModelParameter.pollution_frac,
            mp.ModelParameter.accretion_timescale,
            mp.ModelParameter.pressure,
            mp.ModelParameter.oxygen_fugacity
        ]
        self.assertEqual(mp.get_model_params_in_order(), model_params_in_order)
    
    def test_model_parameters_logic(self):
        mp.model_definitions_dict = {
            'Model_Mock1': {
            },
            'Model_Mock2': {
                mp.ModelParameter.accretion_timescale: False,
            },
            'Model_Mock2': {
                mp.ModelParameter.formation_distance: False,
                mp.ModelParameter.feeding_zone_size: True,
                mp.ModelParameter.fragment_crust_frac: True,
                mp.ModelParameter.pollution_frac: False,
                mp.ModelParameter.accretion_timescale: True,
                mp.ModelParameter.pressure: True
            },
            'Model_24': {
                mp.ModelParameter.metallicity: True,
                mp.ModelParameter.t_sinceaccretion: True,
                mp.ModelParameter.formation_distance: True,
                mp.ModelParameter.feeding_zone_size: True,
                mp.ModelParameter.parent_core_frac: True,
                mp.ModelParameter.parent_crust_frac: True,
                mp.ModelParameter.fragment_core_frac: True,
                mp.ModelParameter.fragment_crust_frac: True,
                mp.ModelParameter.pollution_frac: True,
                mp.ModelParameter.accretion_timescale: True,
                mp.ModelParameter.pressure: False
            },
            'Model_Andy': {
                mp.ModelParameter.metallicity: True,
                mp.ModelParameter.t_sinceaccretion: True,
                mp.ModelParameter.formation_distance: True,
                mp.ModelParameter.feeding_zone_size: True,
                mp.ModelParameter.parent_core_frac: True,
                mp.ModelParameter.parent_crust_frac: True,
                mp.ModelParameter.fragment_core_frac: True,
                mp.ModelParameter.fragment_crust_frac: True,
                mp.ModelParameter.pollution_frac: True,
                mp.ModelParameter.accretion_timescale: True,
                mp.ModelParameter.pressure: True
            }
        }
        
        golden_param_indices = {
            'Model_Mock1': {
            },
            'Model_Mock2': {
            },
            'Model_Mock2': {
                mp.ModelParameter.feeding_zone_size: 0,
                mp.ModelParameter.fragment_crust_frac: 1,
                mp.ModelParameter.accretion_timescale: 2,
                mp.ModelParameter.pressure: 3
            },
            'Model_24': {
                mp.ModelParameter.metallicity: 0,
                mp.ModelParameter.t_sinceaccretion: 1,
                mp.ModelParameter.formation_distance: 2,
                mp.ModelParameter.feeding_zone_size: 3,
                mp.ModelParameter.parent_core_frac: 4,
                mp.ModelParameter.parent_crust_frac: 5,
                mp.ModelParameter.fragment_core_frac: 6,
                mp.ModelParameter.fragment_crust_frac: 7,
                mp.ModelParameter.pollution_frac: 8,
                mp.ModelParameter.accretion_timescale: 9
            },
            'Model_Andy': {
                mp.ModelParameter.metallicity: 0,
                mp.ModelParameter.t_sinceaccretion: 1,
                mp.ModelParameter.formation_distance: 2,
                mp.ModelParameter.feeding_zone_size: 3,
                mp.ModelParameter.parent_core_frac: 4,
                mp.ModelParameter.parent_crust_frac: 5,
                mp.ModelParameter.fragment_core_frac: 6,
                mp.ModelParameter.fragment_crust_frac: 7,
                mp.ModelParameter.pollution_frac: 8,
                mp.ModelParameter.accretion_timescale: 9,
                mp.ModelParameter.pressure: 10
            }
        }
        
        golden_param_uses_pressure = {
            'Model_Mock1': False,
            'Model_Mock2': False,
            'Model_Mock2': True,
            'Model_24': False,
            'Model_Andy': True
        }
        
        for model_name, model in mp.model_definitions_dict.items():
            self.assertEqual(golden_param_uses_pressure[model_name], mp.model_uses_parameter(model_name, mp.ModelParameter.pressure))
            self.assertEqual(golden_param_indices[model_name], mp.parameter_indices(model_name))

class WhiteDwarfModelTests(unittest.TestCase):
    
    def test_white_dwarf_model(self):
        planetesimal_abundance = np.array([0.2, 0.3, 0.4, 0.5])
        non_zero_planetesimal_abundance = np.array([0.3, 0.5])
        wd_timescales = np.array([50000, 60000, 70000, 80000])
        non_zero_wd_timescales = np.array([60000, 80000])
        args = {
            'Early': [0, 200000, planetesimal_abundance, non_zero_planetesimal_abundance, wd_timescales, non_zero_wd_timescales, -6],
			'Medium': [0.1, 200000, planetesimal_abundance, non_zero_planetesimal_abundance, wd_timescales, non_zero_wd_timescales, -6],
			'Late': [0.4, 200000, planetesimal_abundance, non_zero_planetesimal_abundance, wd_timescales, non_zero_wd_timescales, -6]
        }
        expected_results = {
            'Early': [-6.6020599913279625, -6.425968732272281, -6.301029995663981, -6.204119982655925],
            'Medium': [-6.698032879030232, -6.470520647186045, -6.30670962705419, -6.179429567676116],
            'Late': [-7.305484188418363, -6.7684300218686255, -6.379649072648516, -6.081151879145035]
        }
        expected_lifetime_results = {
            'Early': [-6.698970004336019, -6.522878745280337, -6.3979400086720375, -6.301029995663981 ],
            'Medium': [-6.369503928147244, -6.254407893195727, -6.175686043374323, -6.114496075232561],
            'Late': [-4.930950375586853, -5.025940353021858, -5.089882029233725, -5.131808928937756],
        }
        for arg_name, arg_set in args.items():
            result = wdm.process_abundances(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6])
            expected_result = expected_results[arg_name]
            self.assertEqual(len(result), len(expected_result))
            for i, val in enumerate(result):
                self.assertEqual(val, expected_result[i])
            lifetime_result = wdm.process_abundances(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], False)
            expected_lifetime_result = expected_lifetime_results[arg_name]
            self.assertEqual(len(lifetime_result), len(expected_lifetime_result))
            for i, val in enumerate(lifetime_result):
                self.assertEqual(val, expected_lifetime_result[i])

class CompleteModelTests(unittest.TestCase):
    
    def load_generic_float_data_csv(self, input_filename):
        # TODO: Change this to a with clause to prevent ResourceWarnings
        generic_csv = open(get_path_to_data() + input_filename)
        generic_list =  [row for row in csv.reader(generic_csv)]
        generic_array = np.asarray(generic_list)
        return generic_array.astype(np.float)
    
    def test_complete_model(self):
        manager = mn.Manager(
            Namespace(
                wd_data_filename='BlouinConglomNewTimescales.csv',
                stellar_compositions_filename='StellarCompositionsSortFE.csv',
                n_live_points = 0, # This argument shouldn't matter
                pollution_model_names=['Model_24'],
                enhancement_model='NonEarthlike',
                base_dir=get_path_to_test_output()
            )
        )
        manager.publish_live_data(0)
        args = {
            # Args are: fe_star, t_sinceaccretion, d_formation, z_formation, N_c, N_o, f_c, f_o, pollutionfraction, t_disc
           # 'Run1': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           # 'Run2': [145, 0, 2, 0.05, 0.17, 0.01, 0.17, 0, -1, 5.6],
           # 'Run3': [145, 0.2, 2, 0.05, 0.17, 0.01, 0.17, 0, -1, 5.6],
           # 'Run4': [145, 0.5, 2, 0.05, 0.17, 0.01, 0.17, 0, -1, 5.6],
            'SDSSJ0002+3209': [88, 0.195758, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768]  #TODO Add more of these! Try to explore all possible code branches
        }
        for arg_name, arg_set in args.items():
            print('CompleteModelTests ' + arg_name)
            test_result, ignore = cm.complete_model_calculation(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], 10**(arg_set[9]), 54, -2, 'NonEarthlike')
            golden_result = ocm.PWDCodeMultiple(0, 1, False, arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], arg_set[9])
            # We're comparing apples to oranges here because the model is somewhat different in the golden_result vs the test_result
            # But the results should still be in the same ballpark ish
            # This test is more just to flag up if something goes completely crazy
            self.assertEqual(12, len(test_result))
            self.assertEqual(3, len(golden_result))
            print(test_result)
            print(golden_result)
            test_indices = [2, 4, 6]  # Only check Ca, Fe, Mg
            for i, g_result in enumerate(golden_result):
                el = list(test_result.items())[test_indices[i]][0]
                test_val = list(test_result.items())[test_indices[i]][1]
                self.assertTrue(abs(test_val - g_result) < 3, 'Mismatch at index ' + str(i) + ' of ' + arg_name + ' (Element: ' + str(el) + '). Golden: ' + str(g_result) + ', test: ' + str(test_val))
                
    def test_complete_model_earthlike(self):
        manager = mn.Manager(
            Namespace(
                wd_data_filename='BlouinConglomNewTimescales.csv',
                stellar_compositions_filename='StellarCompositionsSortFE.csv',
                n_live_points = 0, # This argument shouldn't matter
                pollution_model_names=['Model_24'],
                enhancement_model='Earthlike',
                base_dir=get_path_to_test_output()
            )
        )
        manager.publish_live_data(0)
        args = {
            # Args are: fe_star, t_sinceaccretion, d_formation, z_formation, N_c, N_o, f_c, f_o, pollutionfraction, t_disc
           # These 4 input sets are variations on a theme. 3 of them change t_sinceaccretion to test all 3 logic paths in white_dwarf_model.
           # One of them decreases d_formation to cause variation in the hsc chemistry outputs (i.e. not all equal to 1)
            'SDSSJ0002+3209': [88, 0, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768],
            'SDSSJ0002+3209': [88, 0.0195758, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768],
            'SDSSJ0002+3209': [88, 0.195758, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768],
            'SDSSJ0002+3209': [88, 0.195758, -0.8, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768]
        }
        for arg_name, arg_set in args.items():
            print('CompleteModelTestsEarthlike ' + arg_name)
            test_result, ignore = cm.complete_model_calculation(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], 10**(arg_set[9]), 54, -2, 'Earthlike', 1.5, False)
            golden_result = ocm.PWDCodeMultiple(0, 1, False, arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], arg_set[9])
            #The Earthlike version should in theory replicate John's original code
            self.assertEqual(12, len(test_result))
            self.assertEqual(3, len(golden_result))
            print(test_result)
            print(golden_result)
            test_indices = [2, 4, 6]  # Only check Ca, Fe, Mg
            for i, g_result in enumerate(golden_result):
                el = list(test_result.items())[test_indices[i]][0]
                test_val = list(test_result.items())[test_indices[i]][1]
                self.assertEqual(test_val, g_result)

class LoglikeTests(unittest.TestCase):
    
    def setUp(self):
        self.test_args = Namespace(
            wd_data_filename='BlouinConglomNewTimescales.csv',
            stellar_compositions_filename='StellarCompositionsSortFE.csv',
            n_live_points = 20,
            pollution_model_names=['Model_Full_No_Crust'],
            enhancement_model='NonEarthlike',
            base_dir=get_path_to_test_output()
        )
    
    def test_bounds(self):
        manager = mn.Manager(self.test_args)
        manager.publish_live_data(226)  # GD61 has various upper bounds
        manager.publish_live_model('Model_Full_No_Crust')
        dummy_cube = [470, 1, 1.5, 0.02, 0.05, -7, 1, 54, -2]
        likelihood = manager.models['Model_Full_No_Crust'].loglike(dummy_cube)
        self.assertNotEqual(likelihood, -9e89)  # This should not have triggered the bounds
        dummy_cube = [470, 1, 1.5, 0.02, 0.05, -4.5, 1, 54, -2]
        likelihood = manager.models['Model_Full_No_Crust'].loglike(dummy_cube)
        self.assertEqual(likelihood, -9e89)  # This should have triggered the Al bound (we just raised the pollution fraction by 2.5 orders of magnitude, putting Al above its upper bound)

#@unittest.skip("Skip for now")
class IntegrationTests(unittest.TestCase):
    
    def setUp(self):
        self.real_wd_file = 'BlouinConglomNewTimescales.csv'
        self.real_stellar_comps_file = 'StellarCompositionsSortFE.csv'
        self.n_live_points = 20
        self.seed = 1
        self.test_args_golden = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points = self.n_live_points,
            pollution_model_names=['Model_24'],
            enhancement_model='NonEarthlike',
            base_dir=get_path_to_test_output()
        )
        self.test_args_golden_earthlike = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points = self.n_live_points,
            pollution_model_names=['Model_24'],
            enhancement_model='Earthlike',
            base_dir=get_path_to_test_output(),
            seed=self.seed
        )
        self.test_args_hierarchy = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points = self.n_live_points,
            pollution_model_names=['Hierarchy_Default'],
            enhancement_model='NonEarthlike',
            base_dir=get_path_to_test_output(),
            seed=self.seed
        )
        self.start_from_scratch = True  # For a real test, this should be True, but to save time switch this to False

    def reset_output(self):
        # Remove previous results
        if self.start_from_scratch:
            if os.path.exists(get_path_to_test_output()) and os.path.isdir(get_path_to_test_output()):
                shutil.rmtree(get_path_to_test_output())
            
            if os.path.exists(get_path_to_chains()) and os.path.isdir(get_path_to_chains()):
                shutil.rmtree(get_path_to_chains())
    
    def assert_equal_within_errors(self, value_1, error_1, value_2, error_2):
        err_msg = str(value_1) + ' not equal to ' + str(value_2) + ' within errors (' + str(error_1) + ', ' + str(error_2) + ')'
        if value_1 >= value_2:
            self.assertGreaterEqual(value_2 + error_2, value_1 - error_1, err_msg)
        else:
            self.assertGreaterEqual(value_1 + error_1, value_2 - error_2, err_msg)
    
    def test_integration(self):
        self.reset_output()
        test_system = 0
        original_result = oc.PWDCode24(test_system, test_system+1, False, self.n_live_points, self.seed)
        
        manager_earthlike = mn.Manager(self.test_args_golden_earthlike)
        manager_earthlike.run([test_system])
        # This  should  match the original model exactly (or at least within machine precision) as long as the seed is consistent
        print('test_integration results 1:')
        print(original_result['logZ'])
        print(original_result['logZerr'])
        print(manager_earthlike.models['Model_24'].result['logZ'])
        print(manager_earthlike.models['Model_24'].result['logZerr'])
        
        manager = mn.Manager(self.test_args_golden)  # FWIW This model doesn't make much sense. For NonEarthlike, the parent core/crust values are just ignored
        manager.run([test_system])
        # Now we're comparing apples to oranges at this point because the models diverged
        # This test is just to flag up anything utterly crazy and/or broken
        # Errors tripled so it almost always passes...
        # Could seed this so that it's 100% guaranteed to pass but I'd like to have some non-seeded test somewhere
        print('test_integration results 2:')
        print(manager.models['Model_24'].result['logZ'])
        print(manager.models['Model_24'].result['logZerr'])
        self.assert_equal_within_errors(
            original_result['logZ'],
            3*original_result['logZerr'],
            manager.models['Model_24'].result['logZ'],
            3*manager.models['Model_24'].result['logZerr']
        )
		
        manager_h = mn.Manager(self.test_args_hierarchy)
        manager_h.run([203]) # This one (PG0843+516) should really prefer to be differentiated, at least
        self.assertEqual(manager_h.final_hierarchy, [0, 1, 3])  # Depending on the seed, 2 might be present!
        self.assertTrue(manager_h.models['Hierarchy_Default_levels_013'].best_model)
        #self.assertTrue('Hierarchy_Default_levels_013_HP' in manager_h.models.keys())  # Could reinstate this we switch HP/LP comparison priors back on
        
        # Moving this to the end because it's the fiddliest test.
        # This should pass, but only if certain parts of the code are changed to mimic the way the original model interprets likelihoods
        # The line log_zero = self.minimum_likelihood needs to be commented out in pollution_model.py
        # In loglike_functions.py, the lines affecting min_likelihood labelled 'For testing purposes only'
        # should be uncommented and the relevant neighbouring line should be commented (there are 8 of these pairs)
        # (These changes will make the LoglikeTests fail)
        # Then remember to revert these changes before committing!!!
        #required_dp = 13
        #self.assertAlmostEqual(original_result['logZ'], manager_earthlike.models['Model_24'].result['logZ'], required_dp)
        #self.assertAlmostEqual(original_result['logZerr'], manager_earthlike.models['Model_24'].result['logZerr'], required_dp)

if __name__ == '__main__':
    unittest.main()
