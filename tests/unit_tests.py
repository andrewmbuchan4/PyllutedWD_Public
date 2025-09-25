#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np
import os
import unittest
import shutil
import sys

import scipy
import scipy.stats as st

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

def get_path_to_src():
    return get_path_to_parent() + '/src/'

def get_path_to_original_src():
    return get_path_to_parent() + '/original_codebase/'

def get_path_to_data():
    return get_path_to_parent() + '/data/'

def get_path_to_feni():
    return get_path_to_parent() + '/feni_src/feni/'

def get_path_to_chains():
    return get_path_to_tests() + '/chains/'

sys.path.append(get_path_to_src())
sys.path.append(get_path_to_original_src())

import abundance_model as am
import atmosphere_model as atm
import chemistry_info as ci
import complete_model as cm
import disc_model as dm
import enhancement_model as em
import excess_oxygen_calculator as eoc
import geology_info as gi
import live_data as ld
import manager as mn
import model_analyser as ma
import model_parameters as mp
import multivariate_tests as mv
import original_complete_model as ocm
import original_enhancement_model as oe
import original_pressure_model as op
import partition_model as pam
import physical_constants as pc
import pwd_utils as pu
import rpy2.robjects as robjects
import solar_abundances as sa
import synthetic_bandpass as sb
import synthetic_configurations as sc
import synthetic_modeller as sm
import synthetic_population as sp
import synthetic_pipeline as spi
import synthetic_observer as so
import thermohaline_interpolator as thi
import timescale_interpolator as ti
import white_dwarf as wd
import white_dwarf_model as wdm

class TimescaleInterpolatorTests(unittest.TestCase):

    def test_init(self):
        test_ti = ti.TimescaleInterpolator()
        self.assertEqual(6.392, test_ti.all_timescale_data[ti.TimescaleType.KoesterOvershoot][ci.Element.He][8][3000][-7][ci.Element.C])
        self.assertEqual(6.32, test_ti.all_timescale_data[ti.TimescaleType.KoesterOvershoot][ci.Element.He][8][3000][-7][ci.Element.N])
        self.assertEqual(-5.533, test_ti.all_timescale_data[ti.TimescaleType.KoesterOvershoot][ci.Element.He][8][3000][-7]['logq'])
        self.assertEqual(6.464, test_ti.all_timescale_data[ti.TimescaleType.KoesterOvershoot][ci.Element.He][8][3000][-7.5][ci.Element.C])
        self.assertEqual(6.389, test_ti.all_timescale_data[ti.TimescaleType.KoesterOvershoot][ci.Element.He][8][3250][-7.5][ci.Element.C])
        self.assertEqual(5.021, test_ti.all_timescale_data[ti.TimescaleType.KoesterOvershoot][ci.Element.He][8.5][3250][-7.5][ci.Element.C])
        self.assertEqual(5.98, test_ti.all_timescale_data[ti.TimescaleType.KoesterNoOvershoot][ci.Element.He][8][3250][-7.5][ci.Element.C])
        with self.assertRaises(KeyError):
            not_real = test_ti.all_timescale_data[ti.TimescaleType.KoesterOvershoot][ci.Element.He][8][3001][-7][ci.Element.C]
        with self.assertRaises(KeyError):
            not_real = test_ti.all_timescale_data[ti.TimescaleType.KoesterOvershoot][ci.Element.H][8][3000][-7][ci.Element.C]  # -7 is meaningless to a H interpolator: CaHe is not a variable
        with self.assertRaises(KeyError):
            not_real = test_ti.all_timescale_data[ti.TimescaleType.KoesterNoOvershoot][ci.Element.He][8][3001][-7][ci.Element.C]
        with self.assertRaises(KeyError):
            not_real = test_ti.all_timescale_data[ti.TimescaleType.KoesterNoOvershoot][ci.Element.H][8][3000][-7][ci.Element.C]  # -7 is meaningless to a H interpolator: CaHe is not a variable
        self.assertEqual(None, test_ti.wd_data)
        self.assertEqual([
            ti.TimescaleType.KoesterOvershoot,
            ti.TimescaleType.KoesterNoOvershoot,
            ti.TimescaleType.BedardNoOvershoot,
            ti.TimescaleType.BedardOvershoot,
            ti.TimescaleType.BedardVariableOvershoot,
            ti.TimescaleType.Bedard3DOvershoot,
            ti.TimescaleType.MWDD,
            ti.TimescaleType.Bedard3DOvershootPatched
        ], list(test_ti.interpolators.keys()))
        self.assertEqual([ci.Element.H, ci.Element.He], list(test_ti.interpolators[ti.TimescaleType.KoesterOvershoot].keys()))
        self.assertEqual([ci.Element.H, ci.Element.He], list(test_ti.interpolators[ti.TimescaleType.KoesterNoOvershoot].keys()))
        self.assertEqual([ci.Element.H, ci.Element.He], list(test_ti.interpolators[ti.TimescaleType.BedardNoOvershoot].keys()))
        self.assertEqual(2, len(test_ti.interpolators[ti.TimescaleType.KoesterOvershoot][ci.Element.H]['logq'].grid))  # H is 2D (Teff, logg)
        self.assertEqual(3, len(test_ti.interpolators[ti.TimescaleType.KoesterOvershoot][ci.Element.He]['logq'].grid)) # He is 3D (Teff, logg, CaHe)
        self.assertEqual(2, len(test_ti.interpolators[ti.TimescaleType.KoesterNoOvershoot][ci.Element.H]['logq'].grid))  # H is 2D (Teff, logg)
        self.assertEqual(3, len(test_ti.interpolators[ti.TimescaleType.KoesterNoOvershoot][ci.Element.He]['logq'].grid)) # He is 3D (Teff, logg, CaHe)
        self.assertEqual(2, len(test_ti.interpolators[ti.TimescaleType.BedardNoOvershoot][ci.Element.H]['logq'].grid))  # Both 2D in this case
        self.assertEqual(2, len(test_ti.interpolators[ti.TimescaleType.BedardNoOvershoot][ci.Element.He]['logq'].grid))
        expected_els_da_koester = [
             'logq',
             ci.Element.He,
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
             ]
        expected_els_db_koester = [
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
             ]
        expected_els_da_bedard = [
             'logq',
             ci.Element.He,
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
             ci.Element.Cu
             ]
        expected_els_db_bedard = [
             'logq',
             ci.Element.H,
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
             ci.Element.Cu
             ]
        self.assertEqual(expected_els_da_koester, list(test_ti.interpolators[ti.TimescaleType.KoesterOvershoot][ci.Element.H].keys()))
        self.assertEqual(expected_els_db_koester, list(test_ti.interpolators[ti.TimescaleType.KoesterOvershoot][ci.Element.He].keys()))
        self.assertEqual(expected_els_da_koester, list(test_ti.interpolators[ti.TimescaleType.KoesterNoOvershoot][ci.Element.H].keys()))
        self.assertEqual(expected_els_db_koester, list(test_ti.interpolators[ti.TimescaleType.KoesterNoOvershoot][ci.Element.He].keys()))
        self.assertEqual(expected_els_da_bedard, list(test_ti.interpolators[ti.TimescaleType.BedardNoOvershoot][ci.Element.H].keys()))
        self.assertEqual(expected_els_db_bedard, list(test_ti.interpolators[ti.TimescaleType.BedardNoOvershoot][ci.Element.He].keys()))

    def test_interpolation(self):
        test_ti = ti.TimescaleInterpolator()
        wd_data = test_ti.load_wd_data()
        self.assertEqual({'CaHe': -9.1, 'Teff': 6466, 'Type': ci.Element.He, 'logg': 8.257}, wd_data['SDSSJ0002+3209'])
        self.assertEqual({'CaHe': -6.11, 'Teff': 20409, 'Type': ci.Element.H, 'logg': 7.95}, wd_data['WD1929+011'])
        expected_extraction1 = {
            ti.TimescaleType.KoesterNoOvershoot: None,
            ti.TimescaleType.KoesterOvershoot: None,
            ti.TimescaleType.BedardNoOvershoot: {
                'logq': -6.232334199999999,
                ci.Element.H: 5.732727200000001,
                ci.Element.Li: 6.257455600000004,
                ci.Element.Be: 6.164903800000006,
                ci.Element.B: 5.847232400000003,
                ci.Element.C: 5.879926399999999,
                ci.Element.N: 6.002531600000001,
                ci.Element.O: 6.034313400000002,
                ci.Element.F: 6.002959399999999,
                ci.Element.Ne: 6.0586920000000015,
                ci.Element.Na: 6.027485599999997,
                ci.Element.Mg: 6.0033434,
                ci.Element.Al: 5.900314600000005,
                ci.Element.Si: 5.780371400000005,
                ci.Element.P:  5.629068800000002,
                ci.Element.S:  5.658020600000001,
                ci.Element.Cl: 5.502982200000005,
                ci.Element.Ar: 5.526023799999999,
                ci.Element.K:  5.580110399999999,
                ci.Element.Ca: 5.647683600000001,
                ci.Element.Sc: 5.650120999999999,
                ci.Element.Ti: 5.663650400000002,
                ci.Element.V:  5.670551999999999,
                ci.Element.Cr: 5.680871800000002,
                ci.Element.Mn: 5.655886800000001,
                ci.Element.Fe: 5.660532800000004,
                ci.Element.Co: 5.616544599999999,
                ci.Element.Ni: 5.552737200000001,
                ci.Element.Cu: 5.445662399999998
            },
            ti.TimescaleType.BedardOvershoot: {
                'logq': -5.6906326,
                ci.Element.H: 6.221067200000002,
                ci.Element.Li: 6.8880693999999965,
                ci.Element.Be: 6.854940000000003,
                ci.Element.B: 6.7592555999999995,
                ci.Element.C: 6.468953200000001,
                ci.Element.N: 6.443493200000002,
                ci.Element.O: 6.5780556,
                ci.Element.F: 6.553332600000002,
                ci.Element.Ne: 6.6093112000000005,
                ci.Element.Na: 6.580157200000001,
                ci.Element.Mg: 6.619160399999998,
                ci.Element.Al: 6.543061000000004,
                ci.Element.Si: 6.517992600000005,
                ci.Element.P:  6.450583000000002,
                ci.Element.S:  6.293319000000003,
                ci.Element.Cl: 6.0812431999999985,
                ci.Element.Ar: 6.110725399999999,
                ci.Element.K:  6.0352624,
                ci.Element.Ca: 6.109430200000002,
                ci.Element.Sc: 6.083895399999999,
                ci.Element.Ti: 6.110914599999998,
                ci.Element.V:  6.146158400000003,
                ci.Element.Cr: 6.190768000000002,
                ci.Element.Mn: 6.193118800000005,
                ci.Element.Fe: 6.215168,
                ci.Element.Co: 6.181930600000001,
                ci.Element.Ni: 6.191247199999997,
                ci.Element.Cu: 6.158520999999998
            },
            ti.TimescaleType.BedardVariableOvershoot: {
                'logq': -5.939016599999998,
                ci.Element.H: 5.972106799999999,
                ci.Element.Li: 6.611198799999998,
                ci.Element.Be: 6.585248,
                ci.Element.B: 6.3012194000000035,
                ci.Element.C: 6.1153874,
                ci.Element.N: 6.293006400000002,
                ci.Element.O: 6.3362929999999995,
                ci.Element.F: 6.307791400000003,
                ci.Element.Ne: 6.359745400000002,
                ci.Element.Na: 6.336093000000002,
                ci.Element.Mg: 6.368794599999999,
                ci.Element.Al: 6.290168999999998,
                ci.Element.Si: 6.159237600000002,
                ci.Element.P: 5.946629200000002,
                ci.Element.S: 5.9568056,
                ci.Element.Cl: 5.872900799999997,
                ci.Element.Ar: 5.749484800000001,
                ci.Element.K: 5.8404953999999965,
                ci.Element.Ca: 5.8849884,
                ci.Element.Sc: 5.9005862,
                ci.Element.Ti: 5.926497400000001,
                ci.Element.V: 5.9347252,
                ci.Element.Cr: 5.963026799999996,
                ci.Element.Mn: 5.964051400000002,
                ci.Element.Fe: 5.983149800000002,
                ci.Element.Co: 5.9499124000000005,
                ci.Element.Ni: 5.9228407999999995,
                ci.Element.Cu: 5.831818400000003
            },
            ti.TimescaleType.MWDD: {
                'logq': -6.232334199999999,
                ci.Element.Li: 6.257455600000004,
                ci.Element.Be: 6.163928400000006,
                ci.Element.B: 5.848032400000001,
                ci.Element.C: 5.879926399999999,
                ci.Element.N: 6.002531600000001,
                ci.Element.O: 6.034313400000002,
                ci.Element.F: 6.003934799999997,
                ci.Element.Ne: 6.0586920000000015,
                ci.Element.Na: 6.027485599999997,
                ci.Element.Mg: 6.0043188,
                ci.Element.Al: 5.900314600000005,
                ci.Element.Si: 5.780371400000005,
                ci.Element.P:  5.630044200000002,
                ci.Element.S:  5.658020600000001,
                ci.Element.Cl: 5.502982200000005,
                ci.Element.Ar: 5.526023799999999,
                ci.Element.K:  5.580110399999999,
                ci.Element.Ca: 5.648659,
                ci.Element.Sc: 5.650120999999999,
                ci.Element.Ti: 5.663650400000002,
                ci.Element.V:  5.670551999999999,
                ci.Element.Cr: 5.677770200000001,
                ci.Element.Mn: 5.655886800000001,
                ci.Element.Fe: 5.660532800000004,
                ci.Element.Co: 5.616544599999999,
                ci.Element.Ni: 5.552737200000001,
                ci.Element.Cu: 5.445662399999998
            }
        }
        test_extraction1 = test_ti.extract_timescales(ci.Element.He, 8.1, 3123, None)
        self.assertTrue(expected_extraction1, test_extraction1)

        test_extraction = test_ti.extract_timescales(ci.Element.He, 8.1, 3123, -7.2)
        expected_test_extraction = {
            ti.TimescaleType.KoesterOvershoot: {
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
            },
            ti.TimescaleType.KoesterNoOvershoot: {
                'logq': -6.261410879999999,
                ci.Element.Li: 5.79920672,
                ci.Element.Be: 5.723501920000001,
                ci.Element.B: 5.681917760000001,
                ci.Element.C: 5.692329920000001,
                ci.Element.N: 5.621113920000001, # I have checked this value by hand from the raw data
                ci.Element.O: 5.5603017600000015,
                ci.Element.F: 5.457495840000001,
                ci.Element.Ne: 5.450973920000002,
                ci.Element.Na: 5.3797979200000015,
                ci.Element.Mg: 5.367157920000001,
                ci.Element.Al: 5.313067040000001,
                ci.Element.Si: 5.30658192,
                ci.Element.P: 5.255251040000001,
                ci.Element.S: 5.248729120000002,
                ci.Element.Cl: 5.194461120000001,
                ci.Element.Ar: 5.129389920000001,
                ci.Element.K: 5.155675040000001,
                ci.Element.Ca: 5.151353760000001,
                ci.Element.Sc: 5.087960000000001,
                ci.Element.Ti: 5.058255200000001,
                ci.Element.V: 5.029333920000002,
                ci.Element.Cr: 5.024747200000001,
                ci.Element.Mn: 4.999665920000001,
                ci.Element.Fe: 4.996843040000001,
                ci.Element.Co: 4.971721120000001,
                ci.Element.Ni: 4.979757920000001,
                ci.Element.Cu: 4.939876000000002,
                ci.Element.Zn: 4.929053120000001
            },
            ti.TimescaleType.BedardNoOvershoot: {
                'logq': -6.232334199999999,
                ci.Element.H: 5.732727200000001,
                ci.Element.Li: 6.257455600000004,
                ci.Element.Be: 6.164903800000006,
                ci.Element.B: 5.847232400000003,
                ci.Element.C: 5.879926399999999,
                ci.Element.N: 6.002531600000001,
                ci.Element.O: 6.034313400000002,
                ci.Element.F: 6.002959399999999,
                ci.Element.Ne: 6.0586920000000015,
                ci.Element.Na: 6.027485599999997,
                ci.Element.Mg: 6.0033434,
                ci.Element.Al: 5.900314600000005,
                ci.Element.Si: 5.780371400000005,
                ci.Element.P:  5.629068800000002,
                ci.Element.S:  5.658020600000001,
                ci.Element.Cl: 5.502982200000005,
                ci.Element.Ar: 5.526023799999999,
                ci.Element.K:  5.580110399999999,
                ci.Element.Ca: 5.647683600000001,
                ci.Element.Sc: 5.650120999999999,
                ci.Element.Ti: 5.663650400000002,
                ci.Element.V:  5.670551999999999,
                ci.Element.Cr: 5.680871800000002,
                ci.Element.Mn: 5.655886800000001,
                ci.Element.Fe: 5.660532800000004,
                ci.Element.Co: 5.616544599999999,
                ci.Element.Ni: 5.552737200000001,
                ci.Element.Cu: 5.445662399999998
            },
            ti.TimescaleType.BedardOvershoot: {
                'logq': -5.6906326,
                ci.Element.H: 6.221067200000002,
                ci.Element.Li: 6.8880693999999965,
                ci.Element.Be: 6.854940000000003,
                ci.Element.B: 6.7592555999999995,
                ci.Element.C: 6.468953200000001,
                ci.Element.N: 6.443493200000002,
                ci.Element.O: 6.5780556,
                ci.Element.F: 6.553332600000002,
                ci.Element.Ne: 6.6093112000000005,
                ci.Element.Na: 6.580157200000001,
                ci.Element.Mg: 6.619160399999998,
                ci.Element.Al: 6.543061000000004,
                ci.Element.Si: 6.517992600000005,
                ci.Element.P: 6.450583000000002,
                ci.Element.S: 6.293319000000003,
                ci.Element.Cl: 6.0812431999999985,
                ci.Element.Ar: 6.110725399999999,
                ci.Element.K: 6.0352624,
                ci.Element.Ca: 6.109430200000002,
                ci.Element.Sc: 6.083895399999999,
                ci.Element.Ti: 6.110914599999998,
                ci.Element.V: 6.146158400000003,
                ci.Element.Cr: 6.190768000000002,
                ci.Element.Mn: 6.193118800000005,
                ci.Element.Fe: 6.215168,
                ci.Element.Co: 6.181930600000001,
                ci.Element.Ni: 6.191247199999997,
                ci.Element.Cu: 6.158520999999998
            },
            ti.TimescaleType.BedardVariableOvershoot: {
                'logq': -5.939016599999998,
                ci.Element.H: 5.972106799999999,
                ci.Element.Li: 6.611198799999998,
                ci.Element.Be: 6.585248,
                ci.Element.B: 6.3012194000000035,
                ci.Element.C: 6.1153874,
                ci.Element.N: 6.293006400000002,
                ci.Element.O: 6.3362929999999995,
                ci.Element.F: 6.307791400000003,
                ci.Element.Ne: 6.359745400000002,
                ci.Element.Na: 6.336093000000002,
                ci.Element.Mg: 6.368794599999999,
                ci.Element.Al: 6.290168999999998,
                ci.Element.Si: 6.159237600000002,
                ci.Element.P: 5.946629200000002,
                ci.Element.S: 5.9568056,
                ci.Element.Cl: 5.872900799999997,
                ci.Element.Ar: 5.749484800000001,
                ci.Element.K: 5.8404953999999965,
                ci.Element.Ca: 5.8849884,
                ci.Element.Sc: 5.9005862,
                ci.Element.Ti: 5.926497400000001,
                ci.Element.V: 5.9347252,
                ci.Element.Cr: 5.963026799999996,
                ci.Element.Mn: 5.964051400000002,
                ci.Element.Fe: 5.983149800000002,
                ci.Element.Co: 5.9499124000000005,
                ci.Element.Ni: 5.9228407999999995,
                ci.Element.Cu: 5.831818400000003
            },
            ti.TimescaleType.Bedard3DOvershoot: None,
            ti.TimescaleType.MWDD: {
                'logq': -6.232334199999999,
                ci.Element.Li: 6.257455600000004,
                ci.Element.Be: 6.163928400000006,
                ci.Element.B: 5.848032400000001,
                ci.Element.C: 5.879926399999999,
                ci.Element.N: 6.002531600000001,
                ci.Element.O: 6.034313400000002,
                ci.Element.F: 6.003934799999997,
                ci.Element.Ne: 6.0586920000000015,
                ci.Element.Na: 6.027485599999997,
                ci.Element.Mg: 6.0043188,
                ci.Element.Al: 5.900314600000005,
                ci.Element.Si: 5.780371400000005,
                ci.Element.P:  5.630044200000002,
                ci.Element.S:  5.658020600000001,
                ci.Element.Cl: 5.502982200000005,
                ci.Element.Ar: 5.526023799999999,
                ci.Element.K:  5.580110399999999,
                ci.Element.Ca: 5.648659,
                ci.Element.Sc: 5.650120999999999,
                ci.Element.Ti: 5.663650400000002,
                ci.Element.V:  5.670551999999999,
                ci.Element.Cr: 5.677770200000001,
                ci.Element.Mn: 5.655886800000001,
                ci.Element.Fe: 5.660532800000004,
                ci.Element.Co: 5.616544599999999,
                ci.Element.Ni: 5.552737200000001,
                ci.Element.Cu: 5.445662399999998
            },
            ti.TimescaleType.Bedard3DOvershootPatched: None
        }
        self.assertEqual(expected_test_extraction, test_extraction)
        test_H_extraction = test_ti.extract_timescales(ci.Element.H, 8.1, 3123)
        expected_test_H_extraction = {
            'logq': -4.338032799999992,
            ci.Element.He: 7.0982055999999965, # I have checked this value by hand from the raw data
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
        self.assertEqual(expected_test_H_extraction, test_H_extraction[ti.TimescaleType.KoesterOvershoot])
        test_extrapolation = test_ti.extract_timescales(ci.Element.He, 9.1, 6123, -6.2)
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
        self.assertEqual(expected_extrapolation, test_extrapolation[ti.TimescaleType.KoesterOvershoot])
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
        self.assertEqual(expected_wd_timescales_dict, wd_timescales_dict[ti.TimescaleType.KoesterOvershoot])
        wd_timescales_list = test_ti.return_wd_timescales_as_list(test_extraction[ti.TimescaleType.KoesterOvershoot])
        self.assertEqual([-5.827288959999999, 521406.1627526086, 286064.1007553645, 359160.2953306345, 238996.48847621836, 248059.2121373417, 264814.30405580055, 593456.5943745527, 516290.6261463279, 606890.3612551007, 931558.7297104555, 1261900.5219890976, 1071292.7356341446], wd_timescales_list)
        self.assertEqual(wd_timescales_dict, test_ti.get_wd_timescales(ci.Element.He, 8.1, 3123, -7.2))
        processed_wd_data = test_ti.process_wd_data()
        self.assertEqual(len(processed_wd_data), len(wd_data))
        test_real_extraction = test_ti.extract_timescales(ci.Element.He, 8.257, 6466, -9.1)
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
        self.assertEqual(expected_real_timescales, test_real_extraction[ti.TimescaleType.KoesterOvershoot])
        self.assertEqual(expected_real_timescales, processed_wd_data['SDSSJ0002+3209'][ti.TimescaleType.KoesterOvershoot])

class AtmosphereModelTests(unittest.TestCase):

    def assertListsAlmostEqual(self, list1, list2):
        # This function is useful across the board - should make it available at a higher level. Should also add useful output for if the entries aren't equal
        self.assertEqual(len(list1), len(list2))
        for i, val1 in enumerate(list1):
            if np.isnan(val1):
                self.assertTrue(np.isnan(list2[i]))
            else:
                self.assertAlmostEqual(val1, list2[i])

    def test_number_abundance_calculation_non_thermohaline(self):
        argument_sets = [
            # Arguments are: t, t_event, Hx, M_cvz, mu_X, M_X, tau_X, expected_N_X
            [100, 90, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([-0.8792852115115177, -0.6845351853425914, -2.1540777613296918])],
            [0, 90, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([-np.inf, -np.inf, -np.inf])],
            [1000, 10000, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([-2.7604502791595324, -2.335358024443877, -3.9421818055678997])],
            [2000, 10000, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([-2.7604502791595253, -2.335358024443877, -3.9421818055679028])], # This is a sanity check - it's practically the same as previous case because it reached steady state
            [1500, 1000, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([-8.765199987276498, -21.076016292773502, -15.715548920369425])],
            [10500, 10000, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([-9.765199987276489, -22.076016292773474, -16.715548920369407])], # This should match the previous one because they're both 500 yr into declining phase, but one order of magnitude less because we spread the same mass over a 10 times longer accretion event, so there was 1/10 the amount of material in atmosphere at start of dec phase
            [100000000000, 10000, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([-np.inf, -np.inf, -np.inf])],
            [100, 90, ci.Element.He, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([-0.28040333557866465, -0.08565330940973795, -1.5551958853968382])], # 4 times higher than the first case bc He is 4 times heavier than H, so there are 4 times less He atoms
            [100, 90, ci.Element.H, 0, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([np.nan, np.nan, np.nan])],
            [100, 90, ci.Element.H, 300, np.array([0, 8000, 10.9]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), np.array([np.nan, -3.6845351853425914, -2.191504259270315])],
            [100, 90, ci.Element.H, 300, np.array([6, 8, 10]), np.array([-200, 10000, 0]), np.array([31, 11, 17]), np.array([np.nan, -0.6845351853425914, -np.inf])],
            [100, 90, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([0, 100000000, 10]), np.array([np.nan, 0.6232490514919188, -2.561221317793099])]
        ]

        for i, arg_set in enumerate(argument_sets):
            print('AtmosphereModelTests: Testing arg set ' + str(i))
            N_X = atm.calculate_abundance_by_number(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6])
            expected_N_X = arg_set[7]
            self.assertListsAlmostEqual(expected_N_X.tolist(), N_X.tolist())

    def test_number_abundance_calculation_thermohaline(self):
        argument_sets = [
            # Arguments are: t, t_event, Hx, M_cvz, mu_X, M_X, tau_X, consider_thermohaline, Teff, logg, expected_N_X
            [100, 90, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1000, 10000, 200]), np.array([31, 11, 17]), False, 10000, 8, np.array([-0.8792852115115177, -0.6845351853425914, -2.1540777613296918])],
            [100, 90, ci.Element.H, 300, np.array([6, 8, 10]), np.array([1, 10, 0.2]), np.array([31, 11, 17]), False, 10000, 8, np.array([-3.8792852115115177, -3.6845351853425914, -5.1540777613296918])],
            [90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), False, 10000, 8, np.array([-4.262068962629518, -4.6876620286480355, -4.500670410361624])],
            [90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([-5.138124625324908, -5.138124625324908, -5.138124625324908])], # Checking that all abundances are equal: differential sinking not active in this case
            [90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 6000, 8, np.array([-4.261330505769083, -4.6869235717876006, -4.499931953501189])], # Now it is active again since this is cold enough to not trigger strong thermohaline mixing
            [90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 30000, 8, np.array([-6.450575499257915, -6.450575499257915, -6.450575499257915])], # Differential sinking not active again. Abundances lower this time (more extreme thermohaline effect)...
            [90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 300000, 8, np.array([-6.450575499257915, -6.450575499257915, -6.450575499257915])], # ...but the effect caps out at high Teff
            [90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8.1, np.array([-5.0717902980865475, -5.0717902980865475, -5.0717902980865475])],
            [100, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8.1, np.array([-5.211885292248887, -5.46660346345314, -5.3272576403825775])],
            [90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 1000, 8, np.array([-4.262068962629518, -4.6876620286480355, -4.500670410361624])],
            [90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 100, 8, np.array([-4.262068962629518, -4.6876620286480355, -4.500670410361624])], # Should cap out at low Teff too
            [80, 100, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-6*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([4.889741275678141, 4.889741275678141, 4.889741275678141])],
            [80, 100, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-5*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([5.889741275678141, 5.889741275678141, 5.889741275678141])],  # Thermohaline effects cap at high Mdot so it just scales linearly
            [80, 100, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-22*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([-9.317408919886569, -9.733480671129392, -9.54806755834558])], # Test the low Mdot interpolation
            [80, 100, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-27*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([-14.317527880791054, -14.733599632033878, -14.548186519250065])],
            [80, 100, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-28*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([-15.317527880791054, -15.733599632033878, -15.548186519250065])], # Should eventually cap out too
            [110, 100, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-24*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 100000, 8, np.array([-11.50948919722806, -12.19662151392834, -11.86938311170987])], # Declining phase, lowMdot interpolation and max Teff cap out, all in one call
            [110, 100, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-4*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 100000, 8, np.array([5.172666137988294, 4.917947966784041, 5.0572937898546035])], # Declining phase, max Mdot cap out and max Teff cap out, all in one call
            [0, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([-np.inf, -np.inf, -np.inf])],
            [0.00000001, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([-14.418485214181391, -14.418485214181391, -14.418485214181391])],
            [90, 90, ci.Element.He, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([-3.663187086696664, -4.0887801527151835, -3.901788534428771])], # No thermohaline mixing for He
            [90, 90, ci.Element.Ru, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 6000, 8, np.array([-2.26090722903824, -2.68650029505676, -2.4995086767703474])], # The logic still works for crazy atmospheres - Ru is 100x heavier than H, so relative number abundances are 100x higher!
            [90, -90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([np.nan, np.nan, np.nan])],
            [-90, 90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([np.nan, np.nan, np.nan])],
            [-90, -90, ci.Element.H, 1E-14, np.array([6, 6, 6]), 1E-17*np.array([1, 1, 1]), np.array([31, 11, 17]), True, 10000, 8, np.array([np.nan, np.nan, np.nan])]
        ]

        for i, arg_set in enumerate(argument_sets):
            print('AtmosphereModelTests (therm): Testing arg set ' + str(i))
            N_X = atm.calculate_abundance_by_number(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], arg_set[9])
            expected_N_X = arg_set[10]
            print('Comparing the following lists')
            print(expected_N_X.tolist())
            print(N_X.tolist())
            self.assertListsAlmostEqual(expected_N_X.tolist(), N_X.tolist())

class ThermohalineInterpolatorTests(unittest.TestCase):

    def test_init(self):
        therm_interpolator = thi.ThermohalineInterpolator()
        self.assertEqual(dict(), therm_interpolator.accretion_rate_limits_dict)
        self.assertEqual(dict(), therm_interpolator.temperature_limits_dict)
        self.assertTrue(therm_interpolator.interpolate)

    def test_call(self):
        therm_interpolator = thi.ThermohalineInterpolator()
        self.assertEqual(dict(), therm_interpolator.accretion_rate_limits_dict)
        self.assertEqual(dict(), therm_interpolator.temperature_limits_dict)
        therm_factor = therm_interpolator((6, 10000, 8))
        self.assertAlmostEqual(np.array(-0.3746319), therm_factor)
        self.assertEqual(dict(), therm_interpolator.accretion_rate_limits_dict)
        self.assertEqual(dict(), therm_interpolator.temperature_limits_dict)
        therm_factor = therm_interpolator((20, 70000, 8.2))
        self.assertAlmostEqual(np.array(-3.07815201), therm_factor)
        self.assertEqual({(20499.840981016227, 8.2): (3.799526635787294, 11.000000000001494)}, therm_interpolator.accretion_rate_limits_dict)
        self.assertEqual({8.2: 20499.840981016227}, therm_interpolator.temperature_limits_dict)
        therm_factor = therm_interpolator((2, 9000, 8))
        self.assertAlmostEqual(np.array(-0.0015027371331392955), therm_factor)
        self.assertEqual({(20499.840981016227, 8.2): (3.799526635787294, 11.000000000001494), (9000, 8): (3.799526635787294, 11.799526635787338)}, therm_interpolator.accretion_rate_limits_dict)
        self.assertEqual({8.2: 20499.840981016227, 8: 20501.4456482059}, therm_interpolator.temperature_limits_dict)
        therm_factor = therm_interpolator((-2, 9000, 8))
        self.assertAlmostEqual(np.array(0), therm_factor)
        self.assertEqual({(20499.840981016227, 8.2): (3.799526635787294, 11.000000000001494), (9000, 8): (3.799526635787294, 11.799526635787338)}, therm_interpolator.accretion_rate_limits_dict)
        self.assertEqual({8.2: 20499.840981016227, 8: 20501.4456482059}, therm_interpolator.temperature_limits_dict)

    def test_find_accretion_rate_limits(self):
        therm_interpolator = thi.ThermohalineInterpolator()
        self.assertEqual(dict(), therm_interpolator.accretion_rate_limits_dict)
        therm_interpolator.find_accretion_rate_limits(12000, 8.1)
        self.assertEqual({(12000, 8.1): (3.799526635787305, 11.799526635787327)}, therm_interpolator.accretion_rate_limits_dict)
        therm_interpolator.find_accretion_rate_limits(70000, 8.1)
        self.assertEqual({(12000, 8.1): (3.799526635787305, 11.799526635787327), (70000, 8.1): (np.nan, np.nan)}, therm_interpolator.accretion_rate_limits_dict)
        therm_interpolator.find_accretion_rate_limits(1000, 8.1)
        self.assertEqual({(12000, 8.1): (3.799526635787305, 11.799526635787327), (70000, 8.1): (np.nan, np.nan), (1000, 8.1): (np.nan, np.nan)}, therm_interpolator.accretion_rate_limits_dict)
        therm_interpolator.find_accretion_rate_limits(70000, 10)
        self.assertEqual({(12000, 8.1): (3.799526635787305, 11.799526635787327), (70000, 8.1): (np.nan, np.nan), (1000, 8.1): (np.nan, np.nan), (70000, 10): (np.nan, np.nan)}, therm_interpolator.accretion_rate_limits_dict)

    def test_find_max_temperature(self):
        therm_interpolator = thi.ThermohalineInterpolator()
        self.assertEqual(dict(), therm_interpolator.temperature_limits_dict)
        therm_interpolator.find_max_temperature(8.1)
        self.assertEqual({8.1: 20500.643322442924}, therm_interpolator.temperature_limits_dict)
        therm_interpolator.find_max_temperature(8.2)
        self.assertEqual({8.1: 20500.643322442924, 8.2: 20499.840981016227}, therm_interpolator.temperature_limits_dict)
        therm_interpolator.find_max_temperature(6)
        self.assertEqual({8.1: 20500.643322442924, 8.2: 20499.840981016227, 6: np.nan}, therm_interpolator.temperature_limits_dict)

class WhiteDwarfDataPointTests(unittest.TestCase):

    def test_init(self):
        test_datapoint1 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.5, 0.2)
        self.assertEqual(test_datapoint1.data_point_type, wd.WhiteDwarfDataPointType.measurement)
        self.assertEqual(test_datapoint1.value, -6.5)
        self.assertEqual(test_datapoint1.upper_error, 0.2)
        self.assertEqual(test_datapoint1.lower_error, 0.2)
        self.assertEqual(test_datapoint1.included, True)

        test_datapoint2 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-6.25', '0.1', False)
        self.assertEqual(test_datapoint2.value, -6.25)
        self.assertEqual(test_datapoint2.upper_error, 0.1)
        self.assertEqual(test_datapoint2.lower_error, 0.1)
        self.assertEqual(test_datapoint2.included, False)

        with self.assertRaises(ValueError):
            test_datapoint3 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7)

        with self.assertRaises(TypeError):
            test_datapoint4 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, None)

        with self.assertRaises(ValueError):
            test_datapoint5 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 'Raises')

        test_datapoint6 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.6, 0.3, False)
        self.assertEqual(test_datapoint6.data_point_type, wd.WhiteDwarfDataPointType.upper_bound)
        self.assertEqual(test_datapoint6.value, -6.6)
        self.assertEqual(test_datapoint6.upper_error, None)
        self.assertEqual(test_datapoint6.lower_error, None)
        self.assertEqual(test_datapoint6.included, False)

        test_datapoint7 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-8')
        self.assertEqual(test_datapoint7.data_point_type, wd.WhiteDwarfDataPointType.lower_bound)
        self.assertEqual(test_datapoint7.value, -8)
        self.assertEqual(test_datapoint7.upper_error, None)
        self.assertEqual(test_datapoint7.lower_error, None)
        self.assertEqual(test_datapoint7.included, True)

        test_datapoint8 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 6.8, 0.15, True, 0.18)
        self.assertEqual(test_datapoint8.data_point_type, wd.WhiteDwarfDataPointType.measurement)
        self.assertEqual(test_datapoint8.value, 6.8)
        self.assertEqual(test_datapoint8.upper_error, 0.15)
        self.assertEqual(test_datapoint8.lower_error, 0.18)
        self.assertEqual(test_datapoint8.included, True)

    def test_str(self):
        test_datapoint1 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-6.5', 0.2)
        self.assertEqual(str(test_datapoint1), '-6.5±0.2')

        test_datapoint2 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.6, 0.3, False)
        self.assertEqual(str(test_datapoint2), '<-6.6')

        test_datapoint3 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-8')
        self.assertEqual(str(test_datapoint3), '>-8.0')

        test_datapoint4 = wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.8, 0.15, True, 0.18)
        self.assertEqual(str(test_datapoint4), '-6.8^{+0.15}_{-0.18}')

class WhiteDwarfAbundanceDataTests(unittest.TestCase):

    def test_get_measurements_as_lists(self):
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Cu: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        el_list, value_list, upper_error_list, lower_error_list = abundance_data.get_measurements_as_lists(False, False)
        self.assertEqual(el_list, [])
        self.assertEqual(value_list, [])
        self.assertEqual(upper_error_list, [])
        self.assertEqual(lower_error_list, [])

        el_list, value_list, upper_error_list, lower_error_list = abundance_data.get_measurements_as_lists(False, True)
        self.assertEqual(el_list, [ci.Element.Cu])
        self.assertEqual(value_list, [-7.77])
        self.assertEqual(upper_error_list, [0.11])
        self.assertEqual(lower_error_list, [0.11])

        el_list, value_list, upper_error_list, lower_error_list = abundance_data.get_measurements_as_lists()
        self.assertEqual(el_list, [ci.Element.Li, ci.Element.Ag, ci.Element.Ds])
        self.assertEqual(value_list, [-4.6, -7.9, -7.6])
        self.assertEqual(upper_error_list, [0.6, 0.31, 0.4])
        self.assertEqual(lower_error_list, [0.6, 0.33, 0.4])

        el_list, value_list, upper_error_list, lower_error_list = abundance_data.get_measurements_as_lists(True, True)
        self.assertEqual(el_list, [ci.Element.Li, ci.Element.Cu, ci.Element.Ag, ci.Element.Ds])
        self.assertEqual(value_list, [-4.6, -7.77, -7.9, -7.6])
        self.assertEqual(upper_error_list, [0.6, 0.11, 0.31, 0.4])
        self.assertEqual(lower_error_list, [0.6, 0.11, 0.33, 0.4])

        el_list, value_list = abundance_data.get_upper_bounds_as_lists(False, False)
        self.assertEqual(el_list, [])
        self.assertEqual(value_list, [])

        el_list, value_list = abundance_data.get_upper_bounds_as_lists(False, True)
        self.assertEqual(el_list, [ci.Element.Fe])
        self.assertEqual(value_list, [-6.54])

        el_list, value_list = abundance_data.get_upper_bounds_as_lists()
        self.assertEqual(el_list, [ci.Element.C])
        self.assertEqual(value_list, [-6.5])

        el_list, value_list = abundance_data.get_upper_bounds_as_lists(True, True)
        self.assertEqual(el_list, [ci.Element.C, ci.Element.Fe])
        self.assertEqual(value_list, [-6.5, -6.54])

    def test_str(self):
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Ts: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        expected_str = "Li: -4.6±0.6\nC: <-6.5\nK: >-12.8 (not included)\nFe: <-6.54 (not included)\nAg: -7.9^{+0.31}_{-0.33}\nPo: >-12.5\nDs: -7.6±0.4\nTs: -7.77±0.11 (not included)\n"
        self.assertEqual(str(abundance_data), expected_str)

    def test_get_relative_value(self):
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Ts: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        self.assertAlmostEqual(abundance_data.get_relative_value(ci.Element.C, ci.Element.Li, False), -1.9)
        self.assertAlmostEqual(abundance_data.get_relative_value(ci.Element.C, ci.Element.Li, True), -1.9 - (8.46-0.96))

class WhiteDwarfTests(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(ValueError):
            white_dwarf1 = wd.WhiteDwarf('TestWD', wd.WhiteDwarfPropertyData({mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.Hf)}), wd.WhiteDwarfAbundanceData({ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')}))

        with self.assertRaises(ValueError):
            white_dwarf1 = wd.WhiteDwarf('TestWD', wd.WhiteDwarfPropertyData({mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'Error')}), wd.WhiteDwarfAbundanceData({ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')}))

        with self.assertRaises(ValueError):
            white_dwarf1 = wd.WhiteDwarf('TestWD', wd.WhiteDwarfPropertyData({mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'DA'), mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.He)}), wd.WhiteDwarfAbundanceData({ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')}))

        white_dwarf1 = wd.WhiteDwarf('TestWD', wd.WhiteDwarfPropertyData({mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.H)}), wd.WhiteDwarfAbundanceData({ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')}))
        self.assertEqual(white_dwarf1.get_spectral_type().value, 'DA')
        self.assertEqual(white_dwarf1.get_atmospheric_type().value, ci.Element.H)
        white_dwarf1 = wd.WhiteDwarf('TestWD', wd.WhiteDwarfPropertyData({mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'DB')}), wd.WhiteDwarfAbundanceData({ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')}))
        self.assertEqual(white_dwarf1.get_spectral_type().value, 'DB')
        self.assertEqual(white_dwarf1.get_atmospheric_type().value, ci.Element.He)

    def test_str(self):
        wd_property_data_raw1 = {
            mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.H),
            mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -1123, 0),
            mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.123, 0)
        }
        property_data1 = wd.WhiteDwarfPropertyData(wd_property_data_raw1)
        wd_property_data_raw2 = {
            mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'DB'),
            mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -1123, 0),
            mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.123, 0),
            mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 1000, 0),
            mp.WDParameter.distance: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 10000, 0)
        }
        property_data2 = wd.WhiteDwarfPropertyData(wd_property_data_raw2)
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Ts: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, 6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        white_dwarf1 = wd.WhiteDwarf('TestWD1', property_data1, abundance_data)
        expected_str1 = "\n--- TestWD1 ---\n\nSpectral Type: DA\nAtmospheric Type: H\nTemperature: -1123.0±0.0 K\nLogg: 8.123±0.0\nMass: ---\nDistance: ---\n\nH: 0.0±0.0 (not included)\nLi: -4.6±0.6\nC: <-6.5\nK: >-12.8 (not included)\nFe: <-6.54 (not included)\nAg: -7.9^{+0.31}_{-0.33}\nPo: >-12.5\nDs: -7.6±0.4\nTs: -7.77±0.11 (not included)\n\n\n"
        self.assertEqual(str(white_dwarf1), expected_str1)

        white_dwarf2 = wd.WhiteDwarf('TestWD2', property_data2, abundance_data)
        expected_str2 = "\n--- TestWD2 ---\n\nSpectral Type: DB\nAtmospheric Type: He\nTemperature: -1123.0±0.0 K\nLogg: 8.123±0.0\nMass: 1000.0±0.0 M_{Solar}\nDistance: 10000.0±0.0 pc\n\nHe: 0.0±0.0 (not included)\nLi: -4.6±0.6\nC: <-6.5\nK: >-12.8 (not included)\nFe: <-6.54 (not included)\nAg: -7.9^{+0.31}_{-0.33}\nPo: >-12.5\nDs: -7.6±0.4\nTs: -7.77±0.11 (not included)\n\n\n"
        self.assertEqual(str(white_dwarf2), expected_str2)

        atm_data = wd.WhiteDwarfAtmosphereDataset({
            ti.TimescaleType.KoesterNoOvershoot: wd.WhiteDwarfAtmosphereData({
                'logq': -14.22,
                ci.Element.Pm: 345.6,
                ci.Element.Sn: 987.5
            }),
            ti.TimescaleType.KoesterOvershoot: wd.WhiteDwarfAtmosphereData({
                mp.WDParameter.logq: -16.22,
                ci.Element.Pm: 3345.6,
                ci.Element.Mc: 0.02
            })
        })
        white_dwarf3 = wd.WhiteDwarf('TestWD1', property_data1, abundance_data, atm_data)
        expected_str3 = "\n--- TestWD1 ---\n\nSpectral Type: DA\nAtmospheric Type: H\nTemperature: -1123.0±0.0 K\nLogg: 8.123±0.0\nMass: ---\nDistance: ---\n\nH: 0.0±0.0 (not included)\nLi: -4.6±0.6\nC: <-6.5\nK: >-12.8 (not included)\nFe: <-6.54 (not included)\nAg: -7.9^{+0.31}_{-0.33}\nPo: >-12.5\nDs: -7.6±0.4\nTs: -7.77±0.11 (not included)\n\nKoesterNoOvershoot Sinking Timescales:\nlog(q): -14.22±0.0\nPm: 345.6 yr\nSn: 987.5 yr\n\nKoesterOvershoot Sinking Timescales:\nlog(q): -16.22±0.0\nPm: 3345.6 yr\nMc: 0.02 yr\n\n\n"
        self.assertEqual(str(white_dwarf3), expected_str3)
        white_dwarf4 = wd.WhiteDwarf('TestWD1', property_data2, abundance_data, atm_data)
        self.assertEqual(6.0255958607435685e-12, white_dwarf4.get_logq_in_solar_masses(ti.TimescaleType.KoesterNoOvershoot)) # A bit out of place...

    def test_getters(self):
        wd_property_data_raw = {
            mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'DA'),
            mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -1123, 0),
            mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.123, 0),
            mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 0.5, 0)
        }
        property_data = wd.WhiteDwarfPropertyData(wd_property_data_raw)
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Ts: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, 6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        white_dwarf1 = wd.WhiteDwarf('TestWD', property_data, abundance_data)
        self.assertEqual(white_dwarf1.get_spectral_type().value, 'DA')
        self.assertEqual(white_dwarf1.get_atmospheric_type().value, ci.Element.H)
        self.assertEqual(white_dwarf1.get_teff().value, -1123)
        self.assertEqual(white_dwarf1.get_logg().value, 8.123)
        self.assertEqual(white_dwarf1.get_mass().value, 0.5)
        self.assertEqual(white_dwarf1.get_distance(), None)
        self.assertEqual(white_dwarf1.get_abundance(ci.Element.Li).value, -4.6)

    def test_get_plottable_points(self):
        els_to_plot = [ci.Element.Ag, ci.Element.Ds, ci.Element.Ts, ci.Element.C, ci.Element.Fe, ci.Element.K, ci.Element.Po, ci.Element.Na] # deliberately attempt to plot Na instead of Li
        els_to_plot_without_C = [ci.Element.Ag, ci.Element.Ds, ci.Element.Ts, ci.Element.Fe, ci.Element.K, ci.Element.Po, ci.Element.Na]
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Ts: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, 6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        wd_property_data_raw = {
            mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.He),
            mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -1123, 0),
            mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.123, 0)
        }
        property_data = wd.WhiteDwarfPropertyData(wd_property_data_raw)
        white_dwarf1 = wd.WhiteDwarf('TestWD', property_data, abundance_data)

        elements_to_plot_inc_refHf, reference_el_inc_refHf, plottable_values_inc_relHf, plottable_ue_inc_relHf, plottable_le_inc_relHf, plottable_ubb_inc_relHf, plottable_lbb_inc_relHf = white_dwarf1.get_plottable_points(els_to_plot)
        self.assertEqual(elements_to_plot_inc_refHf, els_to_plot)
        self.assertEqual(reference_el_inc_refHf, ci.Element.He)
        self.assertEqual(plottable_values_inc_relHf, [-7.9, -7.6, np.nan, -6.5, np.nan, np.nan, -12.5, np.nan])
        self.assertEqual(plottable_ue_inc_relHf, [0.31, 0.4, np.nan, 0, np.nan, np.nan, 0.1, np.nan])
        self.assertEqual(plottable_le_inc_relHf, [0.33, 0.4, np.nan, 0.1, np.nan, np.nan, 0, np.nan])
        self.assertEqual(plottable_ubb_inc_relHf, [False, False, False, True, False, False, False, False])
        self.assertEqual(plottable_lbb_inc_relHf, [False, False, False, False, False, False, True, False])

        elements_to_plot_exc_refHf, reference_el_exc_refHf, plottable_values_exc_relHf, plottable_ue_exc_relHf, plottable_le_exc_relHf, plottable_ubb_exc_relHf, plottable_lbb_exc_relHf = white_dwarf1.get_plottable_points(els_to_plot, None, False)
        self.assertEqual(elements_to_plot_exc_refHf, els_to_plot)
        self.assertEqual(reference_el_exc_refHf, ci.Element.He)
        self.assertEqual(plottable_values_exc_relHf, [np.nan, np.nan, -7.77, np.nan, -6.54, -12.8, np.nan, np.nan])
        self.assertEqual(plottable_ue_exc_relHf, [np.nan, np.nan, 0.11, np.nan, 0, 0.1, np.nan, np.nan])
        self.assertEqual(plottable_le_exc_relHf, [np.nan, np.nan, 0.11, np.nan, 0.1, 0, np.nan, np.nan])
        self.assertEqual(plottable_ubb_exc_relHf, [False, False, False, False, True, False, False, False])
        self.assertEqual(plottable_lbb_exc_relHf, [False, False, False, False, False, True, False, False])

        elements_to_plot_inc_refC, reference_el_inc_refC, plottable_values_inc_relC, plottable_ue_inc_relC, plottable_le_inc_relC, plottable_ubb_inc_relC, plottable_lbb_inc_relC = white_dwarf1.get_plottable_points(els_to_plot, ci.Element.C)
        self.assertEqual(elements_to_plot_inc_refC, els_to_plot_without_C)
        self.assertEqual(reference_el_inc_refC, ci.Element.C)
        self.assertEqual(plottable_values_inc_relC, [6.5-7.9, 6.5-7.6, np.nan, np.nan, np.nan, 6.5-12.5, np.nan])
        self.assertEqual(plottable_ue_inc_relC, [0.31, 0.4, np.nan, np.nan, np.nan, 0.1, np.nan])
        self.assertEqual(plottable_le_inc_relC, [0.33, 0.4, np.nan, np.nan, np.nan, 0, np.nan])
        self.assertEqual(plottable_ubb_inc_relC, [False, False, False, False, False, False, False])
        self.assertEqual(plottable_lbb_inc_relC, [False, False, False, False, False, True, False])

        elements_to_plot_exc_refC, reference_el_exc_refC, plottable_values_exc_relC, plottable_ue_exc_relC, plottable_le_exc_relC, plottable_ubb_exc_relC, plottable_lbb_exc_relC = white_dwarf1.get_plottable_points(els_to_plot, ci.Element.C, False)
        self.assertEqual(elements_to_plot_exc_refC, els_to_plot_without_C)
        self.assertEqual(reference_el_exc_refC, ci.Element.C)
        self.assertEqual(plottable_values_exc_relC, [np.nan, np.nan, 6.5-7.77, 6.5-6.54, 6.5-12.8, np.nan, np.nan])
        self.assertEqual(plottable_ue_exc_relC, [np.nan, np.nan, 0.11, 0, 0.1, np.nan, np.nan])
        self.assertEqual(plottable_le_exc_relC, [np.nan, np.nan, 0.11, 0.1, 0, np.nan, np.nan])
        self.assertEqual(plottable_ubb_exc_relC, [False, False, False, True, False, False, False])
        self.assertEqual(plottable_lbb_exc_relC, [False, False, False, False, True, False, False])

        with self.assertRaises(ValueError):
            elements_to_plot_exc_refC, reference_el_exc_refC, plottable_values_exc_relC, plottable_ue_exc_relC, plottable_le_exc_relC, plottable_ubb_exc_relC, plottable_lbb_exc_relC = white_dwarf1.get_plottable_points(els_to_plot, ci.Element.Hg, False)

    def test_get_all_elements(self):
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Ts: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, 6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        wd_property_data_raw = {
            mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'DA')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        property_data = wd.WhiteDwarfPropertyData(wd_property_data_raw)
        white_dwarf1 = wd.WhiteDwarf('TestWD', property_data, abundance_data)
        self.assertEqual(white_dwarf1.get_all_elements(), [ci.Element.H, ci.Element.Li, ci.Element.C, ci.Element.K, ci.Element.Fe, ci.Element.Ag, ci.Element.Po, ci.Element.Ds, ci.Element.Ts])

    def test_get_all_elements_overlapping_with_fit_dict(self):
        wd_property_data_raw = {
            mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.He),
            mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -1123, 0),
            mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.123, 0)
        }
        property_data = wd.WhiteDwarfPropertyData(wd_property_data_raw)
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Ts: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, 6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        white_dwarf1 = wd.WhiteDwarf('TestWD', property_data, abundance_data)
        fit_dict = {
            'FD1': {
                ci.Element.Ds: -2,
                ci.Element.Ag: -2,
                ci.Element.Ts: -2,
                ci.Element.C: -2,
                ci.Element.Fe: -2,
                ci.Element.Po: -2,
                ci.Element.Li: -2
            },
            'FD2': {
                ci.Element.Ds: -2,
                ci.Element.Ni: -2,
                ci.Element.Ts: -2,
                ci.Element.C: -2,
                ci.Element.Fe: -2,
                ci.Element.Po: -2,
                ci.Element.Li: -2
            }
        }
        expected_overlapping_elements = [
            ci.Element.Li,
            ci.Element.C,
            ci.Element.Fe,
            ci.Element.Po,
            ci.Element.Ds,
            ci.Element.Ts
        ]
        self.assertEqual(white_dwarf1.get_all_elements_overlapping_with_fit_dict(fit_dict), expected_overlapping_elements)

    def test_get_pollution_fraction_estimates(self):
        # Testing both the implied pollution fraction from (included) measurements and the maximum implied pollution fraction. This is not really relevant now
        wd_property_data_raw = {
            mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.H),
            mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -1123, 0),
            mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.123, 0)
        }
        property_data = wd.WhiteDwarfPropertyData(wd_property_data_raw)
        wd_abundance_data_raw = {
            ci.Element.Ds: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.6, 0.4),
            ci.Element.Ag: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-7.9', 0.31, True, '0.33'),
            ci.Element.Ts: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, '-4.77', 0.11, False),
            ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.5, 0.4, True, 'abc'),
            ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, 6.54, None, False),
            ci.Element.K: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '12.8', None, False),
            ci.Element.Po: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.lower_bound, '-12.5'),
            ci.Element.Li: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.6, '0.6')
        }
        abundance_data = wd.WhiteDwarfAbundanceData(wd_abundance_data_raw)
        white_dwarf1 = wd.WhiteDwarf('TestWD', property_data, abundance_data)
        expected_min_pol_frac_all = np.log10(10**(-7.6) + 10**(-7.9) + 10**(-12.5) + 10**(-4.6))
        expected_max_pol_frac_all = 3.3397825088208446  # Not writing this out in full!
        expected_min_pol_frac_ue = -np.inf # Not polluted by any of the usual elements! (At least, there's only an upper bound on C which doesnt help)
        expected_max_pol_frac_ue = 3.2611921814407916 # This is only slightly smaller than expected_max_pol_frac_all, indicating that all non-usual elements are trace
        self.assertEqual(white_dwarf1.estimate_minimum_pollution_fraction(), expected_min_pol_frac_all)
        self.assertEqual(white_dwarf1.estimate_maximum_pollution_fraction(), expected_max_pol_frac_all)
        self.assertEqual(white_dwarf1.estimate_minimum_pollution_fraction(ci.usual_elements), expected_min_pol_frac_ue)
        self.assertEqual(white_dwarf1.estimate_maximum_pollution_fraction(ci.usual_elements), expected_max_pol_frac_ue)

class SolarAbundancesTests(unittest.TestCase):

    def test_get_rel_abundance(self):
        self.assertAlmostEqual(sa.get_solar_relative_abundance(ci.Element.Fe, ci.Element.O), 7.46-8.69)
        self.assertAlmostEqual(sa.get_solar_relative_abundance(ci.Element.Fe, ci.Element.Be), 7.46-1.38)

class ChemistryTests(unittest.TestCase):

    def test_math_operations(self):
        self.assertEqual(ci.Element.C, ci.Element.C)
        self.assertNotEqual(ci.Element.C, ci.Element.N)
        self.assertTrue(ci.Element.C == ci.Element.C)
        self.assertTrue(ci.Element.C != ci.Element.Pt)
        self.assertFalse(ci.Element.C != ci.Element.C)
        self.assertFalse(ci.Element.C == ci.Element.Pt)
        self.assertTrue(ci.Element.C < ci.Element.Pt)
        self.assertFalse(ci.Element.C < ci.Element.C)
        self.assertTrue(ci.Element.C <= ci.Element.C)
        self.assertTrue(ci.Element.Au > ci.Element.Ag)
        self.assertFalse(ci.Element.Au > ci.Element.Au)
        self.assertTrue(ci.Element.Au >= ci.Element.Au)
        self.assertNotEqual(ci.Element.C, None)

class ManagerTests(unittest.TestCase):

    def setUp(self):
        self.real_wd_file = 'WDInputData.csv'
        self.real_stellar_comps_file = 'StellarCompositionsSortFE.csv'
        self.test_args_none = Namespace(
            wd_data_filename=None,
            stellar_compositions_filename=None,
            n_live_points=None,
            enhancement_model=None,
            pollution_model_names=None
        )
        self.test_args_invalid_sc = Namespace(
            wd_data_filename=None,
            stellar_compositions_filename='foo',
            n_live_points=None,
            enhancement_model='NonEarthlike',
            pollution_model_names=['bar', 'baz']
        )
        self.test_args_invalid_pm = Namespace(
            wd_data_filename=None,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=None,
            enhancement_model='NonEarthlike',
            pollution_model_names=['bar', 'baz']
        )
        self.test_args_invalid_em = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=1500,
            enhancement_model='Dummy',
            pollution_model_names=['Model_24']
        )
        self.test_args_golden = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=1500,
            enhancement_model='NonEarthlike',
            pollution_model_names=['Model_24']
        )
        self.test_args_invalid_hierarchy = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=1500,
            enhancement_model='NonEarthlike',
            pollution_model_names=['Hierarchy_NOTREAL']
        )
        self.test_args_hierarchy = Namespace(
            wd_data_filename=self.real_wd_file,
            stellar_compositions_filename=self.real_stellar_comps_file,
            n_live_points=1500,
            enhancement_model='NonEarthlike',
            pollution_model_names=['Hierarchy_Basic']
        )

    def test_init(self):
        test_mn = mn.Manager() # This should default to the settings in configuration.ini - although really this test should probably not be sensitive to the contents of that file!
        # And should also test on a different configuration file which uses specific models rather than a hierarchy and has malformed data/invalid file names etc

        self.assertEqual('WDInputData.csv', test_mn.wd_data_filename)
        self.assertEqual('StellarCompositionsSortFE.csv', test_mn.stellar_compositions_filename)
        self.assertEqual('NonEarthlike', test_mn.enhancement_model)
        self.assertEqual(-1, test_mn.seed)
        self.assertIsNone(test_mn.model_names)
        self.assertTrue(test_mn.use_hierarchy)
        self.assertEqual('Hierarchy_Default', test_mn.hierarchy_name)
        self.assertEqual(
            {
                0: [mp.ModelParameter.metallicity, mp.ModelParameter.t_sinceaccretion, mp.ModelParameter.fragment_mass, mp.ModelParameter.accretion_timescale],
                1: [mp.ModelParameter.formation_distance],
                2: [mp.ModelParameter.feeding_zone_size],
                3: [mp.ModelParameter.fragment_core_frac, mp.ModelParameter.pressure, mp.ModelParameter.oxygen_fugacity]
            }, test_mn.parameter_hierarchy
        )
        self.assertEqual(
            [
                ti.TimescaleType.KoesterNoOvershoot,
                ti.TimescaleType.KoesterOvershoot,
                ti.TimescaleType.BedardNoOvershoot,
                ti.TimescaleType.BedardOvershoot,
                ti.TimescaleType.BedardVariableOvershoot,
                ti.TimescaleType.Bedard3DOvershoot,
                ti.TimescaleType.MWDD,
                ti.TimescaleType.Bedard3DOvershootPatched
            ], test_mn.timescale_types_to_run
        )
        self.assertEqual([False, True], test_mn.thermohaline_regimes_to_run)
        self.assertEqual(False, test_mn.suppress_graphical_output)
        self.assertEqual(737, len(test_mn.white_dwarfs))
        self.assertEqual(8, test_mn.default_logg)
        self.assertEqual(-15, test_mn.default_Ca_value)
        self.assertEqual(2000, test_mn.n_live_points)
        self.assertEqual(958, len(test_mn.stellar_compositions))
        self.assertEqual(dict(), test_mn.models)
        self.assertEqual(dict(), test_mn.executed_models)

    def test_load_values(self):
        test_mn = mn.Manager()

        # For GALEX1931+0117:
        # For purposes of updating John's code, I filled in all the timescales/logqs in the csv file.
        # So they're not actually being calculated by the code any more - shouldn't matter, this entry is deprecated anyway
        expected_wd_222 = wd.WhiteDwarf(
            'GALEX1931+0117',
            wd.WhiteDwarfPropertyData({
                mp.WDParameter.atmospheric_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, ci.Element.H),
                mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 0.596, 0),
                mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 20409, 0),
                mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 7.95, 0)
            }),
            wd.WhiteDwarfAbundanceData({
                ci.Element.Ni: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.7, 0.3),
                ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.5, 0.3),
                ci.Element.O: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.1, 0.3),
                ci.Element.Si: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.75, 0.2),
                ci.Element.Cr: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.1, 0.3),
                ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.8, 0.3, False),
                ci.Element.Al: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.2, 0.2)
            }),
            wd.WhiteDwarfAtmosphereDataset({
                ti.TimescaleType.KoesterNoOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.192102000000006, #This happens to be the same as the overshoot case for this one!
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
                }),
                ti.TimescaleType.KoesterOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.192102000000006,
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
                }),
                ti.TimescaleType.BedardNoOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.371724,  # Again, this happens to be the same as the overshoot case
                    ci.Element.Al: 0.010849249580537668,
                    ci.Element.Ti: 0.0063664187373546465,
                    ci.Element.Ca: 0.007253942354258323,
                    ci.Element.Ni: 0.00486953172940064,
                    ci.Element.Fe: 0.004988893121240738,
                    ci.Element.Cr: 0.005714556092674965,
                    ci.Element.Mg: 0.011832055551331475,
                    ci.Element.Si: 0.009561114078236264,
                    ci.Element.Na: 0.004004110001747341,
                    ci.Element.O: 0.005746163775920785,
                    ci.Element.C: 0.009213301851464161,
                    ci.Element.N: 0.006946819264894303
                }),
                ti.TimescaleType.BedardOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.371724,
                    ci.Element.Al: 0.010849249580537668,
                    ci.Element.Ti: 0.0063664187373546465,
                    ci.Element.Ca: 0.007253942354258323,
                    ci.Element.Ni: 0.00486953172940064,
                    ci.Element.Fe: 0.004988893121240738,
                    ci.Element.Cr: 0.005714556092674965,
                    ci.Element.Mg: 0.011832055551331475,
                    ci.Element.Si: 0.009561114078236264,
                    ci.Element.Na: 0.004004110001747341,
                    ci.Element.O: 0.005746163775920785,
                    ci.Element.C: 0.009213301851464161,
                    ci.Element.N: 0.006946819264894303
                }),
                ti.TimescaleType.BedardVariableOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.371724,
                    ci.Element.Al: 0.010849249580537668,
                    ci.Element.Ti: 0.0063664187373546465,
                    ci.Element.Ca: 0.007253942354258323,
                    ci.Element.Ni: 0.00486953172940064,
                    ci.Element.Fe: 0.004988893121240738,
                    ci.Element.Cr: 0.005714556092674965,
                    ci.Element.Mg: 0.011832055551331475,
                    ci.Element.Si: 0.009561114078236264,
                    ci.Element.Na: 0.004004110001747341,
                    ci.Element.O: 0.005746163775920785,
                    ci.Element.C: 0.009213301851464161,
                    ci.Element.N: 0.006946819264894303
                }),
                ti.TimescaleType.Bedard3DOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.1352631,
                    ci.Element.Al: 0.017029453066988692,
                    ci.Element.Ti: 0.010097398586142536,
                    ci.Element.Ca: 0.010374044342579441,
                    ci.Element.Ni: 0.007123719326677305,
                    ci.Element.Fe: 0.007782660754446889,
                    ci.Element.Cr: 0.008672368034467732,
                    ci.Element.Mg: 0.01699871920436737,
                    ci.Element.Si: 0.014857833936371927,
                    ci.Element.Na: 0.005735496248057374,
                    ci.Element.O: 0.008639760426209783,
                    ci.Element.C: 0.016712218304083728,
                    ci.Element.N: 0.011344286863441407
                }),
                ti.TimescaleType.MWDD: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.3819181,
                    ci.Element.Al: 0.010650029236420904,
                    ci.Element.Ti: 0.006243632031556506,
                    ci.Element.Ca: 0.007135509685949294,
                    ci.Element.Ni: 0.004785519607664693,
                    ci.Element.Fe: 0.00489774533909043,
                    ci.Element.Cr: 0.005615671565486851,
                    ci.Element.Mg: 0.011640461532735517,
                    ci.Element.Si: 0.009386431173243375,
                    ci.Element.Na: 0.0039389013135738525,
                    ci.Element.O: 0.005641711839238161,
                    ci.Element.C: 0.008991937922430205,
                    ci.Element.N: 0.006798094540750175
                }),
                ti.TimescaleType.Bedard3DOvershootPatched: wd.WhiteDwarfAtmosphereData({
                    # This WD happens to sit right inbetween the VO and Patched grids, so need to interpolate. (We have 3 DAs in this test for this reason! One for each regime)
                    mp.WDParameter.logq: -16.2271654,
                    ci.Element.Al: 0.014248757655435654,
                    ci.Element.Ti: 0.008395479393529479,
                    ci.Element.Ca: 0.00903200563770281,
                    ci.Element.Ni: 0.00615039204598485,
                    ci.Element.Fe: 0.00653519574068265,
                    ci.Element.Cr: 0.007341513470696146,
                    ci.Element.Mg: 0.014788533177719006,
                    ci.Element.Si: 0.012524572180999295,
                    ci.Element.Na: 0.004993053275855988,
                    ci.Element.O: 0.007356047077691492,
                    ci.Element.C: 0.013149302002152646,
                    ci.Element.N: 0.00929902253134939
                }),
            })
        )
        expected_wd_236 = wd.WhiteDwarf(
            'SDSSJ1043+0855',
            wd.WhiteDwarfPropertyData({
                mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'DA'),
                mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 0.626, 0),
                mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 18330, 0),
                mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.05, 0)
            }),
            wd.WhiteDwarfAbundanceData({
                ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.15, 0.3, False),
                ci.Element.O: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.9, 0.2),
                ci.Element.Mg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.11, 0.2),
                ci.Element.Al: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.06, 0.3),
                ci.Element.Si: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.33, 0.5),
                ci.Element.Ca: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 5.96, 0.2),

                ci.Element.Ti: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -7.0, 0),
                ci.Element.Cr: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, 6.5, 0),
                ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 6.15, 0.3),
                ci.Element.Ni: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 7.38, 0.3)
            }),
            wd.WhiteDwarfAtmosphereDataset({
                ti.TimescaleType.KoesterOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.555411999999997,
                    ci.Element.Al: 0.005934279535546039,
                    ci.Element.Ti: 0.003671267248872019,
                    ci.Element.Ca: 0.004336706191208475,
                    ci.Element.Ni: 0.0026668586645214724,
                    ci.Element.Fe: 0.002935567020550475,
                    ci.Element.Cr: 0.0032614896149133145,
                    ci.Element.Mg: 0.007130959994285996,
                    ci.Element.Si: 0.005495965527215545,
                    ci.Element.Na: 0.0024623625649367433,
                    ci.Element.O: 0.003469189314151306,
                    ci.Element.C: 0.005206381643052118,
                    ci.Element.N: 0.004071139650067431
                }),
                ti.TimescaleType.KoesterNoOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.555411999999997,
                    ci.Element.Al: 0.005934279535546039,
                    ci.Element.Ti: 0.003671267248872019,
                    ci.Element.Ca: 0.004336706191208475,
                    ci.Element.Ni: 0.0026668586645214724,
                    ci.Element.Fe: 0.002935567020550475,
                    ci.Element.Cr: 0.0032614896149133145,
                    ci.Element.Mg: 0.007130959994285996,
                    ci.Element.Si: 0.005495965527215545,
                    ci.Element.Na: 0.0024623625649367433,
                    ci.Element.O: 0.003469189314151306,
                    ci.Element.C: 0.005206381643052118,
                    ci.Element.N: 0.004071139650067431
                }),
                ti.TimescaleType.BedardNoOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.639186000000002,
                    ci.Element.Al: 0.005873338075387059,
                    ci.Element.Ti: 0.003627028766101857,
                    ci.Element.Ca: 0.0043245109389184145,
                    ci.Element.Ni: 0.0027672920791231063,
                    ci.Element.Fe: 0.002640164577940586,
                    ci.Element.Cr: 0.003354524746220148,
                    ci.Element.Mg: 0.006879802958963415,
                    ci.Element.Si: 0.004886872406605205,
                    ci.Element.Na: 0.0024129093149292415,
                    ci.Element.O: 0.003327943925537073,
                    ci.Element.C: 0.0048434296938909315,
                    ci.Element.N: 0.00396537258750075
                }),
                ti.TimescaleType.BedardOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.639186000000002,
                    ci.Element.Al: 0.005873338075387059,
                    ci.Element.Ti: 0.003627028766101857,
                    ci.Element.Ca: 0.0043245109389184145,
                    ci.Element.Ni: 0.0027672920791231063,
                    ci.Element.Fe: 0.002640164577940586,
                    ci.Element.Cr: 0.003354524746220148,
                    ci.Element.Mg: 0.006879802958963415,
                    ci.Element.Si: 0.004886872406605205,
                    ci.Element.Na: 0.0024129093149292415,
                    ci.Element.O: 0.003327943925537073,
                    ci.Element.C: 0.0048434296938909315,
                    ci.Element.N: 0.00396537258750075
                }),
                ti.TimescaleType.BedardVariableOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.639186000000002,
                    ci.Element.Al: 0.005873338075387059,
                    ci.Element.Ti: 0.003627028766101857,
                    ci.Element.Ca: 0.0043245109389184145,
                    ci.Element.Ni: 0.0027672920791231063,
                    ci.Element.Fe: 0.002640164577940586,
                    ci.Element.Cr: 0.003354524746220148,
                    ci.Element.Mg: 0.006879802958963415,
                    ci.Element.Si: 0.004886872406605205,
                    ci.Element.Na: 0.0024129093149292415,
                    ci.Element.O: 0.003327943925537073,
                    ci.Element.C: 0.0048434296938909315,
                    ci.Element.N: 0.00396537258750075
                }),
                ti.TimescaleType.Bedard3DOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.318444,
                    ci.Element.Al: 0.010802020570960715,
                    ci.Element.Ti: 0.00634319408695934,
                    ci.Element.Ca: 0.0070418861610901376,
                    ci.Element.Ni: 0.004814442453513755,
                    ci.Element.Fe: 0.0049641054680102815,
                    ci.Element.Cr: 0.005706843986710429,
                    ci.Element.Mg: 0.011461699023744169,
                    ci.Element.Si: 0.009487111395555562,
                    ci.Element.Na: 0.003950026415978833,
                    ci.Element.O: 0.0058830291460484395,
                    ci.Element.C: 0.009371691846749701,
                    ci.Element.N: 0.007162835873271362
                }),
                ti.TimescaleType.MWDD: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.644821999999998,
                    ci.Element.Al: 0.005797930431290188,
                    ci.Element.Ti: 0.003585667377489108,
                    ci.Element.Ca: 0.0042847254588627036,
                    ci.Element.Ni: 0.002734683057440141,
                    ci.Element.Fe: 0.0026008721243337006,
                    ci.Element.Cr: 0.0033213909654360727,
                    ci.Element.Mg: 0.006810798002720251,
                    ci.Element.Si: 0.004819722122114213,
                    ci.Element.Na: 0.002390710560716041,
                    ci.Element.O: 0.0032943140613688024,
                    ci.Element.C: 0.004784559943242473,
                    ci.Element.N: 0.00392530131692707
                }),
                ti.TimescaleType.Bedard3DOvershootPatched: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.318444,
                    ci.Element.Al: 0.010802020570960715,
                    ci.Element.Ti: 0.00634319408695934,
                    ci.Element.Ca: 0.0070418861610901376,
                    ci.Element.Ni: 0.004814442453513755,
                    ci.Element.Fe: 0.0049641054680102815,
                    ci.Element.Cr: 0.005706843986710429,
                    ci.Element.Mg: 0.011461699023744169,
                    ci.Element.Si: 0.009487111395555562,
                    ci.Element.Na: 0.003950026415978833,
                    ci.Element.O: 0.0058830291460484395,
                    ci.Element.C: 0.009371691846749701,
                    ci.Element.N: 0.007162835873271362
                })
            })
        )
        expected_wd_301 = wd.WhiteDwarf(
            'GaiaJ0510+2315',
            wd.WhiteDwarfPropertyData({
                mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'DA'),
                mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 0.746, 0),
                mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 21700, 0),
                mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.2, 0)
            }),
            wd.WhiteDwarfAbundanceData({
                ci.Element.O: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -4.2348, 0.1),
                ci.Element.Mg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.2347, 0.1),
                ci.Element.Si: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.119024, 0.1),
                ci.Element.Ca: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 6.3148, 0.1),
            }),
            wd.WhiteDwarfAtmosphereDataset({
                ti.TimescaleType.KoesterOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.50760000000001,
                    ci.Element.C: 0.006980395019858036,
                    ci.Element.N: 0.003971549684604506,
                    ci.Element.O: 0.0030512655099359434,
                    ci.Element.Na: 0.002040609950289374,
                    ci.Element.Mg: 0.0061665179497525004,
                    ci.Element.Al: 0.006177318585915238,
                    ci.Element.Si: 0.005430002130657389,
                    ci.Element.Ca: 0.0036264358550606277,
                    ci.Element.Ti: 0.0034692851727188426,
                    ci.Element.Cr: 0.0030163379624671663,
                    ci.Element.Fe: 0.0027820184125097835,
                    ci.Element.Ni: 0.0025309292209157543
                }),
                ti.TimescaleType.KoesterNoOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.50760000000001,
                    ci.Element.C: 0.006980395019858036,
                    ci.Element.N: 0.003971549684604506,
                    ci.Element.O: 0.0030512655099359434,
                    ci.Element.Na: 0.002040609950289374,
                    ci.Element.Mg: 0.0061665179497525004,
                    ci.Element.Al: 0.006177318585915238,
                    ci.Element.Si: 0.005430002130657389,
                    ci.Element.Ca: 0.0036264358550606277,
                    ci.Element.Ti: 0.0034692851727188426,
                    ci.Element.Cr: 0.0030163379624671663,
                    ci.Element.Fe: 0.0027820184125097835,
                    ci.Element.Ni: 0.0025309292209157543
                }),
                ti.TimescaleType.BedardNoOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.685,
                    ci.Element.C: 0.0043617422692703725,
                    ci.Element.N: 0.003138773664086122,
                    ci.Element.O: 0.0025172131972423053,
                    ci.Element.Na: 0.001727348630497726,
                    ci.Element.Mg: 0.005092605092574638,
                    ci.Element.Al: 0.004829031073766105,
                    ci.Element.Si: 0.004241505163138565,
                    ci.Element.Ca: 0.0031128629469706326,
                    ci.Element.Ti: 0.002839487753868497,
                    ci.Element.Cr: 0.0025019579136939373,
                    ci.Element.Fe: 0.0022115664858298892,
                    ci.Element.Ni: 0.0021090166093702097
                }),
                ti.TimescaleType.BedardOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.685,
                    ci.Element.C: 0.0043617422692703725,
                    ci.Element.N: 0.003138773664086122,
                    ci.Element.O: 0.0025172131972423053,
                    ci.Element.Na: 0.001727348630497726,
                    ci.Element.Mg: 0.005092605092574638,
                    ci.Element.Al: 0.004829031073766105,
                    ci.Element.Si: 0.004241505163138565,
                    ci.Element.Ca: 0.0031128629469706326,
                    ci.Element.Ti: 0.002839487753868497,
                    ci.Element.Cr: 0.0025019579136939373,
                    ci.Element.Fe: 0.0022115664858298892,
                    ci.Element.Ni: 0.0021090166093702097
                }),
                ti.TimescaleType.BedardVariableOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.685,
                    ci.Element.C: 0.0043617422692703725,
                    ci.Element.N: 0.003138773664086122,
                    ci.Element.O: 0.0025172131972423053,
                    ci.Element.Na: 0.001727348630497726,
                    ci.Element.Mg: 0.005092605092574638,
                    ci.Element.Al: 0.004829031073766105,
                    ci.Element.Si: 0.004241505163138565,
                    ci.Element.Ca: 0.0031128629469706326,
                    ci.Element.Ti: 0.002839487753868497,
                    ci.Element.Cr: 0.0025019579136939373,
                    ci.Element.Fe: 0.0022115664858298892,
                    ci.Element.Ni: 0.0021090166093702097
                }),
                ti.TimescaleType.Bedard3DOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.624200000000002,
                    ci.Element.C: 0.005143514873243852,
                    ci.Element.N: 0.0035988189612589237,
                    ci.Element.O: 0.002798981319634372,
                    ci.Element.Na: 0.0018933024635704647,
                    ci.Element.Mg: 0.005584187600997102,
                    ci.Element.Al: 0.005438761317204103,
                    ci.Element.Si: 0.004722587355729608,
                    ci.Element.Ca: 0.003409887100677623,
                    ci.Element.Ti: 0.0032216621887342744,
                    ci.Element.Cr: 0.002802592787424118,
                    ci.Element.Fe: 0.0024774220576332913,
                    ci.Element.Ni: 0.002318996002505059
                }),
                ti.TimescaleType.MWDD: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.70312,
                    ci.Element.C: 0.00417157446120628,
                    ci.Element.N: 0.003023987558303149,
                    ci.Element.O: 0.0024396077256845523,
                    ci.Element.Na: 0.0016791906215408389,
                    ci.Element.Mg: 0.004950624650796981,
                    ci.Element.Al: 0.004667023632109997,
                    ci.Element.Si: 0.004102041029866081,
                    ci.Element.Ca: 0.003026913428101311,
                    ci.Element.Ti: 0.002742963319467846,
                    ci.Element.Cr: 0.0024242644894489043,
                    ci.Element.Fe: 0.0021408192439307974,
                    ci.Element.Ni: 0.002050123387357036
                }),
                ti.TimescaleType.Bedard3DOvershootPatched: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -16.685,
                    ci.Element.C: 0.0043617422692703725,
                    ci.Element.N: 0.003138773664086122,
                    ci.Element.O: 0.0025172131972423053,
                    ci.Element.Na: 0.001727348630497726,
                    ci.Element.Mg: 0.005092605092574638,
                    ci.Element.Al: 0.004829031073766105,
                    ci.Element.Si: 0.004241505163138565,
                    ci.Element.Ca: 0.0031128629469706326,
                    ci.Element.Ti: 0.002839487753868497,
                    ci.Element.Cr: 0.0025019579136939373,
                    ci.Element.Fe: 0.0022115664858298892,
                    ci.Element.Ni: 0.0021090166093702097
                })
            })
        )
        expected_wd_360 = wd.WhiteDwarf(
            'GD61_Corr',
            wd.WhiteDwarfPropertyData({
                mp.WDParameter.spectral_type: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.label, 'DB'), # Should infer He
                mp.WDParameter.mass: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 0.71, 0),
                mp.WDParameter.temperature: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 17280, 0),
                mp.WDParameter.logg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, 8.2, 0)
            }),
            wd.WhiteDwarfAbundanceData({
                ci.Element.Ni: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -8.8, 0),
                ci.Element.Fe: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.6, 0.0667),
                ci.Element.O: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -5.95, 0.0434),
                ci.Element.Si: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.82, 0.0367),
                ci.Element.Na: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -6.8, 0),
                ci.Element.Mg: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -6.69, 0.0467),
                ci.Element.Cr: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -8, 0),
                ci.Element.C: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -9.1, 0, False),
                ci.Element.N: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -8, 0, False),
                ci.Element.Al: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, -7.8, 0),
                ci.Element.Ca: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.measurement, -7.9, 0.0634),
                ci.Element.Ti: wd.WhiteDwarfDataPoint(wd.WhiteDwarfDataPointType.upper_bound, 8.6, 0)  # The minus sign gets added automatically
            }),
            wd.WhiteDwarfAtmosphereDataset({
                ti.TimescaleType.KoesterNoOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -6.675699199999988,
                    ci.Element.Al: 112266.76280212057,
                    ci.Element.Ti: 72910.8229201416,
                    ci.Element.Ca: 88672.14208493519,
                    ci.Element.Ni: 64888.53961297752,
                    ci.Element.Fe: 66958.35876460896,
                    ci.Element.Cr: 68860.15564372398,
                    ci.Element.Mg: 125628.44956781625,
                    ci.Element.Si: 113922.9441342239,
                    ci.Element.Na: 124936.11584994856,
                    ci.Element.O: 167753.65835695906,
                    ci.Element.C: 204158.75094004275,
                    ci.Element.N: 184649.78231741278
                }),
                ti.TimescaleType.KoesterOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -6.247375999999992,
                    ci.Element.Al: 232973.86252231084,
                    ci.Element.Ti: 154620.54509961695,
                    ci.Element.Ca: 189274.1044509627,
                    ci.Element.Ni: 137270.11578430876,
                    ci.Element.Fe: 139050.0662120362,
                    ci.Element.Cr: 147321.86225016386,
                    ci.Element.Mg: 256372.8315223011,
                    ci.Element.Si: 238523.51461649648,
                    ci.Element.Na: 255006.0005828049,
                    ci.Element.O: 335529.87881929334,
                    ci.Element.C: 404414.9483500713,
                    ci.Element.N: 368830.81649282615
                }),
                ti.TimescaleType.BedardNoOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -6.876639999999999,
                    ci.Element.Al: 83964.55681227875,
                    ci.Element.Ti: 60594.331704769946,
                    ci.Element.Ca: 64419.2998025116,
                    ci.Element.Ni: 56290.09738312519,
                    ci.Element.Fe: 59387.630806465495,
                    ci.Element.Cr: 61435.01943511929,
                    ci.Element.Mg: 93473.3913912631,
                    ci.Element.Si: 78376.16331580719,
                    ci.Element.Na: 93735.4783465025,
                    ci.Element.O: 99306.11683163364,
                    ci.Element.C: 104911.71825325783,
                    ci.Element.N: 100374.63984633211
                }),
                ti.TimescaleType.BedardOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -6.3587039999999995,
                    ci.Element.Al: 262325.1923613523,
                    ci.Element.Ti: 159399.88362882257,
                    ci.Element.Ca: 172857.40940931826,
                    ci.Element.Ni: 168562.12592238086,
                    ci.Element.Fe: 168611.81370999877,
                    ci.Element.Cr: 164488.67414412534,
                    ci.Element.Mg: 277955.96593716956,
                    ci.Element.Si: 257660.5918275943,
                    ci.Element.Na: 267167.72970757366,
                    ci.Element.O: 279151.5219556833,
                    ci.Element.C: 308790.5492629715,
                    ci.Element.N: 285522.2780526775
                }),
                ti.TimescaleType.BedardVariableOvershoot: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -6.333536,
                    ci.Element.Al: 274976.76604823413,
                    ci.Element.Ti: 166289.18316346064,
                    ci.Element.Ca: 180810.64690052296,
                    ci.Element.Ni: 176564.74844027977,
                    ci.Element.Fe: 176080.79890136025,
                    ci.Element.Cr: 171550.50466665142,
                    ci.Element.Mg: 289771.72090777167,
                    ci.Element.Si: 270979.2272718993,
                    ci.Element.Na: 278560.79934108024,
                    ci.Element.O: 291018.09935393475,
                    ci.Element.C: 324226.12063769985,
                    ci.Element.N: 298109.626215266
                }),
                ti.TimescaleType.MWDD: wd.WhiteDwarfAtmosphereData({
                    mp.WDParameter.logq: -6.793999999999999,
                    ci.Element.Al: 97884.06573284,
                    ci.Element.Ti: 69282.57171806743,
                    ci.Element.Ca: 73506.8817757163,
                    ci.Element.Ni: 65059.69155790839,
                    ci.Element.Fe: 68773.95497861334,
                    ci.Element.Cr: 70504.36411431662,
                    ci.Element.Mg: 108013.98774569132,
                    ci.Element.Si: 90957.81099372983,
                    ci.Element.Na: 107886.72254261986,
                    ci.Element.O: 113588.92776665377,
                    ci.Element.C: 120153.38212358896,
                    ci.Element.N: 114925.39365030358
                })
            })
        )
        self.assertEqual(test_mn.white_dwarfs[222], expected_wd_222)
        self.assertEqual(test_mn.white_dwarfs[236], expected_wd_236)
        self.assertEqual(test_mn.white_dwarfs[301], expected_wd_301)
        self.assertEqual(test_mn.white_dwarfs[360], expected_wd_360)

    def test_publish_live_data(self):
        test_mn = mn.Manager()
        test_mn.publish_live_data(360, ti.TimescaleType.KoesterOvershoot)
        self.assertEqual(ld._live_white_dwarf, test_mn.white_dwarfs[360])
        expected_result = dict()

        expected_live_all_wd_abundances = np.array([None, None, -7.9, None, -7.6, None, -6.69, -6.82, None, -5.95, None, None])
        expected_live_all_wd_errors = np.array([None, None, 0.0634, None, 0.0667, None, 0.0467, 0.0367, None, 0.0434, None, None])
        expected_live_all_wd_timescales = np.array([
            232973.86252231084,
            154620.54509961695,
            189274.1044509627,
            137270.11578430876,
            139050.0662120362,
            147321.86225016386,
            256372.8315223011,
            238523.51461649648,
            255006.0005828049,
            335529.87881929334,
            404414.9483500713,
            368830.81649282615
        ])
        expected_live_upper_bounds = [-7.8, -8.6, None, -8.8, None, -8.0, None, None, -6.8, None, None, None]
        expected_live_lower_bounds = [None, None, None, None, None, None, None, None, None, None, None, None]
        expected_live_t_mg = expected_live_all_wd_timescales[6]
        expected_live_stellar_compositions_index123 = np.array([
            0.064419323,
            0.002751825,
            0.063013045,
            0.045735057,
            0.831858407,
            0.012310671,
            0.912736277,
            0.034678637,
            18.60872003,
            7.395480328,
            1.862591488
        ])
        expected_live_q = -6.247375999999992
        expected_live_mass = 0.71
        expected_live_type = ci.Element.He
        expected_live_elements_present = [ci.Element.Ca, ci.Element.Fe, ci.Element.Mg, ci.Element.Si, ci.Element.O]
        expected_live_teff = 17280.0
        expected_live_logg = 8.2

        self.assertTrue((ld._live_all_wd_abundances == expected_live_all_wd_abundances).all())
        self.assertTrue((ld._live_all_wd_errors == expected_live_all_wd_errors).all())
        self.assertTrue((ld._live_all_wd_timescales == expected_live_all_wd_timescales).all())
        self.assertTrue((ld._live_upper_bounds == expected_live_upper_bounds).all())
        self.assertTrue((ld._live_lower_bounds == expected_live_lower_bounds).all())
        self.assertEqual(ld._live_t_mg, expected_live_t_mg)
        self.assertTrue((ld._live_stellar_compositions[123] == expected_live_stellar_compositions_index123).all())
        self.assertEqual(ld._live_q, expected_live_q)
        self.assertEqual(ld._live_mass, expected_live_mass)
        self.assertEqual(ld._live_type, expected_live_type)
        self.assertEqual(set(ld._live_elements_present), set(expected_live_elements_present))
        self.assertEqual(ld._live_teff, expected_live_teff)
        self.assertEqual(ld._live_logg, expected_live_logg)

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
        with open(get_path_to_data() + input_filename) as generic_csv:
            generic_list =  [row for row in csv.reader(generic_csv)]
        generic_array = np.asarray(generic_list)
        return generic_array.astype(float)

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
        self.assertAlmostEqual(0.19087291612663068, old_gammas[ci.Element.O])
        self.assertAlmostEqual(5.938117878057908, old_gammas[ci.Element.Si])
        self.assertAlmostEqual(0.8762221554491502, old_gammas[ci.Element.Fe])
        self.assertAlmostEqual(0.6472702432098675, old_gammas[ci.Element.Mn])
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
        self.assertAlmostEqual(0.3906320799707601, gammas[ci.Element.O])
        self.assertAlmostEqual(1.186036978787397, gammas[ci.Element.V])
        self.assertAlmostEqual(1.2185224417219145, gammas_e[ci.Element.V])

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
            ci.Element.Ta: 0.017342227950395302,
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
            ci.Element.Ta: 0.09111458357273382,
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
        expected_abundances = {
            ci.Element.Hf: {
                gi.Layer.bulk: 2.7864515603493478e-08,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 3.457773385643354e-08
            },
            ci.Element.Ta: {
                gi.Layer.bulk: 3.616581827223245e-09,
                gi.Layer.core: 3.2068772969822515e-09,
                gi.Layer.mantle: 3.7152894227473975e-09
            },
            ci.Element.Nb: {
                gi.Layer.bulk: 1.2397067644439534e-07,
                gi.Layer.core: 3.946115224744226e-07,
                gi.Layer.mantle: 5.876693580055409e-08
            },
            ci.Element.Si: {
                gi.Layer.bulk: 0.15013035395386498,
                gi.Layer.core: 0.131781696676525,
                gi.Layer.mantle: 0.15455098277888274
            },
            ci.Element.Mn: {
                gi.Layer.bulk: 0.00038117866818521953,
                gi.Layer.core: 0.001165849993557462,
                gi.Layer.mantle: 0.0001921329047077404
                },
            ci.Element.Zn: {
                gi.Layer.bulk: 1.6014996989478682e-05,
                gi.Layer.core: 7.35533200780222e-06,
                gi.Layer.mantle: 1.8101314571890158e-05
            },
            ci.Element.Ga: {
                gi.Layer.bulk: 1.1263074987867517e-06,
                gi.Layer.core: 1.0971909698821512e-06,
                gi.Layer.mantle: 1.1333223983472394e-06
            },
            ci.Element.V: {
                gi.Layer.bulk: 5.395470924151831e-05,
                gi.Layer.core: 0.00016723772625079367,
                gi.Layer.mantle: 2.6662168047197238e-05
            },
            ci.Element.Cr: {
                gi.Layer.bulk: 0.0023669608105361967,
                gi.Layer.core: 0.007075860701522058,
                gi.Layer.mantle: 0.001232476162161708
            },
            ci.Element.Cu: {
                gi.Layer.bulk: 2.4715808308283354e-05,
                gi.Layer.core: 8.958797025011969e-05,
                gi.Layer.mantle: 9.086579624332024e-06
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.15006893982819816,
                gi.Layer.core: 0.6750980730644544,
                gi.Layer.mantle: 0.023577078497163388
            },
            ci.Element.W: {
                gi.Layer.bulk: 2.420586998738845e-08,
                gi.Layer.core: 1.2000974202323387e-07,
                gi.Layer.mantle: 1.1244647986698433e-09
            },
            ci.Element.P: {
                gi.Layer.bulk: 0.0006042600345224482,
                gi.Layer.core: 0.0030273657178258636,
                gi.Layer.mantle: 2.047690622093159e-05
            },
            ci.Element.Co: {
                gi.Layer.bulk: 0.00039087192643566125,
                gi.Layer.core: 0.0018574823398001305,
                gi.Layer.mantle: 3.7530997605301916e-05
            },
            ci.Element.Ni: {
                gi.Layer.bulk: 0.008120759633256047,
                gi.Layer.core: 0.03848874529973153,
                gi.Layer.mantle: 0.0008043979309279375
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.48615723914543096,
                gi.Layer.core: 0.10797785244937178,
                gi.Layer.mantle: 0.5772695668250886
            },
            ci.Element.C: {
                gi.Layer.bulk: 0.0015917333226106879,
                gi.Layer.core: 0.00819296152013823,
                gi.Layer.mantle: 1.3422474852156032e-06
            },
            ci.Element.S: {
                gi.Layer.bulk: 0.005033944726789019,
                gi.Layer.core: 0.02506831618945311,
                gi.Layer.mantle: 0.00020719361844052199
            },
            ci.Element.N: {
                gi.Layer.bulk: 4.6744203944017475e-05,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 5.800598389388276e-05
            },
            ci.Element.Na: {
                gi.Layer.bulk: 0.0020508290816938467,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0025449221217738738
            },
            ci.Element.Mg: {
                gi.Layer.bulk: 0.1659389539738732,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.20591755724632782
            },
            ci.Element.Al: {
                gi.Layer.bulk: 0.015403870863974397,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.019115026246018785
            },
            ci.Element.Ti: {
                gi.Layer.bulk: 0.00044298713595743367,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0005497131730881219
            },
            ci.Element.Ca: {
                gi.Layer.bulk: 0.011174350504526264,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.013866514791147875
            }
        }

        expected_Ds = {
            ci.Element.Hf: 0.0,
            ci.Element.Ta: 0.8629279364003023,
            ci.Element.Nb: 6.71307582961421,
            ci.Element.Si: 0.8524485194804874,
            ci.Element.Mn: 6.066326250546324,
            ci.Element.Zn: 0.40623468571255655,
            ci.Element.Ga: 0.9678623359106977,
            ci.Element.V:  6.270809823944999,
            ci.Element.Cr: 5.739652291231085,
            ci.Element.Cu: 9.856758082205513,
            ci.Element.Fe: 28.626069153226876,
            ci.Element.W:  106.69780247915111,
            ci.Element.P:  147.80372669960158,
            ci.Element.Co: 49.47883005324342,
            ci.Element.Ni: 47.83520571642847,
            ci.Element.O: 0.1869996786837244,
            ci.Element.C: 6102.294362221525,
            ci.Element.S:  120.95772947244798,
            ci.Element.N: 0,
            ci.Element.Na: 0,
            ci.Element.Mg: 0,
            ci.Element.Al: 0,
            ci.Element.Ti: 0,
            ci.Element.Ca: 0
        }
        expected_pcnf = 0.1941485586953566
        for element, element_data in abundances.items():
            for layer, abundance in element_data.items():
                self.assertAlmostEqual(expected_abundances[element][layer], abundance)
        for element, D in Ds.items():
            self.assertEqual(expected_Ds[element], D)
        self.assertEqual(expected_pcnf, w_met)
        abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(5, -2)
        expected_pcnf = 0.14338933574162993
        self.assertEqual(expected_pcnf, w_met)
        abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(54, -1) # More oxidised = more Fe as FeO in the mantle = smaller core
        expected_pcnf = 0.10245562602329299
        abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(54, -3) # And Vice versa
        expected_pcnf = 0.19639869216833447
        self.assertEqual(expected_pcnf, w_met)
        abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(0, -1)
        expected_pcnf = 0.08746988463747624
        self.assertEqual(expected_pcnf, w_met)
        abundances, w_met, Ds, all_Ds = geology_model.form_a_planet_iteratively(0, -3)
        expected_pcnf = 0.1576935740664529
        self.assertEqual(expected_pcnf, w_met)

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

    def test_cnf_cmf_conversions(self):
        geology_model = gi.GeologyModel()
        crm, mrm = geology_model.get_relative_core_mantle_masses()
        self.assertEqual(50.74093629571959, crm)
        self.assertEqual(20.88252846522537, mrm)
        #Call it twice to check caching
        crm, mrm = geology_model.get_relative_core_mantle_masses()
        self.assertEqual(50.74093629571959, crm)
        self.assertEqual(20.88252846522537, mrm)
        cmf = geology_model.convert_core_number_fraction_to_core_mass_fraction(0.17)
        self.assertEqual(0.33229859947479184, cmf)
        cnf = geology_model.convert_core_mass_fraction_to_core_number_fraction(0.33229859947479184)
        self.assertEqual(0.17, cnf)
        cnf = geology_model.convert_core_mass_fraction_to_core_number_fraction(0.3)
        self.assertEqual(0.1499340993911697, cnf)
        composition_override = { # This makes core density roughly double mantle density
            ci.Element.Al: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.Ti: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.Ca: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.Ni: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.Fe: {
                gi.Layer.core: 1,
                gi.Layer.mantle: 0
            },
            ci.Element.Cr: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.Mg: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.Si: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 1
            },
            ci.Element.Na: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.O: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.C: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            },
            ci.Element.N: {
                gi.Layer.core: 0,
                gi.Layer.mantle: 0
            }
        }
        cmf = geology_model.convert_core_number_fraction_to_core_mass_fraction(0.17, composition_override)
        self.assertEqual(0.28940349101639423, cmf)
        cnf = geology_model.convert_core_mass_fraction_to_core_number_fraction(0.28940349101639423, composition_override)
        self.assertEqual(0.17, cnf)
        # Call this again to check we didn't mess with the cache
        cmf = geology_model.convert_core_number_fraction_to_core_mass_fraction(0.17)
        self.assertEqual(0.33229859947479184, cmf)

    def test_lin_log_conversions(self):
        geology_model = gi.GeologyModel()

        linear_composition = {
            ci.Element.Al: {
                gi.Layer.bulk: 0.01908893026850964
            },
            ci.Element.Ti: {
                gi.Layer.bulk: 0.0005489626996871848
            },
            ci.Element.Ca: {
                gi.Layer.bulk: 0.013847584102512443
            },
            ci.Element.Ni: {
                gi.Layer.bulk: 0.0008103772557189963
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.023816218523713337
            },
            ci.Element.Cr: {
                gi.Layer.bulk: 0.0012537023978132916
            },
            ci.Element.Mg: {
                gi.Layer.bulk: 0.2056364370557375
            },
            ci.Element.Si: {
                gi.Layer.bulk: 0.154912113261916
            },
            ci.Element.Na: {
                gi.Layer.bulk: 0.002541447773895965
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.5774849021972671
            },
            ci.Element.C: {
                gi.Layer.bulk: 1.3976695399380497e-06
            },
            ci.Element.N: {
                gi.Layer.bulk: 5.792679368858128e-05
            }
        }
        log_composition = {
            ci.Element.Al: -5.719438886,
            ci.Element.Ti: -7.260677641,
            ci.Element.Ca: -5.858846466,
            ci.Element.Ni: -7.091533234,
            ci.Element.Fe: -5.623347671,
            ci.Element.Cr: -6.902026021,
            ci.Element.Mg: -4.687120407,
            ci.Element.Si: -4.810135099,
            ci.Element.Na: -6.595139288,
            ci.Element.O: -4.238679843,
            ci.Element.C: -9.854815977,
            ci.Element.N: -8.237340987
        }
        test_log_composition = geology_model.convert_linear_abundances_to_log(linear_composition, -4)
        for el, reference in log_composition.items():
            self.assertAlmostEqual(test_log_composition[el], reference, 3)
        test_linear_composition = geology_model.convert_log_abundances_to_linear(log_composition)
        for el, reference in linear_composition.items():
            self.assertAlmostEqual(test_linear_composition[el][gi.Layer.bulk], reference[gi.Layer.bulk])

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
            mp.ModelParameter.fragment_mass,
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
                mp.ModelParameter.fragment_mass: False,
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
                mp.ModelParameter.fragment_mass: True,
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
                mp.ModelParameter.fragment_mass: True,
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
                mp.ModelParameter.fragment_mass: 8,
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
                mp.ModelParameter.fragment_mass: 8,
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
        wd_timescales = np.array([50000, 60000, 70000, 80000])
        offset = 6.845098040014257 - 6.6020599913279625
        args = {
            'Early': [0, 200000, planetesimal_abundance, wd_timescales, -6],
			'Medium': [0.1, 200000, planetesimal_abundance, wd_timescales, -6],
			'Late': [0.4, 200000, planetesimal_abundance, wd_timescales, -6],
			'Very Late': [1, 200000, planetesimal_abundance, wd_timescales, -6]
        }
        expected_results = {
            'Early': np.array([-6.845098040014257, -6.669006780958576, -6.544068044350276, -6.447158031342219]),
            'Medium': np.array([-6.926929795486369, -6.69941756364218, -6.535606543510326, -6.408326484132252]),
            'Late': np.array([-7.471824349177887, -6.93477018262815, -6.54598923340804, -6.24749203990456]),
            'Very Late': np.array([-9.25412736937022, -7.848484239013979, -6.839282601360653, -6.075469891532259])
        }
        expected_lifetime_results = {
            'Early': np.array([-6.698970004336019, -6.522878745280337, -6.3979400086720375, -6.301029995663981]),
            'Medium': np.array([-6.369503928147244, -6.254407893195727, -6.175686043374323, -6.114496075232561]),
            'Late': np.array([-4.930950375586853, -5.025940353021858, -5.089882029233725, -5.131808928937756]),
            'Very Late': np.array([0.24974190250473516, -0.7322832835445023, -1.4345607295167415, -1.9580342688029324])
        }
        for arg_name, arg_set in args.items():
            result = wdm.process_abundances(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4])
            expected_result = expected_results[arg_name]
            self.assertEqual(len(result), len(expected_result))
            self.assertTrue((result == expected_result).all())
            lifetime_result = wdm.process_abundances(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], False)
            expected_lifetime_result = expected_lifetime_results[arg_name]
            self.assertEqual(len(lifetime_result), len(expected_lifetime_result))
            for i, val in enumerate(lifetime_result):
                self.assertEqual(val, expected_lifetime_result[i])

class ModelAnalyserTests(unittest.TestCase):

    def test_confidence_intervals(self):
        model_analyser = ma.ModelAnalyser(None)
        input_array = np.zeros(shape=(3, 2))  # 3 samples, 2 propertes we want confidence intervals for. The first has values 1, 3, 2 and the second has 2, 6, 7
        input_array[0,:] = [1, 2]
        input_array[1,:] = [3, 6]
        input_array[2,:] = [2, 7]
        arr_low3, arr_low2, arr_low1, arr_median, arr_high1, arr_high2, arr_high3 = model_analyser.confidence_intervals(input_array)
        self.assertEqual(2, arr_median[0]) # this should surely be [2, 6]
        self.assertEqual(6, arr_median[1]) # this should surely be [2, 6]
        self.assertAlmostEqual(1.003, arr_low3[0])
        self.assertAlmostEqual(2.012, arr_low3[1])
        self.assertAlmostEqual(1.046, arr_low2[0])
        self.assertAlmostEqual(2.184, arr_low2[1])
        self.assertAlmostEqual(1.3174, arr_low1[0])
        self.assertAlmostEqual(3.2696, arr_low1[1])
        self.assertAlmostEqual(2.6826, arr_high1[0])
        self.assertAlmostEqual(6.6826, arr_high1[1])
        self.assertAlmostEqual(2.954, arr_high2[0])
        self.assertAlmostEqual(6.954, arr_high2[1])
        self.assertAlmostEqual(2.997, arr_high3[0])
        self.assertAlmostEqual(6.997, arr_high3[1])
        input_array_1d = [4.5, 4.2, np.nan, None, 4.2]
        arr_low3, arr_low2, arr_low1, arr_median, arr_high1, arr_high2, arr_high3 = model_analyser.confidence_intervals(input_array_1d)
        self.assertAlmostEqual(4.2, arr_low3)
        self.assertAlmostEqual(4.2, arr_low2)
        self.assertAlmostEqual(4.2, arr_low1)
        self.assertAlmostEqual(4.2, arr_median)
        self.assertAlmostEqual(4.40478, arr_high1)
        self.assertAlmostEqual(4.4862, arr_high2)
        self.assertAlmostEqual(4.4991, arr_high3)

class SyntheticSystemTests(unittest.TestCase):

    #TODO: This should be split into smaller tests
    def test_synthetic_system(self):
        test_system_none = sp.SyntheticSystem(None, None, None)
        test_system_none_str = "\nWD [No ID]:\n[No WD properties]\nPollution properties:\n[No Pollution properties]\nPollution:\n[No Pollution values]\nObserved:\nNone\nObserved Abundances:\n[No observed abundance values]\nBandpass Magnitudes:\n[No bandpass magnitude values]\nModelled:\n[No Modelled properties]\n"
        self.assertEqual(test_system_none_str, str(test_system_none))
        wd_properties = {
            mp.WDParameter.logg: 8
        }
        pollution_properties = {
            mp.ModelParameter.metallicity: 100
        }
        pollution_abundances = {
            ci.Element.Fe: -9
        }
        bandpass_magnitudes = {
            sb.Bandpass.u: 20.1
        }
        modelled_parameters = {
            mp.ModelParameter.metallicity: 90
        }
        test_system = sp.SyntheticSystem(wd_properties, pollution_properties, pollution_abundances)
        test_system_str = "\nWD [No ID]:\nLogg: 8\nPollution properties:\nMetallicity: 100\nPollution:\nFe: -9\nObserved:\nNone\nObserved Abundances:\n[No observed abundance values]\nBandpass Magnitudes:\n[No bandpass magnitude values]\nModelled:\n[No Modelled properties]\n"
        self.assertEqual(test_system_str, str(test_system))
        self.assertEqual(['None', '8', '-9'], test_system.to_csv_row([mp.WDParameter.logg], [], [ci.Element.Fe], False, [], [], []))
        self.assertEqual(['None', '8', '100', '-9', 'None', 'None', 'None', 'None'], test_system.to_csv_row([mp.WDParameter.logg], [mp.ModelParameter.metallicity], [ci.Element.Fe], True, [ci.Element.Fe], [sb.Bandpass.u], [mp.ModelParameter.metallicity]))
        test_system.set_modelled_properties(modelled_parameters)
        test_system_str = "\nWD [No ID]:\nLogg: 8\nPollution properties:\nMetallicity: 100\nPollution:\nFe: -9\nObserved:\nNone\nObserved Abundances:\n[No observed abundance values]\nBandpass Magnitudes:\n[No bandpass magnitude values]\nModelled:\nMetallicity: 90\n"
        self.assertEqual(['None', 'None', '-9', 'None'], test_system.to_csv_row([], ['Dummy key'], [ci.Element.Fe], True, [], [], []))
        self.assertEqual(['None', '8', '100', '-9', 'None', 'None', '90'], test_system.to_csv_row([mp.WDParameter.logg], [mp.ModelParameter.metallicity], [ci.Element.Fe], False, [ci.Element.Fe], [sb.Bandpass.u], [mp.ModelParameter.metallicity]))
        self.assertEqual(test_system_str, str(test_system))
        self.assertFalse(test_system.check_elements_detected([ci.Element.Fe]))
        test_system.set_bandpass_magnitudes(bandpass_magnitudes)
        test_system_str = "\nWD [No ID]:\nLogg: 8\nPollution properties:\nMetallicity: 100\nPollution:\nFe: -9\nObserved:\nNone\nObserved Abundances:\n[No observed abundance values]\nBandpass Magnitudes:\nu: 20.1\nModelled:\nMetallicity: 90\n"
        self.assertEqual(['None', '8', '100', '-9', 'None', 'None', '20.1', '90'], test_system.to_csv_row([mp.WDParameter.logg], [mp.ModelParameter.metallicity], [ci.Element.Fe], True, [ci.Element.Fe], [sb.Bandpass.u], [mp.ModelParameter.metallicity]))
        self.assertEqual(test_system_str, str(test_system))
        observed = True
        observed_abundances = {
            ci.Element.Fe: -9.1
        }
        test_system.observed = True
        test_system.set_observed_abundances(observed_abundances)
        self.assertTrue(test_system.check_elements_detected([ci.Element.Fe]))
        self.assertFalse(test_system.check_elements_detected([ci.Element.Ca, ci.Element.Zn]))
        test_system_str = "\nWD [No ID]:\nLogg: 8\nPollution properties:\nMetallicity: 100\nPollution:\nFe: -9\nObserved:\nTrue\nObserved Abundances:\nFe: -9.1\nBandpass Magnitudes:\nu: 20.1\nModelled:\nMetallicity: 90\n"
        self.assertEqual(['None', '8', '100', '-9', 'True', '-9.1', '20.1', '90'], test_system.to_csv_row([mp.WDParameter.logg], [mp.ModelParameter.metallicity], [ci.Element.Fe], True, [ci.Element.Fe], [sb.Bandpass.u], [mp.ModelParameter.metallicity]))
        self.assertEqual(test_system_str, str(test_system))
        test_system_2 = sp.SyntheticSystem(wd_properties, pollution_properties, pollution_abundances, observed, observed_abundances, bandpass_magnitudes, None, 1)
        self.assertEqual(True, test_system != test_system_2)
        self.assertEqual(False, test_system == test_system_2)
        test_system_2.set_modelled_properties(modelled_parameters)
        self.assertEqual(True, test_system == test_system_2)
        self.assertEqual(False, test_system != test_system_2)
        self.assertEqual(test_system, test_system_2)
        test_system_2_str = "\nWD 1:\nLogg: 8\nPollution properties:\nMetallicity: 100\nPollution:\nFe: -9\nObserved:\nTrue\nObserved Abundances:\nFe: -9.1\nBandpass Magnitudes:\nu: 20.1\nModelled:\nMetallicity: 90\n"
        self.assertEqual(['1', '8', '100', '-9', 'True', '-9.1', '20.1', '90'], test_system_2.to_csv_row([mp.WDParameter.logg], [mp.ModelParameter.metallicity], [ci.Element.Fe], True, [ci.Element.Fe], [sb.Bandpass.u], [mp.ModelParameter.metallicity]))
        self.assertEqual(test_system_2_str, str(test_system_2))
        test_system.id = 1
        test_system_2.id = 3
        self.assertEqual(True, test_system == test_system_2)
        self.assertEqual(False, test_system != test_system_2)
        self.assertEqual(test_system, test_system_2)

    def test_distance_to(self):
        abundance_set1 = {
            ci.Element.Ca: -6,
            ci.Element.Mg: -6,
            ci.Element.Fe: -6
        }
        abundance_set2 = {
            ci.Element.Ca: -7,
            ci.Element.Mg: -7,
            ci.Element.Fe: -7
        }
        abundance_set3 = {
            ci.Element.Ca: -7.3,
            ci.Element.Mg: -6.6,
            ci.Element.Fe: -7.5
        }
        test_system1 = sp.SyntheticSystem(None, None, None, True, abundance_set1)
        test_system2 = sp.SyntheticSystem(None, None, None, True, abundance_set2)
        test_system3 = sp.SyntheticSystem(None, None, None, True, abundance_set3)
        self.assertEqual(0, test_system1.distance_to(test_system1))
        self.assertEqual(0, test_system1.distance_to(test_system2))
        self.assertEqual(test_system1.distance_to(test_system3), test_system2.distance_to(test_system3))
        self.assertEqual(0.7364260006390226, test_system1.distance_to(test_system3))
        self.assertEqual(0.7364260006390226, test_system3.distance_to(test_system1))

class SyntheticConfigurationsTests(unittest.TestCase):

    def test_predefined_distributions(self):
        #Basically just want to test that all the distributions got loaded directly, but there's so many that I'll just check one index for each key!
        # Want to catch issues such as a wrong column index leading to everything being nan. If one index is correct they probably all are
        expected_bin_index_5_value = {
            'combined': (0.0641025641025641, 0.07692307692307693),
            'combined_coarse': (0.1282051282051282, 0.15384615384615385),
            'Mdot_7': (0.0641025641025641, 0.07692307692307693),
            'Mdot_7_coarse': (0.1282051282051282, 0.15384615384615385),
            'Mdot_8': (0.0641025641025641, 0.07692307692307693),
            'Mdot_8_coarse': (0.1282051282051282, 0.15384615384615385),
            'Mdot_10': (0.0641025641025641, 0.07692307692307693),
            'Mdot_10_coarse': (0.1282051282051282, 0.15384615384615385),
            'CMF_tidal_disruption': (0.04115384615384615, 0.04938461538461539),
            'CMF_tidal_disruption_coarse': (0.1282051282051282, 0.15384615384615385),
            'm_cf_035f6nogas': (0.06329113924050633, 0.0759493670886076),
            'm_cf_035f6nogas_coarse': (0.1282051282051282, 0.15384615384615385),
            'MWDD_DA_Teffs_40pc': (13000.0, 14800.0),
            'MWDD_DA_Loggs_40pc': (7.8, 7.96),
            'MWDD_DB_Teffs_40pc': (25200.0, 30000.0),
            'MWDD_DB_Loggs_40pc': (8.77, 9.0),
            'HollandsTeffs': (6500.0, 7000.0),
            'HollandsLoggs': (8.1, 8.3),
            'NS20_Teffs_DA': (9111.111111111111, 9333.333333333334),
            'NS20_Loggs_DA': (7.101010101010101, 7.121212121212121),
            'NS20_Teffs_DB': (12040.816326530612, 12448.979591836734),
            'NS20_Loggs_DB': (7.204081632653061, 7.244897959183674),
            'GSPCWD_Teffs_DA': (4850.0, 5220.0),
            'GSPCWD_Loggs_DA': (7.1, 7.12),
            'GSPCWD_Teffs_DB': (4850.0, 5220.0),
            'GSPCWD_Loggs_DB': (7.1, 7.12)

        }
        expected_count_index_5_value = {
            'combined': 31293340.540264696,
            'combined_coarse': 570725773.8724912,
            'Mdot_7': 7705140.855275676,
            'Mdot_7_coarse': 176137327.34542602,
            'Mdot_8': 1198582.624243103,
            'Mdot_8_coarse': 16198184.512538336,
            'Mdot_10': 5.889886479893114,
            'Mdot_10_coarse': 288788.32326637954,
            'CMF_tidal_disruption': 495,
            'CMF_tidal_disruption_coarse': 3000,
            'm_cf_035f6nogas': 14,
            'm_cf_035f6nogas_coarse': 381,
            'MWDD_DA_Teffs_40pc': 24,
            'MWDD_DA_Loggs_40pc': 91,
            'MWDD_DB_Teffs_40pc': 1,
            'MWDD_DB_Loggs_40pc': 2,
            'HollandsTeffs': 39,
            'HollandsLoggs': 37,
            'NS20_Teffs_DA': 5,
            'NS20_Loggs_DA': 42,
            'NS20_Teffs_DB': 332,
            'NS20_Loggs_DB': 16,
            'GSPCWD_Teffs_DA': 1598,
            'GSPCWD_Loggs_DA': 147,
            'GSPCWD_Teffs_DB': 87,
            'GSPCWD_Loggs_DB': 12
        }
        for config_name, dist in sc.predefined_distributions.items():
            if config_name in ['MWDD_DB_Teffs_40pc', 'MWDD_DB_Loggs_40pc']:
                # Not enough indices...
                self.assertAlmostEqual(expected_bin_index_5_value[config_name][0], dist[0][4][0])
                self.assertAlmostEqual(expected_bin_index_5_value[config_name][1], dist[0][4][1])
                self.assertAlmostEqual(expected_count_index_5_value[config_name], dist[1][4])
            else:
                self.assertAlmostEqual(expected_bin_index_5_value[config_name][0], dist[0][5][0])
                self.assertAlmostEqual(expected_bin_index_5_value[config_name][1], dist[0][5][1])
                self.assertAlmostEqual(expected_count_index_5_value[config_name], dist[1][5])

class SyntheticPopulationTests(unittest.TestCase):

    def write_golden_csvs(self):
        number_of_systems = 3
        with open('test_synth_dump_v1_golden.csv', 'w', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow(['ID','WD Spectral Type','WD Temperature','WD Log(g)','WD Mass','WD Distance','WD Timescale Type','WD Consider Thermohaline','True log(Al/Hx)','True log(Ti/Hx)','True log(Ca/Hx)','True log(Ni/Hx)','True log(Fe/Hx)','True log(Cr/Hx)','True log(Mg/Hx)','True log(Si/Hx)','True log(Na/Hx)','True log(O/Hx)','True log(C/Hx)','True log(N/Hx)'])
            i = 0
            while i < number_of_systems:
                to_write.writerow([str(i),'DB','5000','8','0.6','40','KoesterOvershoot','False','-7.637195577826703','-9.048070253524186','-7.6402695555270945','-7.826546928534431','-6.558848354928465',
                '-8.353997596544648','-6.369790844454792','-6.375223341025334','-7.688587430031069','-5.763919539645433','-inf','-inf'])
                i += 1
        with open('test_synth_dump_v2_golden.csv', 'w', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow(['ID','WD Spectral Type','WD Temperature','WD Log(g)','WD Mass','WD Distance','WD Timescale Type','WD Consider Thermohaline','Input Metallicity','Input Time Since Accretion','Input Formation Distance','Input Feeding Zone Size','Input Parent Core Number Fr\
action','Input Parent Crust Number Fraction','Input Fragment Core Number Fraction','Input Fragment Crust Number Fraction','In\
put Fragment Mass','Input Accretion Timescale','Input Pressure','Input Oxygen Fugacit\
y','True log(Al/Hx)','True log(Ti/Hx)','True log(Ca/Hx)','True log(Ni/Hx)','True log(Fe/Hx)','True log(Cr/Hx)','True log(Mg/Hx)','True l\
og(Si/Hx)','True log(Na/Hx)','True log(O/Hx)','True log(C/Hx)','True log(N/Hx)','Observed?','Observed log(Al/Hx)','Observed log(Ti/Hx)','Observed log(\
Ca/Hx)','Observed log(Ni/Hx)','Observed log(Fe/Hx)','Observed log(Cr/Hx)','Observed log(Mg/Hx)','Observed l\
og(Si/Hx)','Observed log(Na/Hx)','Observed log(O/Hx)','Observed log(C/Hx)','Observed log(N/Hx)','u','g','r','i','z','Output Metallicity','Output Time Since \
Accretion','Output Formation Distance','Output Feeding Zone Size','Output Parent C\
ore Number Fraction','Output Parent Crust Number Fraction','Output Fragment Core Number Fraction','Output Fragment Crus\
t Number Fraction','Output Fragment Mass','Output Accretion Timescale','Output Pressure','Output O\
xygen Fugacity'])
            i = 0
            while i < number_of_systems:
                to_write.writerow([str(i),'DB','5000','8','0.6','40','KoesterOvershoot','False','400','1.2','0.0','0.03','None','None','0.1','None','20.0','1000000','10','-2',
                '-7.637195577826703','-9.048070253524186','-7.6402695555270945','-7.826546928534431','-6.558848354928465',
                '-8.353997596544648','-6.369790844454792','-6.375223341025334','-7.688587430031069','-5.763919539645433','-inf','-inf',
                'None','None','None','None','None','None','None','None','None','None','None','None','None',
                'None','None','None','None','None','None','None','None','None','None','None','None','None','None','None','None','None'])
                i += 1
        with open('test_synth_dump_v3_golden.csv', 'w', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow(['ID','WD Spectral Type','WD Temperature','WD Log(g)','WD Mass','WD Distance','WD Timescale Type','WD Consider Thermohaline','Input Metallicity',\
            'Input Time Since Accretion','Input Formation Distance','Input Feeding Zone Size','Input Parent Core Number Fr\
action','Input Parent Crust Number Fraction','Input Fragment Core Number Fraction','Input Fragment Crust Number Fraction','In\
put Fragment Mass','Input Accretion Timescale','Input Pressure','Input Oxygen Fugacit\
y','True log(Al/Hx)','True log(Ti/Hx)','True log(Ca/Hx)','True log(Ni/Hx)','True log(Fe/Hx)','True log(Cr/Hx)','True log(Mg/Hx)','True l\
og(Si/Hx)','True log(Na/Hx)','True log(O/Hx)','True log(C/Hx)','True log(N/Hx)','Observed?','Observed log(Al/Hx)','Observed log(Ti/Hx)','Observed log(\
Ca/Hx)','Observed log(Ni/Hx)','Observed log(Fe/Hx)','Observed log(Cr/Hx)','Observed log(Mg/Hx)','Observed l\
og(Si/Hx)','Observed log(Na/Hx)','Observed log(O/Hx)','Observed log(C/Hx)','Observed log(N/Hx)','u','g','r','i','z','Output Metallicity','Output Time S\
ince Accretion','Output Formation Distance','Output Feeding Zone Size','Output Parent C\
ore Number Fraction','Output Parent Crust Number Fraction','Output Fragment Core Number Fraction','Output Fragment Crus\
t Number Fraction','Output Fragment Mass','Output Accretion Timescale','Output Pressure','Output O\
xygen Fugacity'])
            i = 0
            while i < number_of_systems:
                id_no = i if i != 1 else 123
                to_write.writerow([str(id_no),'DB','5000','8','0.6','40','KoesterOvershoot','False','400','1.2','0.0','0.03','None','None','0.1','None','20.0',
                '1000000','10','-2','-7.637195577826703','-9.048070253524186','-7.6402695555270945','-7.826546928534431','-6.558848354928465',
                '-8.353997596544648','-6.369790844454792','-6.375223341025334','-7.688587430031069','-5.763919539645433','-inf','-inf',
                'True','None','None','None','None','-8.7','None','None','None','None','None','None','None','None','None','None','None','None',
                '300','None','None','None','None','None','None','None','None','None','None','None'])
                i += 1

    def test_setup(self):
        for file_which_will_be_affected in ['test_synth_dump_v1.csv', 'test_synth_dump_v1_golden.csv', 'test_synth_dump_v2.csv', 'test_synth_dump_v2_golden.csv', 'test_synth_dump_v3.csv', 'test_synth_dump_v3_golden.csv', 'test_synth_dump_v4.csv']:
            if os.path.exists(file_which_will_be_affected):
                raise IOError('This test deletes file ' + file_which_will_be_affected + ', cannot proceed until it is removed or renamed')
        self.write_golden_csvs()
        with self.assertRaises(ValueError):
            synthetic_pop = sp.SyntheticPopulation(3, 'NonExistantWDConfig', 'TestPollutionConfig')
        with self.assertRaises(ValueError):
            synthetic_pop = sp.SyntheticPopulation(3, 'TestWDConfigDB', 'NonExistantPollutionConfig')
        number_of_systems = 3
        synthetic_pop = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig')
        self.assertEqual(number_of_systems, len(synthetic_pop))
        system_count_manual = 0
        for system in synthetic_pop:
            system_count_manual += 1
        self.assertEqual(number_of_systems, system_count_manual)
        expected_logg = [8, 8, 8]
        expected_metallicity = [400, 400, 400]
        expected_Mo = [None, None, None]
        expected_Fe = [-6.558848354928465, -6.558848354928465, -6.558848354928465]
        expected_CaFe = [-1.08142120059863, -1.08142120059863, -1.08142120059863]
        expected_modelled_metallicity = [None, None, None]
        self.assertEqual(expected_logg, synthetic_pop.wd_values(mp.WDParameter.logg))
        self.assertEqual(expected_metallicity, synthetic_pop.pollution_input_values(mp.ModelParameter.metallicity))
        print(synthetic_pop)
        synthetic_pop.dump_to_csv('test_synth_dump_v1.csv', True)
        synthetic_pop.dump_to_csv('test_synth_dump_v2.csv')
        self.assertEqual(expected_Fe, synthetic_pop.pollution_abundance_values(ci.Element.Fe))
        self.assertEqual(expected_Mo, synthetic_pop.pollution_abundance_values(ci.Element.Mo))
        self.assertEqual(expected_CaFe, synthetic_pop.pollution_abundance_log_ratios(ci.Element.Ca, ci.Element.Fe))
        self.assertEqual(expected_modelled_metallicity, synthetic_pop.modelled_values(mp.ModelParameter.metallicity))
        for system in synthetic_pop:
            system.observed = True
            system.set_observed_abundances({ci.Element.Fe: -8.7})
            system.set_modelled_properties({mp.ModelParameter.metallicity: 300})
        synthetic_pop.population[1].id = 123
        synthetic_pop.dump_to_csv('test_synth_dump_v3.csv')
        expected_modelled_metallicity = [300, 300, 300]
        self.assertEqual(expected_modelled_metallicity, synthetic_pop.modelled_values(mp.ModelParameter.metallicity))
        expected_observed_Fe = [-8.7, -8.7, -8.7]
        self.assertEqual(expected_observed_Fe, synthetic_pop.observed_abundances(ci.Element.Fe))
        synthetic_pop2 = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig')
        self.assertEqual(False, synthetic_pop2 == synthetic_pop)
        for system in synthetic_pop2:
            system.observed = True
            system.set_observed_abundances({ci.Element.Fe: -8.7})
            system.set_modelled_properties({mp.ModelParameter.metallicity: 300})
        self.assertEqual(True, synthetic_pop2 == synthetic_pop)
        with open('test_synth_dump_v1.csv') as f1:
            with open('test_synth_dump_v1_golden.csv') as f2:
                self.assertTrue(f1.read() == f2.read())
        with open('test_synth_dump_v2.csv') as f1:
            with open('test_synth_dump_v2_golden.csv') as f2:
                self.assertTrue(f1.read() == f2.read())
        with open('test_synth_dump_v3.csv') as f1:
            with open('test_synth_dump_v3_golden.csv') as f2:
                self.assertTrue(f1.read() == f2.read())
        self.assertEqual(mp.WDParameter.spectral_type, synthetic_pop.convert_readable_str_to_property('WD Spectral Type'))
        self.assertEqual(mp.WDParameter.logg, synthetic_pop.convert_readable_str_to_property('WD Log(g)'))
        self.assertEqual(mp.ModelParameter.t_sinceaccretion, synthetic_pop.convert_readable_str_to_property('Input Time Since Accretion'))
        self.assertEqual(mp.ModelParameter.fragment_core_frac, synthetic_pop.convert_readable_str_to_property('Output Fragment Core Fraction'))
        self.assertEqual(mp.ModelParameter.formation_distance, synthetic_pop.convert_readable_str_to_property('Output Formation Distance'))
        self.assertEqual(ci.Element.Si, synthetic_pop.convert_readable_str_to_property('True log(Si/Hx)'))
        self.assertEqual(ci.Element.O, synthetic_pop.convert_readable_str_to_property('Observed log(O/Hx)'))
        synthetic_pop3 = sp.SyntheticPopulation(None, None, None, 'test_synth_dump_v3_golden.csv')
        matching_systems = synthetic_pop3.find_systems_with_abundances({ci.Element.Ti: -9.048070253524186, ci.Element.Ca: -7.6402695555270945})
        self.assertEqual(len(synthetic_pop3), len(matching_systems))
        matching_systems = synthetic_pop3.find_systems_with_abundances({ci.Element.Ti: -10.57, ci.Element.Ca: -9.166481121051756})
        self.assertEqual(0, len(matching_systems))
        self.assertEqual(True, synthetic_pop3 == synthetic_pop)
        self.assertEqual(123, synthetic_pop3[1].id) # Throwing in a subscript test
        self.assertEqual(124, synthetic_pop.get_next_id())
        self.assertEqual(3, synthetic_pop2.get_next_id())
        self.assertEqual(124, synthetic_pop3.get_next_id())
        test_new_system = synthetic_pop.create_system()
        self.assertEqual(124, test_new_system.id)
        self.assertEqual(124, synthetic_pop.get_next_id())
        self.assertEqual(None, synthetic_pop.find_systems_with_id(321))
        self.assertEqual(synthetic_pop.population[1], synthetic_pop.find_systems_with_id(123))
        synthetic_pop.population[2].id = 123
        self.assertEqual(synthetic_pop.population[1], synthetic_pop.find_systems_with_id(123))
        test_new_system.observed = False
        synthetic_pop.population.append(test_new_system)
        self.assertEqual(False, synthetic_pop3 == synthetic_pop)
        self.assertEqual(True, synthetic_pop3 == synthetic_pop.get_observed_subset())
        # Now check the population caching logic
        # These 2 pops should be identical despite the random sampling because synthetic_pop4 got dumped
        synthetic_pop4 = sp.SyntheticPopulation(number_of_systems, 'UniformTeff', 'UniformFDT', 'test_synth_dump_v4.csv')
        synthetic_pop4.dump_to_csv()
        synthetic_pop5 = sp.SyntheticPopulation(number_of_systems, 'UniformTeff', 'UniformFDT', 'test_synth_dump_v4.csv')
        self.assertEqual(True, synthetic_pop4 == synthetic_pop5)

        os.remove('test_synth_dump_v4.csv')
        modeller = sm.Modeller(sm.ModellerType.SimpleFcfInterpolation)
        for system in synthetic_pop4:
            system.observed = True
            system.set_observed_abundances(
                {
                    ci.Element.Al: system.pollution_abundances[ci.Element.Al] + 0.01,
                    ci.Element.Ca: system.pollution_abundances[ci.Element.Ca] + 0.015,
                    ci.Element.Mg: system.pollution_abundances[ci.Element.Mg] - 0.01,
                    ci.Element.Fe: system.pollution_abundances[ci.Element.Fe] - 0.02,
                }
            )
        modeller.model_population(synthetic_pop4)  # This automatically dumps the population, given that test_synth_dump_v4.csv is not present

        synthetic_pop6 = sp.SyntheticPopulation(number_of_systems, 'UniformTeff', 'UniformFDT', 'test_synth_dump_v4.csv')
        self.assertEqual(False, synthetic_pop4 == synthetic_pop5)
        self.assertEqual(True, synthetic_pop4 == synthetic_pop6)

        #Test what happens if you request more systems than are present in a csv
        synthetic_pop7 = sp.SyntheticPopulation(number_of_systems + 1, 'UniformTeff', 'UniformFDT', 'test_synth_dump_v4.csv')
        self.assertEqual(number_of_systems + 1, len(synthetic_pop7))

        os.remove('test_synth_dump_v1.csv')
        os.remove('test_synth_dump_v1_golden.csv')
        os.remove('test_synth_dump_v2.csv')
        os.remove('test_synth_dump_v2_golden.csv')
        os.remove('test_synth_dump_v3.csv')
        os.remove('test_synth_dump_v3_golden.csv')
        os.remove('test_synth_dump_v4.csv')

    def test_sample_arbitrary_function(self):
        synthetic_pop = sp.SyntheticPopulation(0, 'TestWDConfigDB', 'TestPollutionConfig') #This doesn't really matter
        def function_to_sample(x):
            if x > 10 or x < 2:
                return 0
            if x > 7:
                return 20
            else:
                return 4*(x-2)
        synthetic_pop.cache_inverse_cdf('test_func', function_to_sample, 2, 10)
        self.assertAlmostEqual(6, synthetic_pop.inverse_cdf_dict['test_func'](16/55))
        self.assertAlmostEqual(8.9, synthetic_pop.inverse_cdf_dict['test_func'](0.8))
        sample_size = 100
        i = 0
        while i < sample_size:
            ran_val = synthetic_pop.draw_variable_from_distribution((sc.Distribution.CustomFunction, ['test_func', function_to_sample, 2, 10]))
            self.assertTrue(ran_val >= 2)
            self.assertTrue(ran_val <= 10)
            i += 1

    def test_sample_arbitrary_distribution(self):
        synthetic_pop = sp.SyntheticPopulation(0, 'TestWDConfigDB', 'TestPollutionConfig') #This doesn't really matter
        bins_to_sample = [(3, 4), (4, 5), (5, 7)]
        counts_to_sample = [40, 60, 100]
        synthetic_pop.cache_inverse_cdf_from_table('test_dist', bins_to_sample, counts_to_sample)
        self.assertEqual(3, synthetic_pop.inverse_cdf_dict['test_dist'](0))
        self.assertEqual(7, synthetic_pop.inverse_cdf_dict['test_dist'](1))
        self.assertEqual(5, synthetic_pop.inverse_cdf_dict['test_dist'](0.5))
        self.assertEqual(3.5, synthetic_pop.inverse_cdf_dict['test_dist'](0.1))
        self.assertEqual(6, synthetic_pop.inverse_cdf_dict['test_dist'](0.75))
        sample_size = 100
        i = 0
        while i < sample_size:
            ran_val = synthetic_pop.draw_variable_from_distribution((sc.Distribution.CustomDistribution, ['test_dist', bins_to_sample, counts_to_sample]))
            self.assertTrue(ran_val >= 3)
            self.assertTrue(ran_val <= 7)
            i += 1

    def test_draw_variable_from_distribution(self):
        synthetic_pop = sp.SyntheticPopulation(1, 'TestWDConfigDB', 'TestPollutionConfig')
        with self.assertRaises(ValueError):
            value = synthetic_pop.draw_variable_from_distribution('NonExistant')
        uniform_distribution = (sc.Distribution.Uniform, [0, 1])
        normal_distribution = (sc.Distribution.Normal, [0.5, 0.001])
        delta_distribution = (sc.Distribution.Delta, [0, 1])
        triangle_distribution = (sc.Distribution.Triangle, [0, 0.7, 1])
        slope_distribution = (sc.Distribution.Slope, [0, 0.7, 1])
        coll_casc_distribution = (sc.Distribution.CustomFunction, ['CollisionalCascade', sc.predefined_functions['CollisionalCascade'], 0.1, 1])
        custom_distribution = (sc.Distribution.CustomDistribution, ['CD', [(0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0)], np.array([1, 4, 6, 5, 3])])
        sample_size = 1000
        distribution_samples = {
            sc.Distribution.Uniform: list(),
            sc.Distribution.Normal: list(),
            sc.Distribution.Delta: list(),
            sc.Distribution.Triangle: list(),
            sc.Distribution.Slope: list(),
            sc.Distribution.CustomFunction: list(),
            sc.Distribution.CustomDistribution: list()
        }
        i = 0
        while i < sample_size:
            distribution_samples[sc.Distribution.Uniform].append(synthetic_pop.draw_variable_from_distribution(uniform_distribution))
            distribution_samples[sc.Distribution.Normal].append(synthetic_pop.draw_variable_from_distribution(normal_distribution))
            distribution_samples[sc.Distribution.Delta].append(synthetic_pop.draw_variable_from_distribution(delta_distribution))
            distribution_samples[sc.Distribution.Triangle].append(synthetic_pop.draw_variable_from_distribution(triangle_distribution))
            distribution_samples[sc.Distribution.Slope].append(synthetic_pop.draw_variable_from_distribution(slope_distribution))
            distribution_samples[sc.Distribution.CustomFunction].append(synthetic_pop.draw_variable_from_distribution(coll_casc_distribution))
            distribution_samples[sc.Distribution.CustomDistribution].append(synthetic_pop.draw_variable_from_distribution(custom_distribution))
            i += 1
        for dist, vals in distribution_samples.items():
            self.assertEqual(sample_size, len(vals))
            for val in vals:
                self.assertTrue(0 <= val <= 1)

    def test_draw_variable_from_custom_distribution(self):
        synthetic_pop = sp.SyntheticPopulation(1, 'TestWDConfigDB', 'TestPollutionConfig')
        a = synthetic_pop.draw_variable_from_distribution([sc.Distribution.CustomDistribution, ['Combined', sc.predefined_distributions['combined'][0], sc.predefined_distributions['combined'][1]]])
        self.assertTrue(isinstance(a, float))
        self.assertFalse(isinstance(a, np.ndarray))
        self.assertTrue(0 <= a <= 1)
        b = synthetic_pop.draw_variable_from_distribution([sc.Distribution.CustomDistribution, ['MWDD_DA_Loggs_40pc', sc.predefined_distributions['MWDD_DA_Loggs_40pc'][0], sc.predefined_distributions['MWDD_DA_Loggs_40pc'][1]]])
        self.assertTrue(isinstance(b, float))
        self.assertFalse(isinstance(b, np.ndarray))
        self.assertTrue(7 <= b <= 9)

    def test_cast_string_to_float_or_int(self):
        synthetic_pop = sp.SyntheticPopulation(1, 'TestWDConfigDB', 'TestPollutionConfig')
        self.assertEqual(0.65, synthetic_pop.cast_string_to_float_or_int('0.65'))
        self.assertEqual(0.65, synthetic_pop.cast_string_to_float_or_int("0.65"))
        self.assertEqual(6, synthetic_pop.cast_string_to_float_or_int('6'))
        self.assertEqual(6, synthetic_pop.cast_string_to_float_or_int("6"))
        self.assertEqual(None, synthetic_pop.cast_string_to_float_or_int('None'))
        self.assertEqual(None, synthetic_pop.cast_string_to_float_or_int("None"))
        self.assertTrue(np.isnan(synthetic_pop.cast_string_to_float_or_int('nan')))
        self.assertTrue(np.isnan(synthetic_pop.cast_string_to_float_or_int("nan")))
        self.assertTrue(np.isnan(synthetic_pop.cast_string_to_float_or_int('np.nan')))
        self.assertTrue(np.isnan(synthetic_pop.cast_string_to_float_or_int("np.nan")))
        self.assertEqual(float('-inf'), synthetic_pop.cast_string_to_float_or_int('-inf'))
        self.assertEqual(float('-inf'), synthetic_pop.cast_string_to_float_or_int("-inf"))
        self.assertEqual([6, 0.65], synthetic_pop.cast_string_to_float_or_int('[6, 0.65]'))
        self.assertEqual([6, 0.65], synthetic_pop.cast_string_to_float_or_int("[6, 0.65]"))
        self.assertEqual([0.65, np.nan], synthetic_pop.cast_string_to_float_or_int('[0.65, nan]'))
        self.assertEqual([0.65, np.nan], synthetic_pop.cast_string_to_float_or_int("[0.65, nan]"))

    def test_random_sample(self):
        number_of_systems = 3
        synthetic_pop = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig')
        self.assertEqual(number_of_systems, len(synthetic_pop))
        replica_synthetic_pop = synthetic_pop.get_random_subset()  # Providing no argument should default it to a sample size equal to number_of_systems
        self.assertEqual(True, synthetic_pop == replica_synthetic_pop)
        smaller_synthetic_pop = synthetic_pop.get_random_subset(number_of_systems - 1)  # Now we're randomly drawing 2 of the (identical) systems out
        expected_synthetic_pop = sp.SyntheticPopulation(number_of_systems - 1, 'TestWDConfigDB', 'TestPollutionConfig')
        self.assertEqual(number_of_systems - 1, len(smaller_synthetic_pop))
        self.assertEqual(number_of_systems - 1, len(expected_synthetic_pop))
        self.assertEqual(True, smaller_synthetic_pop == expected_synthetic_pop)
        no_pop = synthetic_pop.get_random_subset(1, True)  # We've sampled all systems already so if we prevent reuse it will fail
        self.assertIsNone(no_pop)
        # Start with a fresh pop:
        synthetic_pop2 = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig')
        smaller_synthetic_pop = synthetic_pop2.get_random_subset(number_of_systems - 1, True)
        remaining_pop = synthetic_pop2.get_random_subset(1, True) # Because we prevent reuse, the final system should always be the one we didn't select yet
        all_ids = list()
        for system in smaller_synthetic_pop:
            all_ids.append(system.id)
        for system in remaining_pop:
            all_ids.append(system.id)
        self.assertEqual([0, 1, 2], sorted(all_ids))

    def test_most_polluted_sample(self):
        synthetic_pop_a = sp.SyntheticPopulation(1, 'TestWDConfigDB', 'TestPollutionConfig') # fragment mass = 20
        synthetic_pop_b = sp.SyntheticPopulation(1, 'TestWDConfigDB', 'TestPollutionConfig4') # fragment mass = 19.5
        synthetic_pop_c = sp.SyntheticPopulation(1, 'TestWDConfigDB', 'TestPollutionConfig3') # fragment mass = 20.5
        synthetic_pop = sp.SyntheticPopulation(None, None, None, None, [synthetic_pop_a.population[0], synthetic_pop_b.population[0], synthetic_pop_c.population[0]])

        subset = synthetic_pop.get_most_polluted_subset(2)

        self.assertEqual(2, len(subset))
        self.assertEqual(20.5, subset.population[0].pollution_properties[mp.ModelParameter.fragment_mass])
        self.assertEqual(20, subset.population[1].pollution_properties[mp.ModelParameter.fragment_mass])
        subset2 = synthetic_pop.get_most_polluted_subset()  # This will effectively just reorder the whole population
        self.assertEqual(3, len(subset2))
        self.assertEqual(20.5, subset2.population[0].pollution_properties[mp.ModelParameter.fragment_mass])
        self.assertEqual(20, subset2.population[1].pollution_properties[mp.ModelParameter.fragment_mass])
        self.assertEqual(19.5, subset2.population[2].pollution_properties[mp.ModelParameter.fragment_mass])
        observer = so.Observer(so.ObservationType.IndividualElementCutoff)
        observer.observe_population(synthetic_pop)
        synthetic_pop.population[0].observed_abundances[ci.Element.Mg] = -7.4
        synthetic_pop.population[1].observed_abundances[ci.Element.Mg] = -7.2
        synthetic_pop.population[2].observed_abundances[ci.Element.Mg] = -7.1
        subset3 = synthetic_pop.get_most_polluted_subset(2, ci.Element.Mg)
        self.assertEqual(2, len(subset))
        self.assertEqual(-7.1, subset3.population[0].observed_abundances[ci.Element.Mg])
        self.assertEqual(-7.2, subset3.population[1].observed_abundances[ci.Element.Mg])

    def test_get_sample_with_CaFeMg(self):
        number_of_systems = 4
        synthetic_pop = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig')
        synthetic_pop.population[0].observed = True
        synthetic_pop.population[1].observed = True
        synthetic_pop.population[2].observed = True
        synthetic_pop.population[3].observed = False
        synthetic_pop.population[0].observed_abundances = {ci.Element.Ca: -7.4, ci.Element.Mg: -7.3}
        synthetic_pop.population[1].observed_abundances = {ci.Element.Ca: -7.5, ci.Element.Mg: -7.1, ci.Element.Fe: -6.5}
        synthetic_pop.population[2].observed_abundances = {ci.Element.Ca: -7.8, ci.Element.Mg: -7.9, ci.Element.Fe: -6.2, ci.Element.Zn: -5}
        synthetic_pop.population[3].observed_abundances = {ci.Element.Ca: -4.8, ci.Element.Mg: -4.9, ci.Element.Fe: -4.2}

        subset = synthetic_pop.get_subset_with_detected_elements([ci.Element.Ca, ci.Element.Fe, ci.Element.Mg])
        self.assertEqual(2, len(subset))
        self.assertEqual(-7.5, subset.population[0].observed_abundances[ci.Element.Ca])
        self.assertEqual(-7.1, subset.population[0].observed_abundances[ci.Element.Mg])
        self.assertEqual(-6.5, subset.population[0].observed_abundances[ci.Element.Fe])
        self.assertEqual(-7.8, subset.population[1].observed_abundances[ci.Element.Ca])
        self.assertEqual(-7.9, subset.population[1].observed_abundances[ci.Element.Mg])
        self.assertEqual(-6.2, subset.population[1].observed_abundances[ci.Element.Fe])

    def test_get_knn_p_value(self):
        self.assertEqual(0.74609375, sp.SyntheticPopulation.get_knn_p_value(4, 5))
        self.assertEqual(0.5, sp.SyntheticPopulation.get_knn_p_value(5, 4))
        self.assertEqual(None, sp.SyntheticPopulation.get_knn_p_value(0, 0))
        self.assertEqual(1, sp.SyntheticPopulation.get_knn_p_value(0, 2))
        self.assertEqual(0.25, sp.SyntheticPopulation.get_knn_p_value(2, 0))
        self.assertEqual(None, sp.SyntheticPopulation.get_knn_p_value(-1, -1))
        self.assertEqual(None, sp.SyntheticPopulation.get_knn_p_value(0, -1))
        self.assertEqual(None, sp.SyntheticPopulation.get_knn_p_value(-1, 0))
        self.assertEqual(0.0042095398971260245, sp.SyntheticPopulation.get_knn_p_value(490, 410))
        self.assertEqual(0.623046875, sp.SyntheticPopulation.get_knn_p_value(5, 5))
        self.assertEqual(0.63671875, sp.SyntheticPopulation.get_knn_p_value(4, 4))
        self.assertEqual(0.01953125, sp.SyntheticPopulation.get_knn_p_value(8, 1))
        self.assertEqual(0.998046875, sp.SyntheticPopulation.get_knn_p_value(1, 8))

    def tearDown(self):
        pass


class ObserverTests(unittest.TestCase):

    def assertListsAlmostEqual(self, list1, list2):
        # This function is useful across the board - should make it available at a higher level. Should also add useful output for if the entries aren't equal
        self.assertEqual(len(list1), len(list2))
        for i, val1 in enumerate(list1):
            if np.isnan(val1):
                self.assertTrue(np.isnan(list2[i]))
            else:
                self.assertAlmostEqual(val1, list2[i])

    def test_init(self):
        with self.assertRaises(ValueError):
            observer = so.Observer('Not an observation type')
        observer = so.Observer(so.ObservationType.CaMgFeCutoff)
        self.assertEqual(observer.observation_function, observer.apply_CaMgFeCutoff)

    def test_noise(self):
        observer = so.Observer(so.ObservationType.NoCut, {ci.Element.Fe: 0.1})
        abundances = {
            ci.Element.Fe: -9,
            ci.Element.Ca: -9
        }
        noisy_abundances = observer.apply_errors(abundances)
        self.assertEqual(noisy_abundances[ci.Element.Ca], abundances[ci.Element.Ca])
        self.assertNotEqual(noisy_abundances[ci.Element.Fe], abundances[ci.Element.Fe])
        self.assertTrue(noisy_abundances[ci.Element.Fe] > abundances[ci.Element.Fe] - 1)
        self.assertTrue(noisy_abundances[ci.Element.Fe] < abundances[ci.Element.Fe] + 1)

    def test_NoCut(self):
        observer = so.Observer(so.ObservationType.NoCut)
        abundances = {
            ci.Element.Pt: -9
        }
        system = sp.SyntheticSystem(None, None, abundances)
        self.assertEqual(abundances, observer.apply_NoCut(system))

    def test_CaMgFeCutoff(self):
        observer = so.Observer(so.ObservationType.CaMgFeCutoff)
        abundances_detectable = {
            ci.Element.Fe: -9,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9
        }
        abundances_undetectable = {
            ci.Element.Fe: -9,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9.6
        }
        abundances_na = {
            ci.Element.Fe: -9,
            ci.Element.U: -9
        }
        system_detectable = sp.SyntheticSystem(None, None, abundances_detectable)
        system_undetectable = sp.SyntheticSystem(None, None, abundances_undetectable)
        system_na = sp.SyntheticSystem(None, None, abundances_na)
        self.assertEqual(abundances_detectable, observer.apply_CaMgFeCutoff(system_detectable))
        self.assertEqual(dict(), observer.apply_CaMgFeCutoff(system_undetectable))
        self.assertEqual(dict(), observer.apply_CaMgFeCutoff(system_na))

    def test_IndividualElementCutoff(self):
        observer = so.Observer(so.ObservationType.IndividualElementCutoff)
        abundances = {
            ci.Element.Al: -8.4,
            ci.Element.Ti: -9.6,
            ci.Element.Ca: -7.6,
            ci.Element.Ni: -7.4,
            ci.Element.Fe: -8.4,
            ci.Element.Cr: -7.7,
            ci.Element.Mg: -7.7,
            ci.Element.Si: -6.9,
            ci.Element.Na: -7.1,
            ci.Element.O: -8,
            ci.Element.C: -8,
            ci.Element.N: -8,
            ci.Element.Xe: -4
        }
        expected_observation = {
            ci.Element.Ca: -7.6,
            ci.Element.Mg: -7.7,
            ci.Element.Si: -6.9
        }
        system = sp.SyntheticSystem(None, None, abundances)
        observer.observe_system(system)
        self.assertEqual(True, system.observed)
        self.assertEqual(expected_observation, system.observed_abundances)
        undetectable_abundances = {
            ci.Element.Al: -8.4,
            ci.Element.Ti: -9.6,
            ci.Element.Ca: -9.6,
            ci.Element.Ni: -7.4,
            ci.Element.Fe: -8.4,
            ci.Element.Cr: -7.7,
            ci.Element.Mg: -9.7,
            ci.Element.Si: -9.9,
            ci.Element.Na: -7.1,
            ci.Element.O: -8,
            ci.Element.C: -8,
            ci.Element.N: -8,
            ci.Element.Xe: -4
        }
        undetectable_system = sp.SyntheticSystem(None, None, undetectable_abundances)
        observer.observe_system(undetectable_system)
        self.assertEqual(False, undetectable_system.observed)
        self.assertEqual(None, undetectable_system.observed_abundances)

    def test_TeffIndividualElementCutoff(self):
        observer = so.Observer(so.ObservationType.TeffIndividualElementCutoff)
        abundances = {
            ci.Element.Al: -10,
            ci.Element.Ti: -10,
            ci.Element.Ca: -7.3,
            ci.Element.Ni: -10.3,
            ci.Element.Fe: -7,
            ci.Element.Cr: -10.7,
            ci.Element.Mg: -9.7,
            ci.Element.Si: -7.9,
            ci.Element.Na: -6,
            ci.Element.O: -9,
            ci.Element.C: -15,
            ci.Element.N: -15,
            ci.Element.Xe: -4
        }
        expected_observation_1 = {
            ci.Element.Ca: -7.3,
            ci.Element.Fe: -7,
            ci.Element.Si: -7.9,
            ci.Element.Na: -6
        }
        expected_observation_2 = {
            ci.Element.Ca: -7.3, # We lose Fe in this one: can't observe Fe that low when Teff = 20000
            ci.Element.Si: -7.9,
            ci.Element.Na: -6
        }
        CaFe_ratio = abundances[ci.Element.Ca] - abundances[ci.Element.Fe]
        system = sp.SyntheticSystem({mp.WDParameter.spectral_type: 'DA', mp.WDParameter.temperature: 12345}, None, abundances)
        observer.observe_system(system)
        self.assertEqual(True, system.observed)
        self.assertEqual(expected_observation_1, system.observed_abundances)
        self.assertEqual(CaFe_ratio, system.pollution_abundance_log_ratio(ci.Element.Ca, ci.Element.Fe))
        self.assertEqual(CaFe_ratio, system.observed_abundance_log_ratio(ci.Element.Ca, ci.Element.Fe))
        system2 = sp.SyntheticSystem({mp.WDParameter.spectral_type: 'DA', mp.WDParameter.temperature: 14000}, None, abundances)
        observer.observe_system(system2)
        self.assertEqual(True, system2.observed)
        self.assertEqual(expected_observation_2, system2.observed_abundances)
        self.assertEqual(CaFe_ratio, system2.pollution_abundance_log_ratio(ci.Element.Ca, ci.Element.Fe))
        self.assertEqual(None, system2.observed_abundance_log_ratio(ci.Element.Ca, ci.Element.Fe))

        observer2 = so.Observer(so.ObservationType.TeffIndividualElementCutoff, dict(), -1)
        observer2.observe_system(system)
        self.assertEqual(True, system.observed)
        expected_observation_3 = {
            ci.Element.Ca: -7.3,
            ci.Element.Fe: -7,
            ci.Element.Si: -7.9,
            ci.Element.Na: -6,
            ci.Element.Mg: -9.7
        }
        self.assertEqual(expected_observation_3, system.observed_abundances) # Should now observe additional elements (Mg)

    def test_HollandsColourCut(self):
        observer = so.Observer(so.ObservationType.HollandsColourCut)
        abundances = {
            ci.Element.Bi: -9.5
        }
        bandpasses_detectable = {
            sb.Bandpass.u: 15,
            sb.Bandpass.g: 12.5,
            sb.Bandpass.r: 12
        }
        bandpasses_y_undetectable = {
            sb.Bandpass.u: 18,
            sb.Bandpass.g: 12.5,
            sb.Bandpass.r: 12
        }
        bandpasses_x_undetectable = {
            sb.Bandpass.u: 15,
            sb.Bandpass.g: 12.5,
            sb.Bandpass.r: 12.4
        }
        bandpasses_na = {
            sb.Bandpass.u: 15,
            sb.Bandpass.g: 12.5,
            sb.Bandpass.z: 12
        }
        system_detectable = sp.SyntheticSystem(None, None, abundances, None, None, bandpasses_detectable)
        system_y_undetectable = sp.SyntheticSystem(None, None, abundances, None, None, bandpasses_y_undetectable)
        system_x_undetectable = sp.SyntheticSystem(None, None, abundances, None, None, bandpasses_x_undetectable)
        system_na = sp.SyntheticSystem(None, None, abundances, None, None, bandpasses_na)
        observer.observe_system(system_detectable)
        observer.observe_system(system_y_undetectable)
        observer.observe_system(system_x_undetectable)
        observer.observe_system(system_na)
        self.assertEqual(True, system_detectable.observed)
        self.assertEqual(False, system_y_undetectable.observed)
        self.assertEqual(False, system_x_undetectable.observed)
        self.assertEqual(False, system_na.observed)
        self.assertEqual(abundances, system_detectable.observed_abundances)
        self.assertEqual(None, system_y_undetectable.observed_abundances)
        self.assertEqual(None, system_x_undetectable.observed_abundances)
        self.assertEqual(None, system_na.observed_abundances)

    def test_observe_system(self):
        observer = so.Observer(so.ObservationType.CaMgFeCutoff)
        abundances_detectable = {
            ci.Element.Fe: -9,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9
        }
        abundances_undetectable = {
            ci.Element.Fe: -9,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9.6
        }
        abundances_na = {
            ci.Element.Fe: -9,
            ci.Element.U: -9
        }
        system_detectable = sp.SyntheticSystem(None, None, abundances_detectable)
        system_undetectable = sp.SyntheticSystem(None, None, abundances_undetectable)
        system_na = sp.SyntheticSystem(None, None, abundances_na)
        observer.observe_system(system_detectable)
        self.assertEqual(True, system_detectable.observed)
        observer.observe_system(system_undetectable)
        self.assertEqual(False, system_undetectable.observed)
        observer.observe_system(system_na)
        self.assertEqual(False, system_na.observed)

    def test_observe_populations(self):
        observer = so.Observer(so.ObservationType.CaMgFeCutoff)
        number_of_systems = 2
        pop1 = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig5')
        pop2 = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig2')
        populations = {'Pop1': pop1, 'Pop2': pop2}
        observer.observe_populations(populations)
        empty_pop = sp.SyntheticPopulation(0, 'TestWDConfigDB', 'TestPollutionConfig5')
        expected_observed_populations = {'Pop1': empty_pop, 'Pop2': pop2}
        self.assertEqual(expected_observed_populations['Pop1'], populations['Pop1'].get_observed_subset())
        self.assertEqual(expected_observed_populations['Pop2'], populations['Pop2'].get_observed_subset())
        self.assertListsAlmostEqual([0.05139185220436815, 0.05139185220436815], populations['Pop2'].pollution_abundance_log_ratios(ci.Element.Al, ci.Element.Na))
        self.assertListsAlmostEqual([0.05139185220436815, 0.05139185220436815], populations['Pop2'].observed_abundance_log_ratios(ci.Element.Al, ci.Element.Na))

    def test_reload(self):
        dump_file = 'test_obs_dump.csv'
        for file_which_will_be_affected in [dump_file]:
            if os.path.exists(file_which_will_be_affected):
                raise IOError('This test deletes file ' + file_which_will_be_affected + ', cannot proceed until it is removed or renamed')
        observer1 = so.Observer(so.ObservationType.NoCut)
        observer2 = so.Observer(so.ObservationType.CaMgFeCutoff)
        abundances_detectable = {
            ci.Element.Ca: -9,
            ci.Element.Fe: -9,
            ci.Element.Mg: -9
        }
        abundances_undetectable = {
            ci.Element.Ca: -9,
            ci.Element.Fe: -9,
            ci.Element.Mg: -9.6
        }
        abundances_na = {
            ci.Element.Fe: -9,
            ci.Element.U: -9
        }
        wd_properties = {
            mp.WDParameter.spectral_type: 'DB',
            mp.WDParameter.temperature: 9876,
            mp.WDParameter.logg: 8,
            mp.WDParameter.mass: 0.7,
            mp.WDParameter.distance: 30
        }
        pollution_properties = {
            mp.ModelParameter.metallicity: None,
            mp.ModelParameter.t_sinceaccretion: None,
            mp.ModelParameter.formation_distance: None,
            mp.ModelParameter.feeding_zone_size: None,
            mp.ModelParameter.parent_core_frac: None,
            mp.ModelParameter.parent_crust_frac: None,
            mp.ModelParameter.fragment_core_frac: 0.5,
            mp.ModelParameter.fragment_crust_frac: None,
            mp.ModelParameter.fragment_mass: None,
            mp.ModelParameter.accretion_timescale: None,
            mp.ModelParameter.pressure: None,
            mp.ModelParameter.oxygen_fugacity: None
        }
        system_detectable = sp.SyntheticSystem(wd_properties, pollution_properties, abundances_detectable, None, None, None, None, 1)
        system_undetectable = sp.SyntheticSystem(wd_properties, pollution_properties, abundances_undetectable, None, None, None, None, 2)
        pop1 = sp.SyntheticPopulation(None, 'TestWDConfigDB', 'TestPollutionConfig', dump_file, [system_detectable, system_undetectable])
        observer1.observe_population(pop1)
        pop1.dump_to_csv(dump_file)
        pop2 = sp.SyntheticPopulation(None, 'TestWDConfigDB', 'TestPollutionConfig', dump_file)
        observer2.observe_population(pop2)  # Shouldn't do anything - should realise the pop was already observed
        self.assertEqual(pop1.get_observed_subset(), pop2.get_observed_subset())
        observer2.observe_population(pop2, True)  # Now tell it to overwrite the previous results
        self.assertEqual(2, len(pop1.get_observed_subset()))
        self.assertEqual(1, len(pop2.get_observed_subset()))
        os.remove(dump_file)

class MultiVariateTests(unittest.TestCase):

    def test_convert_data_to_r_format(self):
        data_set1 = [
            [1, 2, 3, 4, 5],
            [10, 20, 30, 40, 50],
            [11, 21, 31, 41, 51]
        ]
        data_set2 = [
            [2, 1, 6, 1.7],
            [11, 28, 20, 50],
            [11, 22, 31, 42]
        ]
        data_set1_r = np.asarray(mv.convert_data_to_r_format(data_set1))
        data_set2_r = np.asarray(mv.convert_data_to_r_format(data_set2))
        expected_data_set1_r = np.asarray([[1, 10, 11], [2, 20, 21], [3, 30, 31], [4, 40, 41], [5, 50, 51]])
        expected_data_set2_r = np.asarray([[2, 11, 11], [1, 28, 22], [6, 20, 31], [1.7, 50, 42]])
        self.assertTrue(np.allclose(expected_data_set1_r, data_set1_r))
        self.assertTrue(np.allclose(expected_data_set2_r, data_set2_r))

    def test_convert_rlistvector_to_dict(self):
        expected_converted_vector = {'foo': 3.12, 'bar': 'baz'}
        rlistvector = robjects.ListVector(expected_converted_vector)
        converted_vector = mv.convert_rlistvector_to_dict(rlistvector)
        self.assertEqual(expected_converted_vector, converted_vector)

    def test_cramer_test_p_value(self):
        if mv.imported_correctly:
            data_set1 = [
                [1, 2, 3, 4, 5],
                [10, 20, 30, 40, 50],
                [11, 21, 31, 41, 51]
            ]
            data_set2 = [
                [10, 20, 30, 40, 50],
                [1, 2, 3, 4, 5],
                [11, 21, 31, 41, 51]
            ]
            data_set1_r = mv.convert_data_to_r_format(data_set1)
            data_set2_r = mv.convert_data_to_r_format(data_set2)
            p_value_match = mv.cramer_test_p_value(data_set1_r, data_set1_r)
            self.assertTrue(p_value_match > 0.8) # This ought to be high (0.8 is quite conservative)
            p_value_mismatch = mv.cramer_test_p_value(data_set1_r, data_set2_r)
            self.assertTrue(p_value_mismatch < 0.1) # This ought to be low (0.1 is conservative)

    # No longer using this test
    #def test_fasano_franceschini_test_p_value(self):
    #    data_set1 = [
    #        [1, 2, 3, 4, 5],
    #        [10, 20, 30, 40, 50],
    #        [11, 21, 31, 41, 51]
    #    ]
    #    data_set2 = [
    #        [10, 20, 30, 40, 50],
    #        [1, 2, 3, 4, 5],
    #        [11, 21, 31, 41, 51]
    #    ]
    #    data_set1_r = mv.convert_data_to_r_format(data_set1)
    #    data_set2_r = mv.convert_data_to_r_format(data_set2)
    #    p_value_match = mv.fasano_franceschini_test_p_value(data_set1_r, data_set1_r)
    #    self.assertTrue(p_value_match > 0.8) # This ought to be high (0.8 is quite conservative)
    #    p_value_mismatch = mv.cramer_test_p_value(data_set1_r, data_set2_r)
    #    self.assertTrue(p_value_mismatch < 0.1) # This ought to be low (0.1 is conservative)

class ExcessOxygenTests(unittest.TestCase):

    def test_init(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        expected_priority = [gi.Layer.mantle, gi.Layer.crust, gi.Layer.bulk, gi.Layer.core]
        self.assertEqual(expected_priority, eo_calc.layer_priority)

    def test_get_favoured_layer(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        abundances = {ci.Element.O: {gi.Layer.core: 0.1, gi.Layer.crust: 0.1, gi.Layer.mantle: None}}
        self.assertEqual(gi.Layer.crust, eo_calc.get_favoured_layer(abundances))

    def test_normalise_abundances(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        abundances = {
            ci.Element.O: {gi.Layer.core: 0.1, gi.Layer.crust: 1},
            ci.Element.Fe: {gi.Layer.core: 0.4, gi.Layer.crust: 1}
        }
        expected_abundances = {
            ci.Element.O: {gi.Layer.core: 0.2, gi.Layer.crust: 0.5},
            ci.Element.Fe: {gi.Layer.core: 0.8, gi.Layer.crust: 0.5}
        }
        self.assertEqual(expected_abundances, eo_calc.normalise_abundances(abundances))

    def test_assign_oxygen(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        abundances = {
            ci.Element.O: {gi.Layer.mantle: 0.2},
            ci.Element.Fe: {gi.Layer.mantle: 0.2},
            ci.Element.Mg: {gi.Layer.mantle: 0.2},
            ci.Element.Si: {gi.Layer.mantle: 0.2},
            ci.Element.Al: {gi.Layer.mantle: 0.1},
            ci.Element.Ca: {gi.Layer.mantle: 0.1}
        }
        layer = gi.Layer.mantle
        ox_strat = eoc.OxidationStrategy.default
        expected_assignation = {
            ci.Element.Fe: 0.2,
            ci.Element.Mg: 0.2,
            ci.Element.Si: 0.4,
            ci.Element.Al: 0.15,
            ci.Element.Ca: 0.1
        }
        for element, ea in expected_assignation.items():
            self.assertAlmostEqual(ea, eo_calc.assign_oxygen(abundances, layer, ox_strat)[element])

    def test_calculate_excess_oxygen(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        abundances = {
            ci.Element.O: {gi.Layer.mantle: 0.6},
            ci.Element.Fe: {gi.Layer.mantle: 0.1},
            ci.Element.Mg: {gi.Layer.mantle: 0.1},
            ci.Element.Si: {gi.Layer.mantle: 0.1},
            ci.Element.Al: {gi.Layer.mantle: 0.05},
            ci.Element.Ca: {gi.Layer.mantle: 0.05}
        }
        layer = gi.Layer.mantle
        ox_strat = eoc.OxidationStrategy.default
        expected_assignation = {
            ci.Element.Fe: 0.1,
            ci.Element.Mg: 0.1,
            ci.Element.Si: 0.2,
            ci.Element.Al: 0.075,
            ci.Element.Ca: 0.05
        }
        excess_oxygen, fractional_excess_oxygen, oxygen_assignations = eo_calc.calculate_excess_oxygen(abundances, layer, ox_strat)
        self.assertAlmostEqual(0.075, excess_oxygen)
        self.assertAlmostEqual(0.125, fractional_excess_oxygen)
        for element, ea in expected_assignation.items():
            self.assertAlmostEqual(ea, oxygen_assignations[element])

    def test_calculate_excess_oxygen_mass(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        abundances = {
            ci.Element.O: {gi.Layer.mantle: 0.6},
            ci.Element.Fe: {gi.Layer.mantle: 0.1},
            ci.Element.Mg: {gi.Layer.mantle: 0.1},
            ci.Element.Si: {gi.Layer.mantle: 0.1},
            ci.Element.Al: {gi.Layer.mantle: 0.05},
            ci.Element.Ca: {gi.Layer.mantle: 0.05}
        }
        layer = gi.Layer.mantle
        fractional_excess_oxygen = 0.125
        excess_o_mass = eo_calc.calculate_excess_oxygen_mass(abundances, layer, fractional_excess_oxygen, 10)
        self.assertAlmostEqual(0.5046817007160045, excess_o_mass)

    def test_calculate_water_abundance(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        abundances = {
            ci.Element.O: {gi.Layer.mantle: 0.6},
            ci.Element.Fe: {gi.Layer.mantle: 0.1},
            ci.Element.Mg: {gi.Layer.mantle: 0.1},
            ci.Element.Si: {gi.Layer.mantle: 0.1},
            ci.Element.Al: {gi.Layer.mantle: 0.05},
            ci.Element.Ca: {gi.Layer.mantle: 0.05}
        }
        layer = gi.Layer.mantle
        fractional_excess_oxygen = 0.125
        water_mass, water_mass_fraction, excess_o_mass = eo_calc.calculate_water_abundance(abundances, layer, fractional_excess_oxygen, 10)
        self.assertAlmostEqual(0.5682755696230278, water_mass)
        self.assertAlmostEqual(0.05682755696230278, water_mass_fraction)
        self.assertAlmostEqual(0.5046817007160045, excess_o_mass)

    def test_reformat_abundances(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        initial_dict = {ci.Element.Fe: 0.2}
        expected_dict1 = {ci.Element.Fe: {gi.Layer.bulk: 0.2}}
        expected_dict2 = {ci.Element.Fe: {gi.Layer.bulk: 0.2}}
        reformatted_dict1 = eo_calc.reformat_abundances(initial_dict)
        reformatted_dict2 = eo_calc.reformat_abundances(reformatted_dict1)
        self.assertEqual(expected_dict1, reformatted_dict1)
        self.assertEqual(expected_dict2, reformatted_dict2)

    def test_sigma_calculation(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        abundances = {
            ci.Element.O: {gi.Layer.mantle: 0.6},
            ci.Element.Fe: {gi.Layer.mantle: 0.1},
            ci.Element.Mg: {gi.Layer.mantle: 0.1},
            ci.Element.Si: {gi.Layer.mantle: 0.1},
            ci.Element.Al: {gi.Layer.mantle: 0.05},
            ci.Element.Ca: {gi.Layer.mantle: 0.05}
        }
        favoured_layer = eo_calc.get_favoured_layer(abundances)
        ox_strat = eoc.OxidationStrategy.default
        o_error = 0.1
        sigma_significance = eo_calc.calculate_sigma_significance(abundances, favoured_layer, ox_strat, o_error)
        expected_sigma = 0.024133645601883087/o_error
        self.assertAlmostEqual(expected_sigma, sigma_significance)

    def test_calculate_conservative_composition(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        #calculate_conservative_composition(self, nominal_model, observations, reference_element, reference_layer, elements_to_alter=None, force_data_fit=False)
        simple_model = {
            ci.Element.Mg: {
                gi.Layer.bulk: 0.1
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.85 # Should try to move to observed abundance
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.05 # Should try to move to solar abundance
            }
        }
        simple_observations = {
            ci.Element.Mg: -8.3,
            ci.Element.O: -7.4,
            ci.Element.Fe: -8.4
        }
        simple_cc = eo_calc.calculate_conservative_composition(simple_model, simple_observations, ci.Element.Mg, gi.Layer.bulk)
        expected_cc = {
            ci.Element.Mg: {
                gi.Layer.bulk: 0.10229180723840627
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.8125327067043979
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.08517548605719581
            }
        }
        self.assertEqual(expected_cc, simple_cc)
        with self.assertRaises(AssertionError):
            noref_cc = eo_calc.calculate_conservative_composition(simple_model, simple_observations, ci.Element.U, gi.Layer.bulk)
        with self.assertRaises(KeyError):
            noref_cc = eo_calc.calculate_conservative_composition(simple_model, simple_observations, ci.Element.Mg, gi.Layer.crust)
        expected_cc_only_o = {
            ci.Element.Mg: {
                gi.Layer.bulk: 0.10589538290062592
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.8411569256490611
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.05294769145031296
            }
        }
        simple_cc_only_o = eo_calc.calculate_conservative_composition(simple_model, simple_observations, ci.Element.Mg, gi.Layer.bulk, [ci.Element.O])
        self.assertEqual(expected_cc_only_o, simple_cc_only_o)
        expected_cc_data = {
            ci.Element.Mg: {
                gi.Layer.bulk: 0.10269459756912862
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.8157321840280648
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.08157321840280649
            }
        }
        simple_cc_data = eo_calc.calculate_conservative_composition(simple_model, simple_observations, ci.Element.Mg, gi.Layer.bulk, None, True)
        self.assertEqual(expected_cc_data, simple_cc_data)
        realistic_observations = {
            ci.Element.Al: 0.0,
            ci.Element.Ti: 0.0,
            ci.Element.Ca: -10.1,
            ci.Element.Ni: 0.0,
            ci.Element.Fe: -8.5,
            ci.Element.Cr: 0.0,
            ci.Element.Mg: -8.5,
            ci.Element.Si: 0.0,
            ci.Element.Na: 0.0,
            ci.Element.O: 0.0,
            ci.Element.C: 0.0,
            ci.Element.N: 0.0
        }
        realistic_model = {
            ci.Element.Hf: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Ta: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Nb: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Si: {
                gi.Layer.bulk: 0.11651439896972245,
                gi.Layer.core: 0.03684790596954492,
                gi.Layer.mantle: 0.13202739688278847
            },
            ci.Element.Mn: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Zn: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Ga: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.V: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Cr: {
                gi.Layer.bulk: 0.0023244276062221506,
                gi.Layer.core: 0.005309869081576006,
                gi.Layer.mantle: 0.0017430897610048158
            },
            ci.Element.Cu: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.1682998378407617,
                gi.Layer.core: 0.8807358616602161,
                gi.Layer.mantle: 0.029571269097568523
            },
            ci.Element.W: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.P: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Co: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Ni: {
                gi.Layer.bulk: 0.008836565622895802,
                gi.Layer.core: 0.0524941854839732,
                gi.Layer.mantle: 0.0003353684408783415
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.5241667429990097,
                gi.Layer.core: 0.02461217780468976,
                gi.Layer.mantle: 0.6214421298709385
            },
            ci.Element.C: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.S: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.N: {
                gi.Layer.bulk: 0.0,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0
            },
            ci.Element.Na: {
                gi.Layer.bulk: 0.007016108650306535,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.008382315128729115
            },
            ci.Element.Mg: {
                gi.Layer.bulk: 0.1499803239444996,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.17918512968819916
            },
            ci.Element.Al: {
                gi.Layer.bulk: 0.010840528321201531,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.012951442042770243
            },
            ci.Element.Ti: {
                gi.Layer.bulk: 0.0003941440918770745,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.0004708935036369562
                },
            ci.Element.Ca: {
                gi.Layer.bulk: 0.011626921953503611,
                gi.Layer.core: 0.0,
                gi.Layer.mantle: 0.013890965583485931
            }
        }
        realistic_cc_data = eo_calc.calculate_conservative_composition(realistic_model, realistic_observations, ci.Element.Mg, gi.Layer.bulk)
        expected_realistic_cc = {
            ci.Element.Hf: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Ta: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Nb: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Si: {
                gi.Layer.bulk: 0.06027418094294684
            },
            ci.Element.Mn: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Zn: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Ga: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.V: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Cr: {
                gi.Layer.bulk: 0.0009776748084421546
            },
            ci.Element.Cu: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.07078840024157185
            },
            ci.Element.W: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.P: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Co: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Ni: {
                gi.Layer.bulk: 0.0037167376516803875
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.22046916784223025
            },
            ci.Element.C: {
                gi.Layer.bulk: 0.4559083375185133
            },
            ci.Element.S: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.N: {
                gi.Layer.bulk: 0.1122156241881471
            },
            ci.Element.Na: {
                gi.Layer.bulk: 0.002951037348866442
            },
            ci.Element.Mg: {
                gi.Layer.bulk: 0.06308305067880741
            },
            ci.Element.Al: {
                gi.Layer.bulk: 0.004559622085657477
            },
            ci.Element.Ti: {
                gi.Layer.bulk: 0.00016578049085848686
            },
            ci.Element.Ca: {
                gi.Layer.bulk: 0.004890386202278275
            }
        }
        self.assertEqual(expected_realistic_cc, realistic_cc_data)
        solar_Fe_comp = eo_calc.calculate_conservative_composition(realistic_model, None, ci.Element.Mg, gi.Layer.bulk, [ci.Element.Fe])
        expected_solar_Fe_comp = dict()
        for element, abundances in realistic_model.items():
            self.assertAlmostEqual(abundances[gi.Layer.bulk], solar_Fe_comp[element][gi.Layer.bulk])
        data_O_comp = eo_calc.calculate_conservative_composition(realistic_model, realistic_observations, ci.Element.Mg, gi.Layer.bulk, [ci.Element.O], True)
        for element, abundances in realistic_model.items():
            self.assertAlmostEqual(abundances[gi.Layer.bulk], data_O_comp[element][gi.Layer.bulk])
        data_comp = eo_calc.calculate_conservative_composition(realistic_model, realistic_observations, ci.Element.Mg, gi.Layer.bulk, None, True)
        expected_data_comp = {
            ci.Element.Hf: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Ta: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Nb: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Si: {
                gi.Layer.bulk: 0.11964664038601591
            },
            ci.Element.Mn: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Zn: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Ga: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.V: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Cr: {
                gi.Layer.bulk: 0.002386914890899102
            },
            ci.Element.Cu: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Fe: {
                gi.Layer.bulk: 0.15401222546432933
            },
            ci.Element.W: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.P: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Co: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Ni: {
                gi.Layer.bulk: 0.009074117865936782
            },
            ci.Element.O: {
                gi.Layer.bulk: 0.5382578492998866
            },
            ci.Element.C: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.S: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.N: {
                gi.Layer.bulk: 0.0
            },
            ci.Element.Na: {
                gi.Layer.bulk: 0.007204721785593065
            },
            ci.Element.Mg: {
                gi.Layer.bulk: 0.15401222546432933
            },
            ci.Element.Al: {
                gi.Layer.bulk: 0.011131952832527324
            },
            ci.Element.Ti: {
                gi.Layer.bulk: 0.0004047398161779445
            },
            ci.Element.Ca: {
                gi.Layer.bulk: 0.003868612194304434
            }
        }
        self.assertAlmostEqual(expected_data_comp, data_comp)

    def test_check_for_nans(self):
        abundances1 = {ci.Element.O: 0.1}
        abundances2 = {ci.Element.O: None}
        abundances3 = {ci.Element.O: np.nan}
        abundances4 = {ci.Element.O: 0.3, ci.Element.Mg: {gi.Layer.bulk: 0.1, gi.Layer.core: np.nan}}
        abundances5 = {ci.Element.O: {gi.Layer.bulk: 0.1, gi.Layer.core: 0.2}, ci.Element.Mg: {gi.Layer.bulk: 0.1, gi.Layer.core: np.nan}}
        eo_calc = eoc.ExcessOxygenCalculator()
        self.assertFalse(eo_calc.check_for_nans(abundances1, gi.Layer.bulk))
        self.assertTrue(eo_calc.check_for_nans(abundances2, gi.Layer.bulk))
        self.assertTrue(eo_calc.check_for_nans(abundances3, gi.Layer.bulk))
        self.assertFalse(eo_calc.check_for_nans(abundances4, gi.Layer.bulk))
        self.assertTrue(eo_calc.check_for_nans(abundances4, gi.Layer.core))
        self.assertTrue(eo_calc.check_for_nans(abundances4, gi.Layer.crust))
        self.assertFalse(eo_calc.check_for_nans(abundances5, gi.Layer.bulk))
        self.assertTrue(eo_calc.check_for_nans(abundances5, gi.Layer.core))
        self.assertTrue(eo_calc.check_for_nans(abundances5, gi.Layer.crust))

    def test_full_calculation(self):
        eo_calc = eoc.ExcessOxygenCalculator()
        abundances = {
            'Test': {
                ci.Element.O: {gi.Layer.mantle: 0.3},
                ci.Element.Fe: {gi.Layer.mantle: 0.05},
                ci.Element.Mg: {gi.Layer.mantle: 0.05},
                ci.Element.Si: {gi.Layer.mantle: 0.05},
                ci.Element.Al: {gi.Layer.mantle: 0.025},
                ci.Element.Ca: {gi.Layer.mantle: 0.025}
            }
        }

        eoc_stat_dict = eo_calc.run_full_calculation(abundances, 10)
        excess_oxygen, fractional_excess_oxygen, oxygen_assignations, water_mass, water_mass_fraction, excess_o_mass, normalised_abundances, species_names, sigma, layer = eoc_stat_dict['Test'][eoc.OxidationStrategy.default]
        self.assertEqual(gi.Layer.mantle, layer)
        self.assertAlmostEqual(0.075, excess_oxygen)
        self.assertAlmostEqual(0.125, fractional_excess_oxygen)
        expected_assignation = {
            ci.Element.Fe: 0.1,
            ci.Element.Mg: 0.1,
            ci.Element.Si: 0.2,
            ci.Element.Al: 0.075,
            ci.Element.Ca: 0.05
        }
        for element, ea in expected_assignation.items():
            self.assertAlmostEqual(ea, oxygen_assignations[element])
        self.assertAlmostEqual(0.5682755696230278, water_mass)
        self.assertAlmostEqual(0.05682755696230278, water_mass_fraction)
        self.assertAlmostEqual(0.5046817007160045, excess_o_mass)
        expected_normalised_abundances = {
            ci.Element.O: {gi.Layer.mantle: 0.6},
            ci.Element.Fe: {gi.Layer.mantle: 0.1},
            ci.Element.Mg: {gi.Layer.mantle: 0.1},
            ci.Element.Si: {gi.Layer.mantle: 0.1},
            ci.Element.Al: {gi.Layer.mantle: 0.05},
            ci.Element.Ca: {gi.Layer.mantle: 0.05}
        }
        expected_species_names = {
            ci.Element.Mg: 'MgO',
            ci.Element.Al: 'Al$_2$O$_3$',
            ci.Element.Si: 'SiO$_2$',
            ci.Element.Ca: 'CaO',
            ci.Element.Fe: 'FeO'
        }
        self.assertEqual(expected_normalised_abundances, normalised_abundances)
        self.assertEqual(expected_species_names, species_names)
        self.assertIsNone(sigma)

#@unittest.skip("Skip for now")
class ModellerTests(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(ValueError):
            modeller = sm.Modeller('Not a modeller type')
        modeller = sm.Modeller(sm.ModellerType.SimpleFcfInterpolation)
        self.assertEqual(modeller.model_system, modeller.apply_simple_fcf_interpolation)

    def test_null_modeller(self):
        modeller = sm.Modeller(sm.ModellerType.Null)
        system1 = sp.SyntheticSystem(None, None, None, None, None, None, None)
        system2 = sp.SyntheticSystem(None, None, None, None, None, None, {mp.ModelParameter.fragment_crust_frac: -0.1})
        modeller.model_system(system1)
        self.assertEqual({
            mp.ModelParameter.metallicity: [None],
            mp.ModelParameter.t_sinceaccretion: [None],
            mp.ModelParameter.formation_distance: [None],
            mp.ModelParameter.feeding_zone_size: [None],
            mp.ModelParameter.parent_core_frac: [None],
            mp.ModelParameter.parent_crust_frac: [None],
            mp.ModelParameter.fragment_core_frac: [None],
            mp.ModelParameter.fragment_crust_frac: [None],
            mp.ModelParameter.fragment_mass: [None],
            mp.ModelParameter.accretion_timescale: [None],
            mp.ModelParameter.pressure: [None],
            mp.ModelParameter.oxygen_fugacity: [None]
            }, system1.modelled_properties)
        modeller.model_system(system2)
        self.assertEqual({
            mp.ModelParameter.metallicity: [None],
            mp.ModelParameter.t_sinceaccretion: [None],
            mp.ModelParameter.formation_distance: [None],
            mp.ModelParameter.feeding_zone_size: [None],
            mp.ModelParameter.parent_core_frac: [None],
            mp.ModelParameter.parent_crust_frac: [None],
            mp.ModelParameter.fragment_core_frac: [None],
            mp.ModelParameter.fragment_crust_frac: -0.1,
            mp.ModelParameter.fragment_mass: [None],
            mp.ModelParameter.accretion_timescale: [None],
            mp.ModelParameter.pressure: [None],
            mp.ModelParameter.oxygen_fugacity: [None]
            }, system2.modelled_properties)

    def test_simple_fcf_interpolation(self):
        modeller = sm.Modeller(sm.ModellerType.SimpleFcfInterpolation)
        abundances_detectable = {
            ci.Element.Fe: -9,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9
        }
        abundances_undetectable = {
            ci.Element.Fe: -9,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9.5
        }
        abundances_na = {
            ci.Element.Fe: -9,
            ci.Element.U: -9
        }
        system_detectable = sp.SyntheticSystem(None, None, abundances_detectable)
        system_detectable.set_observed_abundances(abundances_detectable)
        system_undetectable = sp.SyntheticSystem(None, None, abundances_undetectable)
        system_undetectable.set_observed_abundances(abundances_undetectable)
        system_na = sp.SyntheticSystem(None, None, abundances_na)
        system_na.set_observed_abundances(abundances_na)
        expected_values_detectable = {
            mp.ModelParameter.fragment_core_frac: [1/3]
        }
        expected_values_undetectable = {
            mp.ModelParameter.fragment_core_frac: [0.5]
        }
        self.assertEqual(None, system_detectable.modelled_properties)
        modeller.model_system(system_detectable)
        self.assertEqual(expected_values_detectable, system_detectable.modelled_properties)

        self.assertEqual(None, system_undetectable.modelled_properties)
        modeller.model_system(system_undetectable)
        self.assertEqual(expected_values_undetectable, system_undetectable.modelled_properties)

        self.assertEqual(None, system_na.modelled_properties)
        modeller.model_system(system_na)
        self.assertEqual(None, system_na.modelled_properties)

    def test_model_populations(self):
        modeller = sm.Modeller(sm.ModellerType.SimpleFcfInterpolation)
        number_of_systems = 2
        pop1 = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig')
        pop2 = sp.SyntheticPopulation(number_of_systems, 'TestWDConfigDB', 'TestPollutionConfig2')
        populations = {'Pop1': pop1, 'Pop2': pop2}
        for system in populations['Pop1']:
            self.assertEqual(None, system.modelled_properties)
            system.observed = True
            system.set_observed_abundances(system.pollution_abundances)
        for system in populations['Pop2']:
            self.assertEqual(None, system.modelled_properties)
            system.observed = True
            system.set_observed_abundances(system.pollution_abundances)
        modeller.model_populations(populations)
        expected_pop_dict = {
            mp.ModelParameter.fragment_core_frac: [0.2703141631754426]
        }
        for system in populations['Pop1']:
            self.assertEqual(expected_pop_dict, system.modelled_properties)
        expected_pop_dict_2 = {
            mp.ModelParameter.fragment_core_frac: [0.27031416317544227] # This is the same as above because the only difference is the fragment mass, which doesn't affect relative abundances
        }
        for system in populations['Pop2']:
            self.assertEqual(expected_pop_dict_2, system.modelled_properties)

        pop3 = sp.SyntheticPopulation(None, 'TestWDConfigDB', 'TestPollutionConfig', None, [sp.SyntheticSystem(None, None, None)])
        population_none = {'Pop3': pop3}
        modeller.model_populations(population_none)
        self.assertEqual(None, population_none['Pop3'].population[0].modelled_properties)
        population_none['Pop3'].population[0].pollution_abundances = {ci.Element.Yb: -3}
        population_none['Pop3'].population[0].observed = True
        population_none['Pop3'].population[0].set_observed_abundances(population_none['Pop3'].population[0].pollution_abundances)
        modeller.model_populations(population_none)
        self.assertEqual(None, population_none['Pop3'].population[0].modelled_properties)
        population_none['Pop3'].population[0].pollution_abundances = {ci.Element.Mg: -3, ci.Element.Fe: -12}
        population_none['Pop3'].population[0].observed = True
        population_none['Pop3'].population[0].set_observed_abundances(population_none['Pop3'].population[0].pollution_abundances)
        modeller.model_populations(population_none)
        expected_pop_none_dict = {
            mp.ModelParameter.fragment_core_frac: [0]
        }
        self.assertEqual(expected_pop_none_dict, population_none['Pop3'].population[0].modelled_properties)
        pop4 = sp.SyntheticPopulation(1, 'DAstar', 'DApollution')
        DApopulations = {'pop4': pop4}
        for system in DApopulations['pop4']:
            system.observed = True
            system.set_observed_abundances(system.pollution_abundances)
        modeller2 = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        modeller2.model_populations(DApopulations)
        expected_pop_dict_4 = {
            mp.ModelParameter.metallicity: [388, 478],
            mp.ModelParameter.formation_distance: [-0.4302811371772311, -0.4378175613917902],
            mp.ModelParameter.fragment_core_frac: [0.1789431804147672, 0.17939637253110072],
            mp.ModelParameter.accretion_timescale: [34836.13746783229, 34836.13746783229],
            mp.ModelParameter.t_sinceaccretion: [0.017418068733916145, 0.017418068733916145]
        }
        for system in DApopulations['pop4']:
            print()
            print(system.modelled_properties)
            self.assertEqual(expected_pop_dict_4, system.modelled_properties)

    def test_analytic_sinking(self):
        modeller = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        element1 = ci.Element.Al
        element2 = ci.Element.Ca
        element_ratio = 2
        t_1 = 1000000
        t_2 = 800000
        t_disc, t_since_accretion = modeller.calculate_sinking_from_element_pair(element1, element2, element_ratio, t_1, t_2)
        expected_t_disc = 3*max(t_1, t_2)
        self.assertEqual(expected_t_disc, t_disc[478])
        self.assertEqual(389480.5908203125, t_since_accretion[478])

    def test_analytic_fcf(self):
        modeller = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        element1 = ci.Element.Ca
        element2 = ci.Element.Fe
        element_ratio = 0.4
        fcf = modeller.calculate_fcf_from_element_pair(element1, element2, element_ratio)
        self.assertEqual(0.019090892197292972, fcf)

    def test_analytic_heating(self):
        manager = mn.Manager()  #This is just to load the stellar compositions
        manager.publish_live_data(0, ti.TimescaleType.KoesterOvershoot)
        modeller = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        element1 = ci.Element.Na
        element2 = ci.Element.Ca
        target_ratio = 0.6
        t_disc = 100000
        t_sinceaccretion = 0.01
        t_1 = 1000000
        t_2 = 800000
        d_formation = modeller.find_d_formation(element1, element2, target_ratio, [t_disc], [t_sinceaccretion], t_1, t_2, [478])
        self.assertEqual([0.34086889028549194], d_formation)

    def test_analytic_approximation(self):
        manager = mn.Manager()  #This is just to load the stellar compositions
        manager.publish_live_data(0, ti.TimescaleType.KoesterOvershoot)
        modeller = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        wd_properties = {
            mp.WDParameter.spectral_type: 'DB',
            mp.WDParameter.logg: 8,
            mp.WDParameter.temperature: 9876,
            mp.WDParameter.consider_thermohaline: False
        }
        abundances = {
            ci.Element.Al: -8.5,
            ci.Element.Fe: -7,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9,
            ci.Element.Na: -9.2
        }
        system = sp.SyntheticSystem(wd_properties, None, abundances)
        system.observed = True
        system.set_observed_abundances(system.pollution_abundances)
        modeller.apply_analytic_approximation(system)
        expected_modelled_properties_for_index_478 = {
            mp.ModelParameter.formation_distance: -0.5126873560288674,
            mp.ModelParameter.fragment_core_frac: 0.8101129703838985,
            mp.ModelParameter.accretion_timescale: 13262007.058409285,
            mp.ModelParameter.t_sinceaccretion: 18.172781439835305
        }
        self.assertEqual([None], system.modelled_properties[mp.ModelParameter.metallicity])
        self.assertEqual(expected_modelled_properties_for_index_478[mp.ModelParameter.formation_distance], system.modelled_properties[mp.ModelParameter.formation_distance][478])
        self.assertEqual(expected_modelled_properties_for_index_478[mp.ModelParameter.fragment_core_frac], system.modelled_properties[mp.ModelParameter.fragment_core_frac][478])
        self.assertEqual(expected_modelled_properties_for_index_478[mp.ModelParameter.accretion_timescale], system.modelled_properties[mp.ModelParameter.accretion_timescale][478])
        self.assertEqual(expected_modelled_properties_for_index_478[mp.ModelParameter.t_sinceaccretion], system.modelled_properties[mp.ModelParameter.t_sinceaccretion][478])

    def test_collapse_list_of_repeats(self):
        modeller = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        test_list1 = [1]
        test_list2 = [1, 1]
        test_list3 = [1, 2]
        test_list4 = [None, None]
        self.assertEqual([1], modeller.collapse_list_of_repeats(test_list1))
        self.assertEqual([1], modeller.collapse_list_of_repeats(test_list2))
        self.assertEqual([1, 2], modeller.collapse_list_of_repeats(test_list3))
        self.assertEqual([None], modeller.collapse_list_of_repeats(test_list4))

    def test_modelling_with_missing_elements(self):
        wd_properties = {
            mp.WDParameter.spectral_type: 'DA',
            mp.WDParameter.logg: 8,
            mp.WDParameter.temperature: 9876,
            mp.WDParameter.consider_thermohaline: False
        }
        abundances_set_1 = {
            ci.Element.Al: -8.5,
            ci.Element.Fe: -7,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9,
            ci.Element.Na: -9.2
        }
        system = sp.SyntheticSystem(wd_properties, None, abundances_set_1)
        observer = so.Observer(so.ObservationType.IndividualElementCutoff)
        observer.observe_system(system)
        self.assertEqual(True, system.observed)
        self.assertEqual({ci.Element.Fe: -7}, system.observed_abundances)
        modeller = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        modeller.apply_analytic_approximation(system)
        expected_modelled_properties = {
            mp.ModelParameter.formation_distance: [None],
            mp.ModelParameter.fragment_core_frac: [None],
            mp.ModelParameter.accretion_timescale: [6743.8825121139835],
            mp.ModelParameter.t_sinceaccretion: [0.003371941256056992],
            mp.ModelParameter.metallicity: [None]
        }
        self.assertEqual(expected_modelled_properties, system.modelled_properties)
        wd_properties_2 = {
            mp.WDParameter.spectral_type: 'DB',
            mp.WDParameter.logg: 8,
            mp.WDParameter.temperature: 9876,
            mp.WDParameter.consider_thermohaline: False
        }
        abundances_set_1 = {
            ci.Element.Al: -8.5,
            ci.Element.Fe: -7,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9,
            ci.Element.Na: -9.2
        }
        system2 = sp.SyntheticSystem(wd_properties_2, None, abundances_set_1)
        observer.observe_system(system2)
        self.assertEqual(True, system2.observed)
        self.assertEqual({ci.Element.Fe: -7}, system2.observed_abundances)
        modeller.apply_analytic_approximation(system2)
        expected_modelled_properties2 = {
            mp.ModelParameter.formation_distance: [None],
            mp.ModelParameter.fragment_core_frac: [None],
            mp.ModelParameter.accretion_timescale: [None],
            mp.ModelParameter.t_sinceaccretion: [None],
            mp.ModelParameter.metallicity: [None]
        }
        self.assertEqual(expected_modelled_properties2, system2.modelled_properties)
        abundances_set_2 = {
            ci.Element.Al: -8.5,
            ci.Element.Fe: -7,
            ci.Element.Ca: -8.1,
            ci.Element.Mg: -7.9,
            ci.Element.Na: -9.2
        }
        system3 = sp.SyntheticSystem(wd_properties, None, abundances_set_2)
        observer.observe_system(system3)
        self.assertEqual(True, system3.observed)
        self.assertEqual({ci.Element.Fe: -7, ci.Element.Ca: -8.1, ci.Element.Mg: -7.9}, system3.observed_abundances)
        modeller.apply_analytic_approximation(system3)
        expected_modelled_properties3 = {
            mp.ModelParameter.formation_distance: [-0.8070552938588017],
            mp.ModelParameter.fragment_core_frac: [0.7640731816983493],
            mp.ModelParameter.accretion_timescale: [6743.8825121139835],
            mp.ModelParameter.t_sinceaccretion: [0.003371941256056992],
            mp.ModelParameter.metallicity: [749]
        }
        self.assertEqual(expected_modelled_properties3, system3.modelled_properties)
        abundances_set_3 = {
            ci.Element.Al: -8.5,
            ci.Element.Fe: -7,
            ci.Element.Ca: -8.1,
            ci.Element.Mg: -7.9,
            ci.Element.Na: -6.9
        }
        system4 = sp.SyntheticSystem(wd_properties, None, abundances_set_3)
        observer.observe_system(system4)
        self.assertEqual(True, system4.observed)
        self.assertEqual({ci.Element.Fe: -7, ci.Element.Ca: -8.1, ci.Element.Mg: -7.9, ci.Element.Na: -6.9}, system4.observed_abundances)
        modeller.apply_analytic_approximation(system4)
        expected_modelled_properties4 = {
            mp.ModelParameter.formation_distance: [0.47712125471966244],
            mp.ModelParameter.fragment_core_frac: [0.20736772075259863],
            mp.ModelParameter.accretion_timescale: [6743.8825121139835],
            mp.ModelParameter.t_sinceaccretion: [0.003371941256056992],
            mp.ModelParameter.metallicity: [749]
        }
        self.assertEqual(expected_modelled_properties4, system4.modelled_properties)
        abundances_set_4 = {
            ci.Element.Al: -6.9,
            ci.Element.Fe: -9,
            ci.Element.Ca: -8.1,
            ci.Element.Mg: -7.9,
            ci.Element.Na: -6.9
        }
        system5 = sp.SyntheticSystem(wd_properties, None, abundances_set_4)
        observer.observe_system(system5)
        self.assertEqual(True, system5.observed)
        self.assertEqual({ci.Element.Al: -6.9, ci.Element.Ca: -8.1, ci.Element.Mg: -7.9, ci.Element.Na: -6.9}, system5.observed_abundances)
        modeller.apply_analytic_approximation(system5)
        expected_modelled_properties5 = {
            mp.ModelParameter.formation_distance: [0.47712125471966244],
            mp.ModelParameter.fragment_core_frac: [None],
            mp.ModelParameter.accretion_timescale: [6743.8825121139835],
            mp.ModelParameter.t_sinceaccretion: [0.003371941256056992],
            mp.ModelParameter.metallicity: [17]
        }
        self.assertEqual(expected_modelled_properties5, system5.modelled_properties)
        abundances_set_5 = {
            ci.Element.Al: -6.9,
            ci.Element.Fe: -7,
            ci.Element.Ca: -8.1,
            ci.Element.Mg: -7.9,
            ci.Element.Na: -6.9
        }
        system6 = sp.SyntheticSystem(wd_properties, None, abundances_set_5)
        observer.observe_system(system6)
        self.assertEqual(True, system6.observed)
        self.assertEqual({ci.Element.Al: -6.9, ci.Element.Fe: -7, ci.Element.Ca: -8.1, ci.Element.Mg: -7.9, ci.Element.Na: -6.9}, system6.observed_abundances)
        modeller.apply_analytic_approximation(system6)
        expected_modelled_properties6 = {
            mp.ModelParameter.formation_distance: [0.47712125471966244],
            mp.ModelParameter.fragment_core_frac: [0.20736772075259863],
            mp.ModelParameter.accretion_timescale: [6743.8825121139835],
            mp.ModelParameter.t_sinceaccretion: [0.003371941256056992],
            mp.ModelParameter.metallicity: [17]
        }
        self.assertEqual(expected_modelled_properties6, system6.modelled_properties)

        abundances_set_6 = {
            ci.Element.Al: -10,
            ci.Element.Fe: -7.8,
            ci.Element.Ca: -10.1,
            ci.Element.Mg: -7.9,
            ci.Element.Na: -10.9
        }
        system7 = sp.SyntheticSystem(wd_properties, None, abundances_set_6)
        observer.observe_system(system7)
        self.assertEqual(True, system7.observed)
        self.assertEqual({ci.Element.Fe: -7.8, ci.Element.Mg: -7.9}, system7.observed_abundances)
        modeller.apply_analytic_approximation(system7)
        expected_modelled_properties7 = {
            mp.ModelParameter.formation_distance: [None],
            mp.ModelParameter.fragment_core_frac: [0.36879615312359404],
            mp.ModelParameter.accretion_timescale: [6743.8825121139835],
            mp.ModelParameter.t_sinceaccretion: [0.003371941256056992],
            mp.ModelParameter.metallicity: [None]
        }
        self.assertEqual(expected_modelled_properties7, system7.modelled_properties)

    def test_reload(self):
        dump_file = 'test_mod_dump.csv'
        for file_which_will_be_affected in [dump_file]:
            if os.path.exists(file_which_will_be_affected):
                raise IOError('This test deletes file ' + file_which_will_be_affected + ', cannot proceed until it is removed or renamed')
        modeller1 = sm.Modeller(sm.ModellerType.SimpleFcfInterpolation)
        modeller2 = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        abundances_detectable = {
            ci.Element.Al: -8.8,
            ci.Element.Fe: -8,
            ci.Element.Ca: -9,
            ci.Element.Mg: -9,
            ci.Element.Na: -8.5
        }

        expected_values_detectable1 = {
            mp.ModelParameter.fragment_core_frac: [2/3]
        }
        expected_values_detectable2 = {
            mp.ModelParameter.metallicity: [17],
            mp.ModelParameter.formation_distance: [0.47712125471966244],
            mp.ModelParameter.fragment_core_frac: [0.1688303282450301],
            mp.ModelParameter.t_sinceaccretion: [0.003371941256056992],
            mp.ModelParameter.accretion_timescale: [6743.8825121139835]
        }
        expected_values_detectable4 = {
            mp.ModelParameter.metallicity: [None],
            mp.ModelParameter.formation_distance: [0.47712125471966244, 0.2],
            mp.ModelParameter.fragment_core_frac: [0.16991075071445785, 0.3],
            mp.ModelParameter.t_sinceaccretion: [13.271849026779673, 0.4],
            mp.ModelParameter.accretion_timescale: [13262007.058409285, 0.5]
        }
        wd_properties = {
            mp.WDParameter.spectral_type: 'DA',
            mp.WDParameter.logg: 8,
            mp.WDParameter.temperature: 9876,
            mp.WDParameter.consider_thermohaline: False
        }
        system_detectable = sp.SyntheticSystem(wd_properties, None, abundances_detectable, None, None, None, None, 1)
        system_detectable.observed = True
        system_detectable.set_observed_abundances(abundances_detectable)
        pop1 = sp.SyntheticPopulation(None, 'TestWDConfigDB', 'TestPollutionConfig', dump_file, [system_detectable])
        modeller1.model_population(pop1)

        pop2 = sp.SyntheticPopulation(None, 'TestWDConfigDB', 'TestPollutionConfig', dump_file)

        self.assertEqual(expected_values_detectable1, pop1.population[0].modelled_properties)
        self.assertEqual(expected_values_detectable1, pop2.population[0].modelled_properties)

        modeller2.model_population(pop2)  # This shouldn't do anything
        self.assertEqual(expected_values_detectable1, pop2.population[0].modelled_properties)

        os.remove(dump_file) # If it's not removed, it will abort the dump
        modeller2.model_population(pop2, True)  # This should overwrite

        self.assertEqual(expected_values_detectable2, pop2.population[0].modelled_properties)
        pop3 = sp.SyntheticPopulation(None, 'TestWDConfigDB', 'TestPollutionConfig', dump_file)
        self.assertEqual(expected_values_detectable2, pop3.population[0].modelled_properties)

        pop3.population[0].set_modelled_properties(expected_values_detectable4)
        os.remove(dump_file) # If it's not removed, it will abort the dump
        pop3.dump_to_csv()

        pop4 = sp.SyntheticPopulation(None, 'TestWDConfigDB', 'TestPollutionConfig', dump_file)
        self.assertEqual(expected_values_detectable4, pop4.population[0].modelled_properties)
        os.remove(dump_file)

class PipelineTests(unittest.TestCase):

    def test_init(self):
        population_parameter_dict = {
            'TestPop': {
                sp.PopulationParameter.size: 3,
                sp.PopulationParameter.wd_config: 'TestWDConfig',
                sp.PopulationParameter.pollution_config: 'TestPollutionConfig'
            }
        }
        observer_dict = {
            'TestObs': [
                so.ObservationType.NoCut,
                0
            ]
        }
        modeller_dict = {
            'TestMod': [
                sm.ModellerType.AnalyticApproximation,
                [False, 'synthetic_grid_dummy.csv']
            ]
        }
        pipeline = spi.Pipeline('TestPipeline', population_parameter_dict, observer_dict, modeller_dict)

        self.assertEqual(pipeline.name, 'TestPipeline')
        self.assertEqual(pipeline.population_parameter_dict, population_parameter_dict)
        self.assertEqual(pipeline.population_dict, dict())
        self.assertEqual(pipeline.observer_dict, observer_dict)
        self.assertEqual(pipeline.modeller_dict, modeller_dict)
        self.assertEqual(pipeline.base_dir, pu.get_path_to_pipeline_base_dir())
        self.assertEqual(pipeline.results_dict, dict())

        self.assertEqual(pipeline.get_popdump_dir(), pu.get_path_to_pipeline_base_dir() + 'popdumps/')
        self.assertEqual(pipeline.get_pipeline_name('a', 'b', 'c'), 'a_b_c')
        self.assertEqual(pipeline.get_base_pipeline_dir(), pu.get_path_to_pipeline_base_dir() + 'TestPipeline/')
        self.assertEqual(pipeline.get_pipeline_dir('a', 'b', 'c'), pu.get_path_to_pipeline_base_dir() + 'TestPipeline/a_b_c/')
        self.assertEqual(pipeline.get_variable_type_list(), ['Input', 'Pollution', 'Observed', 'Modelled'])
        mock_plot_description = (('a', 0, 1, 0.025), ('b', -20, -2, 0.1))
        self.assertEqual(pipeline.get_vfilter('Input', mock_plot_description), ('a', 0, 1, 0.025))
        self.assertEqual(pipeline.get_vfilter('Pollution', mock_plot_description), ('b', -np.inf, np.inf, 0.1))
        self.assertEqual(pipeline.get_vfilter('Observed', mock_plot_description), ('b', -20, -2, 0.1))
        self.assertEqual(pipeline.get_vfilter('Modelled', mock_plot_description), ('a', 0, 1, 0.025))

    def test_ks_test(self):
        pipeline = spi.Pipeline('TestPipeline', dict(), dict(), dict())

        data1 = [3.4, 5.6, 9.9]
        data2 = [4.4, 5.5, 6.6]
        min_bin_edge = 2
        max_bin_edge = 12
        half_bin_size = 1
        ks_statistic, p_value = pipeline.run_ks_test_on_samples(data1, data2, min_bin_edge, max_bin_edge, half_bin_size)
        golden_result = st.ks_2samp(data1, data2, "two-sided", "asymp")
        scipy_version = scipy.__version__
        version_parse = scipy_version.split('.')
        major_version = int(version_parse[0])
        minor_version = int(version_parse[1])
        if not (major_version >=1 and minor_version >=3):
            # Skip this test for versions >1.3 - see note in run_ks_test_on_cdf in synthetic_pipeline.py
            self.assertEqual(golden_result[0], ks_statistic)
            self.assertEqual(golden_result[1], p_value)

        mu1 = 7
        sigma1 = 1
        mu2 = 7.1
        sigma2 = 0.9
        sample_size = 100000
        half_small_bin_size = 0.1

        high_pval_failure_count = 0
        high_pval_threshold = 0.5
        high_pval_allowed_failures = 10
        low_pval_failure_count = 0
        low_pval_threshold = 0.0000000001
        low_pval_allowed_failures = 0
        trials = 50
        trial_no = 0
        while trial_no < trials:
            large_data1 = np.random.normal(mu1, sigma1, sample_size)
            large_data2 = np.random.normal(mu1, sigma1, sample_size)
            large_data3 = np.random.normal(mu2, sigma2, sample_size)
            ks_statistic1, p_value_high = pipeline.run_ks_test_on_samples(large_data1, large_data2, min_bin_edge, max_bin_edge, half_bin_size)
            ks_statistic2, p_value_low = pipeline.run_ks_test_on_samples(large_data1, large_data3, min_bin_edge, max_bin_edge, half_bin_size)
            if p_value_high < high_pval_threshold:
                # For identical underlying dists, the pvalue ought to be pretty close to 1. For some reason, this doesn't always seem to be true for ks_2samp
                high_pval_failure_count += 1
            if p_value_low > low_pval_threshold:
                # For such a large sample, and different underlying dists, the pvalue ought to be pretty close to 0. For some reason, this doesn't always seem to be true for ks_2samp
                low_pval_failure_count += 1
            trial_no += 1
        self.assertTrue(high_pval_failure_count <= high_pval_allowed_failures) # Allowing 10 failures as this test is inherantly slightly unreliable
        self.assertTrue(low_pval_failure_count <= low_pval_allowed_failures)

        ks_statistic, p_value = pipeline.run_ks_test_on_samples(data1, data2, min_bin_edge, max_bin_edge, half_bin_size, [1, 2, 3], [0.3, 0.2, 0.5])
        self.assertEqual(0.5, ks_statistic)
        self.assertEqual(0.6823104680778247, p_value)

        data2 = [4.4, 5.5, [6.6, 10.3]]
        ks_statistic, p_value = pipeline.run_ks_test_on_samples(data1, data2, min_bin_edge, max_bin_edge, half_bin_size)
        self.assertEqual(0.3333333333333333, ks_statistic)
        self.assertEqual(0.9762126488644778, p_value)

        data2 = [[2.1, 4.4], [3.1, 5.5, 9.8], [6.6, 10.3]]
        ks_statistic, p_value = pipeline.run_ks_test_on_samples(data1, data2, min_bin_edge, max_bin_edge, half_bin_size)
        self.assertEqual(0.16666666666666663, ks_statistic)
        self.assertEqual(0.9999999955543378, p_value)

        ks_statistic, p_value = pipeline.run_ks_test_on_samples(data1, data1, min_bin_edge, max_bin_edge, half_bin_size)
        self.assertEqual(0, ks_statistic)
        self.assertEqual(1, p_value)

    def test_different_timescales(self):
        # Going to check that we can make multiple populations with different timescales used
        population_parameter_dict = {
            'TestKOPop': {
                sp.PopulationParameter.size: 1,
                sp.PopulationParameter.wd_config: 'TestWDConfigDA', #KO
                sp.PopulationParameter.pollution_config: 'TestPollutionConfig6'
            },
            'TestKNPop': {
                sp.PopulationParameter.size: 1,
                sp.PopulationParameter.wd_config: 'TestWDConfigDA2', #KN
                sp.PopulationParameter.pollution_config: 'TestPollutionConfig6'
            }
        }
        observer_dict = {
            'TestObs': [
                so.ObservationType.NoCut,
                0
            ]
        }
        modeller_dict = {
            'TestKOMod': [
                sm.ModellerType.AnalyticApproximation,
                [ti.TimescaleType.KoesterOvershoot, False]
            ],
            'TestKNMod': [
                sm.ModellerType.AnalyticApproximation,
                [ti.TimescaleType.KoesterNoOvershoot, False]
            ]
        }
        pipeline = spi.Pipeline('TestPipeline', population_parameter_dict, observer_dict, modeller_dict)
        file_to_dump_to1 = pipeline.get_popdump_dir() + 'popdump_' + pipeline.get_pipeline_name('TestKOPop', 'TestObs', 'TestKOMod') + '.csv'
        file_to_dump_to2 = pipeline.get_popdump_dir() + 'popdump_' + pipeline.get_pipeline_name('TestKOPop', 'TestObs', 'TestKNMod') + '.csv'
        file_to_dump_to3 = pipeline.get_popdump_dir() + 'popdump_' + pipeline.get_pipeline_name('TestKNPop', 'TestObs', 'TestKOMod') + '.csv'
        file_to_dump_to4 = pipeline.get_popdump_dir() + 'popdump_' + pipeline.get_pipeline_name('TestKNPop', 'TestObs', 'TestKNMod') + '.csv'

        pipeline.run_all_combinations()

        # I've picked configurations that should avoid [None] metallicity results (which then lead to sampling every single one of the 958 stars...)
        self.assertEqual(1, len(pipeline.results_dict['TestKOMod']['TestObs']['TestKOPop'][0].modelled_properties[mp.ModelParameter.fragment_core_frac]))
        self.assertEqual(1, len(pipeline.results_dict['TestKNMod']['TestObs']['TestKOPop'][0].modelled_properties[mp.ModelParameter.fragment_core_frac]))
        self.assertEqual(1, len(pipeline.results_dict['TestKOMod']['TestObs']['TestKNPop'][0].modelled_properties[mp.ModelParameter.fragment_core_frac]))
        self.assertEqual(1, len(pipeline.results_dict['TestKNMod']['TestObs']['TestKNPop'][0].modelled_properties[mp.ModelParameter.fragment_core_frac]))

        os.remove(file_to_dump_to1)
        os.remove(file_to_dump_to2)
        os.remove(file_to_dump_to3)
        os.remove(file_to_dump_to4)

        self.assertNotEqual(pipeline.results_dict['TestKOMod']['TestObs']['TestKOPop'][0].pollution_abundances[ci.Element.Fe], pipeline.results_dict['TestKOMod']['TestObs']['TestKNPop'][0].pollution_abundances[ci.Element.Fe])
        self.assertNotEqual(pipeline.results_dict['TestKOMod']['TestObs']['TestKNPop'][0].pollution_abundances[ci.Element.Fe], pipeline.results_dict['TestKNMod']['TestObs']['TestKOPop'][0].pollution_abundances[ci.Element.Fe])
        self.assertEqual(pipeline.results_dict['TestKOMod']['TestObs']['TestKNPop'][0].pollution_abundances[ci.Element.Fe], pipeline.results_dict['TestKNMod']['TestObs']['TestKNPop'][0].pollution_abundances[ci.Element.Fe])
        self.assertEqual(pipeline.results_dict['TestKOMod']['TestObs']['TestKOPop'][0].pollution_abundances[ci.Element.Fe], pipeline.results_dict['TestKNMod']['TestObs']['TestKOPop'][0].pollution_abundances[ci.Element.Fe])
        # Should get back the same favoured metallicity as long as the timescales used in Mod and Pop are matched
        self.assertEqual(pipeline.results_dict['TestKOMod']['TestObs']['TestKOPop'][0].modelled_properties[mp.ModelParameter.metallicity], pipeline.results_dict['TestKNMod']['TestObs']['TestKNPop'][0].modelled_properties[mp.ModelParameter.metallicity])

    def test_thermohaline_usage(self):
        # Going to check that we can make multiple populations with different timescales used
        population_parameter_dict = {
            'TestThermPop1': {
                sp.PopulationParameter.size: 3,
                sp.PopulationParameter.wd_config: 'TestDAsThermohalineOff',
                sp.PopulationParameter.pollution_config: 'TestPollutionConfig6'
            },
            'TestThermPop2': {
                sp.PopulationParameter.size: 3,
                sp.PopulationParameter.wd_config: 'TestDAsThermohalineOn',
                sp.PopulationParameter.pollution_config: 'TestPollutionConfig6'
            }
        }
        observer_dict = {
            'TestObs': [
                so.ObservationType.NoCut,
                0
            ]
        }
        modeller_dict = {
            'TestMod': [
                sm.ModellerType.Null,
                None
            ]
        }
        pipeline = spi.Pipeline('TestPipeline', population_parameter_dict, observer_dict, modeller_dict)
        file_to_dump_to1 = pipeline.get_popdump_dir() + 'popdump_' + pipeline.get_pipeline_name('TestThermPop1', 'TestObs', 'TestMod') + '.csv'
        file_to_dump_to2 = pipeline.get_popdump_dir() + 'popdump_' + pipeline.get_pipeline_name('TestThermPop2', 'TestObs', 'TestMod') + '.csv'
        pipeline.run_all_combinations()
        #The thermohaline one should have less Fe (less of everything!)
        print('Removing the following files:')
        print(file_to_dump_to1)
        print(file_to_dump_to2)
        os.remove(file_to_dump_to1)
        os.remove(file_to_dump_to2)
        self.assertTrue(pipeline.results_dict['TestMod']['TestObs']['TestThermPop1'][0].pollution_abundances[ci.Element.Fe] > pipeline.results_dict['TestMod']['TestObs']['TestThermPop2'][0].pollution_abundances[ci.Element.Fe])

    def test_typical_use(self): # This is slightly not the point of unit tests - should ideally split into smaller tests
        population_parameter_dict = {
            'TestPop': {
                sp.PopulationParameter.size: 3,
                sp.PopulationParameter.wd_config: 'TestWDConfigDB',
                sp.PopulationParameter.pollution_config: 'TestPollutionConfig'
            }
        }
        observer_dict = {
            'TestObs': [
                so.ObservationType.NoCut,
                0
            ]
        }
        modeller_dict = {
            'TestMod': [
                sm.ModellerType.AnalyticApproximation,
                [ti.TimescaleType.KoesterOvershoot, False]
            ]
        }
        pipeline = spi.Pipeline('TestPipeline', population_parameter_dict, observer_dict, modeller_dict)
        file_to_dump_to = pipeline.get_popdump_dir() + 'popdump_' + pipeline.get_pipeline_name('TestPop', 'TestObs', 'TestMod') + '.csv'
        pipeline.run_all_combinations()
        os.remove(file_to_dump_to)
        expected_population = sp.SyntheticPopulation(3, 'TestWDConfigDB', 'TestPollutionConfig')
        observer = so.Observer(so.ObservationType.NoCut)
        observer.observe_population(expected_population)
        modeller = sm.Modeller(sm.ModellerType.AnalyticApproximation, [ti.TimescaleType.KoesterOvershoot, False])
        modeller.model_population(expected_population)
        expected_results_dict = {
            'TestMod': {
                'TestObs': {
                    'TestPop': expected_population
                }
            }
        }
        for i, system in enumerate(expected_population):
            self.assertEqual(system.modelled_properties[mp.ModelParameter.fragment_core_frac], pipeline.results_dict['TestMod']['TestObs']['TestPop'].population[0].modelled_properties[mp.ModelParameter.fragment_core_frac])
            self.assertEqual(system.modelled_properties[mp.ModelParameter.formation_distance], pipeline.results_dict['TestMod']['TestObs']['TestPop'].population[0].modelled_properties[mp.ModelParameter.formation_distance])
            self.assertEqual(system.modelled_properties[mp.ModelParameter.metallicity], pipeline.results_dict['TestMod']['TestObs']['TestPop'].population[0].modelled_properties[mp.ModelParameter.metallicity])
            self.assertEqual(system.modelled_properties[mp.ModelParameter.accretion_timescale], pipeline.results_dict['TestMod']['TestObs']['TestPop'].population[0].modelled_properties[mp.ModelParameter.accretion_timescale])
            self.assertEqual(system.modelled_properties[mp.ModelParameter.t_sinceaccretion], pipeline.results_dict['TestMod']['TestObs']['TestPop'].population[0].modelled_properties[mp.ModelParameter.t_sinceaccretion])
            #self.assertEqual(system, pipeline.results_dict['TestMod']['TestObs']['TestPop'].population[0])
        self.assertEqual(pipeline.results_dict, expected_results_dict)
        # Now let's manipulate the results a little to test some features
        expected_population.population[0].modelled_properties[mp.ModelParameter.fragment_core_frac] = [0.1, 0.2, 0.5, 0.8]
        expected_population.population[1].modelled_properties[mp.ModelParameter.fragment_core_frac] = [np.nan]
        expected_population.population[2].modelled_properties[mp.ModelParameter.fragment_core_frac] = [0.55]
        plot_descriptions = {
             # input description, output description (each is variable, min bin value, max bin value, half bin size)
            'fcf': ((mp.ModelParameter.fragment_core_frac, 0, 1, 0.025), (ci.Element.Fe, -20, -2, 0.1))
        }
        variable_dict, weights_dict, variable_below_threshold_dict, variable_above_threshold_dict = pipeline.extract_plottables(
            'TestPop',
            expected_population,
            plot_descriptions
        )
        self.assertEqual(variable_dict['Modelled']['fcf']['TestPop'], [0.1, 0.2, 0.5, 0.8, 0.55])
        self.assertEqual(weights_dict['Modelled']['fcf']['TestPop'], [0.25, 0.25, 0.25, 0.25, 1])
        self.assertEqual(variable_below_threshold_dict['Modelled']['fcf']['TestPop'], 0)
        self.assertEqual(variable_above_threshold_dict['Modelled']['fcf']['TestPop'], 0)

class CompleteModelTests(unittest.TestCase):

    def load_generic_float_data_csv(self, input_filename):
        generic_list = None
        with open(get_path_to_data() + input_filename, encoding='utf-8') as generic_csv:
            generic_list = [row for row in csv.reader(generic_csv)]
        generic_array = np.asarray(generic_list)
        return generic_array.astype(float)

    def test_complete_model(self):
        manager = mn.Manager()
        manager.publish_live_data(0, ti.TimescaleType.KoesterOvershoot)
        args = {
           # Args are: fe_star, t_sinceaccretion, d_formation, z_formation, N_c, N_o, f_c, f_o, fragment_mass, t_disc, pressure, fO2
           # (Skipping the following:          enhancement_model (defaults to 'NonEarthlike'), consider_thermohaline (False), t_formation (1.5), normalise_abundances (True) )
            'CMTEST1': [88, 0, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, 10, 7, 54, -2], # At t = 0, there is no pollution yet so all outputs are -inf
            'CMTEST2': [88, 0.0195758, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, 22, 7, 54, -2],
            'CMTEST3': [88, 5, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, 22, 7, 54, -2],
            'CMTEST4': [88, 15, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, 22, 7, 54, -2],
            'CMTEST5': [88, 0.195758, -0.8, 0.05, 0.17, 0.01, 0.128178, 0.01, np.log10(5E22), 4.773768, 54, -2],
            'CMTEST6': [199, 3.11, -0.3, 0.05, 0.17, 0.01, 0.8, 0.1, np.log10(5E22), np.log10(3000000), 54, -2],
            'CMTEST7': [199, 3.11, -0.3, 0.05, 0.17, 0.01, 0.8, 0, np.log10(5E22), np.log10(3000000), 54, -2],
             # If we set the fcf to the pcnf at a given P/fO2...
            'CMTEST8': [199, 3.11, -0.3, 0.05, 0.17, 0.01, 0.15147328034891958, 0, np.log10(5E22), np.log10(3000000), 54, -2],
             # ...then Changing pressure and/or fO2 should barely, if at all, change the result.
             #(But if you hold fcf constant, the result will change which may seem counterintuitive at first - but remember P/fO2 changes the pcnf)
            'CMTEST9': [199, 3.11, -0.3, 0.05, 0.17, 0.01, 0.11566386300964879, 0, np.log10(5E22), np.log10(3000000), 5, -2]
        }

        expected_results = {
            'CMTEST1': {
                ci.Element.Al: -np.inf,
                ci.Element.Ti: -np.inf,
                ci.Element.Ca: -np.inf,
                ci.Element.Ni: -np.inf,
                ci.Element.Fe: -np.inf,
                ci.Element.Cr: -np.inf,
                ci.Element.Mg: -np.inf,
                ci.Element.Si: -np.inf,
                ci.Element.Na: -np.inf,
                ci.Element.O:  -np.inf,
                ci.Element.C:  -np.inf,
                ci.Element.N:  -np.inf
            },
            'CMTEST2': {
                ci.Element.Al: -8.184854999968964,
                ci.Element.Ti: -9.658907891176113,
                ci.Element.Ca: -8.326851723061791,
                ci.Element.Ni: -8.420100882496675,
                ci.Element.Fe: -7.199948466256751,
                ci.Element.Cr: -9.00892211168204,
                ci.Element.Mg: -7.053117362414797,
                ci.Element.Si: -7.083478217193907,
                ci.Element.Na: -8.5530279995805,
                ci.Element.O:  -5.8521682162849595,
                ci.Element.C:  -6.209698053914535,
                ci.Element.N:  -7.030672490535127
            },
            'CMTEST3': {
                ci.Element.Al: -6.665381093519236,
                ci.Element.Ti: -8.366008559810357,
                ci.Element.Ca: -6.939614404279943,
                ci.Element.Ni: -7.190577473374786,
                ci.Element.Fe: -5.958916190758224,
                ci.Element.Cr: -7.744251552231134,
                ci.Element.Mg: -5.482245243817824,
                ci.Element.Si: -5.56107534186259,
                ci.Element.Na: -6.980172894759156,
                ci.Element.O:  -4.118666039571906,
                ci.Element.C:  -4.374779113819992,
                ci.Element.N:  -5.249715865643942
            },
            'CMTEST4': {
                ci.Element.Al: -10.069581713347697,
                ci.Element.Ti: -14.166162477536453,
                ci.Element.Ca: -11.58310968618205,
                ci.Element.Ni: -13.930236520337578,
                ci.Element.Fe: -12.517002434681643,
                ci.Element.Cr: -13.945025689320138,
                ci.Element.Mg: -8.49926912638722,
                ci.Element.Si: -8.941957052838013,
                ci.Element.Na: -9.98315296648056,
                ci.Element.O:  -6.164298451538868,
                ci.Element.C:  -5.955929082520835,
                ci.Element.N:  -7.067511585456074
            },
            'CMTEST5': {
                ci.Element.Al: -3.2730032119496038,
                ci.Element.Ti: -4.8700830140290465,
                ci.Element.Ca: -3.779861302717375,
                ci.Element.Ni: -4.49442480201954,
                ci.Element.Fe: -3.4032277233846835,
                ci.Element.Cr: -5.5424772482367635,
                ci.Element.Mg: -3.4472389727260344,
                ci.Element.Si: -3.432904040649441,
                ci.Element.Na: -6.901760807925317,
                ci.Element.O:  -2.642549658386914,
                ci.Element.C:  -np.inf,
                ci.Element.N:  -np.inf
            },
            'CMTEST6': {
                ci.Element.Al: -6.140598546767587,
                ci.Element.Ti: -7.780292631076679,
                ci.Element.Ca: -6.259010482970321,
                ci.Element.Ni: -5.3067881506584325,
                ci.Element.Fe: -4.059066372592598,
                ci.Element.Cr: -6.083118838735554,
                ci.Element.Mg: -4.841520932877515,
                ci.Element.Si: -4.327840342625156,
                ci.Element.Na: -6.232083116338056,
                ci.Element.O:  -4.025562898714041,
                ci.Element.C:  -np.inf,
                ci.Element.N:  -np.inf
            },
            'CMTEST7': {
                ci.Element.Al: -6.140598546767587,
                ci.Element.Ti: -7.780292631076679,
                ci.Element.Ca: -6.259010482970321,
                ci.Element.Ni: -5.3067881506584325,
                ci.Element.Fe: -4.059066372592598,
                ci.Element.Cr: -6.083118838735554,
                ci.Element.Mg: -4.841520932877515,
                ci.Element.Si: -4.327840342625156,
                ci.Element.Na: -6.232083116338056,
                ci.Element.O:  -4.025562898714041,
                ci.Element.C:  -np.inf,
                ci.Element.N:  -np.inf
            },
            'CMTEST8': {
                ci.Element.Al: -5.081151884957815,
                ci.Element.Ti: -6.720845969266905,
                ci.Element.Ca: -5.199563821160548,
                ci.Element.Ni: -5.55598049212516,
                ci.Element.Fe: -4.280336174911599,
                ci.Element.Cr: -6.057944193636121,
                ci.Element.Mg: -3.7820742710677404,
                ci.Element.Si: -3.856456320225422,
                ci.Element.Na: -5.172636454528283,
                ci.Element.O:  -3.1203581961138376,
                ci.Element.C:  -np.inf,
                ci.Element.N:  -np.inf
            },
            'CMTEST9': {
                ci.Element.Al: -5.081151884957815,
                ci.Element.Ti: -6.720845969266905,
                ci.Element.Ca: -5.199563821160548,
                ci.Element.Ni: -5.55598049212516,
                ci.Element.Fe: -4.280336174911599,
                ci.Element.Cr: -6.057944193636121,
                ci.Element.Mg: -3.7820742710677404,
                ci.Element.Si: -3.856456320225422,
                ci.Element.Na: -5.172636454528283,
                ci.Element.O:  -3.1203581961138376,
                ci.Element.C:  -np.inf,
                ci.Element.N:  -np.inf
            }
        }

        for arg_name, arg_set in args.items():
            print('CompleteModelTests: Testing arg set ' + arg_name)
            result, diagnostics = cm.complete_model_calculation(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], 10**(arg_set[9]), arg_set[10], arg_set[11])
            expected = expected_results[arg_name]
            print('Comparing the following results for ' + arg_name)
            print(result)
            print(expected)
            print('PCNF = ' + str(diagnostics['Enhancements']['ParentCoreNumberFraction']))
            for el, val in expected.items():
                if np.isnan(val):
                    self.assertTrue(np.isnan(result[el]))
                else:
                    self.assertAlmostEqual(val, result[el])

    def test_complete_model_earthlike(self):
        manager = mn.Manager()
        manager.publish_live_data(0, ti.TimescaleType.KoesterOvershoot)
        args = {
            # Args are: fe_star, t_sinceaccretion, d_formation, z_formation, N_c, N_o, f_c, f_o, pollutionfraction, t_disc, fragment_mass
           # One of them decreases d_formation to cause variation in the hsc chemistry outputs (i.e. not all equal to 1)
            'SDSSJ0002+3209_1': [88, 0.0195758, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768, np.log10(5E22)],
            'SDSSJ0002+3209_2': [88, 0.195758, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768, np.log10(5E22)],
            'SDSSJ0002+3209_3': [88, 0.195758, -0.8, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768, np.log10(5E22)]
        }
        for arg_name, arg_set in args.items():
            test_result, ignore = cm.complete_model_calculation(arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[10], 10**(arg_set[9]), 54, -2, 'Earthlike', False, 1.5, False)
            golden_result = ocm.PWDCodeMultiple(0, 1, False, arg_set[0], arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], arg_set[9])
            #The Earthlike version should in theory replicate John's original code, except for an offset -> ratios should be consistent
            self.assertEqual(12, len(test_result))
            self.assertEqual(3, len(golden_result))
            # We'll check for consistency in the Ca/Fe ratio
            test_cafe = test_result[ci.Element.Ca] - test_result[ci.Element.Fe]
            golden_cafe = golden_result[0] - golden_result[1]
            self.assertAlmostEqual(test_cafe, golden_cafe)

class LoglikeTests(unittest.TestCase):

    def test_bounds(self):
        manager = mn.Manager()
        manager.model_names = ['Model_Full_No_Crust']
        manager.use_hierarchy = False
        manager.parameter_hierarchy = None
        manager.hierarchy_name = None
        manager.load_models()

        timescale_type = ti.TimescaleType.KoesterOvershoot
        consider_thermohaline = True

        manager.publish_live_data(226, timescale_type)  # GD61 has various upper bounds
        manager.publish_live_model('Model_Full_No_Crust', consider_thermohaline)

        dummy_cube = [470, 1, 1.5, 0.02, 0.05, 17, 1, 54, -2]
        likelihood = manager.models[timescale_type][consider_thermohaline]['Model_Full_No_Crust'].loglike(dummy_cube)
        self.assertNotEqual(likelihood, -9e89)  # This should not have triggered the bounds (but it will if C is mistakenly not excluded)

        dummy_cube = [470, 1, 1.5, 0.02, 0.05, 21.2, 1, 54, -2]
        likelihood = manager.models[timescale_type][consider_thermohaline]['Model_Full_No_Crust'].loglike(dummy_cube)
        self.assertEqual(likelihood, -9e89)  # This should have triggered the Al bound (we just raised the fragment mass, putting Al above its upper bound)

        dummy_cube = [470, 1, -1.5, 0.02, 0.05, 21.2, 1, 54, -2] # Trying to trigger a situation where all element values are just nan
        likelihood = manager.models[timescale_type][consider_thermohaline]['Model_Full_No_Crust'].loglike(dummy_cube)
        self.assertEqual(likelihood, -1.1e90)

#@unittest.skip("Can skip tests by uncommenting these decorators")
class IntegrationTests(unittest.TestCase):

    def remove_output_dir(self, output_dir):
        if os.path.exists(output_dir) and os.path.isdir(output_dir):
            shutil.rmtree(output_dir)

    def assert_equal_within_errors(self, value_1, error_1, value_2, error_2):
        # no longer used
        err_msg = str(value_1) + ' not equal to ' + str(value_2) + ' within errors (' + str(error_1) + ', ' + str(error_2) + ')'
        if value_1 >= value_2:
            self.assertGreaterEqual(value_2 + error_2, value_1 - error_1, err_msg)
        else:
            self.assertGreaterEqual(value_1 + error_1, value_2 - error_2, err_msg)

    def test_integration(self):
        test_system_id = 203 # PG0843+516
        manager = mn.Manager()
        manager.n_live_points = 20
        manager.seed = 13022024
        manager.timescale_types_to_run = [ti.TimescaleType.KoesterOvershoot] # For testing purposes, I don't think there's much value in trying different timescales
        manager.thermohaline_regimes_to_run = [False, True] #...but there is value in trying thermohaline mixing on/off
        manager.resume = False

        manager.run([test_system_id])

        self.assertEqual([ti.TimescaleType.KoesterOvershoot], list(manager.models.keys()))
        self.assertTrue(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].best_model)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].result['logZ'], -10.153117770994657)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].result['logZerr'], 0.5445478817689606)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].comparison['HD0']['ln_Z_model'], -14.02636861323811)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].comparison['HD0']['ln_Z_base'], -17.389452481763307)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].comparison['HD0']['chi_model'], 11.692651467896253)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].comparison['HD0']['chi_model_per_data_point'], 1.9487752446493756)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].comparison['HD0']['Bayes_factor_model_base'], 28.87810999582329)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].comparison['HD0']['n_sigma_model_base'], 3.081503086678803)
        self.assertTrue(manager.models[ti.TimescaleType.KoesterOvershoot][False]['HD013'].comparison['HD0']['good_fit'])

        self.assertTrue(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD0'].best_model)
        self.assertTrue(not manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].best_model)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].result['logZ'], -19.023845878911327)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].result['logZerr'], 0.7955094690326314)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].comparison['HD0']['ln_Z_model'], -20.803229019887183)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].comparison['HD0']['ln_Z_base'], -19.079007364602248)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].comparison['HD0']['chi_model'], 15.959903381709513)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].comparison['HD0']['chi_model_per_data_point'], 2.6599838969515854)
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].comparison['HD0']['Bayes_factor_model_base'], 0.17831178581056126)
        self.assertTrue(np.isnan(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].comparison['HD0']['n_sigma_model_base']))
        self.assertTrue(not manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD01'].comparison['HD0']['good_fit'])
        self.assertAlmostEqual(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD0'].comparison['HD0']['Bayes_factor_model_base'], 1.0)
        self.assertTrue(np.isnan(manager.models[ti.TimescaleType.KoesterOvershoot][True]['HD0'].comparison['HD0']['n_sigma_model_base']))

if __name__ == '__main__':
    unittest.main()
