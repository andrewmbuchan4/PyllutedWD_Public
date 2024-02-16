#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
from enum import Enum
import numpy as np

import chemistry_info as ci
import partition_model as pam
import pwd_utils as pu

class Layer(Enum):
    bulk = 0
    core = 1
    mantle = 2
    crust = 3

    def __str__(self):
        return self.name

    def __hash__(self):
        # For performance purposes:
        # Normally python will try to hash an Enum by hashing its name (which guaranteess uniqueness)
        # In this case I know that the values are also unique, and should satisfy properties of a hash
        # So I do this, which is much quicker:
        return self.value

class GeologyModel():

    # normalise_abundances is only ever set to False for testing purposes.
    def __init__(self, bulk_abundance_overrides=None, normalise_abundances=True):
        self.earth_layer_number_fractions = {
            Layer.core: 0.17,
            Layer.mantle: 0.829,
            Layer.crust: 0.001
        }

        self.apply_bulk_abundances(bulk_abundance_overrides)
        if normalise_abundances:
            self.normalise_non_mantle_abundances()
        self.fill_in_mantle_abundances()
        self.initialise_stellar_info_for_graphs()

        # These filenames should probably be arguments to GeologyModel.__init__()
        self.pamela = pam.PartitionModel(
            pu.get_path_to_feni() + 'data/part_param_fischer_blanchard_epsilon_update.dat',
            pu.get_path_to_feni() + 'data/int_param_fischer_blanchard_update.dat',
            pu.get_path_to_feni() + 'data/composition.dat',
            pu.get_path_to_feni() + 'data/e_param_fischer_epsilon_update.dat'
        )
        self.config_name = 'sisi'
        self.G = 6.67408E-11  # Gravitational constant (SI units)
        self.M_Earth = 5.9736E24  # Mass of Earth in kg (McDonough 2003)
        self.M_Sun = 1.9885E30 # Mass of Sun in kg. Maybe solar/stellar info could have its own class
        self.mars_abundances = { # Taking this from table 10 of Yoshizaki 2020
            ci.Element.Al: {Layer.bulk: 0.013, Layer.mantle: 0.015, Layer.core: 0},
            ci.Element.Ti: {Layer.bulk: 0, Layer.mantle: 0, Layer.core: 0},
            ci.Element.Ca: {Layer.bulk: 0.01, Layer.mantle: 0.011, Layer.core: 0},
            ci.Element.Ni: {Layer.bulk: 0.005, Layer.mantle: 0.0001, Layer.core: 0.042},
            ci.Element.Fe: {Layer.bulk: 0.1, Layer.mantle: 0.044, Layer.core: 0.48},
            ci.Element.Cr: {Layer.bulk: 0, Layer.mantle: 0, Layer.core: 0},
            ci.Element.Mg: {Layer.bulk: 0.15, Layer.mantle: 0.17, Layer.core: 0},
            ci.Element.Si: {Layer.bulk: 0.14, Layer.mantle: 0.16, Layer.core: 0},
            ci.Element.Na: {Layer.bulk: 0, Layer.mantle: 0, Layer.core: 0},
            ci.Element.O: {Layer.bulk: 0.53, Layer.mantle: 0.59, Layer.core: 0.11},
            ci.Element.C: {Layer.bulk: 0, Layer.mantle: 0, Layer.core: 0},
            ci.Element.N: {Layer.bulk: 0, Layer.mantle: 0, Layer.core: 0}
        }
        self.core_relative_mass = None
        self.mantle_relative_mass = None

    def reinit(self, bulk_abundance_overrides=None):
        self.apply_bulk_abundances(bulk_abundance_overrides)
        self.normalise_non_mantle_abundances()
        self.fill_in_mantle_abundances()

    def initialise_stellar_info_for_graphs(self):
        # Solar info - this could potentially be split into its own class
        local_stellar_abundances = {
            # These are log number abundances relative to Mg
            ci.Element.C: 0.8507,
            ci.Element.N: 0.3113,
            ci.Element.O: 1.2284,
            ci.Element.Na: -1.346,
            ci.Element.Mg: 0,
            ci.Element.Al: -1.145,
            ci.Element.Si: -0.011,
            ci.Element.Ca: -1.169,
            ci.Element.Ti: -2.587,
            ci.Element.Cr: -1.868,
            ci.Element.Fe: -0.047,
            ci.Element.Ni: -1.289
        }
        # solar_ratiod_to_stellar =  log(  (X/Mg)solar / (X/Mg)stellar  )     = log(X/Mg)solar - log(X/Mg)stellar
        # => solar_abundance = log(X/Mg)solar = solar_ratiod_to_stellar + log(X/Mg)stellar
        solar_ratiod_to_stellar = {
            ci.Element.Al: -0.015854489,
            ci.Element.Ti: -0.042910874,
            ci.Element.Ca: -0.05177587,
            ci.Element.Ni: -0.011069397,
            ci.Element.Fe: -0.032526231,
            ci.Element.Cr: -0.021908295,
            ci.Element.Mg: 0,  # I guess this is technically undefined (0/0), but it's 0 for our purposes
            ci.Element.Si: -0.00878137,
            ci.Element.Na: -0.013652197,
            ci.Element.O: -0.098731306,
            ci.Element.C: 0.008264847,
            ci.Element.N: -0.061159358
        }
        # upper_ratiod_to_stellar = log(  (X/Mg)upper / (X/Mg)stellar  ) = log(X/Mg)upper - log(X/Mg)stellar
        # We want log(X/Mg)upper - log(X/Mg)solar
        # So add log(X/Mg)stellar - log(X/Mg)solar
        # or equivalently, subtract solar_ratiod_to_stellar
        # ditto for lower
        upper_ratiod_to_stellar = {
            ci.Element.Al: 0.114145511,
            ci.Element.Ti: 0.131417374,
            ci.Element.Ca: 0.152552378,
            ci.Element.Ni: 0.138930603,
            ci.Element.Fe: 0.151802017,
            ci.Element.Cr: 0.132419952,
            ci.Element.Mg: 0,
            ci.Element.Si: 0.095546878,
            ci.Element.Na: 0.230676051,
            ci.Element.O: 0.281268694,
            ci.Element.C: 0.172593095,
            ci.Element.N: 0.318840642
        }
        lower_ratiod_to_stellar = {
            ci.Element.Al: -0.300126301,
            ci.Element.Ti: -0.062910874,
            ci.Element.Ca: -0.13177587,
            ci.Element.Ni: -0.161069397,
            ci.Element.Fe: -0.252526231,
            ci.Element.Cr: -0.251908295,
            ci.Element.Mg: 0,
            ci.Element.Si: -0.12878137,
            ci.Element.Na: -0.203652197,
            ci.Element.O: -0.228731306,
            ci.Element.C: -0.234382055,
            ci.Element.N: -0.321159357
        }
        # These are the raw [X/H] ratios from https://docs.google.com/spreadsheets/d/1B2WXfFtx3KhP4Ucl2kmJGmhq73ldJqqZiofG8m10n48
        self.upper_XH_ratiod_to_solar = {
            ci.Element.C: 0.3745,
            ci.Element.N: 0.5445,
            ci.Element.O: 0.45,
            ci.Element.Na: 0.57,
            ci.Element.Mg: 0.34,
            ci.Element.Al: 0.44,
            ci.Element.Si: 0.38,
            ci.Element.Ca: 0.4,
            ci.Element.Ti: 0.3741,
            #ci.Element.V: 0.3545,
            ci.Element.Cr: 0.3841,
            #ci.Element.Mn: 0.4841,
            ci.Element.Fe: 0.4,
            ci.Element.Ni: 0.45,
            #ci.Element.Y: 0.5141
        }
        self.lower_XH_ratiod_to_solar = {
            ci.Element.C: -0.42,
            ci.Element.N: -0.54,
            ci.Element.O: -0.2145,
            ci.Element.Na: -0.5141,
            ci.Element.Mg: -0.39,
            ci.Element.Al: -0.45,
            ci.Element.Si: -0.4041,
            ci.Element.Ca: -0.42,
            ci.Element.Ti: -0.3341,
            #ci.Element.V: -0.4045,
            ci.Element.Cr: -0.5541,
            #ci.Element.Mn: -0.7582,
            ci.Element.Fe: -0.5182,
            ci.Element.Ni: -0.4982,
            #ci.Element.Y: -0.5482
        }
        self.solar_abundances = dict()
        self.upper_ratiod_to_solar = dict()
        self.lower_ratiod_to_solar = dict()
        self.solar_ratiod_to_H = {  # Asplund 2021
            ci.Element.H: 12,
            ci.Element.He: 10.914,
            ci.Element.Al: 6.43,
            ci.Element.Ti: 4.97,
            ci.Element.Ca: 6.3,
            ci.Element.Ni: 6.2,
            ci.Element.Fe: 7.46,
            ci.Element.Cr: 5.62,
            ci.Element.Mg: 7.55,
            ci.Element.Si: 7.51,
            ci.Element.Na: 6.22,
            ci.Element.O: 8.69,
            ci.Element.C: 8.46,
            ci.Element.N: 7.83
        }
        for el in self.solar_ratiod_to_H:
            self.solar_ratiod_to_H[el] -= 12  # These are given with an offset of 12, which we need to remove
        for element in ci.usual_elements:
            self.solar_abundances[element] = local_stellar_abundances[element] + solar_ratiod_to_stellar[element]
            self.upper_ratiod_to_solar[element] = upper_ratiod_to_stellar[element] - solar_ratiod_to_stellar[element]
            self.lower_ratiod_to_solar[element] = lower_ratiod_to_stellar[element] - solar_ratiod_to_stellar[element]

    def apply_bulk_abundances(self, bulk_abundance_overrides=None):
        if bulk_abundance_overrides is None:
            self.element_info = {
                # McDonough 2003
                # NB: These are not normalised! Use the loaded normalised version for reference
                # NB: The Fe abundance in Table 5 of McDonough 2003 was written as 0.49 but it should actually be 0.149
                ci.Element.Al: {Layer.bulk: 0.0153, Layer.crust: 0.0641, Layer.core: 0},
                ci.Element.Ti: {Layer.bulk: 0.00044, Layer.crust: 0.0042, Layer.core: 0},
                ci.Element.Ca: {Layer.bulk: 0.011099, Layer.crust: 0.04452, Layer.core: 0},
                ci.Element.Ni: {Layer.bulk: 0.008066, Layer.crust: 0.0000371, Layer.core: 0.0444},
                ci.Element.Fe: {Layer.bulk: 0.149057, Layer.crust: 0.0314, Layer.core: 0.7676},
                ci.Element.Cr: {Layer.bulk: 0.002351, Layer.crust: 0.000139, Layer.core: 0.00868},
                ci.Element.Mg: {Layer.bulk: 0.16482, Layer.crust: 0.04167, Layer.core: 0},
                ci.Element.Si: {Layer.bulk: 0.149118, Layer.crust: 0.181509, Layer.core: 0.1071},
                ci.Element.Na: {Layer.bulk: 0.002037, Layer.crust: 0.01771, Layer.core: 0},
                ci.Element.O: {Layer.bulk: 0.482879, Layer.crust: 0.6011, Layer.core: 0},
                ci.Element.C: {Layer.bulk: 0.001581, Layer.crust: 0, Layer.core: 0.0083482},
                ci.Element.N: {Layer.bulk: 0.000046429, Layer.crust: 0, Layer.core: 0.0002684}, # NB: This implies a negative mantle abundance! Luckily it doesn't matter in actual runs
                ci.Element.S: {Layer.bulk: 0.005, Layer.core: 0.03},
                ci.Element.Hf: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Hf, 0.00000019)},
                ci.Element.U: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.U, 0.000000015)},
                ci.Element.Ta: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Ta, 0.000000025)},
                ci.Element.Nb: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Nb, 0.00000044)},
                ci.Element.Pb: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Pb, 0.00000023)},
                ci.Element.Zn: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Zn, 0.00004)},
                ci.Element.Mn: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Mn, 0.0008)},
                ci.Element.Ga: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Ga, 0.000003)},
                ci.Element.V: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.V, 0.000105)},
                ci.Element.Cu: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Cu, 0.00006)},
                ci.Element.W: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.W, 0.00000017)},
                ci.Element.P: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.P, 0.000715)},
                ci.Element.Co: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Co, 0.00088)}

                #ci.Element.Pt: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Pt, 0.0000019)},
                #ci.Element.Pd: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Pd, 0.000001)},
                #ci.Element.Rh: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Rh, 0.00000024)},
                #ci.Element.Re: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Re, 0.000000075)},
                #ci.Element.Ru: {Layer.bulk: self.convert_bulk_abundance_mass_to_number(ci.Element.Ru, 0.0000013)},
            }
        else:
            self.element_info = dict()
            for el, override in bulk_abundance_overrides.items():
                self.element_info[el] = {Layer.bulk: override}
        self.convergence_elements = list()
        for element, abundances in self.element_info.items():
            if abundances[Layer.bulk] > 0:
                self.convergence_elements.append(element)

    def normalise_non_mantle_abundances(self):
        total_bulk_abundance = 0
        total_core_abundance = 0
        total_crust_abundance = 0
        for element in self.element_info.keys():
            total_bulk_abundance += self.element_info[element][Layer.bulk]
            try:
                total_core_abundance += self.element_info[element][Layer.core]
            except KeyError:
                pass
            try:
                total_crust_abundance += self.element_info[element][Layer.crust]
            except KeyError:
                pass
        for element in self.element_info.keys():
            self.element_info[element][Layer.bulk] /= total_bulk_abundance
            try:
                self.element_info[element][Layer.core] /= total_core_abundance
            except KeyError:
                pass
            try:
                self.element_info[element][Layer.crust] /= total_crust_abundance
            except KeyError:
                pass

    def convert_bulk_abundance_mass_to_number(self, element, mass_abundance):
        return (mass_abundance*self.get_earth_mean_molecular_weight())/ci.get_element_mass(element)

    def get_earth_differentiation_pressure(self):
        return 54 #GPa

    def get_earth_oxygen_fugacity(self):
        return -2 # i.e. 2 log units below IW buffer

    def get_mars_differentiation_pressure(self):
        return 13 #GPa. See Rai/van Westrenen 2013

    def get_mars_oxygen_fugacity(self):
        return -1 # i.e. 1 log unit below IW buffer. See Rai/van Westrenen 2013 (although they also suggest -1.3, but then just use -1 in their calcs. and conclusion)

    def get_earth_mean_molecular_weight(self):
        # Estimate based on Table 5 in McDonough 2003
        # Comparing number to weight abundances gives an estimate of mmw: mass_abundance = (element_mass * number_abundance) / mmw
        return 26

    def get_earth_layer_number_fraction(self, layer):
        return self.earth_layer_number_fractions[layer]

    def get_bulk_abundance(self, element):
        try:
            return self.element_info[element][Layer.bulk]
        except KeyError:
            # Presumably the element is not present
            return None

    def get_core_abundance(self, element):
        try:
            return self.element_info[element][Layer.core]
        except KeyError:
            # Presumably the element is not present
            return None

    def get_mantle_abundance(self, element):
        try:
            return self.element_info[element][Layer.mantle]
        except KeyError:
            # Presumably the element is not present
            return None

    def get_crust_abundance(self, element):
        try:
            return self.element_info[element][Layer.crust]
        except KeyError:
            # Presumably the element is not present
            return None

    def fill_in_mantle_abundances(self):
        for element in self.element_info:
            try:
                self.element_info[element][Layer.mantle] = (self.element_info[element][Layer.bulk] - ((self.element_info[element][Layer.core] * self.get_earth_layer_number_fraction(Layer.core)) + (self.element_info[element][Layer.crust] * self.get_earth_layer_number_fraction(Layer.crust))))/self.get_earth_layer_number_fraction(Layer.mantle)
            except KeyError as e:
                # Then I assume we don't have enough information to calculate an observed mantle abundance (eg Hf, U etc)
                pass

    def get_mu(self, element):
        # At some point this function could become an element-dependent function saying what proportion of it gets a chance to react with the metal
        return 1

    def find_system_specific_abundances(self, abundances, core_number_fraction):
        # Input: A set of core/mantle abundances and a core_number_fraction
        # Output: The bulk abundances of a body which forms with the specified cnf and core/mantle composition
        output = dict()
        mantle_number_fraction = 1 - core_number_fraction
        for element, abundance in abundances.items():
            output[element] = core_number_fraction*abundance[Layer.core] + mantle_number_fraction*abundance[Layer.mantle]
        return output

    def find_Ds_deviation(self, test_Ds, new_Ds, important_elements):
        tot_err_sq = 0
        for element in important_elements:
            tot_err_sq += ((new_Ds[element] - test_Ds[element])/test_Ds[element])**2
        return tot_err_sq

    def grid_method(self, Ds, core_number_fraction, pressure, fO2, abundances, temp, nbot):
        important_elements = [ci.Element.Fe, ci.Element.Si, ci.Element.Ni, ci.Element.O, ci.Element.Cr, ci.Element.C]  # If we do this for every element, it'll take way too long
        grid_points = [0.9, 0.95, 1, 1.05, 1.1]
        grid_points_dict = dict()
        for element in important_elements:
            grid_points_dict[element] = [g*Ds[element] for g in grid_points]
        min_deviation = np.inf
        solution_Ds = None
        solution_abundances = None
        solution_cnf = None
        for c_val in grid_points_dict[ci.Element.C]:
            for cr_val in grid_points_dict[ci.Element.Cr]:
                for o_val in grid_points_dict[ci.Element.O]:
                    for ni_val in grid_points_dict[ci.Element.Ni]:
                        for si_val in grid_points_dict[ci.Element.Si]:
                            for fe_val in grid_points_dict[ci.Element.Fe]:
                                test_Ds = dict()
                                for el in Ds.keys():
                                    test_Ds[el] = Ds[el]
                                test_Ds[ci.Element.Fe] = fe_val
                                test_Ds[ci.Element.Si] = si_val
                                test_Ds[ci.Element.Ni] = ni_val
                                test_Ds[ci.Element.O] = o_val
                                test_Ds[ci.Element.Cr] = cr_val
                                test_Ds[ci.Element.C] = c_val
                                abundances, cnf = self.calculate_abundances(test_Ds, core_number_fraction, True)
                                new_Ds = self.pamela.get_all_partition_coefficients(pressure, fO2, abundances, temp, nbot)
                                deviation = self.find_Ds_deviation(test_Ds, new_Ds, important_elements)
                                if deviation < min_deviation:
                                    # I assume there will be no ties
                                    min_deviation = deviation
                                    solution_abundances = abundances
                                    solution_cnf = cnf
                                    solution_Ds = test_Ds
        return solution_abundances, solution_cnf, solution_Ds

    def calculate_abundances(self, Ds, core_number_fraction, calculate_mantle=False):
        mantle_number_fraction = 1 - core_number_fraction
        abundances = dict()
        for element, D in Ds.items():
            abundances[element] = dict()
            bulk_abundance_by_number = self.get_bulk_abundance(element)
            if bulk_abundance_by_number is None:
                # TODO: It's pointless to calculate any Ds for which the bulk abundance is 0. Could be an easy performance gain
                bulk_abundance_by_number = 0
            abundances[element][Layer.bulk] = bulk_abundance_by_number
            mu = self.get_mu(element)
            if element == ci.Element.Placeholder:  # This used to be special logic for Oxygen
                # Then "D" is not actually D
                # This calculation assumes mu = 1...
                abundances[element][Layer.core] = D
                if calculate_mantle:
                    abundances[element][Layer.mantle] = (bulk_abundance_by_number - (D*core_number_fraction))/(1 - core_number_fraction)
            else:
                metal_abundance = mu * bulk_abundance_by_number * ((D * core_number_fraction) / (mantle_number_fraction + (D * core_number_fraction)))
                abundances[element][Layer.core] = metal_abundance
                if calculate_mantle:
                    # We only need to calculate mantle stuff on the final iteration, otherwise it just wastes time (unless S is present...)
                    sil_abundance = bulk_abundance_by_number * (1 - (mu * (1 - (mantle_number_fraction / (mantle_number_fraction + (D * core_number_fraction))))))
                    abundances[element][Layer.mantle] = sil_abundance
        total_mantle_abundance = 0
        total_core_abundance = 0

        # Badro/Siebert O
        for e in abundances.keys():
            total_core_abundance += abundances[e][Layer.core]
            if calculate_mantle:
                total_mantle_abundance += abundances[e][Layer.mantle]
        inv_core_abundance = 1/total_core_abundance
        if calculate_mantle:
            inv_mantle_abundance = 1/total_mantle_abundance
        for e in abundances.keys():
            abundances[e][Layer.core] *= inv_core_abundance
            if calculate_mantle:
                abundances[e][Layer.mantle] *= inv_mantle_abundance

        # Fischer O
        # This logic keeps O the same, normalises everything else. Reinstate if using old O calculation
        #for e in abundances.keys():
        #    if e != ci.Element.O:
        #        total_core_abundance += abundances[e][Layer.core]
        #        if calculate_mantle:
        #            total_mantle_abundance += abundances[e][Layer.mantle]
        #inv_core_abundance = 1/total_core_abundance
        #if calculate_mantle:
        #    inv_mantle_abundance = 1/total_mantle_abundance
        #for e in abundances.keys():
        #    if e != ci.Element.O:
        #        abundances[e][Layer.core] *= (1 - abundances[ci.Element.O][Layer.core])*inv_core_abundance
        #        if calculate_mantle:
        #            abundances[e][Layer.mantle] *= (1 - abundances[ci.Element.O][Layer.mantle])*inv_mantle_abundance
        #total_core_abundance /= (1 - abundances[ci.Element.O][Layer.core])
        return abundances, total_core_abundance

    def check_convergence(self, Ds, prev_Ds, tolerance=0.01):
        #important_elements = None # This will check all elements
        important_elements = [el for el in self.convergence_elements if el in Ds]
        if Ds is None or prev_Ds is None:
            return False
        if important_elements is None:
            important_elements = Ds.keys()
        for element in important_elements:
            D = Ds[element]
            prev_D = prev_Ds[element]
            if prev_D == 0:
                if D == 0:
                    pass
                else:
                    return False
            else:
                rel_diff = abs((D - prev_D)/prev_D)
                if rel_diff >= tolerance:
                    return False
        return True

    def form_a_planet_iteratively(self, pressure, fO2, temp=None, nbot=None, initial_w_met=None, initial_Ds=None, return_all_Ds=False):  # ..with no crust
        # Setup: Don't touch!
        w_met = initial_w_met if initial_w_met is not None else 0.5
        w_sil = 1 - w_met
        abundances = None
        converged = False
        prev_Ds = None
        iteration_count = 0
        average_next_iteration = False
        all_Ds = list()
        o_increased_last_time = False
        o_cap_set = False

        # Internal parameters: Do touch
        nudge_iterations = 1000  # Every so often, nudge the algorithm in the right direction to help convergence
        grid_iterations = 1001 # Invoke grid solution after this many iterations (if O increases 2 times in a row)
        max_iterations = 1000 # Give up eventually
        double_o_detection = False
        cap_o = True
        cap_si = True
        o_cap = 0.3 # Justification: O and Si should never be highly siderophilic (i.e. D much greater than 1) otherwise no silicate is left
        si_cap = 1 #               and the whole exercise loses its meaning! To be conservative, set D_Si <= 1 and D_O <= 1. But actually D_O = 1 is still a bit silly and leads to sudden X_O spikes so I tuned it to 0.3

        #NB: if using a parametrisation with direct calculation of X_O rather than D_O, the o_cap should be dropped! (to something in the region of 0.2 ish)

        movement_fraction = 1  # Setting this to 0.5 is the same as setting nudge_iterations to 1, in principle
        while not converged:
            # If using the old Oxygen logic: the Ds value for oxygen is NOT actually D_O, it's just skipped straight to abundances[O][Layer.core] and that value is stored in Ds
            Ds = self.pamela.get_all_partition_coefficients(pressure, fO2, abundances, temp, nbot)
            if movement_fraction != 1 and prev_Ds is not None:
                for element in Ds.keys():
                    Ds[element] = prev_Ds[element] + movement_fraction*(Ds[element] - prev_Ds[element])
            if initial_Ds is not None and iteration_count == 0:
                if set(initial_Ds.keys()) != set(self.pamela.ele_set):
                    print('Warning! Initial Ds do not have the same set of elements as PAMELA')
                for ele in self.pamela.ele_set:
                    try:
                        Ds[ele] = initial_Ds[ele]
                    except KeyError:
                        pass
            if double_o_detection and not o_cap_set and o_increased_last_time and Ds[ci.Element.O] > prev_Ds[ci.Element.O]:
                cap_o = True
                o_cap = (prev_Ds[ci.Element.O]**2)/Ds[ci.Element.O]  # i.e. reverse the direction in log space
                print('Runaway detected - setting O cap to ' + str(o_cap))
                o_cap_set = True
            if average_next_iteration and prev_Ds is not None:
                for element in Ds.keys():
                    Ds[element] = (Ds[element] + prev_Ds[element])/2
            if cap_si:
                Ds[ci.Element.Si] = min(Ds[ci.Element.Si], si_cap)
            if cap_o:
                Ds[ci.Element.O] = min(Ds[ci.Element.O], o_cap)
            if return_all_Ds:
                all_Ds.append(Ds)

            calculate_mantle = False
            if abundances is None:
                calculate_mantle = True
            elif abundances[ci.Element.S].get(Layer.core, 0) != 0 or abundances[ci.Element.S].get(Layer.mantle, 0) != 0:
                calculate_mantle = True
            abundances, w_met = self.calculate_abundances(Ds, w_met, calculate_mantle)

            if not average_next_iteration:
                # To prevent 'cheating', we are only allowed to declare convergence if we didn't do a nudge this iteration
                # This way convergence is only declared if the Ds are genuinely a solution (i.e. they didn't significantly change for a whole, legit, iteration)
                converged = self.check_convergence(Ds, prev_Ds)
            else:
                average_next_iteration = False

            w_sil = 1 - w_met
            iteration_count += 1
            if iteration_count % nudge_iterations == 0:
                # Next iteration, we'll give the algorithm a nudge
                average_next_iteration = True
            if iteration_count > grid_iterations and o_increased_last_time:
                print('Warning! Failed to converge normally within ' + str(grid_iterations) + ' iterations, using grid method')
                abundances, w_met, Ds = self.grid_method(prev_Ds, w_met, pressure, fO2, abundances, temp, nbot)
                if return_all_Ds:
                    all_Ds.append(Ds)
                return abundances, w_met, Ds, all_Ds if return_all_Ds else None
            if iteration_count > max_iterations:
                print('Warning! Failed to converge!')
                return None, None, None, None
            if prev_Ds is not None:
                o_increased_last_time = Ds[ci.Element.O] > prev_Ds[ci.Element.O]
            prev_Ds = Ds
        abundances, w_met = self.calculate_abundances(Ds, w_met, True)
        return abundances, w_met, Ds, all_Ds if return_all_Ds else None

    def tabulate_output(self, abundances, w_met, Ds, outfile=None):
        print('Core Number Fraction     |       ' + str(w_met))
        print('')
        print('Element | D')
        for element, D in Ds.items():
            print(str(element).ljust(7) + ' | ' + str(D))
        print('')
        print('Element | Mantle Number Fraction  | Core Number Fraction')
        for element, el_abundances in abundances.items():
            print(str(element).ljust(7) + ' | ' + str(el_abundances[Layer.mantle]).ljust(22) + '  | ' + str(el_abundances[Layer.core]))
        print('')
        if outfile is not None:
            with open(outfile + '.csv', 'w', newline='', encoding='utf-8') as f:
                to_write = csv.writer(f)
                to_write.writerow(['Core Number Fraction', w_met])
                to_write.writerow([])
                to_write.writerow(['Element', 'D'])
                for element, D in Ds.items():
                    to_write.writerow([element, D])
                to_write.writerow([])
                to_write.writerow(['Element', 'Mantle Number Fraction', 'Core Number Fraction'])
                for element, el_abundances in abundances.items():
                    to_write.writerow([element, el_abundances[Layer.mantle], el_abundances[Layer.core]])
                to_write.writerow([])

    #Make a @staticmethod?
    def convert_number_abundances_to_mass(self, number_abundance_dict):
        # Take a number abundance dict and return the equivalent mass abundance dict
        mass_abundance_dict = dict()
        for element in number_abundance_dict.keys():
            mass_abundance_dict[element] = dict()
        for layer in Layer:
            total_mass = 0
            for element in number_abundance_dict.keys():
                X_E = number_abundance_dict[element].get(layer, 0)  # Assume the abundance of any absent element is 0
                M_E = ci.get_element_mass(element)
                total_mass_contribution = X_E * M_E
                total_mass += total_mass_contribution
            # Now loop through again to build the output dict
            if total_mass > 0:  # If total_mass was zero, then just ignore that layer (Layer was probably not present in the input)
                for element in number_abundance_dict.keys():
                    X_E = number_abundance_dict[element].get(layer, 0)
                    M_E = ci.get_element_mass(element)
                    mass_abundance_dict[element][layer] = (X_E * M_E)/total_mass
        return mass_abundance_dict

    #Make a @staticmethod?
    def convert_mass_abundances_to_number(self, mass_abundance_dict):
        # Take a mass abundance dict and return the equivalent number abundance dict
        number_abundance_dict = dict()
        for element in mass_abundance_dict.keys():
            number_abundance_dict[element] = dict()
        for layer in Layer:
            total_number = 0
            for element in mass_abundance_dict.keys():
                M_E = mass_abundance_dict[element].get(layer, 0)  # Assume the abundance of any absent element is 0
                X_E = M_E/ci.get_element_mass(element)
                total_number += X_E
            # Now loop through again to build the output dict
            if total_number > 0:  # If total_number was zero, then just ignore that layer (Layer was probably not present in the input)
                for element in mass_abundance_dict.keys():
                    M_E = mass_abundance_dict[element].get(layer, 0)  # Assume the abundance of any absent element is 0
                    X_E = M_E/ci.get_element_mass(element)
                    number_abundance_dict[element][layer] = X_E/total_number
        return number_abundance_dict

    def get_relative_core_mantle_masses(self, composition_override=None):
        if composition_override is None and self.core_relative_mass is not None and self.mantle_relative_mass is not None:
            return self.core_relative_mass, self.mantle_relative_mass
        core_relative_mass = 0
        mantle_relative_mass = 0
        for element in ci.usual_elements:  # Only considering major elements for now. Also may need to exclude N, currently its mantle content is negative...
            if composition_override is None:
                info = self.element_info[element]
            else:
                info = composition_override[element]
            core_relative_mass += info[Layer.core]*ci.element_masses[element]
            mantle_relative_mass += info[Layer.mantle]*ci.element_masses[element]
        if composition_override is None:
            # Cache the default, non-override case if we get the chance:
            self.core_relative_mass = core_relative_mass
            self.mantle_relative_mass = mantle_relative_mass
        return core_relative_mass, mantle_relative_mass

    def convert_core_number_fraction_to_core_mass_fraction(self, core_number_fraction, composition_override=None):
        #Need to calculate relative mass of core and mantle like material
        core_relative_mass, mantle_relative_mass = self.get_relative_core_mantle_masses(composition_override)
        mantle_number_fraction = 1 - core_number_fraction
        cmf = (core_number_fraction*core_relative_mass)/((core_number_fraction*core_relative_mass) + (mantle_number_fraction*mantle_relative_mass))
        return cmf

    def convert_core_mass_fraction_to_core_number_fraction(self, core_mass_fraction, composition_override=None):
        #Need to calculate relative mass of core and mantle like material
        core_relative_mass, mantle_relative_mass = self.get_relative_core_mantle_masses(composition_override)
        cnf = (mantle_relative_mass*core_mass_fraction)/(((mantle_relative_mass*core_mass_fraction) + core_relative_mass) - (core_relative_mass*core_mass_fraction))
        return cnf

    def get_core_fraction_from_abundance_dict(self, abundance_dict):
        # abundance_dict could be by number, or by mass. The output will then be either number or mass accordingly
        # also the mantle frac is then 1 - core_fraction
        # This does require mantle abundances to be provided
        tol = 0.00000001
        implied_cfs = list()
        for element, abundance_by_layer in abundance_dict.items():
            if abundance_by_layer[Layer.bulk] > 0:
                trial_cf = 0.5
                diff = tol + 1
                increased = False
                decreased = False
                step = 0.1
                while abs(diff) >= tol:
                    implied_bulk_abundance = trial_cf*abundance_by_layer[Layer.core] + (1-trial_cf)*abundance_by_layer[Layer.mantle]
                    diff = implied_bulk_abundance - abundance_by_layer[Layer.bulk]
                    if abs(diff) < tol:
                        step = 0
                    if (diff > 0 and (abundance_by_layer[Layer.core] > abundance_by_layer[Layer.mantle])) or (diff < 0 and (abundance_by_layer[Layer.core] < abundance_by_layer[Layer.mantle])):
                        trial_cf -= step
                        increased = True
                    else:
                        trial_cf += step
                        decreased = True
                    if increased and decreased:
                        step /= 2
                        increased = False
                        decreased = False
                implied_cfs.append(trial_cf)
        return np.mean(implied_cfs)

    def calculate_pressure_and_mass(self, planet_radius, iron_mass_fraction, lower_limit=False, alpha=1067.44, beta=0.329, gamma=0.31, lmbda=7008.42, mu=18.29, nu=0.313):
        # An implementation of Lena's new parametrisation
        # Planet radius in km, iron mass fraction on a scale from 0 to 1. NB NOT NUMBER FRACTION!
        planet_mass = (planet_radius/(lmbda - (100*mu*iron_mass_fraction)))**(1/nu)  # In Earth masses
        planet_mass_kg = planet_mass*self.M_Earth # In kg
        core_radius = alpha*((100*iron_mass_fraction)**beta)*(planet_mass**gamma) # In km
        planet_radius_metres = 1000*planet_radius
        core_radius_metres = 1000*core_radius
        g_surface = (self.G*planet_mass_kg)/(planet_radius_metres**2)  # In m/s
        g_core = (self.G*planet_mass_kg*iron_mass_fraction)/(core_radius_metres**2) # In m/s
        density_mantle = (3*(1-iron_mass_fraction)*planet_mass_kg)/(4*np.pi*((planet_radius_metres**3) - (core_radius_metres**3))) # In kg/m^3
        if lower_limit:
            # Then we are going to find a lower limit on radius/mass by returning the CMB pressure in this function instead of the (more realistic) mid-mantle pressure
            g_mantle = (g_surface/2) + (g_core/2) # In m/s
            depth_CMB_metres = planet_radius_metres - core_radius_metres
            pressure_CMB = (g_mantle * density_mantle * depth_CMB_metres)/1000000000 # In GPa
            return pressure_CMB, planet_mass
        else:
            g_uppermantle = ((3*g_surface)/4) + (g_core/4) # In m/s
            density_surface = 3100
            density_uppermantle = 0.5*(density_mantle + density_surface)
            depth_midmantle_metres = planet_radius_metres - (0.5*(planet_radius_metres + core_radius_metres))
            pressure_midmantle = (g_uppermantle * density_uppermantle * depth_midmantle_metres)/1000000000 # In GPa
            return pressure_midmantle, planet_mass

    def calculate_radius_and_mass(self, pressure_in_GPa, abundance_dict, lower_limit=False):
        if pressure_in_GPa <= 0:
            print('Warning! Pressure was <= 0, returning None, None')
            return None, None
        lower_radius_bound = 0 # Strictly speaking, this should be a bit bigger than 0 because Lena's model breaks down for very small objects
        upper_radius_bound = 12000  # This is roughly where Lena's model will start breaking down
        mass_abundance_dict = self.convert_number_abundances_to_mass(abundance_dict)
        try:
            iron_mass_fraction = mass_abundance_dict[ci.Element.Fe][Layer.bulk]
        except KeyError:
            # Then there's not much we can do
            print('Warning! Iron wt% was not supplied, returning None, None')
            return None, None
        if not (0 < iron_mass_fraction < 1):
            print('Warning! Iron wt% was non-physical, returning None, None')
            return None, None
        max_pressure, max_mass = self.calculate_pressure_and_mass(upper_radius_bound, iron_mass_fraction, lower_limit)
        if pressure_in_GPa > max_pressure:
            # Then we will never replicate the pressure and will get stuck in an infinite loop
            print('Warning! Pressure was too high to replicate. Considering raising the upper_radius_bound in GeologyModel.calculate_radius. Returning None, None.')
            return None, None
        tolerance = 0.0000000001
        iteration_count = 0
        max_iterations = 1000000
        while True:
            # Do a simple binary search
            test_radius = 0.5*(lower_radius_bound + upper_radius_bound)
            test_pressure, test_mass = self.calculate_pressure_and_mass(test_radius, iron_mass_fraction, lower_limit)
            rel_diff = abs((test_pressure - pressure_in_GPa)/pressure_in_GPa)
            if rel_diff < tolerance:
                return test_radius, test_mass
            else:
                # This logic relies on the fact that P increases with R
                if test_pressure > pressure_in_GPa:
                    # We guessed too high
                    upper_radius_bound = test_radius
                else:
                    # We guessed too low
                    lower_radius_bound = test_radius
            iteration_count += 1
            if iteration_count > max_iterations:
                # To prevent infinite loops
                print('Warning! Radius/Mass calculation exceeded iteration limit (' + str(max_iterations) + '), returning None, None')
                return None, None

    def get_earth_mix(self, fragment_core_fraction, abundance_dict=None, elements_to_mix=ci.usual_elements):
        # Function to calculate the composition of a body with an Earth-like mantle and Earth-like core, but arbitrary fragment core fraction
        # bulk = fcf*core_abundance + (1-fcf)*mantle_abundance
        # Abundances can also be overridden
        mix = dict()
        for element in elements_to_mix:  # For now we only care about these ones
            if abundance_dict is None:
                core_abundance = self.get_core_abundance(element)
                mantle_abundance = self.get_mantle_abundance(element)
            else:
                core_abundance = abundance_dict[element][Layer.core]
                mantle_abundance = abundance_dict[element][Layer.mantle]
            bulk_comp = (fragment_core_fraction*core_abundance) + ((1-fragment_core_fraction)*mantle_abundance)
            mix[element] = bulk_comp
        return mix

    def get_critical_fragment_core_fraction(self, D, parent_core_frac, dD_dV, dpcf_dV):
        # dD_dV is the derivative of the D w.r.t our variable of interest V (could be pressure, fO2...anything that affects the core/mantle compositions of the parent body of fixed bulk composition)
        # dpcf_dV is the derivative of the parent core fraction w.r.t V
        numerator = (parent_core_frac*dD_dV) + ((D-1)*dpcf_dV)
        denominator = dD_dV - ((D-1)*(D-1)*dpcf_dV)
        try:
            toret = numerator/denominator
        except ZeroDivisionError:
            print('Warning: hit a ZeroDivisionError. Returning None.')
            return None
        return toret

    def get_dX_dV(self, D, parent_core_frac, fragment_core_frac, dD_dV, dpcf_dV, parent_bulk_abundance=1):
        # This calculates the rate at which we change the final (normalised) number abundance of an element as we change some variable in the parent like pressure
        # It scales linearly with parent_bulk_abundance so if we don't care about the actual numbers we can just take that to be 1 and work in relative terms
        # dD_dV is the derivative of the D w.r.t our variable of interest V (could be pressure, fO2...anything that affects the core/mantle compositions of the parent body of fixed bulk composition)
        # dpcf_dV is the derivative of the parent core fraction w.r.t V
        numerator = (parent_bulk_abundance*(fragment_core_frac-parent_core_frac)*dD_dV) - (parent_bulk_abundance*dpcf_dV*(D-1)*(1+(fragment_core_frac*(D-1))))
        denominator = (1+(parent_core_frac*(D-1)))**2
        try:
            toret = numerator/denominator
        except ZeroDivisionError:
            print('Warning: hit a ZeroDivisionError. Returning None.')
            return None
        return toret

    def get_dlogX_dV(self, D, parent_core_frac, fragment_core_frac, dD_dV, dpcf_dV):
        # Like get_dX_dV, except for log (base 10) X instead of X since that is actually more interesting for observations (it also eliminates parent_bulk_abundance)
        numerator = ((fragment_core_frac-parent_core_frac)*dD_dV) - (((fragment_core_frac*(D-1))+1)*(D-1)*dpcf_dV)
        denominator = np.log(10)*((fragment_core_frac*(D-1))+1)*(1+((D-1)*parent_core_frac))
        try:
            toret = numerator/denominator
        except ZeroDivisionError:
            print('Warning: hit a ZeroDivisionError. Returning None.')
            return None
        return toret

    def convert_log_abundances_to_linear(self, abundance_dict):
        # Assuming input is a {element: log(X/Hx)} dict, as used in WD data
        unnormalised = dict()
        total = 0
        for el, log_abundance in abundance_dict.items():
            linear_abundance = 10**log_abundance
            unnormalised[el] = linear_abundance
            total += linear_abundance
        toret = dict()
        for el, unnormalised_abundance in unnormalised.items():
            toret[el] = dict()
            toret[el][Layer.bulk] = unnormalised_abundance/total
        return toret

    def convert_linear_abundances_to_log(self, abundance_dict, pollution_fraction=-6, layer=Layer.bulk):
        # Assuming input is a normalised {element: {layer: abundance}} dict, as used elsewhere in geology_info
        toret = dict()
        for el, sub_dict in abundance_dict.items():
            toret[el] = np.log10(sub_dict[layer]) + pollution_fraction
        return toret

