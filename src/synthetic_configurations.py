#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum

import csv
import geology_info as gi
import model_parameters as mp
import numpy as np

class Distribution(Enum):
    Uniform = 0
    Normal = 1
    Delta = 2
    Triangle = 3
    Slope = 4
    CustomFunction = 5
    CustomDistribution = 6

    def __str__(self):
        return self.name

predefined_distributions = dict()

# These distributions are deprecated
def load_marc_distributions():
    geology_model = gi.GeologyModel()
    file_names = ['combined', 'Mdot_7', 'Mdot_8', 'Mdot_10']
    for file_name in file_names:
        data = np.load('../data/' + file_name + '.npz', allow_pickle=True)
        CMF = data['bins_CMF'] # core mass fraction
        N_CMF = data['N_CMF'] # weight of this core mass fraction
        bins = np.linspace(0, 0.5, 40)
        coarse_bins = np.linspace(0, 1, 40)  # For plotting purposes
        CNF = [geology_model.convert_core_mass_fraction_to_core_number_fraction(cmf) for cmf in CMF]
        marc_counts, edges = np.histogram(CNF, bins=bins, weights=N_CMF)
        marc_bins_tuples = [(edges[i], edges[i+1]) for i in range(0, len(edges) - 1)]
        predefined_distributions[file_name] = (marc_bins_tuples, marc_counts)
        coarse_counts, coarse_edges = np.histogram(CNF, bins=coarse_bins, weights=N_CMF)
        coarse_bins_tuples = [(coarse_edges[i], coarse_edges[i+1]) for i in range(0, len(coarse_edges) - 1)]
        predefined_distributions[file_name + '_coarse'] = (coarse_bins_tuples, coarse_counts)

def load_tidal_distributions():
    geology_model = gi.GeologyModel()
    file_names = ['CMF_tidal_disruption'] # This is a list of core mass fractions from Marc Brouwers - need to be converted to number fractions and binned
    # WARNING: The caching (synthetic_population.py, cache_inverse_cdf_from_table) seems to have a strange bug where it's possible to sample above the maximum input value
    # if you supply it with bins that go above that value (even if there's nothing in those bins). Workaround:
    # Make sure that your bins extend as far as the values and no further! Hence the precisely chosen upper limit of 0.321 here
    # This doesn't affect the collisional equivalent because that extends to 1
    bins = np.linspace(0, 0.321, 40)
    coarse_bins = np.linspace(0, 1, 40)  # For plotting purposes only (the bug above isn't relevant here)
    for file_name in file_names:
        cnfs = list()
        with open('../data/' + file_name + '.csv', encoding='utf-8') as cmf_csv:
            for row in csv.reader(cmf_csv):
                cnf = geology_model.convert_core_mass_fraction_to_core_number_fraction(float(row[0]))
                cnfs.append(cnf)
        tidal_counts, edges = np.histogram(cnfs, bins=bins)
        tidal_bins_tuples = [(edges[i], edges[i+1]) for i in range(0, len(edges) - 1)]
        predefined_distributions[file_name] = (tidal_bins_tuples, tidal_counts)
        coarse_counts, coarse_edges = np.histogram(cnfs, bins=coarse_bins)
        coarse_bins_tuples = [(coarse_edges[i], coarse_edges[i+1]) for i in range(0, len(coarse_edges) - 1)]
        predefined_distributions[file_name + '_coarse'] = (coarse_bins_tuples, coarse_counts)

def load_amy_distributions():
    geology_model = gi.GeologyModel()
    amy_filenames = ['m_cf_035f6nogas']
    mass_cutoff = 0.1 # Ignore everything heavier than 0.1 Earth mass (see Bonsor+ 2020)
    bins = np.linspace(0, 1, 80)
    coarse_bins = np.linspace(0, 1, 40)
    for file_name in amy_filenames:
        cnfs = list()
        with open('../data/' + file_name + '.dat', encoding='utf-8') as config_csv:
            for row in csv.reader(config_csv, delimiter=' '):
                mass = float(row[0])
                if mass < mass_cutoff:
                    cnf = geology_model.convert_core_mass_fraction_to_core_number_fraction(float(row[1]))
                    cnfs.append(cnf)
        amy_counts, edges = np.histogram(cnfs, bins=bins)
        amy_bins_tuples = [(edges[i], edges[i+1]) for i in range(0, len(edges) - 1)]
        predefined_distributions[file_name] = (amy_bins_tuples, amy_counts)
        coarse_counts, coarse_edges = np.histogram(cnfs, bins=coarse_bins)
        coarse_bins_tuples = [(coarse_edges[i], coarse_edges[i+1]) for i in range(0, len(coarse_edges) - 1)]
        predefined_distributions[file_name + '_coarse'] = (coarse_bins_tuples, coarse_counts)

def load_mwdd_distributions():
    DA_DB_40pc_base_filename = 'MWDD-export-40pc'

    dadb_teff_bins = {
        'DA': np.linspace(4000, 40000, 21),
        'DB': np.linspace(6000, 30000, 6)
    }
    dadb_logg_bins = {
        'DA': np.linspace(7, 9.4, 16),
        'DB': np.linspace(7.85, 9, 6)
    }

    for stellar_suffix in ['DA', 'DB']:
        teffs = list()
        loggs = list()
        teff_dist_name = 'MWDD_' + stellar_suffix + '_Teffs_40pc'
        logg_dist_name = 'MWDD_' + stellar_suffix + '_Loggs_40pc'
        with open('../data/' + DA_DB_40pc_base_filename + stellar_suffix + 's.csv', encoding='utf-8') as config_csv:
            for row in csv.reader(config_csv, delimiter=','):
                try:
                    teff = float(row[5])
                except ValueError:
                    teff = None
                try:
                    logg = float(row[6])
                except ValueError:
                    logg = None
                if teff is not None:
                    teffs.append(teff)
                if logg is not None:
                    loggs.append(logg)
        teff_counts, teff_edges = np.histogram(teffs, bins=dadb_teff_bins[stellar_suffix])
        teff_tuples = [(teff_edges[i], teff_edges[i+1]) for i in range(0, len(teff_edges) - 1)]
        predefined_distributions[teff_dist_name] = (teff_tuples, teff_counts)
        logg_counts, logg_edges = np.histogram(loggs, bins=dadb_logg_bins[stellar_suffix])
        logg_tuples = [(logg_edges[i], logg_edges[i+1]) for i in range(0, len(logg_edges) - 1)]
        predefined_distributions[logg_dist_name] = (logg_tuples, logg_counts)

def load_hollands_distributions():
    teffs = list()
    loggs = list()
    teff_bins = np.linspace(4000, 9000, 11)
    logg_bins = np.linspace(7.1, 8.7, 9) # important to have a bin centred on 8: a lot of log(g)s are equal to 8 exactly
    file_name = 'WDInputData'
    min_row = 1
    max_row = 201
    with open('../data/' + file_name + '.csv', encoding='utf-8') as config_csv:
        row_count = 0
        for row in csv.reader(config_csv, delimiter=','):
            if min_row <= row_count <= max_row:
                try:
                    teff = int(row[4])
                except ValueError:
                    teff = None
                if teff is not None:
                    teffs.append(teff)
                try:
                    logg = float(row[6])
                except ValueError:
                    logg = None
                if logg is not None:
                    loggs.append(logg)
            row_count += 1
        teff_counts, teff_edges = np.histogram(teffs, bins=teff_bins)
        logg_counts, logg_edges = np.histogram(loggs, bins=logg_bins)
        teff_tuples = [(teff_edges[i], teff_edges[i+1]) for i in range(0, len(teff_edges) - 1)]
        logg_tuples = [(logg_edges[i], logg_edges[i+1]) for i in range(0, len(logg_edges) - 1)]
    predefined_distributions['HollandsTeffs'] = (teff_tuples, teff_counts)
    predefined_distributions['HollandsLoggs'] = (logg_tuples, logg_counts)

def load_distributions():
    print('Loading Distributions')
    load_marc_distributions()
    load_tidal_distributions()
    load_amy_distributions()
    load_mwdd_distributions()
    load_hollands_distributions()

load_distributions()

wd_configurations = {
    'TestWDConfig': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Delta, [5000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'DummyWDConfig': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Uniform, [3000, 9000]),
        mp.WDParameter.logg: (Distribution.Normal, [8, 0.2]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.05]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'TestA': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Slope, [5000, 5, 50000]),
        mp.WDParameter.logg: (Distribution.Normal, [8, 0.2]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'TestB': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA', 'DB']),
        mp.WDParameter.temperature: (Distribution.Slope, [5000, 0.2, 60000]),
        mp.WDParameter.logg: (Distribution.Normal, [8, 0.2]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'TestC': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Slope, [25000, 0.1, 40000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'TestD': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Slope, [0, 4, 40000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'Homogeneous': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Delta, [9000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'UniformTeff': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Uniform, [5000, 15000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'DAstar': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.Delta, [9000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'RealisticDAs': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['MWDD_DA_Teffs_40pc', predefined_distributions['MWDD_DA_Teffs_40pc'][0], predefined_distributions['MWDD_DA_Teffs_40pc'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['MWDD_DA_Loggs_40pc', predefined_distributions['MWDD_DA_Loggs_40pc'][0], predefined_distributions['MWDD_DA_Loggs_40pc'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'RealisticDBs': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['MWDD_DB_Teffs_40pc', predefined_distributions['MWDD_DB_Teffs_40pc'][0], predefined_distributions['MWDD_DB_Teffs_40pc'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['MWDD_DB_Loggs_40pc', predefined_distributions['MWDD_DB_Loggs_40pc'][0], predefined_distributions['MWDD_DB_Loggs_40pc'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'HollandsDBs': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['HollandsTeffs', predefined_distributions['HollandsTeffs'][0], predefined_distributions['HollandsTeffs'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['HollandsLoggs', predefined_distributions['HollandsLoggs'][0], predefined_distributions['HollandsLoggs'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    }
}

pollution_configurations = {
    'TestPollutionConfig': {
        mp.ModelParameter.metallicity: (Distribution.Delta, [400]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [1.2]),  # Myr
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),   # years
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'DummyPollutionConfig': {
        mp.ModelParameter.metallicity: (Distribution.Triangle, [50, 450, 950]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [2, 2.0005]),
        #mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [2.0005]),
        mp.ModelParameter.formation_distance: (Distribution.Slope, [-0.3, 0.3, 0.8]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Uniform, [0.01, 0.99]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 2.5]),
        #mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [999999, 1000001]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [2000000]),
        mp.ModelParameter.pressure: (Distribution.Slope, [0, 5, 60]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Triangle, [-3, -2, -1])
    },
    'TestPollutionConfig2': {
        mp.ModelParameter.metallicity: (Distribution.Delta, [400]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [1.2]),  # Myr
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.pollution_frac: (Distribution.Delta, [-6]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),   # years
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'Test1': {
        mp.ModelParameter.metallicity: (Distribution.Triangle, [100, 400, 600]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [1.2]),  # Myr
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.01, 0.99]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),   # years
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'Test2': {
        mp.ModelParameter.metallicity: (Distribution.Triangle, [100, 400, 600]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [1.0005]),
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Slope, [0, 2, 0.017]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'Test3': {
        mp.ModelParameter.metallicity: (Distribution.Triangle, [100, 400, 600]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [1.0005]),
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Normal, [0.5, 0.1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'UniformFcf': {
        #mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [10]),
        mp.ModelParameter.formation_distance: (Distribution.Delta, [-0.5]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Uniform, [0, 1]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.5]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        #mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [100000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'UniformDistance': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        #mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [0]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.6, -0.25]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [0.1]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [100000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'UniformTime': {
        #mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0.5]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [100000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'UniformFDT': {
        #mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
        #mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [15]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.6, -0.25]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [-0.3]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Uniform, [0, 1]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        #mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [100000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'GaussianFUniformDT': {
        #mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.6, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Normal, [0.17, 0.05]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.5]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        #mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [100000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DeltaFUniformDT': {
        #mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.6, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.5]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-7, 0.5]),
        #mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [100000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'Amy035f6nogas': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.5]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Normal, [-6, 1]), # This distribution is roughly based on output from amy's sample from the Al-26 paper
        #mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DApollution': {
        #mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
        #mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [15]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
        mp.ModelParameter.formation_distance: (Distribution.Delta, [-0.3]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [-0.3]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        #mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DAUniformMDeltaFDT': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        #mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [15]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
        mp.ModelParameter.formation_distance: (Distribution.Delta, [-0.3]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [-0.3]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        #mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DAUniformMFDDeltaT': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        #mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [15]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.6, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Uniform, [0, 1]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        #mp.ModelParameter.pollution_frac: (Distribution.Delta, [-7]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DAControl': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Uniform, [0, 1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Triangle, [-12, -12, -3]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DACollisional': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Triangle, [-12, -12, -3]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBCollisional': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Triangle, [-12, -12, -5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DATidal': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['CMF_tidal_disruption', predefined_distributions['CMF_tidal_disruption'][0], predefined_distributions['CMF_tidal_disruption'][1]]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Triangle, [-12, -12, -3]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBTidal': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['CMF_tidal_disruption', predefined_distributions['CMF_tidal_disruption'][0], predefined_distributions['CMF_tidal_disruption'][1]]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Triangle, [-12, -12, -5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBControl': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]), #Myr
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Uniform, [0, 1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Triangle, [-12, -12, -3]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]), #yr
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    #'DAForSynthBayesCompMantle': {
    #    mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
    #    mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
    #    mp.ModelParameter.formation_distance: (Distribution.Delta, [-0.5]),
    #    mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
    #    mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
    #    mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
    #    mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.05]),
    #    mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
    #    mp.ModelParameter.pollution_frac: (Distribution.Delta, [-6]),
    #    mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
    #    mp.ModelParameter.pressure: (Distribution.Delta, [45]),
    #    mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    #},
    #'DAForSynthBayesCompCore': {
    #    mp.ModelParameter.metallicity: (Distribution.Delta, [478]),
    #    mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
    #    mp.ModelParameter.formation_distance: (Distribution.Delta, [-0.5]),
    #    mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
    #    mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
    #    mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
    #    mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.75]),
    #    mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
    #    mp.ModelParameter.pollution_frac: (Distribution.Delta, [-6]),
    #    mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
    #    mp.ModelParameter.pressure: (Distribution.Delta, [45]),
    #    mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    #},
    'DADeltaFcf': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [4]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Triangle, [-12, -12, -3]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [5000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBDeltaFcf': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]), #Myr
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.pollution_frac: (Distribution.Triangle, [-12, -12, -5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]), #yr
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    }
}

def plot_predefined_distributions(reference_dist_for_bins, dist_names, dist_label_names, file_prefix, file_suffix, x_label, bin_width, additional_line_dict=None):
    import graph_factory as gf
    graph_fac = gf.GraphFactory()
    xbar = [(edges[0]+edges[1])/2 for edges in predefined_distributions[reference_dist_for_bins][0]]
    all_heights = list()
    for dn in dist_names:
        norm = sum(predefined_distributions[dn][1])
        #if dn == 'm_cf_035f6nogas':
        #    all_heights.append(predefined_distributions[dn][1][0:39]/norm)
        #else:
        all_heights.append(predefined_distributions[dn][1]/norm)

    graph_fac.make_histogram(
        xbar,
        all_heights,
        dist_label_names,
        file_prefix,
        bin_width,
        1,
        x_label,
        file_suffix,
        dict(),
        {'Delta': {'x_start': 0.17}}
    )

def main():
    plot_predefined_distributions(
        'CMF_tidal_disruption_coarse',
        ['CMF_tidal_disruption_coarse', 'm_cf_035f6nogas_coarse'],
        ['Orbit-by-orbit', 'Collisional'],
        'tidal',
        'dists',
        'Fragment Core Number Fraction',
        1/19,
        {'Delta': {'x_start': 0.17}}
    )
    plot_predefined_distributions(
        'MWDD_DA_Teffs_40pc',
        ['MWDD_DA_Teffs_40pc'],
        ['DAs in 40pc sample'],
        'DATeff',
        'dists',
        'Teff /K',
        1800
    )
    plot_predefined_distributions(
        'MWDD_DA_Loggs_40pc',
        ['MWDD_DA_Loggs_40pc'],
        ['DAs in 40pc sample'],
        'DALogg',
        'dists',
        'log(g)',
        0.16
    )
    plot_predefined_distributions(
        'HollandsTeffs',
        ['HollandsTeffs'],
        ['Cool DZs'],
        'DBTeff',
        'dists',
        'Teff /K',
        500
    )
    plot_predefined_distributions(
        'HollandsLoggs',
        ['HollandsLoggs'],
        ['Cool DZs'],
        'DBLogg',
        'dists',
        'log(g)',
        0.2
    )

if __name__ == '__main__':
    main()

