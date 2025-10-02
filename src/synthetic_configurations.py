#!/usr/bin/env python

from enum import Enum

import csv
import geology_info as gi
import model_parameters as mp
import numpy as np
import physical_constants as pc
import timescale_interpolator as ti

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
predefined_functions = dict()

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

def load_gspcwd_distributions(): # Vincent+ 2024 Table 3
    gspcwd_filename = 'GSPCWD_catalogue_101224.csv'
    max_teff = 40000
    min_teff = 3000
    max_logg = 9
    min_logg = 7
    dadb_teff_bins = {
        'DA': np.linspace(min_teff, max_teff, 101),
        'DB': np.linspace(min_teff, max_teff, 101)
    }
    dadb_logg_bins = {
        'DA': np.linspace(min_logg, max_logg, 101),
        'DB': np.linspace(min_logg, max_logg, 101)
    }
    teffs_da = list()
    teffs_db = list()
    loggs_da = list()
    loggs_db = list()
    teff_dist_name = 'GSPCWD_Teffs'
    logg_dist_name = 'GSPCWD_Loggs'
    with open('../data/' + gspcwd_filename, encoding='utf-8') as config_csv:
        next(config_csv, None)
        for row in csv.reader(config_csv, delimiter=','):
            spt = row[1]
            teff = float(row[18])
            logg = float(row[19])
            #comp = row[26]
            if min_teff <= teff <= max_teff and min_logg <= logg <= max_logg:
                if spt in ['DA', 'DA:']: # Luckily for us, the DAs and DBs fall neatly into H and He dominated here
                    teffs_da.append(teff)
                    loggs_da.append(logg)
                elif spt in ['DB', 'DB:']:
                    teffs_db.append(teff)
                    loggs_db.append(logg)
                elif spt in ['DZ', 'DZ:']:
                    p_da = float(row[12])
                    #p_db = float(row[13])
                    #p_dc = float(row[14])
                    #if p_da > 0.5 and p_db < 0.05 and p_dc < 0.05:
                    if p_da > 0.5:
                        # We get 142 of these...
                        teffs_da.append(teff)
                        loggs_da.append(logg)
                    else:
                        # ... and 1686 of these
                        teffs_db.append(teff)
                        loggs_db.append(logg)
                elif spt in ['DC:', 'DO', 'DC', 'DQ:', 'DQ', 'DO:']:
                    pass
                else:
                    pass
    #print(min(teffs_da)) #8358
    #print(max(teffs_da)) #556127
    #print(len(teffs_da)) #33226
    #print(min(teffs_db)) #10413
    #print(max(teffs_db)) #53975
    #print(len(teffs_db)) #5499
    #print(min(loggs_da)) #6.138
    #print(max(loggs_da)) #11.059
    #print(len(loggs_da)) #33226
    #print(min(loggs_db)) #5.481
    #print(max(loggs_db)) #9.702
    #print(len(loggs_db)) #5499
    teff_counts_da, teff_edges_da = np.histogram(teffs_da, bins=dadb_teff_bins['DA'])
    teff_tuples_da = [(teff_edges_da[i], teff_edges_da[i+1]) for i in range(0, len(teff_edges_da) - 1)]
    predefined_distributions[teff_dist_name + '_DA'] = (teff_tuples_da, teff_counts_da)
    logg_counts_da, logg_edges_da = np.histogram(loggs_da, bins=dadb_logg_bins['DA'])
    logg_tuples_da = [(logg_edges_da[i], logg_edges_da[i+1]) for i in range(0, len(logg_edges_da) - 1)]
    predefined_distributions[logg_dist_name + '_DA'] = (logg_tuples_da, logg_counts_da)

    teff_counts_db, teff_edges_db = np.histogram(teffs_db, bins=dadb_teff_bins['DB'])
    teff_tuples_db = [(teff_edges_db[i], teff_edges_db[i+1]) for i in range(0, len(teff_edges_db) - 1)]
    predefined_distributions[teff_dist_name + '_DB'] = (teff_tuples_db, teff_counts_db)
    logg_counts_db, logg_edges_db = np.histogram(loggs_db, bins=dadb_logg_bins['DB'])
    logg_tuples_db = [(logg_edges_db[i], logg_edges_db[i+1]) for i in range(0, len(logg_edges_db) - 1)]
    predefined_distributions[logg_dist_name + '_DB'] = (logg_tuples_db, logg_counts_db)
    #print(teff_edges_da[1] - teff_edges_da[0])
    #print(teff_edges_db[1] - teff_edges_db[0])
    #print(logg_edges_da[1] - logg_edges_da[0])
    #print(logg_edges_db[1] - logg_edges_db[0])
    #print(predefined_distributions['GSPCWD_Teffs_DA'])
    #print(predefined_distributions['GSPCWD_Teffs_DB'])
    #print(predefined_distributions['GSPCWD_Loggs_DA'])
    #print(predefined_distributions['GSPCWD_Loggs_DB'])

def load_hollands_distributions():
    teffs = list()
    loggs = list()
    teff_bins = np.linspace(4000, 9000, 11)
    logg_bins = np.linspace(7.1, 8.7, 9) # important to have a bin centred on 8: a lot of log(g)s are equal to 8 exactly
    file_name = 'WDInputData'
    #min_row = 0
    #max_row = 200
    hollands_ids = list(range(0, 201))
    hollands_ids.remove(85)
    hollands_ids.extend([249, 250])
    with open('../data/' + file_name + '.csv', encoding='utf-8') as config_csv:
        row_count = 0
        for row in csv.DictReader(config_csv):
            if row_count in hollands_ids:
                try:
                    teff = int(row['T_eff'])
                except ValueError:
                    teff = None
                if teff is not None:
                    teffs.append(teff)
                try:
                    logg = float(row['logg'])
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

def collisional_cascade_distribution(mass):
    alpha = 11/6
    # For a power law where dn = m^-alpha dm
    # (Don't need to normalise here, it will be normalised later)
    # The sanity check here is that for a uniform distribution we would have alpha = 0 - and indeed, setting alpha = 0 gives a uniform distribution
    return mass**(-alpha)

def log_collisional_cascade_distribution(log_mass):
    alpha = 11/6
    base = 10
    # For a power law where dn = m^-alpha dm, let l = log_b(m) and you get this:
    # (Don't need to normalise here, it will be normalised later) - so technically we might not actually need the base conversion prefactor
    # The sanity check here is that if you set alpha to zero, sample from this, and take 10**(answer), that should be uniformly distributed, and indeed it is
    mass = base**log_mass
    return np.log(base) * (mass**(1 - alpha))

def log_collisional_cascade_distribution_wyatt(log_mass):
    alpha = 1.57
    base = 10
    mass = base**log_mass
    return np.log(base) * (mass**(1 - alpha))

def load_distributions():
    print('Loading Distributions')
    #TODO: This should be called during main execution so that the file name can be set dynamically from the configuration file

    load_marc_distributions()
    load_tidal_distributions()
    load_amy_distributions()
    load_mwdd_distributions()
    load_hollands_distributions()
    #load_ns20_distributions()
    load_gspcwd_distributions()

def load_functions():
    predefined_functions['CollisionalCascade'] = collisional_cascade_distribution
    predefined_functions['LogCollisionalCascade'] = log_collisional_cascade_distribution
    predefined_functions['LogCollisionalCascadeWyatt'] = log_collisional_cascade_distribution_wyatt

load_distributions()
load_functions()

wd_configurations = {
    'TestWDConfigDB': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Delta, [5000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'TestWDConfigDB2': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.Delta, [5000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterNoOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'TestWDConfigDA': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.Delta, [10000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'TestWDConfigDA2': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.Delta, [10000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterNoOvershoot,
        mp.WDParameter.consider_thermohaline: False,
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
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'DAstar': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.Delta, [9000]),
        mp.WDParameter.logg: (Distribution.Delta, [8]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40])
    },
    'RealisticDAs': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['MWDD_DA_Teffs_40pc', predefined_distributions['MWDD_DA_Teffs_40pc'][0], predefined_distributions['MWDD_DA_Teffs_40pc'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['MWDD_DA_Loggs_40pc', predefined_distributions['MWDD_DA_Loggs_40pc'][0], predefined_distributions['MWDD_DA_Loggs_40pc'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.Bedard3DOvershootPatched,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'RealisticDBs': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['MWDD_DB_Teffs_40pc', predefined_distributions['MWDD_DB_Teffs_40pc'][0], predefined_distributions['MWDD_DB_Teffs_40pc'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['MWDD_DB_Loggs_40pc', predefined_distributions['MWDD_DB_Loggs_40pc'][0], predefined_distributions['MWDD_DB_Loggs_40pc'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.BedardVariableOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'HollandsDBs': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['HollandsTeffs', predefined_distributions['HollandsTeffs'][0], predefined_distributions['HollandsTeffs'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['HollandsLoggs', predefined_distributions['HollandsLoggs'][0], predefined_distributions['HollandsLoggs'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'RealisticDAsATT': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['MWDD_DA_Teffs_40pc', predefined_distributions['MWDD_DA_Teffs_40pc'][0], predefined_distributions['MWDD_DA_Teffs_40pc'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['MWDD_DA_Loggs_40pc', predefined_distributions['MWDD_DA_Loggs_40pc'][0], predefined_distributions['MWDD_DA_Loggs_40pc'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'HollandsDBsKN': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['HollandsTeffs', predefined_distributions['HollandsTeffs'][0], predefined_distributions['HollandsTeffs'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['HollandsLoggs', predefined_distributions['HollandsLoggs'][0], predefined_distributions['HollandsLoggs'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterNoOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'HollandsDBsBN': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['HollandsTeffs', predefined_distributions['HollandsTeffs'][0], predefined_distributions['HollandsTeffs'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['HollandsLoggs', predefined_distributions['HollandsLoggs'][0], predefined_distributions['HollandsLoggs'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.BedardNoOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'HollandsDBsKO': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['HollandsTeffs', predefined_distributions['HollandsTeffs'][0], predefined_distributions['HollandsTeffs'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['HollandsLoggs', predefined_distributions['HollandsLoggs'][0], predefined_distributions['HollandsLoggs'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'HollandsDBsBO': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['HollandsTeffs', predefined_distributions['HollandsTeffs'][0], predefined_distributions['HollandsTeffs'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['HollandsLoggs', predefined_distributions['HollandsLoggs'][0], predefined_distributions['HollandsLoggs'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.BedardOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    #'HWODAsKO': {
    #    mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
    #    mp.WDParameter.temperature: (Distribution.CustomDistribution, ['NS20_Teffs_DA', predefined_distributions['NS20_Teffs_DA'][0], predefined_distributions['NS20_Teffs_DA'][1]]),
    #    mp.WDParameter.logg: (Distribution.CustomDistribution, ['NS20_Loggs_DA', predefined_distributions['NS20_Loggs_DA'][0], predefined_distributions['NS20_Loggs_DA'][1]]),
    #    mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
    #    mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
    #    mp.WDParameter.consider_thermohaline: False,
    #    mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    #},
    #'HWODBsKO': {
    #    mp.WDParameter.spectral_type: (Distribution.Delta, ['DB']),
    #    mp.WDParameter.temperature: (Distribution.CustomDistribution, ['NS20_Teffs_DB', predefined_distributions['NS20_Teffs_DB'][0], predefined_distributions['NS20_Teffs_DB'][1]]),
    #    mp.WDParameter.logg: (Distribution.CustomDistribution, ['NS20_Loggs_DB', predefined_distributions['NS20_Loggs_DB'][0], predefined_distributions['NS20_Loggs_DB'][1]]),
    #    mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
    #    mp.WDParameter.timescale_type: ti.TimescaleType.KoesterOvershoot,
    #    mp.WDParameter.consider_thermohaline: False,
    #    mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    #},
    'DAsThermohalineOn': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['MWDD_DA_Teffs_40pc', predefined_distributions['MWDD_DA_Teffs_40pc'][0], predefined_distributions['MWDD_DA_Teffs_40pc'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['MWDD_DA_Loggs_40pc', predefined_distributions['MWDD_DA_Loggs_40pc'][0], predefined_distributions['MWDD_DA_Loggs_40pc'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.Bedard3DOvershoot,
        mp.WDParameter.consider_thermohaline: True,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'DAsThermohalineOff': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.CustomDistribution, ['MWDD_DA_Teffs_40pc', predefined_distributions['MWDD_DA_Teffs_40pc'][0], predefined_distributions['MWDD_DA_Teffs_40pc'][1]]),
        mp.WDParameter.logg: (Distribution.CustomDistribution, ['MWDD_DA_Loggs_40pc', predefined_distributions['MWDD_DA_Loggs_40pc'][0], predefined_distributions['MWDD_DA_Loggs_40pc'][1]]),
        #mp.WDParameter.temperature: (Distribution.CustomDistribution, ['GSPCWD_Teffs_DA', predefined_distributions['GSPCWD_Teffs_DA'][0], predefined_distributions['GSPCWD_Teffs_DA'][1]]),
        #mp.WDParameter.logg: (Distribution.CustomDistribution, ['GSPCWD_Loggs_DA', predefined_distributions['GSPCWD_Loggs_DA'][0], predefined_distributions['GSPCWD_Loggs_DA'][1]]),
        mp.WDParameter.mass: (Distribution.Normal, [0.6, 0.02]), # This doesn't actually matter UPDATE/WARNING: It actually does now! The new atmosphere model takes M_cvz as an input, which depends on M_WD via logq
        mp.WDParameter.timescale_type: ti.TimescaleType.Bedard3DOvershoot,
        mp.WDParameter.consider_thermohaline: False,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'TestDAsThermohalineOn': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.Delta, [17000]),
        mp.WDParameter.logg: (Distribution.Delta, [8.0]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.timescale_type: ti.TimescaleType.Bedard3DOvershootPatched,
        mp.WDParameter.consider_thermohaline: True,
        mp.WDParameter.distance: (Distribution.Delta, [40]) # This doesn't actually matter
    },
    'TestDAsThermohalineOff': {
        mp.WDParameter.spectral_type: (Distribution.Delta, ['DA']),
        mp.WDParameter.temperature: (Distribution.Delta, [17000]),
        mp.WDParameter.logg: (Distribution.Delta, [8.0]),
        mp.WDParameter.mass: (Distribution.Delta, [0.6]),
        mp.WDParameter.timescale_type: ti.TimescaleType.Bedard3DOvershootPatched,
        mp.WDParameter.consider_thermohaline: False,
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
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [20]), # kg, 10^20kg corresponds to an asteroid ish type mass (Vesta-esque, so a large asteroid)
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [21]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),   # years
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'TestPollutionConfig3': {
        mp.ModelParameter.metallicity: (Distribution.Delta, [400]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [1.2]),  # Myr
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [20.5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),   # years
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'TestPollutionConfig4': {
        mp.ModelParameter.metallicity: (Distribution.Delta, [400]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [1.2]),  # Myr
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [19.5]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),   # years
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'TestPollutionConfig5': {
        mp.ModelParameter.metallicity: (Distribution.Delta, [400]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [1.2]),  # Myr
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [18]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [1000000]),   # years
        mp.ModelParameter.pressure: (Distribution.Delta, [10]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-2])
    },
    'TestPollutionConfig6': {
        mp.ModelParameter.metallicity: (Distribution.Delta, [400]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [0.5]),  # Myr
        mp.ModelParameter.formation_distance: (Distribution.Delta, [0]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.03]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.1]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [None]),
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [18]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [20]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Normal, [20, 1]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [20]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [20]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Delta, [20]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Uniform, [10, 25]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Uniform, [10, 25]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Uniform, [10, 25]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Uniform, [10, 25]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Uniform, [10, 25]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Uniform, [10, 25]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Uniform, [10, 25]),
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
        mp.ModelParameter.fragment_mass: (Distribution.Uniform, [10, 25]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]), #yr
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'Wyatt2014': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 5]), #Myr - I reduced the rage but this is still such a long timescale that this is only applicable to DBs in its current form, really
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascadeWyatt', predefined_functions['LogCollisionalCascadeWyatt'], 10, 21.50515]), # Bottom end is still arbitrary - could constrain from their mu, sigma parameters perhaps?
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [20]), #yr
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'TestCollCasc': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]), #Myr
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15.7, np.log10(pc.M_Earth)]), # Maybe trymmax = 3.2 × 1024 g, μ = 8.0, σ = 1.3, q = 1.57
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]), #yr
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBDeltaCC': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]), #Myr
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 13.5, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]), #yr
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBCollisionalCC': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 13.5, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBTidalCC': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-0.55, -0.25]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['CMF_tidal_disruption', predefined_distributions['CMF_tidal_disruption'][0], predefined_distributions['CMF_tidal_disruption'][1]]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 13.5, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DADevoWet': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [5]),
        mp.ModelParameter.formation_distance: (Distribution.Slope, [-1, 0.1, 1]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [2]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        #mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DADevoDry': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [5]),
        mp.ModelParameter.formation_distance: (Distribution.Slope, [-1, 10, 1]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [2]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        #mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBDevoWet': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Slope, [-1, 0.1, 1]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [2]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        #mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBDevoDry': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Slope, [-1, 10, 1]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [2]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        #mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
   'DALightEl': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [5]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-1, 1]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [2]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        #mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBLightEl': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-1, 1]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [2]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        #mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.CustomDistribution, ['m_cf_035f6nogas', predefined_distributions['m_cf_035f6nogas'][0], predefined_distributions['m_cf_035f6nogas'][1]]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DARef': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Delta, [5]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-1, 1]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [2]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        #mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Delta, [10000000]),
        mp.ModelParameter.pressure: (Distribution.Delta, [45]),
        mp.ModelParameter.oxygen_fugacity: (Distribution.Delta, [-1.3])
    },
    'DBRef': {
        mp.ModelParameter.metallicity: (Distribution.Uniform, [0, 958]),
        mp.ModelParameter.t_sinceaccretion: (Distribution.Uniform, [0, 20]),
        mp.ModelParameter.formation_distance: (Distribution.Uniform, [-1, 1]),
        #mp.ModelParameter.formation_distance: (Distribution.Delta, [2]),
        mp.ModelParameter.feeding_zone_size: (Distribution.Uniform, [0, 0.15]),
        #mp.ModelParameter.feeding_zone_size: (Distribution.Delta, [0.05]),
        mp.ModelParameter.parent_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.parent_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        #mp.ModelParameter.fragment_core_frac: (Distribution.Delta, [0.17]),
        mp.ModelParameter.fragment_crust_frac: (Distribution.Delta, [0]),
        mp.ModelParameter.fragment_mass: (Distribution.CustomFunction, ['LogCollisionalCascade', predefined_functions['LogCollisionalCascade'], 15, np.log10(pc.M_Earth)]),
        mp.ModelParameter.accretion_timescale: (Distribution.Uniform, [0, 10000000]),
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
    plot_predefined_distributions(
        'NS20_Teffs_DA',
        ['NS20_Teffs_DA'],
        ['DAs in NS20'],
        'NS20DATeff',
        'dists',
        'Teff',
        222.22222222222263
    )
    plot_predefined_distributions(
        'NS20_Teffs_DB',
        ['NS20_Teffs_DB'],
        ['DBs in NS20'],
        'NS20DBTeff',
        'dists',
        'Teff',
        408.1632653061224
    )
    plot_predefined_distributions(
        'NS20_Loggs_DA',
        ['NS20_Loggs_DA'],
        ['DAs in NS20'],
        'NS20DALogg',
        'dists',
        'log(g)',
        0.020202020202019888
    )
    plot_predefined_distributions(
        'NS20_Loggs_DB',
        ['NS20_Loggs_DB'],
        ['DBs in NS20'],
        'NS20DBLogg',
        'dists',
        'log(g)',
        0.040816326530611846
    )
    plot_predefined_distributions(
        'GSPCWD_Teffs_DA',
        ['GSPCWD_Teffs_DA'],
        ['DAs in GSPCWD'],
        'GSPCWDDATeff',
        'dists',
        'Teff',
        370
    )
    plot_predefined_distributions(
        'GSPCWD_Teffs_DB',
        ['GSPCWD_Teffs_DB'],
        ['DBs in GSPCWD'],
        'GSPCWDDBTeff',
        'dists',
        'Teff',
        370
    )
    plot_predefined_distributions(
        'GSPCWD_Loggs_DA',
        ['GSPCWD_Loggs_DA'],
        ['DAs in GSPCWD'],
        'GSPCWDDALogg',
        'dists',
        'log(g)',
        0.02
    )
    plot_predefined_distributions(
        'GSPCWD_Loggs_DB',
        ['GSPCWD_Loggs_DB'],
        ['DBs in GSPCWD'],
        'GSPCWDDBLogg',
        'dists',
        'log(g)',
        0.02
    )

if __name__ == '__main__':
    main()
