#!/usr/bin/env python
# -*- coding: utf-8 -*-

# A script to produce some new plots/results for my thesis, based
# on having a second go at the analysis of the sample of oxygen-bearing WDs,
# Based on Brouwers+ 2023

from argparse import Namespace
import csv
import graph_factory as gf
import manager as mn
import pwd_utils as pu

def get_sample_dict():
    # This should be a superset of Brouwers+ 2022b and Trierweiler+ 2023
    # From Trierweiler, we exclude WD1226+110 (alias of SDSSJ1228+1040) WD1248+1004 (alias of 1248+1005) Ton345 (alias of 0845+2257). Should note these aliases in the write-up! They are actually not excluded at all!

    # The IDs of the systems which have O, Mg, Fe, Si detections or upper bounds (i.e. the main rock-forming elements)
    # Exclude obsolete datasets
    # Tie breaker for dupes: O detections rather than upper bounds, then lowest O/Si
    #EC22211-2525 = 530 NEW
    #G241-6 = 359 OLD
    #G29-38 = 366 OLD
    #GaiaJ0218+3625 = 524 NEW
    #GaiaJ0510+2315Spec = 465 (REMOVE: No Fe in later version) OBSOLETE
    #GaiaJ0644-0352 = 302 REMOVE: obsolete
    #GaiaJ0644-0352Photo = 486 REMOVE: obsolete
    #GaiaJ0644-0352Spec = 457 REMOVE: obsolete
    #GaiaJ0644-op-P = 553 NEW
    #GaiaJ0644-op-S = 552 REMOVE: O/Si higher than for op-P OBSOLETE
    #GaiaJ1922+4709 = 529 NEW
    #GALEX1931+0117 = 352 NEW DATASET, OLD SYSTEM
    #GALEX1931+0117 = 353 REMOVE: O/Si higher than for 352 NB: Marc used this set of data instead of 352!!! He cited the wrong source (i.e., he cited Vennes but used the Melis data) OBSOLETE
    #GALEXJ2339 = 260 OLD
    #GD133 = 522 NEW
    #GD378 = 261 OLD
    #GD40 = 358 OLD
    #GD424 = 253 OLD
    #GD61 = 360 OLD
    #HS2253+8023 = 362 OLD
    #PG0843+516 = 349 OLD
    #PG1015+161 = 356 OLD
    #SDSSJ0738+1835 = 363 OLD
    #SDSSJ0845+2257 = 368 OLD
    #SDSSJ0956+5912 = 462 (note: not GTC, which would have lower O) (should probably use GTC then!) REMOVE OBSOLETE
    #SDSSJ0956+5912 = 479 (GTC) OLD
    #SDSSJ1043+0855 = 370 OLD
    #SDSSJ1228+1040 = 365 This is the Opt one OLD
    #SDSSJ1228+1040 = 366 REMOVE: higher O/Si
    #SDSSJ1242+5226 = 367 OLD
    #SDSSJ1248+1005 = 526 NEW
    #SDSSJ1734+6052 = 528 NEW
    #SDSSJ2047-1259 = 354 OLD
    #SDSSJ2248+2632 = 531 NEW
    #WD0446-255 = 347 NEW
    #WD1145+017 = 375 (higher O/Si) REMOVE
    #WD1145+017 = 472 (Budaj, would be the more conservative choice!) NEW DATASET, OLD SYSTEM
    #WD1232+563 = 373 OLD
    #WD1244+498 = 525 NEW
    #WD1350-162 = 377 OLD
    #WD1415+234 = 527 NEW
    #WD1425+540 = 371 OLD
    #WD1536+520 = 369 OLD
    #WD1551+175 = 355 OLD
    #WD1622-uv-P = 547 NEW
    #WD1622-uv-S = 545 REMOVE: higher O/Si
    #WD2207+121 = 374 OLD

    #Also adding ones that only have an O upper bound:
    #PG1225-079 = 380 OLD
    #WD2115-560 = 372 OLD
    #WD0449-259 = 376 (No Na) NEW
    #WD2157-574 = 240 OLD
    #WD0122-227 = 345 NEW
    #WD2216-657 = 378 NEW
    #WD2230-125 = 346 NEW
    #GD362 = 357 OLD

    #Now adding ones from the newest version of Laura's data. (Which also removes GaiaJ0611-op-S = 536, and WD0611-6931 by proxy)
    #GaiaJ0006-op-P EXCLUDE (has higher O/Si)
    #GaiaJ0006-op-S = 532
    #SDSSJ0006+2858 EXCLUDE (obsolete)
    #SDSSJ0006+2858Spec EXCLUDE (obsolete)
    #GaiaJ0510-op-P = 542
    #GaiaJ0510-op-S EXCLUDE (has higher O/Si)
    #GaiaJ2100-op-P EXCLUDE (has higher O/Si)
    #GaiaJ2100-op-S = 554
    #GaiaJ2100+2122  EXCLUDE (obsolete)
    #GaiaJ2100+2122Photo  EXCLUDE (obsolete)
    #GaiaJ2100+2122Spec  EXCLUDE (obsolete)
    #WD1622-op-P EXCLUDE (no O detection)
    #WD1622-op-S EXCLUDE (no O detection)
    #WD1622+587Spec EXCLUDE (obsolete)

    # Note to self, need to graft in 532 542 and 554 manually

    sample_dict = {
        'EC22211-2525': {'id': 530},
        'G241-6': {'id': 359, 'csv_suffix': 'Corr'},
        'G29-38': {'id': 366, 'csv_suffix': 'Corr'},
        'GaiaJ0218+3625': {'id': 524},
        'GaiaJ0644-op-P': {'id': 553},
        'GaiaJ1922+4709': {'id': 529},
        'GALEX1931+0117': {'id': 352, 'csv_suffix': 'VCorr'},
        'GALEXJ2339': {'id': 260},
        'GD133': {'id': 522},
        'GD378': {'id': 261},
        'GD40': {'id': 358, 'csv_suffix': 'Corr'},
        'GD424': {'id': 253},
        'GD61': {'id': 360, 'csv_suffix': 'Corr'},
        'HS2253+8023': {'id': 362, 'csv_suffix': 'Corr'},
        'PG0843+516': {'id': 349, 'csv_suffix': 'GCorr'},
        'PG1015+161': {'id': 356, 'csv_suffix': 'Corr'},
        'SDSSJ0738+1835': {'id': 363, 'csv_suffix': 'Corr'},
        'SDSSJ0845+2257': {'id': 368, 'csv_suffix': 'Corr'},
        'SDSSJ0956+5912': {'id': 479, 'csv_suffix': 'H21GTC'},
        'SDSSJ1043+0855': {'id': 370, 'csv_suffix': 'Corr'},
        'SDSSJ1228+1040': {'id': 365, 'csv_suffix': 'Corr'},
        'SDSSJ1242+5226': {'id': 367, 'csv_suffix': 'Corr'},
        'SDSSJ1248+1005': {'id': 526},
        'SDSSJ1734+6052': {'id': 528},
        'SDSSJ2047-1259': {'id': 354, 'csv_suffix': 'Corr'},
        'SDSSJ2248+2632': {'id': 531},
        'WD0446-255': {'id': 347},
        'WD1145+017': {'id': 472, 'csv_suffix': 'Budaj'},
        'WD1232+563': {'id': 373, 'csv_suffix': 'Corr'},
        'WD1244+498': {'id': 525},
        'WD1350-162': {'id': 377, 'csv_suffix': 'NoNaCorr'},
        'WD1415+234': {'id': 527},
        'WD1425+540': {'id': 371, 'csv_suffix': 'Corr'},
        'WD1536+520': {'id': 369, 'csv_suffix': 'Corr'},
        'WD1551+175': {'id': 355, 'csv_suffix': 'Corr'},
        'WD1622-uv-P': {'id': 547},
        'WD2207+121': {'id': 374, 'csv_suffix': 'Corr'},
        'PG1225-079': {'id': 380, 'csv_suffix': 'Corr'},
        'WD2115-560': {'id': 372, 'csv_suffix': 'Corr'},
        'WD0449-259': {'id': 376, 'csv_suffix': 'NoNaCorr'},
        'WD2157-574': {'id': 240},
        'WD0122-227': {'id': 345, 'csv_suffix': 'SiO'},
        'WD2216-657': {'id': 378, 'csv_suffix': 'SiCorr'},
        'WD2230-125': {'id': 346},
        'GD362': {'id': 357, 'csv_suffix': 'Corr'},
        'GaiaJ0006-op-S': {'id': 532},
        'GaiaJ0510-op-P': {'id': 542},
        'GaiaJ2100-op-S': {'id': 554}
    }
    return sample_dict

def get_sample_names(id_list):
    manager = mn.Manager(  #This is just to load the stellar compositions
        Namespace(
            wd_data_filename='WDInputData.csv',
            stellar_compositions_filename='StellarCompositionsSortFE.csv',
            n_live_points = 0,
            pollution_model_names=['Model_24'],
            enhancement_model='Earthlike',
            base_dir=pu.get_path_to_data()
        )
    )
    toret = list()
    i = 0
    for name in manager.wd_names:
        if i in id_list:
            toret.append(name)
        i += 1
    return toret

def get_output_csvs(sample_dict):
    print('Warning! This script has not yet been updated to use the new output directory structure/file names!')
    print('If using on newly created outputs, you should update get_output_csvs to give the correct names (sorry)')
    for name in sample_dict:
        sample_dict[name]['csv'] = 'stats_lp2000_obs' + str(sample_dict[name]['id']) + '_' + name + sample_dict[name].get('csv_suffix', '') + '_Hierarchy_Default_NEL.csv'
    return sample_dict

def compile_oxygen_stats(sample_dict):
    # We should also read:
    #Oxygen error!
    oxygen_strat = 'default'
    for name in sample_dict:
        median_fractional_excess_oxygen = None
        median_fractional_excess_oxygen_upper_error = None
        median_fractional_excess_oxygen_lower_error = None
        oxygen_sigma_significance = None
        p_excess = None
        p_deficit = None
        buildup_pc = None
        ss_pc = None
        dec_pc = None
        oxygen_error = None
        csv_name = pu.get_path_to_output_statfiles_dir() + sample_dict[name]['csv']
        with open(csv_name, encoding='utf-8') as output_csv:
            print()
            print('Reading ' + csv_name)
            best_heated_model_name = None
            good_fit = None
            reading_from_correct_section = False
            for row in csv.reader(output_csv):
                if len(row) > 0:
                    if row[0] == 'Best heated model name:':
                        best_heated_model_name = row[1]
                        print('Best model with heating was ' + best_heated_model_name)
            output_csv.seek(0)
            for row in csv.reader(output_csv):
                if len(row) > 0:
                    if row[0] == best_heated_model_name:
                        good_fit = row[12] == 'True'
                        if not good_fit:
                            sample_dict[name]['Good Fit'] = False
                            print('WARNING! Fit was not good')
                        else:
                            sample_dict[name]['Good Fit'] = True
                            print('(Good fit!)')
                    if row[0] == 'Results from model:':
                        if row[1] == best_heated_model_name:
                            reading_from_correct_section = True
                        else:
                            reading_from_correct_section = False
                    if row[0] == 'Median excess (default)' and reading_from_correct_section:
                        median_fractional_excess_oxygen = float(row[1])
                    if row[0] == 'Upper error excess (default)' and reading_from_correct_section:
                        median_fractional_excess_oxygen_upper_error = float(row[1])
                    if row[0] == 'Lower error excess (default)' and reading_from_correct_section:
                        median_fractional_excess_oxygen_lower_error = float(row[1])
                    if row[0] == 'P(excess) (default)' and reading_from_correct_section:
                        p_excess = float(row[1])
                    if row[0] == 'P(deficit) (default)' and reading_from_correct_section:
                        p_deficit = float(row[1])
                    if row[0] == 'Sigma excess (default)' and reading_from_correct_section:
                        oxygen_sigma_significance = float(row[1])
                    if row[0] == 'Build Up % (sampled):' and reading_from_correct_section:
                        buildup_pc = float(row[1])
                    if row[0] == 'Steady State % (sampled):' and reading_from_correct_section:
                        ss_pc = float(row[1])
                    if row[0] == 'Declining % (sampled):' and reading_from_correct_section:
                        dec_pc = float(row[1])
                    if row[0] == 'Input Errors:' and reading_from_correct_section:
                        oxygen_error = float(row[10])
        sample_dict[name]['Median Fractional Excess Oxygen'] = median_fractional_excess_oxygen
        sample_dict[name]['Median Fractional Excess Oxygen Upper Error'] = median_fractional_excess_oxygen_upper_error
        sample_dict[name]['Median Fractional Excess Oxygen Lower Error'] = median_fractional_excess_oxygen_lower_error
        sample_dict[name]['Oxygen Sigma Significance'] = oxygen_sigma_significance
        sample_dict[name]['P(excess)'] = p_excess
        sample_dict[name]['P(deficit)'] = p_deficit
        sample_dict[name]['Build Up %'] = buildup_pc
        sample_dict[name]['Steady State %'] = ss_pc
        sample_dict[name]['Declining %'] = dec_pc
        sample_dict[name]['Oxygen Error'] = oxygen_error
        sample_dict[name]['ScaledTimePlot'] = name + sample_dict[name].get('csv_suffix', '') + '_' + best_heated_model_name + '_lp2000_obs' + str(sample_dict[name]['id']) + '_NEL_D__scaled_times_thesis.pdf.txt'
        sample_dict[name]['SemiSampledEOPlot'] = name + sample_dict[name].get('csv_suffix', '')  + '_' + best_heated_model_name + '_lp2000_obs' + str(sample_dict[name]['id']) + '_NEL_D__semisampled_oxygen_excess_thesis.pdf.txt'
    return sample_dict

def write_sample_dict(sample_dict):
    with open('oxygen_excess_stats.csv', 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            'System',
            'Oxygen Error',
            'Good Fit',
            'Build Up %',
            'Steady State %',
            'Declining %',
            'Median Fractional Excess Oxygen',
            'Median Fractional Excess Oxygen Upper Error',
            'Median Fractional Excess Oxygen Lower Error',
            'Oxygen Sigma Significance',
            'P(excess)',
            'P(deficit)'
        ])
        for sys, vals in sample_dict.items():
            to_write.writerow([
                sys,
                vals['Oxygen Error'],
                'TRUE' if vals['Good Fit'] else 'FALSE',
                vals['Build Up %'],
                vals['Steady State %'],
                vals['Declining %'],
                vals['Median Fractional Excess Oxygen'],
                vals['Median Fractional Excess Oxygen Upper Error'],
                vals['Median Fractional Excess Oxygen Lower Error'],
                vals['Oxygen Sigma Significance'],
                vals['P(excess)'],
                vals['P(deficit)']
            ])

def multipanelise_from_dumps(sample_dict, plot_key):
    y_dimension = 24
    x_dimension = 2
    file_name_root = 'oxygen_multipanel_' + plot_key
    extensions = ['pdf', 'png']
    filenames = [file_name_root + '.' + extension for extension in extensions]
    fig_height = 30
    fig_width = 15
    gridspec_wspace = 0.1
    gridspec_hspace = 0
    sharey_axes = False
    sharex_axes = True
    list_of_plot_dumps = list()
    for sys, vals in sorted(sample_dict.items()):
        list_of_plot_dumps.append(pu.get_path_to_output_graphs_dir() + vals[plot_key])
    print(list_of_plot_dumps)
    graph_fac = gf.GraphFactory()
    graph_fac.multipanelise_from_dumps(
        list_of_plot_dumps,
        y_dimension,
        x_dimension,
        filenames,
        fig_height,
        fig_width,
        gridspec_wspace,
        gridspec_hspace,
        sharey_axes,
        sharex_axes
    )

def main():
    sample_dict = get_sample_dict()
    sample_dict = get_output_csvs(sample_dict)
    sample_dict = compile_oxygen_stats(sample_dict)
    #write_sample_dict(sample_dict)
    multipanelise_from_dumps(sample_dict, 'ScaledTimePlot')
    multipanelise_from_dumps(sample_dict, 'SemiSampledEOPlot')
    print(sample_dict['WD1425+540'])
    print(len(sample_dict))

if __name__ == '__main__':
    main()
