#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np
import os
import xlrd

# This script is designed to be used on xlsx outputs from the old version of the code (in original_codebase) - it's essentially deprecated now and included just for completeness

old_model_params_dict = {
    'M1': ["Stellar metallicity indices", "log(Pollution Fraction)"],
    'M2': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M3': ["Stellar metallicity indices", "log(Formation Distance/AU)", "log(Pollution Fraction)"],
    'M4': ["Stellar metallicity indices", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "log(Pollution Fraction)"],
    'M5': ["Stellar metallicity indices", "log(Formation Distance/AU)", "Fragment Core Fraction", "log(Pollution Fraction)"],
    'M6': ["Stellar metallicity indices", "log(Formation Distance/AU)", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)"],
    'M7': ["Stellar metallicity indices", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Fragment Core Fraction", "log(Pollution Fraction)"],
    'M8': ["Stellar metallicity indices", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)"],
    #'M9': ?,
    'M10': ["Stellar metallicity indices", "log(Formation Distance/AU)", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)"],
    'M11': ["Stellar metallicity indices", "log(Formation Distance/AU)", "Parent Core Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)"],
    'M12': ["Stellar metallicity indices", "log(Formation Distance/AU)", "Parent Core Fraction", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)"],
    #'M13': ?,
    #'M14': ?,
    #'M15': ?,
    #'M16': ?,
    'M17': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M18': ["Stellar metallicity indices","Time since Accretion/Myrs", "log(Formation Distance/AU)", "Parent Core Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M19': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M20': ["Stellar metallicity indices","Time since Accretion/Myrs", "log(Formation Distance/AU)", "Parent Core Fraction", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M21': ["Stellar metallicity indices","Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M22': ["Stellar metallicity indices","Time since Accretion/Myrs", "log(Formation Distance/AU)","Feeding Zone Size/AU", "Parent Core Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M23': ["Stellar metallicity indices","Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M24': ["Stellar metallicity indices","Time since Accretion/Myrs", "log(Formation Distance/AU)","Feeding Zone Size/AU", "Parent Core Fraction", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M25': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M26': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)","Feeding Zone Size/AU", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M27': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Fragment Core Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"],
    'M28': ["Stellar metallicity indices","Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Fragment Core Fraction", "log(Pollution Fraction)","log(Accretion Event Timescale/Yrs)"]
}

model_params_dict = {
    'M1 = (Stellar Index, Time since Accretion, Accretion Event Lifetime, Pollution Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M2 = (M1 + Formation Distance)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M3 = (M1 + Formation Distance + Feeding Zone)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M4 = (M1 + Fragment Core Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "Fragment Core Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M5 = (M1 + Fragment Core Fraction + Fragment Crust Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M6 = (M2 + Fragment Core Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Fragment Core Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M7 = (M2 + Fragment Core Fraction + Fragment Crust Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M8 = (M3 + Fragment Core Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Fragment Core Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M9 = (M3 + Fragment Core Fraction + Fragment Crust Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M10 = (M7 + Parent Core Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Parent Core Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M11 = (M7 + Parent Crust Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M12 = (M7 + Parent Core Fraction + Parent Crust Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Parent Core Fraction", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M13 = (M9 + Parent Core Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Parent Core Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M14 = (M9 + Parent Crust Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"],
    'M15 = (M9 + Parent Core Fraction + Parent Crust Fraction)': ["Stellar metallicity indices", "Time since Accretion/Myrs", "log(Formation Distance/AU)", "Feeding Zone Size/AU", "Parent Core Fraction", "Parent Crust Fraction", "Fragment Core Fraction", "Fragment Crust Fraction", "log(Pollution Fraction)", "log(Accretion Event Timescale/Yrs)"]
}

elements_present_dict = {
    'SDSSJ0002+3209': 3,
    'SDSSJ0004+0819': 3,
    'SDSSJ0006+0520': 3,
    'SDSSJ0010-0430': 3,
    'SDSSJ0013+1109': 3,
    'SDSSJ0019+2209': 4,
    'SDSSJ0044+0418': 3,
    'SDSSJ0046+2717': 3,
    'SDSSJ0047+1628': 4,
    'SDSSJ0052+1846': 3,
    'SDSSJ0053+3115': 3,
    'SDSSJ0056+2453': 3,
    'SDSSJ0108-0537': 3,
    'SDSSJ0114+3505': 3,
    'SDSSJ0116+2050': 5,
    'SDSSJ0117+0021': 3,
    'SDSSJ0126+2534': 3,
    'SDSSJ0135+1302': 3,
    'SDSSJ0143+0113': 4,
    'SDSSJ0144+1920': 4,
    'SDSSJ0144+0305': 3,
    'SDSSJ0148-0112': 3,
    'SDSSJ0150+1354': 4,
    'SDSSJ0158-0942': 3,
    'SDSSJ0201+2015': 3,
    'SDSSJ0208-0542': 3,
    'SDSSJ0234-0510': 3,
    'SDSSJ0252-0401': 4,
    'SDSSJ0252+0054': 4,
    'SDSSJ0447+1124': 3,
    'SDSSJ0721+3928': 3,
    'SDSSJ0736+4118': 5,
    'SDSSJ0739+3112': 3,
    'SDSSJ0741+3146': 5,
    'SDSSJ0744+4649': 5,
    'SDSSJ0744+4408': 3,
    'SDSSJ0744+2701': 3,
    'SDSSJ0744+1640': 3,
    'SDSSJ0758+1013': 4,
    'SDSSJ0800+2242': 3,
    'SDSSJ0806+3055': 4,
    'SDSSJ0807+4930': 4,
    'SDSSJ0816+2330': 3,
    'SDSSJ0818+1247': 3,
    'SDSSJ0830-0319': 3,
    'SDSSJ0838+2322': 3,
    'SDSSJ0842+1406': 3,
    'SDSSJ0842+1536': 3,
    'SDSSJ0843+5614': 4,
    'SDSSJ0851+1543': 3,
    'SDSSJ0852+3402': 3,
    'SDSSJ0901+0752': 6,
    'SDSSJ0906+1141': 3,
    'SDSSJ0908+5136': 3,
    'SDSSJ0908+4119': 3,
    'SDSSJ0913+2627': 3,
    'SDSSJ0913+4127': 3,
    'SDSSJ0916+2540': 6,
    'SDSSJ0924+4301': 3,
    'SDSSJ0925+3130': 3,
    'SDSSJ0929+4247': 4,
    'SDSSJ0933+6334': 3,
    'SDSSJ0937+5228': 3,
    'SDSSJ0939+4136': 5,
    'SDSSJ0939+5019': 4,
    'SDSSJ0946+2024': 3,
    'SDSSJ0948+3008': 3,
    'SDSSJ0956+5912': 4,
    'SDSSJ1005+2244': 3,
    'SDSSJ1006+1752': 3,
    'SDSSJ1014+2827': 4,
    'SDSSJ1017+3447': 3,
    'SDSSJ1017+2419': 4,
    'SDSSJ1019+3535': 3,
    'SDSSJ1019+2045': 3,
    'SDSSJ1024+4531': 3,
    'SDSSJ1024+1014': 5,
    'SDSSJ1032+1338': 4,
    'SDSSJ1033+1809': 3,
    'SDSSJ1038-0036': 4,
    'SDSSJ1038+0432': 3,
    'SDSSJ1040+2407': 6,
    'SDSSJ1041+3432': 3,
    'SDSSJ1043+3516': 5,
    'SDSSJ1046+1329': 3,
    'SDSSJ1055+3725': 4,
    'SDSSJ1058+3143': 3,
    'SDSSJ1102+2827': 4,
    'SDSSJ1102+0214': 3,
    'SDSSJ1103+4144': 4,
    'SDSSJ1105+0228': 3,
    'SDSSJ1112+0700': 3,
    'SDSSJ1132+3323': 3,
    'SDSSJ1134+1542': 3,
    'SDSSJ1144+3720': 3,
    'SDSSJ1144+1218': 4,
    'SDSSJ1147+5429': 3,
    'SDSSJ1149+0519': 3,
    'SDSSJ1150+4928': 3,
    'SDSSJ1152+5101': 3,
    'SDSSJ1157+6138': 3,
    'SDSSJ1158+0454': 4,
    'SDSSJ1158+1845': 3,
    'SDSSJ1158+4712': 3,
    'SDSSJ1158+5448': 3,
    'SDSSJ1158+5942': 3,
    'SDSSJ1205+3536': 4,
    'SDSSJ1211+2326': 4,
    'SDSSJ1217+1157': 3,
    'SDSSJ1218+0023': 3,
    'SDSSJ1220+0929': 4,
    'SDSSJ1224+2838': 3,
    'SDSSJ1229+0743': 5,
    'SDSSJ1230+3143': 3,
    'SDSSJ1234+5208': 6,
    'SDSSJ1238+2149': 4,
    'SDSSJ1245+0822': 4,
    'SDSSJ1254+3551': 3,
    'SDSSJ1257+3238': 3,
    'SDSSJ1257-0310': 3,
    'SDSSJ1259+3112': 3,
    'SDSSJ1259+4729': 3,
    'SDSSJ1303+4055': 3,
    'SDSSJ1308+0957': 3,
    'SDSSJ1308+0258': 3,
    'SDSSJ1314+3748': 4,
    'SDSSJ1316+1918': 3,
    'SDSSJ1319+3641': 3,
    'SDSSJ1320+0204': 3,
    'SDSSJ1321-0237': 5,
    'SDSSJ1329+1301': 3,
    'SDSSJ1336+3547': 5,
    'SDSSJ1339+2643': 3,
    'SDSSJ1340+2702': 4,
    'SDSSJ1342+1813': 3,
    'SDSSJ1345+1153': 4,
    'SDSSJ1347+1415': 3,
    'SDSSJ1350+1058': 3,
    'SDSSJ1351+2645': 3,
    'SDSSJ1356+2416': 3,
    'SDSSJ1356+0236': 3,
    'SDSSJ1401+3659': 3,
    'SDSSJ1404+3620': 3,
    'SDSSJ1405+2542': 3,
    'SDSSJ1405+1549': 4,
    'SDSSJ1411+3410': 5,
    'SDSSJ1421+1843': 3,
    'SDSSJ1428+4403': 3,
    'SDSSJ1429+3841': 3,
    'SDSSJ1430-0151': 6,
    'SDSSJ1443+5833': 3,
    'SDSSJ1443+3014': 3,
    'SDSSJ1445+0913': 4,
    'SDSSJ1448+1047': 3,
    'SDSSJ1500+2315': 3,
    'SDSSJ1502+3744': 3,
    'SDSSJ1507+4034': 3,
    'SDSSJ1518+0506': 3,
    'SDSSJ1524+4049': 5,
    'SDSSJ1534+1242': 3,
    'SDSSJ1535+1247': 7,
    'SDSSJ1537+3608': 3,
    'SDSSJ1540+5352': 3,
    'SDSSJ1542+4650': 4,
    'SDSSJ1543+2024': 3,
    'SDSSJ1545+5236': 3,
    'SDSSJ1549+2633': 3,
    'SDSSJ1549+1906': 3,
    'SDSSJ1554+1735': 4,
    'SDSSJ1604+1830': 3,
    'SDSSJ1610+4006': 3,
    'SDSSJ1612+3534': 3,
    'SDSSJ1616+3303': 3,
    'SDSSJ1624+3310': 3,
    'SDSSJ1626+3303': 3,
    'SDSSJ1627+4646': 3,
    'SDSSJ1636+1619': 3,
    'SDSSJ1641+1856': 3,
    'SDSSJ1649+2238': 4,
    'SDSSJ1706+2541': 3,
    'SDSSJ2109-0039': 3,
    'SDSSJ2110+0512': 3,
    'SDSSJ2123+0016': 4,
    'SDSSJ2157+1206': 4,
    'SDSSJ2225+2338': 3,
    'SDSSJ2230+1905': 5,
    'SDSSJ2231+0906': 3,
    'SDSSJ2235-0056': 3,
    'SDSSJ2238+0213': 3,
    'SDSSJ2238-0113': 3,
    'SDSSJ2304+2415': 3,
    'SDSSJ2319+3018': 4,
    'SDSSJ2328+0830': 3,
    'SDSSJ2330+2805': 3,
    'SDSSJ2333+1058': 4,
    'SDSSJ2340+0124': 3,
    'SDSSJ2340+0817': 3,
    'SDSSJ2343-0010': 3,
    'SDSSJ2352+1922': 3,
    'SDSSJ2352+3344': 3,
    'SDSSJ2357+2348': 3
}

def find_xlsx_files(path=None):
    xlsx_files = list()
    for filename in os.listdir(path):
        if filename.endswith('.xlsx'):
            xlsx_files.append(filename)
    xlsx_files.sort()
    return xlsx_files

def collect_best_fits():
    best_fit_dict = dict()
    na_string = 'NA'
    volatile_rich_temp = 1000  # Below this temp, we're volatile rich
    volatile_poor_temp = 1400 # Above this temp, we're volatile poor (including moderate volatiles) (loosely basing this on Lodders 2003 table 8)
    for xlsx_file in find_xlsx_files():
        system = xlsx_file.split('PWDOutputs')[0]
        if system.startswith('Gaia') or system.startswith('SDSSJ2047') or system.startswith('PG') or system.startswith('WD'):
            continue # Not from the Hollands data
        print(system)
        workbook = xlrd.open_workbook(xlsx_file)
        sheet = workbook.sheet_by_index(0)
        
        model_list = sheet.row_values(0)
        evidences = sheet.row_values(1)
        chi_sq_list = sheet.row_values(2)
        num_elements = elements_present_dict[system]

        core_diff_sigma = sheet.cell_value(5, 15)  #P6
        crust_diff_sigma = sheet.cell_value(6, 15)  #P7
        temperature = sheet.cell_value(10, 15)  #P11
        
        steady_state_sigma = sheet.cell_value(3, 15)  #P4
        declining_phase_sigma = sheet.cell_value(8, 15)  #P9
        
        if temperature == '':
            temperature = 0
        
        
        highest_evidences = [np.NINF]
        highest_evidence_indices = [0]
        
        for index, evidence in enumerate(evidences):
            try:
                if evidence > highest_evidences[0]:
                    highest_evidences = [evidence]
                    highest_evidence_indices = [index]
                elif evidence == highest_evidences[0]:
                    highest_evidences.append(evidence)
                    highest_evidence_indices.append(index)
            except TypeError:
                # Then it was blank, or otherwise ignorable text
                pass
                
        # If this isn't true then we need to know
        assert len(highest_evidences) == 1
        assert len(highest_evidence_indices) == 1
        
        best_model = model_list[highest_evidence_indices[0]]
        chi_sq = chi_sq_list[highest_evidence_indices[0]]
        chi_sq_per_data_point = chi_sq/num_elements
        best_model_evidence = highest_evidences[0]
        best_model_description = model_params_dict[best_model]
        good_fit = chi_sq_per_data_point < 1
        
        primitive = core_diff_sigma == na_string and crust_diff_sigma == na_string
        
        fragment_core_frac = None
        fragment_crust_frac = None
        fragment_mantle_frac = None
        core_rich = False
        crust_rich = False
        mantle_rich = False
        
        # True: we will determine whether something is core/crust/mantle rich based on whether the fragment is enriched in those layers relative to its parent
        # False: we will determine whether something is core/crust/mantle rich based on whether the fragment is primarily core/crust/mantle
        measure_layers_relative_to_parent = True
        
        if not primitive:
            fragment_core_frac = sheet.cell_value(4, 8)  #I5
            if fragment_core_frac == '':
                fragment_core_frac = 0
            fragment_crust_frac = sheet.cell_value(4, 9)  #J5
            if fragment_crust_frac == '':
                fragment_crust_frac = 0
            fragment_mantle_frac = 1 - (fragment_core_frac + fragment_crust_frac)
            
            if measure_layers_relative_to_parent:
                
                delta_core = sheet.cell_value(11, 15)  #P12
                if delta_core == '':
                    delta_core = 0
                delta_crust = sheet.cell_value(12, 15)  #P13
                if delta_crust == '':
                    delta_crust = 0
                else:
                    # This fixes a bug in the original code which assumed an earth value of 0.01, not 0.001.
                    # So it subtracts 0.01 in the default case, ie the cases where a parent crust fragment is not specified
                    # This corrects these cases by adding 0.009
                    if "Parent Crust Fraction" not in best_model_description:
                        delta_crust += 0.009
                delta_mantle = 0 - (delta_core + delta_crust)
                
                core_rich = delta_core > delta_crust and delta_core > delta_mantle
                crust_rich = delta_crust > delta_core and delta_crust > delta_mantle
                mantle_rich = delta_mantle > delta_crust and delta_mantle > delta_core
            else:
                core_rich = fragment_core_frac > fragment_crust_frac and fragment_core_frac > fragment_mantle_frac
                crust_rich = fragment_crust_frac > fragment_core_frac and fragment_crust_frac > fragment_mantle_frac
                mantle_rich = fragment_mantle_frac > fragment_crust_frac and fragment_mantle_frac > fragment_core_frac
        
        volatile_rich = temperature < volatile_rich_temp
        volatile_depleted = temperature >= volatile_rich_temp and temperature <= volatile_poor_temp
        moderate_volatile_depleted = temperature > volatile_poor_temp

        buildup_phase = steady_state_sigma == na_string and declining_phase_sigma == na_string
        steady_state_phase = False
        declining_phase = False
        if not buildup_phase:
            if declining_phase_sigma == na_string:
                steady_state_phase = True
            else:
                declining_phase = True

        best_fit_dict[system] = {
            'Model': best_model,
            'Evidence': best_model_evidence,
            'ChiSq': chi_sq,
            'NumElements': num_elements,
            'ChiSqPerDataPoint': chi_sq_per_data_point,
            'GoodFit': good_fit,
            'Primitive': primitive,
            'CoreRich': core_rich,
            'MantleRich': mantle_rich,
            'CrustRich': crust_rich,
            'VolatileRich': volatile_rich,
            'VolatileDepleted': volatile_depleted,
            'ModerateVolatileDepleted': moderate_volatile_depleted,
            'BuildUpPhase': buildup_phase,
            'SteadyStatePhase': steady_state_phase,
            'DecliningPhase': declining_phase
        }
        
    return best_fit_dict
    
def write_best_fits(best_fit_dict):
    with open('best_fits_hb20.csv', 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            'System',
            'Best fit model',
            'Evidence',
            'Chi Squared',
            'Number of Elements',
            'Chi Squared per data point',
            'Good Fit?',
            'Primitive?',
            'Core Rich?',
            'Mantle Rich?',
            'Crust Rich?',
            'Volatile Rich?',
            'Volatile Depleted?',
            'Moderate Volatile Depleted?',
            'Build Up Phase?',
            'Steady State Phase?',
            'Declining Phase?',
            'Model parameters'
        ])
        for sys, vals in best_fit_dict.items():
            to_write.writerow([
                sys,
                vals['Model'],
                vals['Evidence'],
                vals['ChiSq'],
                vals['NumElements'],
                vals['ChiSqPerDataPoint'],
                vals['GoodFit'],
                vals['Primitive'],
                vals['CoreRich'],
                vals['MantleRich'],
                vals['CrustRich'],
                vals['VolatileRich'],
                vals['VolatileDepleted'],
                vals['ModerateVolatileDepleted'],
                vals['BuildUpPhase'],
                vals['SteadyStatePhase'],
                vals['DecliningPhase'],
                model_params_dict[vals['Model']]
            ])

def write_summary_table(best_fit_dict):
    unclassified = 0
    classifications = {
        'VolatileRich': {
            'Primitive': 0,
            'CoreRich': 0,
            'MantleRich': 0,
            'CrustRich': 0
        },
        'VolatileDepleted': {
            'Primitive': 0,
            'CoreRich': 0,
            'MantleRich': 0,
            'CrustRich': 0
        },
        'ModerateVolatileDepleted': {
            'Primitive': 0,
            'CoreRich': 0,
            'MantleRich': 0,
            'CrustRich': 0
        }
    }
    phases = {
        'BuildUpPhase': 0,
        'SteadyStatePhase': 0,
        'DecliningPhase': 0
    }
    for sys, vals in best_fit_dict.items():
        if not vals['GoodFit']:
            unclassified += 1
        else:
            if vals['VolatileRich']:
                t_classification = 'VolatileRich'
            elif vals['VolatileDepleted']:
                t_classification = 'VolatileDepleted'
            elif vals['ModerateVolatileDepleted']:
                t_classification = 'ModerateVolatileDepleted'
            
            if vals['Primitive']:
                c_classification = 'Primitive'
            elif vals['CoreRich']:
                c_classification = 'CoreRich'
            elif vals['MantleRich']:
                c_classification = 'MantleRich'
            elif vals['CrustRich']:
                c_classification = 'CrustRich'
                
            classifications[t_classification][c_classification] += 1
            if vals['BuildUpPhase']:
                phases['BuildUpPhase'] += 1
            elif vals['SteadyStatePhase']:
                phases['SteadyStatePhase'] += 1
            elif vals['DecliningPhase']:
                phases['DecliningPhase'] += 1

    with open('summary_table_hb20.csv', 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            'Unclassified:',
            unclassified
        ])
        to_write.writerow([])
        to_write.writerow([
            '',
            'Primitive',
            'Core Rich',
            'Mantle Rich',
            'Crust Rich'
        ])
        for t_classification, vals in classifications.items():
            to_write.writerow([
                t_classification,
                vals['Primitive'],
                vals['CoreRich'],
                vals['MantleRich'],
                vals['CrustRich']
            ])
        to_write.writerow([])
        to_write.writerow([
            'Build Up Phase',
            'Steady State Phase',
            'Declining Phase'
        ])
        to_write.writerow([
            phases['BuildUpPhase'],
            phases['SteadyStatePhase'],
            phases['DecliningPhase']
        ])

def main():
    best_fit_dict = collect_best_fits()
    write_best_fits(best_fit_dict)
    write_summary_table(best_fit_dict)
    
if __name__ == '__main__':
    main()
