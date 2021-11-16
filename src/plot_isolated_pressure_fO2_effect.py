#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import csv
import numpy as np
import random
import sys

import abundance_model as am
import chemistry_info as ci
import complete_model as cm
import geology_info as gi
import graph_factory as gf
import live_data as ld
import manager as mn
import pwd_utils as pu
import white_dwarf_model as wdm

from argparse import Namespace

sys.path.append(pu.get_path_to_utils())
sys.path.append(pu.get_path_to_original_src())

import dict_plotter as dp
import original_complete_model as ocm

colour_dict = {
    ci.Element.Ca: '#191970',
    ci.Element.Al: '#66CDAA',
    ci.Element.Fe: '#C0C0C0',
    ci.Element.Si: '#DF2E20',
    ci.Element.Mg: '#87CEFA',
    ci.Element.O: '#FFD700',
    ci.Element.C: '#7FFFD4',
    ci.Element.Cr: '#808000',
    ci.Element.Ni: '#22CD22',
    ci.Element.Placeholder: '#000000'
}

fise_run_args = {
    #These were run with the sisi (Siebert silicon) version
    #'GD61': [226, 437.37508960421, 5.567072577431, 0.296106531106386, 0.077322581088442, None, None, 0.026231629101023, None, -5.82699432826509, 7.10770252512147, 42.9965877524006, -2.47677571461788],
    #'WD0611-6931Full': [305, 411.9126192, 0.005003706711, 0.01447163126, 0.05, None, None, 0.7683208858, None, -3.751726203, 4.145535017, 29.65027528, -1.904223833],
    #'PG0843+516Gaensicke': [219, 418.349121227523, 0.004910999454207, -0.045130378685706, 0.05, None, None, 0.639718848123253, None, -4.26116833054434, 4.10364972393969, 27.0819223781586, -2.56976334306055],
    #'PG0843+516Xu': [220, 506.707928771645, 0.009118408135688, -0.008845612755553, 0.05, None, None, 0.728745771883592, None, -3.77012751244553, 4.37608752540848, 23.0904826689188, -1.76765362383884],
    #'NLTT43806': [227, 791.861144367747, 0.652759311705154, -0.058444624506319, 0.05, None, None, 0.027875631165411, None, -6.56447747052922, 6.27227175172384, 46.2001486427106, -2.67509568637878],
    #'SDSSJ0738+1835': [230, 402.243066054597, 3.81446489839799, 0.078757704335966, 0.05, None, None, 0.43722029937676, None, -3.7875350361852, 2.72629139196373, 31.2042441423048, -2.54048938334322],
    #'SDSSJ0738+1835Corr': [363, 375.900885063282, 0.49247702176291, 0.099327565963983, 0.05, None, None, 0.440526840524824, None, -3.78919557507339, 2.38503940644603, 32.3298145531553, -2.53861020691926],
    #'SDSSJ0845+2257': [234, 546.739784661982, 0.67204910652165, -0.263042850580013, 0.086961789532605, None, None, 0.331040611889366, None, -3.78522283127578, 6.39639198601126, 29.5502718185378, -2.01081246790329],
    #'GD 378': [238, 822.745317332263, 1.14996669381006, -0.581638291228524, 0.132955598673431, None, None, 0.054500320061958, None, -5.6162722931722, 4.6532747836274, 37.0105226526938, -2.2182777080587],
    #'WD0449-259': [246, 904.781458244381, 16.8028612492612, 0.157516477093874, 0.05, None, None, 0.66877586362797, None, -6.88101618978015, 3.27787535095666, 19.157755678069, -1.36618477730749],
    #'WD1350-162': [247, 849.651227784859, 7.59645704009145, -0.053608382038046, 0.05, None, None, 0.39945856964945, None, -5.99350642223607, 2.99142162619961, 20.0942317714276, -1.95588390243148],
    #'WD0122-227': [245, 485.616928387098, 7.8156937259736, 0.03441691784735, 0.05, None, None, 0.567305307981525, None, -8.2243016670192, 3.89834523189547, 29.4379454699773, -1.99238401733305],
    #'SDSSJ0512-0505': [249, 501.475293813679, 3.65911964115896, -1.03730199796092, 0.130668392751986, None, None, 0.545002949374465, None, -7.38179071174634, 3.18501234066977, 25.2549663467436, -1.77609214975635],
    #'SDSSJ0823+0546': [250, 430.984857988105, 17.6224780818674, 0.019933303211355, 0.05, None, None, 0.956795923975906, None, -7.23008858378296, 3.2368803133493, 29.0070252532306, -2.10081082736568],
    #'WD0611-6931Full': [305, 409.186098187532, 0.004823887570456, 0.014358832466395, 0.05, None, None, 0.768404825227812, None, -3.75258757440603, 4.13863172530002, 29.661637997752, -1.91571596482581],
    #'LHS2534': [257, 926.238972508352, 4.70497840145297, 0.550173774272923, 0.05, None, None, 0.518208753763491, None, -8.37966354619442, 3.61489926028609, 35.415573421677, -1.28532198095239],
    #'WD2105-820': [248, 210.87845678345, 0.0000961997252700482, 0.071917807890479, 0.05, None, None, 0.838348492316552, None, -5.95926441189598, 1.53448835011035,18.7277178945314,-1.6750110526307],
    #'GALEX1931+0117Melis': [306,562.466911343358,0.00179454738721,0.346798759498556,0.072993612911824,None,None,0.302430342279287,None,-3.66186108173482,3.67868342589617,15.9961525686632,-2.63547009336713],
    #'IdealHPMantle1': [329,256.253137538389,1.25772345779633,-0.641462521190879,0.125490238005588,None,None,0.002505105355787,None,-4.01460383401613,6.4876927301867,56.4533456889211,-2.88191070564759],
    #'IdealHPCore1': [330,616.996788331321,0.185897746718775,-0.505758096744753,0.05,None,None,0.993936868451227,None,-3.83941711670438,4.54786319355762,50.2971106091221,-1.34753838938756],
    #'IdealLPMantle1': [331,749.500098581405,0.651878222558579,-0.578159080494695,0.110474020121104,None,None,1.19943771354248E-06,None,-4.00533442170302,2.68203528834887,0.024131980675524,-2.98008745573218],
    #'GALEX1931+0117': [222,471.786463177388,0.00366409904558,2,0.05,None,None,None,None,-3.73231616812409,3.97212248224808,None,None],
    ##'GaiaJ0510+2315': [301,593.431875941104,0.003758664511902,2,0.05,None,None,None,None,-4.10796214773616,4.02014846228164,None,None]  # The median values
    #'GaiaJ0510+2315': [301,581.3,0.000522,2,0.05,None,None,None,None,-4.111,2.727,None,None], # The Max Likelihood values (they fit better!)
    ##'GaiaJ0644-0352': [302,288.619992153966,0.285278451983339,-0.640580889530974,0.05,None,None,None,None,-5.02224885319909,3.75510039449635,None,None], # Highest log Z model, median values
    #'GaiaJ0644-0352': [302,876.1,59.4638320940943288,-0.576849933829693073,0.05,None,None,0.0708865391912797821,None,-4.99875561537941060,7.93747927503662343,44.3974402943637756,-2.37815863222664259], # Best fit model, max likelihood values
    #'SynthEarthfcf0': [307,384.560579847084,0.16726184373474,-0.550590766715912,0.092326209761797,None,None,0.013857052365627,None,-3.99157181454789,3.83744247471605,23.4507397312912,-2.35147228174733],
    #'WD0449-259NoNa': [342,557.317452003682,6.97618833434471,-0.985765500119387,0.103870359979021,None,None,0.637580010112239,None,-7.45781133397279,4.01805779444807,20.1788399789155,-1.44238021207705],
    #'WD1350-162NoNa': [343,171.818212941235,7.31541433533503,0.025577571044642,0.05,None,None,0.414381955635034,None,-6.03824452251335,3.14227638271032,23.0031345846313,-2.08319750488267],
    #'IdealHPMantle4': [340,616.187306876256,0.101288647261848,-0.509276468517034,0.05,None,None,0.012473279169699,None,-4.0068215773501,3.92043001595668,46.6028000969223,-2.12808997103755],
    #'IdealLPMantle2': [337,320.542433461403,0.647239361906416,-0.750639925657234,0.141993949232361,None,None,4.16153157149695E-06,None,-3.98617291312008,2.7204957285566,0.087532511419774,-2.95460193289795],
    #'IdealHPCfcf75err5': [397,356.91492883401,11.8351315556828,-0.584288132998062,0.097038753008786,None,None,0.824414901894526,None,-4.02152647155219,7.4333379374253,11.7852127666372,-2.40027437365142],
    #'SynthEarthfcf0.7': [314,596.938597489607,0.150314806568319,-0.510798429826005,0.05,None,None,0.66802004674132,None,-4.01984252881437,4.25676799308911,46.1640507207147,-2.75832533857166],
    #'SynthEarthfcf0.8': [315,383.01136690618,4.99677630667471,-0.584687934079343,0.099466147558644,None,None,0.855818734724399,None,-4.00660297203982,7.14425290494588,24.2043203306747,-2.77307189857887],
    #'WD1350-162NoNaCorr': [377,162.019231367967,6.83503111540107,0.01033097266606,0.05,None,None,0.414644299207587,None,-6.03802276078376,3.12952116775951,22.4471729043176,-2.05890354283997],
    #'GALEX1931+0117GCorr': [351,459.436688813984,0.002795356013718,2,0.05,None,None,None,None,-3.74182233851499,3.85652188295597,54,-2],
    #'GALEX1931+0117VCorr': [352,137.895272219274,0.003483767472115,0.552328519115545,0.05,None,None,0.158660318708795,None,-3.46325777513603,3.93522350145396,45.0079465185255,-2.62218607860118],
    #'GALEX1931+0117MCorr': [353,561.534557730359,0.000702881452596,0.433069991893536,0.070534707971139,None,None,0.314143312681015,None,-3.58125092627202,3.06579389404442,51.5925458254742,-2.73855772669102]
    #'WD2105-820': [248,209.665823846459,9.58421401303261E-05,0.070881042833234,0.05,None,None,0.837580130963948,None,-5.95833726628792,1.53608665508434,18.8654573745211,-1.67716051960656],
    #'WD1232+563Corr': [373,559.52104594338,13.2447022443142,2,0.05,None,None,None,None,-5.10106592894453,7.48555150602173,54,-2]

    # Systems to include in bow tie:
    'G166-58': [413,None,None,None,None,None,None,None,None,None,None,None,None],
    'G241-6Corr': [359,None,None,None,None,None,None,None,None,None,None,None,None],
    'G29-38Corr': [366,None,None,None,None,None,None,None,None,None,None,None,None],
    'GALEX1931+0117GCorr': [351,None,None,None,None,None,None,None,None,None,None,None,None],
    'GALEXJ2339': [260,None,None,None,None,None,None,None,None,None,None,None,None],
    'GD362Corr': [357,None,None,None,None,None,None,None,None,None,None,None,None],
    'GD378': [261,None,None,None,None,None,None,None,None,None,None,None,None],
    'GD40Corr': [358,None,None,None,None,None,None,None,None,None,None,None,None],
    'GD424': [253,None,None,None,None,None,None,None,None,None,None,None,None],
    'GD56': [419,None,None,None,None,None,None,None,None,None,None,None,None],
    'GD61Corr': [360,None,None,None,None,None,None,None,None,None,None,None,None],
    'HE0106-3253': [421,None,None,None,None,None,None,None,None,None,None,None,None],
    'HS2253+8023Corr': [362,None,None,None,None,None,None,None,None,None,None,None,None],
    'LHS2534': [257,None,None,None,None,None,None,None,None,None,None,None,None],
    'NLTT43806Corr': [361,None,None,None,None,None,None,None,None,None,None,None,None],
    'PG0843+516XCorr': [350,None,None,None,None,None,None,None,None,None,None,None,None],
    'PG1015+161Xu': [422,None,None,None,None,None,None,None,None,None,None,None,None],
    'PG1225-079Corr': [380,None,None,None,None,None,None,None,None,None,None,None,None],
    'SDSSJ0512-0505': [249,None,None,None,None,None,None,None,None,None,None,None,None],
    'SDSSJ0738+1835Corr': [363,None,None,None,None,None,None,None,None,None,None,None,None],
    'SDSSJ0823+0546': [250,None,None,None,None,None,None,None,None,None,None,None,None],
    'SDSSJ0845+2257Corr': [368,None,None,None,None,None,None,None,None,None,None,None,None],
    'SDSSJ1043+0855Corr': [370,None,None,None,None,None,None,None,None,None,None,None,None],
    'SDSSJ1228+1040Corr': [365,None,None,None,None,None,None,None,None,None,None,None,None],
    'SDSSJ1242+5226Corr': [367,None,None,None,None,None,None,None,None,None,None,None,None],
    'SDSSJ2047-1259Corr': [354,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD0122-227SiO': [345,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD0446-255': [347,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD0449-259NoNaCorr': [376,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD1145+017Corr': [375,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD1145+288': [417,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD1232+563Corr': [373,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD1350-162NoNaCorr': [377,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD1425+540Corr': [371,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD1536+520Corr': [369,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD1551+175Corr': [355,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD2105-820': [248,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD2115-560Corr': [372,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD2157-574': [240,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD2207+121Corr': [374,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD2216-657SiCorr': [378,None,None,None,None,None,None,None,None,None,None,None,None],
    'WD2230-125': [346,None,None,None,None,None,None,None,None,None,None,None,None],
    'WDJ1814-7354': [348,None,None,None,None,None,None,None,None,None,None,None,None]
}

limit_to_sample = True
systems_in_sample = [
    'G166-58', # No evidence of differentiation
    'G241-6Corr',
    'G29-38Corr',
    'GALEX1931+0117MCorr', # Poor fit
    'GALEXJ2339',
    'GD362Corr',
    'GD378',
    'GD40Corr',
    'GD424', # High pressure mantle
    'GD56',
    'GD61Corr',
    'HE0106-3253', # Pressure degenerate
    'HS2253+8023Corr',
    'LHS2534',
    'NLTT43806Corr',
    'PG0843+516XCorr',
    'PG1015+161Xu', # Pressure unconstrained
    'PG1225-079Corr',
    'SDSSJ0512-0505',
    'SDSSJ0738+1835Corr',
    'SDSSJ0823+0546',
    'SDSSJ0845+2257Corr',
    'SDSSJ1043+0855Corr',
    'SDSSJ1228+1040Corr',
    'SDSSJ1242+5226Corr',
    'SDSSJ2047-1259Corr',
    'WD0122-227SiO',
    'WD0446-255',
    'WD0449-259NoNaCorr',
    'WD1145+017Corr',
    'WD1145+288',
    'WD1232+563Corr',
    'WD1350-162NoNaCorr',
    'WD1425+540Corr',
    'WD1536+520Corr',
    'WD1551+175Corr',
    'WD2105-820',
    'WD2115-560Corr',
    'WD2157-574',
    'WD2207+121Corr',
    'WD2216-657SiCorr',
    'WDJ1814-7354'
]


system_categories = {
    'G166-58': 'NED', # No evidence of differentiation
    'G241-6Corr': 'NED',
    'G29-38Corr': 'NED',
    'GALEX1931+0117MCorr': 'NED',
    'GALEXJ2339': 'NED',
    'GD362Corr': 'NED',
    'GD378': 'NED',
    'GD40Corr': 'NED',
    'GD424': 'NED', 
    'GD56': 'NED',
    'GD61Corr': 'HPM', # High pressure mantle
    'HE0106-3253': 'Pdegen', # Pressure degenerate
    'HS2253+8023Corr': 'NED',
    'LHS2534': 'Unphysical',
    'NLTT43806Corr': 'HPM',
    'PG0843+516XCorr': 'Pdegen',
    'PG1015+161Xu': 'Pdegen', 
    'PG1225-079Corr': 'NED',
    'SDSSJ0512-0505': 'Pdegen',
    'SDSSJ0738+1835Corr': 'Puncon', # Pressure unconstrained
    'SDSSJ0823+0546': 'Puncon',
    'SDSSJ0845+2257Corr': 'Puncon',
    'SDSSJ1043+0855Corr': 'NED',
    'SDSSJ1228+1040Corr': 'NED',
    'SDSSJ1242+5226Corr': 'NED',
    'SDSSJ2047-1259Corr': 'NED',
    'WD0122-227SiO': 'Puncon',
    'WD0446-255': 'HPM',
    'WD0449-259NoNaCorr': 'LPC',  # Low pressure core
    'WD1145+017Corr': 'NED',
    'WD1145+288': 'Puncon',
    'WD1232+563Corr': 'NED',
    'WD1350-162NoNaCorr': 'LPC',
    'WD1425+540Corr': 'NED',
    'WD1536+520Corr': 'NED',
    'WD1551+175Corr': 'NED',
    'WD2105-820': 'LPC',
    'WD2115-560Corr': 'NED',
    'WD2157-574': 'NED',
    'WD2207+121Corr': 'NED',
    'WD2216-657SiCorr': 'NED',
    'WDJ1814-7354': 'NED'
}

def generate_pressure_vals():
    return list()
    print('Warning: using only 2 pressure values!')
    return [1, 60]
    return list(range(0, 10, 1)) + list(range(10, 55, 5)) + [54] + list(range(55, 105, 5))
    #return list(range(0, 101, 1))
    
def generate_fO2_vals():
    return list()
    return list(range(-3, 0, 1))
    
def generate_fcf_vals():
    return list()
    return [0.1, 0.2]
    return [0, 0.0001, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 0.9999, 1]

def generate_synthetic_wd_data(wd_abundances, wd_errors, N_synth=10000):
    # Structure of synthetic data to be the same as structure of real data ie a list of N lists of 11 (X/Mg) ratios (in a dict)
    synth_data = dict()
    for wd_name, abundances in wd_abundances.items():
        toret = list()
        i = 0
        while i < N_synth:
            synthetic_wd = list()
            for element in ci.usual_elements:
                XH_scatter = np.random.normal(0, wd_errors[wd_name][element])
                synthetic_value = abundances[element] + XH_scatter
                synthetic_wd.append(synthetic_value)
            to_append = list()
            mg_index = ci.usual_elements.index(ci.Element.Mg)
            for el_index, element_abundance in enumerate(synthetic_wd):
                if el_index != mg_index:
                    to_append.append(10**(element_abundance-synthetic_wd[mg_index]))
            toret.append(to_append)
            i += 1
        synth_data[wd_name] = toret
    return synth_data

def generate_synthetic_stellar_data(stellar_data, N_synth=10000):
    toret = list()  # Structure of synthetic data to be the same as structure of real data ie a list of N lists of 11 (X/Mg) ratios
    i = 0
    sigma_dict = {
        ci.Element.Al: 0.175,  # These values are from my Hollands2017WhiteDwarfObservationalDataCompositionsExpanded sheet
        ci.Element.Ti: 0.145,  # I just did a mean of the relevant WD data non-zero errors
        ci.Element.Ca: 0.227,  # Exception: C and N have no observations so assuming standard 0.1 error
        ci.Element.Ni: 0.157,
        ci.Element.Fe: 0.291,
        ci.Element.Cr: 0.176,
        ci.Element.Mg: 0.24,
        ci.Element.Si: 0.123,
        ci.Element.Na: 0.172,
        ci.Element.O: 0.153,
        ci.Element.C: 0.1,
        ci.Element.N: 0.1
    }

    while i < N_synth:
        random_index = random.choice(list(range(0, len(stellar_data))))
        base_star = stellar_data[random_index]
        synthetic_star = list()
        el_index = 0
        MgH_scatter = np.random.normal(0, sigma_dict[ci.Element.Mg])   # Make sure Mg scatter is consistent for each star
        for element in ci.usual_elements:
            if element is not ci.Element.Mg:
                # We want to put scatter on the raw (logarithmic) X/H data
                # But the base_star data is expressed as X/Mg (non-logarithmic)
                # So base_star[X] = 10^(log(X/H) - log(Mg/H))
                # We need to recover log(X/H) and log(Mg/H), but we don't know H
                # So assume a value for the (pre-scatter) log(Mg/H)
                # then recover log(X/H) as log(Mg/H) + log10(base_star[X])
                # then put scatter on log(X/H) and log(Mg/H) and recalculate X/Mg
                # But actually the assumed log(Mg/H) shouldn't matter in the end
                # It all cancels down to
                # new value = base_star[X] * 10 ^ (XH_scatter - MgH_scatter)
                
                XH_scatter = np.random.normal(0, sigma_dict[element])
                synthetic_value = base_star[el_index] * (10**(XH_scatter - MgH_scatter))
                synthetic_star.append(synthetic_value)
                
                el_index += 1
        toret.append(synthetic_star)
        i += 1
    return toret

def load_generic_float_data_csv(input_filename):
    # TODO: Change this to a with clause to prevent ResourceWarnings
    generic_csv = open(pu.get_path_to_data() + input_filename)
    generic_list =  [row for row in csv.reader(generic_csv)]
    if input_filename == 'Hollands2017WhiteDwarfObservationalDataCompositions.csv':
        generic_list.append(['-6.5', '0.2', '0', '0', '0', '0', '-6.3', '0.3', '-4.6', '0.2', '-5.8', '0.3', '-5', '0.2', '-5.2', '0.2', '0', '0', '-5', '0.3', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'])  # PG0843+516
        generic_list.append(['-6.99', '0.15', '-8.68', '0.11', '-6.93', '0.07', '0', '0', '-6.6', '0.1', '-8.25', '0.07', '-6.29', '0.05', '-6.33', '0.1', '0', '0', '-5.48', '0.15', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'])  # WD1551+175    
    if input_filename == 'WhiteDwarfObservationalDataTimescales.csv':
        # The timescales for PG... are super short, but the model seems to fit the t_sinceaccretion to vaguely normal values?! ==> lots of nan abundances (ie 0)
        generic_list.append(['0.006', '0.0035', '0.0036', '0.0025', '0.0027', '0.003', '0.0059', '0.0052', '0.002', '0.003', '0.006', '0.004'])  # PG0843+516
        generic_list.append(['837529', '523600', '568853', '576766', '559758', '540754', '866962', '843335', '833681', '862979', '961612', '881049'])  # WD1551+175    
    generic_array = np.asarray(generic_list)
    return generic_array.astype(np.float)

def run_complete_model(trial_fits, extended_fits):
    # This section is basically replicating some set up functionality in the Manager which I want to avoid using directly
    observations = load_generic_float_data_csv('Hollands2017HB20CompositionsSpuriousRemoved.csv')#'Hollands2017HB20CompositionsSpuriousRemoved.csv')
    timescales = load_generic_float_data_csv('Hollands2017HB20Timescales.csv')#'Hollands2017HB20Timescales.csv')
    
    wd_data = dict()
    wd_abundances = dict()
    wd_timescales = dict()
    for i in range(0, len(observations)):
        wd_data[i] = observations[i][:]
    for i in range(0, len(timescales)):
        wd_timescales[i] = timescales[i][:]
    for i in range(0, len(observations)):
        wd_abundances[i] = np.transpose([
            wd_data[i][0],
            wd_data[i][2],
            wd_data[i][4],
            wd_data[i][6],
            wd_data[i][8],
            wd_data[i][10],
            wd_data[i][12],
            wd_data[i][14],
            wd_data[i][16],
            wd_data[i][18],
            wd_data[i][20],
            wd_data[i][22]
        ])
    
    ld._live_stellar_compositions = load_generic_float_data_csv('StellarCompositionsSortFE.csv')
    ld._live_model = 'Model_Full_No_Crust'
    
    data_dump = dict()
    abundances_dict = dict()
    abundance_lower_bounds_dict = dict()
    abundance_upper_bounds_dict = dict()
    errors_dict = dict()
    stellar_data = list()
    
    args = {
        # Args are: wd observation index, fe_star, t_sinceaccretion, d_formation, z_formation, N_c, N_o, f_c, f_o, pollutionfraction, t_disc
        #'Run1': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        #'Run2': [145, 0, 2, 0.05, 0.17, 0.01, 0.17, 0, -1, 5.6],
        #'Run3': [145, 0.2, 2, 0.05, 0.17, 0.01, 0.17, 0, -1, 5.6],
        #'Run4': [145, 0.5, 2, 0.05, 0.17, 0.01, 0.17, 0, -1, 5.6]
        #'Run5': [0, 479, 0, 2, 0.05, 0.17, 0, 0.7, 0, 0, 5.6],
        #'Run6': [0, 479, 0, 2, 0.05, 0.17, 0, 0, 0, 0, 5.6],
        # The following are for an earlier model (fise)
        #'PG0843+516': [208, 443, 0.003064, 0.091101, 0.05, 0.17, 0.01, 0.705035, 0.00454, -4.242601, 3.915285],
        #'WD1551+175': [209, 610, 0.301838, -0.620723, 0.05, 0.17, 0.01, 0.103149, 0.037161, -5.428173, 4.920082],
        #'SDSSJ1535+1247': [0, 88, 0.195758, 1.248983, 0.05, 0.17, 0.01, 0.128178, 0.01, -7.05483, 4.773768],
        #'SDSSJ1430-0151': [10, 257, 0.20625, 1.027286, 0.05, 0.17, 0.01, 0.075559, 0.01, -6.278862, 4.447416],
        #'SDSSJ0150+1354': [20, 445, 4.545403, -1.175611, 0.1, 0.17, 0.01, 0.128163, 0.638003, -6.481734, 3.1132],
        #'SDSSJ0512-0505': [21, 707, 0.330673, -0.500687, 0.05, 0.17, 0.01, 0.185142, 0.009168, -7.579329, 5.046344],
        #'SDSSJ0736+4118': [30, 502, 0.178597, -0.563, 0.05, 0.17, 0.01, 0.158098, 0.24458, -7.180743, 4.563162],
        #'SDSSJ1144+1218': [51, 430, 1.138125, 1.195405, 0.05, 0.17, 0.01, 0.105049, 0.099352, -7.958425, 5.674139],
        #'SDSSJ1055+3725': [91, 545, 0.142434, 1.102149, 0.07848, 0.17, 0.01, 0.144674, 0.739604, -7.408492, 4.786214],
        #'SDSSJ1627+4646': [95, 530, 0.988052, 0.22005, 0.05, 0.102405, 0.063971, 0.628905, 0.005141, -7.090077, 3.943588],
        #'SDSSJ0046+2717': [97, 499, 1.581918, 1.161042, 0.05, 0.17, 0.01, 0.071396, 0.01, -5.9609, 5.073561]
    }
    full_run_args = {
        # Args are: wd observation index, fe_star, t_sinceaccretion, d_formation, z_formation, N_c, N_o, f_c, f_o, pollutionfraction, t_disc, pressure, fO2
        #'SDSSJ1535+1247': [0, 155.904566092763957, 0.143913612069693386, -0.170827464668557916, 0.121024653548637312, 0.0453075222533427177, 0, 0.0841105119136084101, 0, -7.04417530428071448, 0.879914870128519655, 46.3447177213887755, -2.87986220761787415],
        #'SDSSJ0046+2717': [97, 750.503679021099174, 1.09638869752662549, 0.520542804882032861, 0.128782086823985925, 0.149062660713768724, 0, 0.00290528146157534523, 0, -5.90598794770677493, 1.41890778600791401, 12.4388694844717964, -0.939350467005668044],
        #'PG0843+516': [208, 281.829778044233080, 0.905002653021461656, 0.611410066567565913, 0.117038090370890691, 0.169146281116199920, 0, 0.680358779789567159, 0, -4.17912696501628567, 6.50979661721710290, 6.92682107866911156, -4.43378748980489412]
        # Next 3 runs are using mean parameters rather than max likelihood
        #'SDSSJ1535+1247': [0, 162.558754743260693, 0.254128427443081062, 1.29064651004301534, 0.0756003743453170091, 0.103843962865981534, 0, 0.0923759475835635446, 0, -7.06503772615663106, 4.35733021265158715, 41.6306315670475229, -1.89128230060722702],
        #'SDSSJ0046+2717': [97, 497.566649374994824, 4.44983251793118928, 0.986849496686652672, 0.0800765049082564490, 0.101173441512605947, 0, 0.190268783122635343, 0, -5.93627694792752969, 4.72730656336368149, 50.9515380519172751, -2.55232983692679438],
        #'PG0843+516': [208, 506.384908064871354, 1.22015164643037965, 1.65900087151922726, 0.0765287567312275541, 0.119782946405400131, 0, 0.700244502198948293, 0, -4.08092789522224120, 3.93531788746016087, 45.6418294690524249, -2.58645754299788688]
    }
    manager = mn.Manager(
        Namespace(
            wd_data_filename='BlouinConglomNewTimescales.csv',
            stellar_compositions_filename='StellarCompositionsSortFE.csv',
            n_live_points = 0, # This argument shouldn't matter, in fact we only use the manager to access observational data so nothing else matters
            enhancement_model = 'NonEarthlike',
            base_dir = '.',
            pollution_model_names=['Model_24']
        )
    )
    for arg_name, arg_set in fise_run_args.items():
        if limit_to_sample and arg_name not in systems_in_sample:
            continue
        extra_fits = trial_fits.get(arg_name, list())
        extend_fits = extended_fits.get(arg_name, dict())
        print('Running fise results for ' + arg_name)
        observation_to_test_on = arg_set[0]
        #non_zero_wd_abundances = list()
        #non_zero_wd_timescales = list()
        data_dump[arg_name] = dict()
        #num_elements = 12
        #for i in range(0, num_elements):
        #    wd_abundance = wd_abundances[observation_to_test_on][i]
        #    if wd_abundance != 0:
        #        non_zero_wd_abundances.append(wd_abundance)
        #        non_zero_wd_timescales.append(wd_timescales[observation_to_test_on][i])
        #ld._live_planetesimal_abundance = np.zeros(num_elements)
        #ld._live_all_wd_abundances = np.transpose(wd_abundances[observation_to_test_on])
        #ld._live_non_zero_wd_timescales = np.transpose(non_zero_wd_timescales)

        print('Observation:')
        print(observation_to_test_on)
        manager.publish_live_data(observation_to_test_on)
        abundances_dict[arg_name] = manager.wd_abundances[arg_name]
        abundance_lower_bounds_dict[arg_name] = manager.wd_abundance_lower_bounds[arg_name]
        abundance_upper_bounds_dict[arg_name] = manager.wd_abundance_upper_bounds[arg_name]
        print(arg_name)
        print(abundances_dict[arg_name])
        #abundances_dict[arg_name][abundances_dict[arg_name] == 0] = None  # ...why is this here?
        errors_dict[arg_name] = manager.wd_errors[arg_name]
        #errors_dict[arg_name][errors_dict[arg_name] == 0] = None  # ...why is this here?
        for element in abundances_dict[arg_name]:
            if abundances_dict[arg_name][element] == 0:
                abundances_dict[arg_name][element] = np.nan
                errors_dict[arg_name][element] = np.nan
        we_care_about_actual_result = False        
        we_care_about_extra_fits = False
        if we_care_about_actual_result:
            print()
            print(arg_set)
            test_result = cm.complete_model_calculation(arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], arg_set[9], 10**(arg_set[10]), arg_set[11], arg_set[12], 'NonEarthlike')
            print('Ran test fit on  ' + arg_name + ' Diagnostics were: ')
            print('DiscAbundances')
            print(test_result[1]['DiscAbundances'])
            try:
                ld._geo_model.tabulate_output(test_result[1]['Enhancements']['Abundances'], test_result[1]['Enhancements']['ParentCoreNumberFraction'], test_result[1]['Enhancements']['Ds'])
            except TypeError:  # If no differentiation, then this output will be None
                pass
            print('FinalAbundances')
            print(test_result[0])
        else:
            test_result = (None, None)
        print('Result was: ')
        print(test_result[0])
        data_dump[arg_name][('LP run')] = test_result[0]
        if we_care_about_actual_result:
            og_arg_set = args.get(arg_name)
            try:
                data_dump[arg_name][('Model24')] = ocm.PWDCodeMultiple(og_arg_set[0], og_arg_set[0]+1, False, og_arg_set[1], og_arg_set[2], og_arg_set[3], og_arg_set[4], og_arg_set[5], og_arg_set[6], og_arg_set[7], og_arg_set[8], og_arg_set[9], og_arg_set[10])
            except (TypeError, IndexError) as e:
                pass # Occurs for PG... and WD... because they don't appear in original dataset. Also if no args exist
        for ef in extra_fits:
            if we_care_about_extra_fits:
                result = cm.complete_model_calculation(arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], ef[0], arg_set[8], arg_set[9], 10**(arg_set[10]), ef[1], ef[2], 'NonEarthlike')
                data_dump[arg_name][ef] = result[0]
                print('Ran fit ' + str(ef) + ' on system ' + arg_name + ' Diagnostics were: ')
                print('DiscAbundances')
                print(result[1]['DiscAbundances'])
                ld._geo_model.tabulate_output(result[1]['Enhancements']['Abundances'], result[1]['Enhancements']['ParentCoreNumberFraction'], result[1]['Enhancements']['Ds'])
                print('FinalAbundances')
                print(result[0])
            else:
                data_dump[arg_name][ef] = None
        for ef_name, ef in extend_fits.items():
            if we_care_about_extra_fits:
                result = cm.complete_model_calculation(ef[0], ef[1], ef[2], ef[3], ef[4], ef[5], ef[6], ef[7], ef[8], 10**(ef[9]), ef[10], ef[11], 'NonEarthlike')
                data_dump[arg_name][ef_name] = result[0]
                print('Ran fit ' + str(ef) + ' on system ' + arg_name + ' Diagnostics were: ')
                print('DiscAbundances')
                print(result[1]['DiscAbundances'])
                try:
                    ld._geo_model.tabulate_output(result[1]['Enhancements']['Abundances'], result[1]['Enhancements']['ParentCoreNumberFraction'], result[1]['Enhancements']['Ds'])
                except:
                    print('Could not tabulate output (likely some values were None)')
                print('FinalAbundances')
                print(result[0])
            else:
                data_dump[arg_name][ef_name] = None
        if we_care_about_actual_result:
            for p in generate_pressure_vals():
                data_dump[arg_name][(arg_set[7], p, arg_set[12])] = cm.complete_model_calculation(arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], arg_set[7], arg_set[8], arg_set[9], 10**(arg_set[10]), p, arg_set[12], 'NonEarthlike')[0]
            for fcf in generate_fcf_vals():
                data_dump[arg_name][(fcf, arg_set[11], arg_set[12])] = cm.complete_model_calculation(arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], fcf, arg_set[8], arg_set[9], 10**(arg_set[10]), arg_set[11], arg_set[12], 'NonEarthlike')[0]
                for p in generate_pressure_vals():
                    data_dump[arg_name][(fcf, p, arg_set[12])] = cm.complete_model_calculation(arg_set[1], arg_set[2], arg_set[3], arg_set[4], arg_set[5], arg_set[6], fcf, arg_set[8], arg_set[9], 10**(arg_set[10]), p, arg_set[12], 'NonEarthlike')[0]
    
    stellar_data = manager.stellar_compositions

    return data_dump, abundances_dict, errors_dict, stellar_data, abundance_lower_bounds_dict, abundance_upper_bounds_dict

def collect_model_output(system, all_observations, all_model_data, trial_fits):
    toret = dict()
    #fits = [(1, -3), (60, -3), (100, -3), (1, -1), (60, -1), (100, -1), (54, -2), 'Golden']
    fits = ['Golden', 'Sulf run'] + trial_fits.get(system, list())
    observations = all_observations[system]
    for f in fits:
        fit_output = dict()
        for element, value in observations.items():
            if value is np.nan:
                fit_output[element] = np.nan
            else:
                try:
                    fit_output[element] = all_model_data[system][f][element]
                except KeyError as e:
                    pass
        toret[str(f)] = fit_output
    return toret

def collect_model_output_mk2(system, wd_type, all_model_data, trial_fits, extended_fits, keys_as_tuple_hack=False, raw=False):
    toret = dict()
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
    #fits = [(1, -3), (60, -3), (100, -3), (1, -1), (60, -1), (100, -1), (54, -2), 'Golden']
    fits = ['LP run'] + trial_fits.get(system, list()) + [k for k in extended_fits.get(system, dict()).keys()]
    for f in fits:
        if all_model_data[system][f] is not None:
            output = list()
            for element in elements:
                if raw or element != ci.Element.Mg:
                    try:
                        if raw:
                            output.append(all_model_data[system][f][element] - (ld._geo_model.solar_ratiod_to_H[element] - ld._geo_model.solar_ratiod_to_H[wd_type]))
                        else:
                            output.append(all_model_data[system][f][element] - all_model_data[system][f][ci.Element.Mg] - ld._geo_model.solar_abundances[element])
                    except TypeError:
                        output.append(None)
            if keys_as_tuple_hack:
                toret[f] = output
            else:
                toret[str(f)] = output
    return toret
    
def collect_xfe_data(system, all_model_data, trial_fits):
    fcf_vals = list()
    crfe_v_fcf = list()
    nife_v_fcf = list()
    femg_v_fcf = list()
    p_vals = list()
    crfe_v_p = list()
    nife_v_p = list()
    femg_v_p = list()
    crfe_array = np.zeros((len(generate_fcf_vals()), len(generate_pressure_vals())))
    nife_array = np.zeros((len(generate_fcf_vals()), len(generate_pressure_vals())))
    f_index = 0
    for fcf in generate_fcf_vals():
        p_index = 0
        fcf_vals.append(fcf)
        fit = (fcf, fise_run_args[system][11], fise_run_args[system][12])
        crfe_v_fcf.append(all_model_data[system][fit][ci.Element.Cr] - all_model_data[system][fit][ci.Element.Fe])
        nife_v_fcf.append(all_model_data[system][fit][ci.Element.Ni] - all_model_data[system][fit][ci.Element.Fe])
        femg_v_fcf.append(all_model_data[system][fit][ci.Element.Fe] - all_model_data[system][fit][ci.Element.Mg])
        for p in generate_pressure_vals():
            system_run_data = all_model_data[system].get((fcf, p, fise_run_args[system][12]))
            if system_run_data is None:
                cr_value = None
                fe_value = None
                ni_value = None
            else:
                cr_value = system_run_data.get(ci.Element.Cr)
                fe_value = system_run_data.get(ci.Element.Fe)
                ni_value = system_run_data.get(ci.Element.Ni)
            if cr_value is None or fe_value is None:
                crfe_array[f_index][p_index] = None
            else:
                crfe_array[f_index][p_index] = cr_value - fe_value
            if ni_value is None or fe_value is None:
                nife_array[f_index][p_index] = None
            else:
                nife_array[f_index][p_index] = ni_value - fe_value
            p_index += 1
        f_index += 1
    for p in generate_pressure_vals():
        p_vals.append(p)
        fit = (fise_run_args[system][7], p, fise_run_args[system][12])
        fit_data = all_model_data[system].get(fit)
        if fit_data is None:
            cr_value = None
            fe_value = None
            ni_value = None
            mg_value = None
        else:
            cr_value = fit_data.get(ci.Element.Cr)
            fe_value = fit_data.get(ci.Element.Fe)
            ni_value = fit_data.get(ci.Element.Ni)
            mg_value = fit_data.get(ci.Element.Mg)
        if cr_value is None or fe_value is None:
            crfe_v_p.append(None)
        else:
            crfe_v_p.append(cr_value - fe_value)
        if ni_value is None or fe_value is None:
            nife_v_p.append(None)
        else:
            nife_v_p.append(ni_value - fe_value)
        if fe_value is None or mg_value is None:
            femg_v_p.append(None)
        else:
            femg_v_p.append(fe_value - mg_value)
    return fcf_vals, {ci.Element.Cr: crfe_v_fcf, ci.Element.Ni: nife_v_fcf}, femg_v_fcf, p_vals, {ci.Element.Cr: crfe_v_p, ci.Element.Ni: nife_v_p}, femg_v_p, {ci.Element.Cr: crfe_array, ci.Element.Ni: nife_array}

def collect_best_fit_ratios_from_stellar_data(stellar_data):
    stellar_crfe = list()
    stellar_nife = list()
    stellar_mgfe = list()
    stellar_crmg = list()
    stellar_sife = list()
    for star in stellar_data:
        ni_index = 3
        fe_index = 4
        cr_index = 5
        si_index = 6
        crfe = np.log10(star[cr_index]/star[fe_index])
        mgfe = np.log10(1/star[fe_index])
        nife = np.log10(star[ni_index]/star[fe_index])
        crmg = np.log10(star[cr_index])
        sife = np.log10(star[si_index]/star[fe_index])
        stellar_crfe.append(crfe)
        stellar_mgfe.append(mgfe)
        stellar_nife.append(nife)
        stellar_crmg.append(crmg)
        stellar_sife.append(sife)
    return stellar_crfe, stellar_nife, stellar_mgfe, stellar_crmg, stellar_sife

def collect_best_fit_ratios(all_model_data, abundances_dict, errors_dict, stellar_data, synthetic_stellar_data, synthetic_wd_data, abundance_lower_bounds_dict, abundance_upper_bounds_dict):
    systems = list()
    # TODO: Make this function less horrible!
    crfe_by_system_obs = list()
    mgfe_by_system_obs = list()
    nife_by_system_obs = list()
    crmg_by_system_obs = list()
    sife_by_system_obs = list()
    crfe_by_system_obs_err = list()
    mgfe_by_system_obs_err = list()
    nife_by_system_obs_err = list()
    crmg_by_system_obs_err = list()
    sife_by_system_obs_err = list()
    crfe_by_system_obs_lb = list()
    mgfe_by_system_obs_lb = list()
    nife_by_system_obs_lb = list()
    crmg_by_system_obs_lb = list()
    sife_by_system_obs_lb = list()
    crfe_by_system_obs_ub = list()
    mgfe_by_system_obs_ub = list()
    nife_by_system_obs_ub = list()
    crmg_by_system_obs_ub = list()
    sife_by_system_obs_ub = list()
    crfe_by_system = list()
    mgfe_by_system = list()
    nife_by_system = list()
    crmg_by_system = list()
    sife_by_system = list()
    system_ps = list()
    system_fcfs = list()
    text_offset_dict = {
        'crfemgfe': collections.OrderedDict({
            'SinkingEffects': (0, 0),  # Having the Sinking/Heating effects and core % here is a bit of a hack
            'HeatingEffects': (0, 0),
            '2 core': (0.25, -0.3),
            '75 core': (0.05, -0.15),
            '90 core': (0.05, -0.15),
            'Local Stars': (-0.25, -0.05),
            'CI chondrites': (-0.05, -0.23),
            'Earth (bulk)': (-0.15, 0.04),
            'Earth (mantle)': (0.05, 0),
            'Earth (crust)': (0, -0.2),
            'Mars (bulk)': (-0.1, 0.01),
            'Mars (silicate)': (0, 0),
            'HED chondrites': (-0.25, 0),
            'SDSSJ1234+5208': (0.3, 0.03),
            'SDSSJ1336+3547': (0.15, 0),
            'SDSSJ1430-0151': (0.43, -0.1),
            'GD61Corr': (-0.15, -0.07),
            'GD424': (0.1, 0),
            'WD1350-162NoNaCorr': (0.15, -0.25),
            'WD0446-255': (-0.2, 0.08),
            'NLTT43806Corr': (0.1, -0.16)
        }),
        'nifemgfe': collections.OrderedDict({
            'SinkingEffects': (0, 0),
            'HeatingEffects': (0, -0.2),
            '2 core': (-0.17, -0.2),
            '75 core': (0, 0),
            '90 core': (0, 0),
            'Local Stars': (-0.25, -0.05),
            'CI chondrites': (0.02, 0.14),
            'Earth (bulk)': (-0.3, 0.025),
            'Earth (mantle)': (0.1, -0.2),
            'Earth (crust)': (0, 0.03),
            'Mars (bulk)': (0, -0.22),
            'Mars (silicate)': (-0.25, 0.04),
            'WD0449-259NoNaCorr': (0, 0.05),
            'GD61Corr': (0.06, 0.05),
            'GD424': (0.07, 0.1),
            'SDSSJ1043+0855Corr': (0.2, -0.2),
            'SDSSJ2047-1259': (0.2, 0.05),
            'WD0446-255': (0, -0.4),
        }),
        'crfenife': collections.OrderedDict({
            'SinkingEffects': (0, -0.15),
            'HeatingEffects': (0, 0),
            '2 core': (0.8, 0),
            '75 core': (0.3, 0),
            '90 core': (0.3, -0.1),
            'Local Stars': (-0.25, -0.05),
            'CI chondrites': (0.55, -0.08),
            'Earth (bulk)': (-0.6, -0.08),
            'Earth (mantle)': (-0.45, -0.15),
            'Earth (crust)': (0.4, 0),
            'Earth (core)': (0.45, -0.15),
            'Mars (bulk)': (-0.3, 0),
            'Mars (silicate)': (0.2, 0)
        }),
        'crmgmgfe': collections.OrderedDict({
            'SinkingEffects': (0, 0),
            'HeatingEffects': (0, 0),
            '2 core': (0, 0),
            '75 core': (-0.15, -0.15),
            '90 core': (-0.1, -0.2),
            'Local Stars': (-0.25, -0.05),
            'CI chondrites': (-0.3, -0.2),
            'Earth (bulk)': (-0.45, -0.1),
            'Earth (mantle)': (0.02, -0.16),
            'Earth (crust)': (-0.25, -0.2),
            'Mars (bulk)': (-0.05, 0.03),
            'Mars (silicate)': (0.2, 0),
            'HED chondrites': (0.3, 0),
            'SDSSJ1234+5208': (0.2, 0.1),
            'SDSSJ1336+3547': (0.17, 0)
        }),
        'sifemgfe': collections.OrderedDict({
            'SinkingEffects': (0, -0.2),
            'HeatingEffects': (0.2, -0.25),
            '2 core': (0.25, -0.1),
            '75 core': (0.1, -0.2),
            '90 core': (0.1, -0.2),
            'Local Stars': (-0.25, -0.05),
            'CI chondrites': (-0.31, -0),
            'HED chondrites': (-0.3, 0),
            'Earth (bulk)': (-0.35, -0.05),
            'Earth (mantle)': (0.15, 0.02),
            'Earth (crust)': (-0.22, 0.05),
            'Mars (bulk)': (0, -0.2),
            'Mars (silicate)': (-0.15, 0.1),
            'WD1551+175': (0.33, -0.15),
            'WD1929+011': (-0.2, 0),
            'PG0843+516XCorr': (0, 0.2),
            'SDSSJ2047-1259': (-0.2, 0),
            'WD0449-259NoNaCorr': (0.1, 0.05),
            'WD1350-162NoNaCorr': (0.35, -0.25),
            'GD424': (0.15, 0),
            'GD61Corr': (0.05, -0.2),
            'SDSSJ1043+0855Corr': (0.2, -0.2),
            'WD0446-255': (0, 0.07),
            'WD2105-820': (-0.1, 0),
            'NLTT43806Corr': (0, -0.26)
        })
    }
    geo_model = gi.GeologyModel()
    #systems.append('Local Stars')
    #crfe_by_system.append(geo_model.solar_abundances[ci.Element.Cr] - geo_model.solar_abundances[ci.Element.Fe])
    #mgfe_by_system.append(geo_model.solar_abundances[ci.Element.Mg] - geo_model.solar_abundances[ci.Element.Fe])
    #nife_by_system.append(geo_model.solar_abundances[ci.Element.Ni] - geo_model.solar_abundances[ci.Element.Fe])
    #crfe_by_system_obs.append(None)
    #mgfe_by_system_obs.append(None)
    #nife_by_system_obs.append(None)
    #crfe_by_system_obs_err.append(None)
    #mgfe_by_system_obs_err.append(None)
    #nife_by_system_obs_err.append(None)
    systems.append('CI chondrites')
    crfe_by_system.append(5.65 - 7.47)  # From Lodders 2003
    mgfe_by_system.append(7.55 - 7.47)
    nife_by_system.append(6.22 - 7.47)
    crmg_by_system.append(5.65 - 7.55)
    sife_by_system.append(7.54 - 7.47)
    crfe_by_system_obs.append(None)
    mgfe_by_system_obs.append(None)
    nife_by_system_obs.append(None)
    crmg_by_system_obs.append(None)
    sife_by_system_obs.append(None)
    crfe_by_system_obs_err.append(None)
    mgfe_by_system_obs_err.append(None)
    nife_by_system_obs_err.append(None)
    crmg_by_system_obs_err.append(None)
    sife_by_system_obs_err.append(None)
    crfe_by_system_obs_lb.append(None)
    mgfe_by_system_obs_lb.append(None)
    nife_by_system_obs_lb.append(None)
    crmg_by_system_obs_lb.append(None)
    sife_by_system_obs_lb.append(None)
    crfe_by_system_obs_ub.append(None)
    mgfe_by_system_obs_ub.append(None)
    nife_by_system_obs_ub.append(None)
    crmg_by_system_obs_ub.append(None)
    sife_by_system_obs_ub.append(None)
    #systems.append('HED chondrites') # Thought to be a proxy for Vesta (at least the surface)
    #crfe_by_system.append(np.log10(0.02937150561))  # From my spreadsheet analysing a table included with Kimura 2018's NIPR data
    #mgfe_by_system.append(np.log10(1.37867479))
    #nife_by_system.append(None)  # Ni is included in the raw data but for now I skipped it (it would be non-trivial to calculate for Ni because along with Co it's the only one given in ppm rather than wt %)
    #crmg_by_system.append(np.log10(0.02937150561/1.37867479))  # = (cr/fe)/(mg/fe)
    #sife_by_system.append(np.log10(3.523104002))  # Looks like there's 2 x HED. One with Si/Mg ~ 1 and one with Si/Mg >> 1
    #crfe_by_system_obs.append(None)
    #mgfe_by_system_obs.append(None)
    #nife_by_system_obs.append(None)
    #crmg_by_system_obs.append(None)
    #sife_by_system_obs.append(None)
    #crfe_by_system_obs_err.append(None)
    #mgfe_by_system_obs_err.append(None)
    #nife_by_system_obs_err.append(None)
    #crmg_by_system_obs_err.append(None)
    #sife_by_system_obs_err.append(None)
    #crfe_by_system_obs_lb.append(None)
    #mgfe_by_system_obs_lb.append(None)
    #nife_by_system_obs_lb.append(None)
    #crmg_by_system_obs_lb.append(None)
    #sife_by_system_obs_lb.append(None)
    #crfe_by_system_obs_ub.append(None)
    #mgfe_by_system_obs_ub.append(None)
    #nife_by_system_obs_ub.append(None)
    #crmg_by_system_obs_ub.append(None)
    #sife_by_system_obs_ub.append(None)
    systems.append('Mars (silicate)')
    crfe_by_system.append(np.log10((6000/10000)/4.4))  # From yoshizaki 2020
    mgfe_by_system.append(np.log10(17/4.4))
    nife_by_system.append(np.log10(0.01/4.4))
    crmg_by_system.append(np.log10((6000/10000)/17))
    sife_by_system.append(np.log10(16/4.4))
    crfe_by_system_obs.append(None)
    mgfe_by_system_obs.append(None)
    nife_by_system_obs.append(None)
    crmg_by_system_obs.append(None)
    sife_by_system_obs.append(None)
    crfe_by_system_obs_err.append(None)
    mgfe_by_system_obs_err.append(None)
    nife_by_system_obs_err.append(None)
    crmg_by_system_obs_err.append(None)
    sife_by_system_obs_err.append(None)
    crfe_by_system_obs_lb.append(None)
    mgfe_by_system_obs_lb.append(None)
    nife_by_system_obs_lb.append(None)
    crmg_by_system_obs_lb.append(None)
    sife_by_system_obs_lb.append(None)
    crfe_by_system_obs_ub.append(None)
    mgfe_by_system_obs_ub.append(None)
    nife_by_system_obs_ub.append(None)
    crmg_by_system_obs_ub.append(None)
    sife_by_system_obs_ub.append(None)
    systems.append('Mars (bulk)')
    crfe_by_system.append(np.log10((4900/10000)/10))  # From yoshizaki 2020
    mgfe_by_system.append(np.log10(0.15/0.1))
    nife_by_system.append(np.log10(0.5/0.1))
    crmg_by_system.append(np.log10((4900/10000)/15))
    sife_by_system.append(14/10)
    crfe_by_system_obs.append(None)
    mgfe_by_system_obs.append(None)
    nife_by_system_obs.append(None)
    crmg_by_system_obs.append(None)
    sife_by_system_obs.append(None)
    crfe_by_system_obs_err.append(None)
    mgfe_by_system_obs_err.append(None)
    nife_by_system_obs_err.append(None)
    crmg_by_system_obs_err.append(None)
    sife_by_system_obs_err.append(None)
    crfe_by_system_obs_lb.append(None)
    mgfe_by_system_obs_lb.append(None)
    nife_by_system_obs_lb.append(None)
    crmg_by_system_obs_lb.append(None)
    sife_by_system_obs_lb.append(None)
    crfe_by_system_obs_ub.append(None)
    mgfe_by_system_obs_ub.append(None)
    nife_by_system_obs_ub.append(None)
    crmg_by_system_obs_ub.append(None)
    sife_by_system_obs_ub.append(None)
    systems.append('Mars (core)')
    crfe_by_system.append(None)  # From yoshizaki 2020
    mgfe_by_system.append(None)
    nife_by_system.append(np.log10(4.2/48))
    crmg_by_system.append(None)
    sife_by_system.append(None)
    crfe_by_system_obs.append(None)
    mgfe_by_system_obs.append(None)
    nife_by_system_obs.append(None)
    crmg_by_system_obs.append(None)
    sife_by_system_obs.append(None)
    crfe_by_system_obs_err.append(None)
    mgfe_by_system_obs_err.append(None)
    nife_by_system_obs_err.append(None)
    crmg_by_system_obs_err.append(None)
    sife_by_system_obs_err.append(None)
    crfe_by_system_obs_lb.append(None)
    mgfe_by_system_obs_lb.append(None)
    nife_by_system_obs_lb.append(None)
    crmg_by_system_obs_lb.append(None)
    sife_by_system_obs_lb.append(None)
    crfe_by_system_obs_ub.append(None)
    mgfe_by_system_obs_ub.append(None)
    nife_by_system_obs_ub.append(None)
    crmg_by_system_obs_ub.append(None)
    sife_by_system_obs_ub.append(None)
    for layer in [gi.Layer.bulk, gi.Layer.mantle, gi.Layer.core]:
        systems.append('Earth (' + str(layer) + ')')
        crfe_by_system.append(np.log10(geo_model.element_info[ci.Element.Cr][layer]/geo_model.element_info[ci.Element.Fe][layer]))
        mgfe_by_system.append(np.log10(geo_model.element_info[ci.Element.Mg][layer]/geo_model.element_info[ci.Element.Fe][layer]))
        nife_by_system.append(np.log10(geo_model.element_info[ci.Element.Ni][layer]/geo_model.element_info[ci.Element.Fe][layer]))
        try:
            crmg_by_system.append(np.log10(geo_model.element_info[ci.Element.Cr][layer]/geo_model.element_info[ci.Element.Mg][layer]))
        except ZeroDivisionError:
            # Occurs for core (Mg is 0)
            crmg_by_system.append(None)
        sife_by_system.append(np.log10(geo_model.element_info[ci.Element.Si][layer]/geo_model.element_info[ci.Element.Fe][layer]))
        crfe_by_system_obs.append(None)
        mgfe_by_system_obs.append(None)
        nife_by_system_obs.append(None)
        crmg_by_system_obs.append(None)
        sife_by_system_obs.append(None)
        crfe_by_system_obs_err.append(None)
        mgfe_by_system_obs_err.append(None)
        nife_by_system_obs_err.append(None)
        crmg_by_system_obs_err.append(None)
        sife_by_system_obs_err.append(None)
        crfe_by_system_obs_lb.append(None)
        mgfe_by_system_obs_lb.append(None)
        nife_by_system_obs_lb.append(None)
        crmg_by_system_obs_lb.append(None)
        sife_by_system_obs_lb.append(None)
        crfe_by_system_obs_ub.append(None)
        mgfe_by_system_obs_ub.append(None)
        nife_by_system_obs_ub.append(None)
        crmg_by_system_obs_ub.append(None)
        sife_by_system_obs_ub.append(None)
    
    for arg_name, arg_set in fise_run_args.items():
        if limit_to_sample and arg_name not in systems_in_sample:
            continue
        try:
            el_dict = all_model_data[arg_name][('LP run')]
        except KeyError:
            el_dict = dict()
        systems.append(arg_name)
        obs_Ni = abundances_dict[arg_name][ci.Element.Ni]
        obs_Fe = abundances_dict[arg_name][ci.Element.Fe]
        obs_Cr = abundances_dict[arg_name][ci.Element.Cr]
        obs_Mg = abundances_dict[arg_name][ci.Element.Mg]
        obs_Si = abundances_dict[arg_name][ci.Element.Si]
        obs_err_Ni = errors_dict[arg_name][ci.Element.Ni]
        obs_err_Fe = errors_dict[arg_name][ci.Element.Fe]
        obs_err_Cr = errors_dict[arg_name][ci.Element.Cr]
        obs_err_Mg = errors_dict[arg_name][ci.Element.Mg]
        obs_err_Si = errors_dict[arg_name][ci.Element.Si]
        obs_lb_Ni = abundance_lower_bounds_dict[arg_name][ci.Element.Ni]
        obs_lb_Fe = abundance_lower_bounds_dict[arg_name][ci.Element.Fe]
        obs_lb_Cr = abundance_lower_bounds_dict[arg_name][ci.Element.Cr]
        obs_lb_Mg = abundance_lower_bounds_dict[arg_name][ci.Element.Mg]
        obs_lb_Si = abundance_lower_bounds_dict[arg_name][ci.Element.Si]
        obs_ub_Ni = abundance_upper_bounds_dict[arg_name][ci.Element.Ni]
        obs_ub_Fe = abundance_upper_bounds_dict[arg_name][ci.Element.Fe]
        obs_ub_Cr = abundance_upper_bounds_dict[arg_name][ci.Element.Cr]
        obs_ub_Mg = abundance_upper_bounds_dict[arg_name][ci.Element.Mg]
        obs_ub_Si = abundance_upper_bounds_dict[arg_name][ci.Element.Si]
        crfe_by_system.append(None if (np.isnan(obs_Cr) or np.isnan(obs_Fe) or el_dict is None) else el_dict.get(ci.Element.Cr) - el_dict.get(ci.Element.Fe))
        mgfe_by_system.append(None if (np.isnan(obs_Mg) or np.isnan(obs_Fe) or el_dict is None) else el_dict.get(ci.Element.Mg) - el_dict.get(ci.Element.Fe))
        nife_by_system.append(None if (np.isnan(obs_Ni) or np.isnan(obs_Fe) or el_dict is None) else el_dict.get(ci.Element.Ni) - el_dict.get(ci.Element.Fe))
        crmg_by_system.append(None if (np.isnan(obs_Cr) or np.isnan(obs_Mg) or el_dict is None) else el_dict.get(ci.Element.Cr) - el_dict.get(ci.Element.Mg))
        sife_by_system.append(None if (np.isnan(obs_Si) or np.isnan(obs_Fe) or el_dict is None) else el_dict.get(ci.Element.Si) - el_dict.get(ci.Element.Fe))
        crfe_by_system_obs.append(None if (np.isnan(obs_Cr) or np.isnan(obs_Fe)) else obs_Cr - obs_Fe)
        mgfe_by_system_obs.append(None if (np.isnan(obs_Mg) or np.isnan(obs_Fe)) else obs_Mg - obs_Fe)
        nife_by_system_obs.append(None if (np.isnan(obs_Ni) or np.isnan(obs_Fe)) else obs_Ni - obs_Fe)
        crmg_by_system_obs.append(None if (np.isnan(obs_Cr) or np.isnan(obs_Mg)) else obs_Cr - obs_Mg)
        sife_by_system_obs.append(None if (np.isnan(obs_Si) or np.isnan(obs_Fe)) else obs_Si - obs_Fe)
        crfe_by_system_obs_err.append(None if (np.isnan(obs_Cr) or np.isnan(obs_Fe)) else np.sqrt(obs_err_Cr**2 + obs_err_Fe**2))
        mgfe_by_system_obs_err.append(None if (np.isnan(obs_Mg) or np.isnan(obs_Fe)) else np.sqrt(obs_err_Mg**2 + obs_err_Fe**2))
        nife_by_system_obs_err.append(None if (np.isnan(obs_Ni) or np.isnan(obs_Fe)) else np.sqrt(obs_err_Ni**2 + obs_err_Fe**2))
        crmg_by_system_obs_err.append(None if (np.isnan(obs_Cr) or np.isnan(obs_Mg)) else np.sqrt(obs_err_Cr**2 + obs_err_Mg**2))
        sife_by_system_obs_err.append(None if (np.isnan(obs_Si) or np.isnan(obs_Fe)) else np.sqrt(obs_err_Si**2 + obs_err_Fe**2))
    
        crfe_by_system_obs_lb.append(None if (obs_lb_Cr is None or np.isnan(obs_Fe)) else obs_lb_Cr - obs_Fe)
        mgfe_by_system_obs_lb.append(None if (obs_lb_Mg is None or np.isnan(obs_Fe)) else obs_lb_Mg - obs_Fe)
        nife_by_system_obs_lb.append(None if (obs_lb_Ni is None or np.isnan(obs_Fe)) else obs_lb_Ni - obs_Fe)
        crmg_by_system_obs_lb.append(None if (obs_lb_Cr is None or np.isnan(obs_Mg)) else obs_lb_Cr - obs_Mg)
        sife_by_system_obs_lb.append(None if (obs_lb_Si is None or np.isnan(obs_Fe)) else obs_lb_Si - obs_Fe)
    
        crfe_by_system_obs_ub.append(None if (obs_ub_Cr is None or np.isnan(obs_Fe)) else obs_ub_Cr - obs_Fe)
        mgfe_by_system_obs_ub.append(None if (obs_ub_Mg is None or np.isnan(obs_Fe)) else obs_ub_Mg - obs_Fe)
        nife_by_system_obs_ub.append(None if (obs_ub_Ni is None or np.isnan(obs_Fe)) else obs_ub_Ni - obs_Fe)
        crmg_by_system_obs_ub.append(None if (obs_ub_Cr is None or np.isnan(obs_Mg)) else obs_ub_Cr - obs_Mg)
        sife_by_system_obs_ub.append(None if (obs_ub_Si is None or np.isnan(obs_Fe)) else obs_ub_Si - obs_Fe)
        
        system_ps.append(None)  # This used to be taken from arg_set but I'm disabling because the arg_set isn't up to date!
        system_fcfs.append(None)  # This used to be taken from arg_set but I'm disabling because the arg_set isn't up to date!
    
    stellar_crfe, stellar_nife, stellar_mgfe, stellar_crmg, stellar_sife = collect_best_fit_ratios_from_stellar_data(stellar_data)
    synth_stellar_crfe, synth_stellar_nife, synth_stellar_mgfe, synth_stellar_crmg, synth_stellar_sife = collect_best_fit_ratios_from_stellar_data(synthetic_stellar_data)
    synth_wd_crfe_dict = dict()
    synth_wd_nife_dict = dict()
    synth_wd_mgfe_dict = dict()
    synth_wd_crmg_dict = dict()
    synth_wd_sife_dict = dict()
    for wd_name, wd_data in synthetic_wd_data.items():
        synth_wd_crfe, synth_wd_nife, synth_wd_mgfe, synth_wd_crmg, synth_wd_sife = collect_best_fit_ratios_from_stellar_data(wd_data)
        synth_wd_crfe_dict[wd_name] = synth_wd_crfe
        synth_wd_nife_dict[wd_name] = synth_wd_nife
        synth_wd_mgfe_dict[wd_name] = synth_wd_mgfe
        synth_wd_crmg_dict[wd_name] = synth_wd_crmg
        synth_wd_sife_dict[wd_name] = synth_wd_sife
    system_obs = {
        'crfe': crfe_by_system_obs,
        'mgfe': mgfe_by_system_obs,
        'nife': nife_by_system_obs,
        'crmg': crmg_by_system_obs,
        'sife': sife_by_system_obs
    }
    system_obs_err = {
        'crfe': crfe_by_system_obs_err,
        'mgfe': mgfe_by_system_obs_err,
        'nife': nife_by_system_obs_err,
        'crmg': crmg_by_system_obs_err,
        'sife': sife_by_system_obs_err
    }
    system_obs_lb = {
        'crfe': crfe_by_system_obs_lb,
        'mgfe': mgfe_by_system_obs_lb,
        'nife': nife_by_system_obs_lb,
        'crmg': crmg_by_system_obs_lb,
        'sife': sife_by_system_obs_lb
    }
    system_obs_ub = {
        'crfe': crfe_by_system_obs_ub,
        'mgfe': mgfe_by_system_obs_ub,
        'nife': nife_by_system_obs_ub,
        'crmg': crmg_by_system_obs_ub,
        'sife': sife_by_system_obs_ub
    }
    elel_by_system = {
        'crfe': crfe_by_system,
        'mgfe': mgfe_by_system,
        'nife': nife_by_system,
        'crmg': crmg_by_system,
        'sife': sife_by_system
    }
    stellar_dict = {
        'crfe': stellar_crfe,
        'mgfe': stellar_mgfe,
        'nife': stellar_nife,
        'crmg': stellar_crmg,
        'sife': stellar_sife
    }
    synth_stellar_dict = {
        'crfe': synth_stellar_crfe,
        'mgfe': synth_stellar_mgfe,
        'nife': synth_stellar_nife,
        'crmg': synth_stellar_crmg,
        'sife': synth_stellar_sife
    }
    synth_wd_dict = {
        'crfe': synth_wd_crfe_dict,
        'mgfe': synth_wd_mgfe_dict,
        'nife': synth_wd_nife_dict,
        'crmg': synth_wd_crmg_dict,
        'sife': synth_wd_sife_dict
    }
    systems_to_show_ellipse_dict = {
        'crfemgfe': [
            'PG0843+516XCorr',
            #'WD1551+175',
            #'SDSSJ1535+1247',
            #'SDSSJ1430-0151',
            #'SDSSJ0150+1354',
            #'SDSSJ0512-0505',
            #'SDSSJ0736+4118',
            #'SDSSJ0956+5912',
            #'SDSSJ1144+1218',
            #'SDSSJ0252-0401',
            #'SDSSJ1055+3725',
            #'SDSSJ1627+4646',
            #'SDSSJ0046+2717',
            #'SDSSJ2047-1259',
            #'WD1551+175',
            #'WD1929+011',
            #'SDSSJ0116+2050',
            #'SDSSJ0741+3146',
            #'SDSSJ0806+3055',
            #'SDSSJ0901+0752',
            #'SDSSJ0916+2540',
            #'SDSSJ0929+4247',
            #'SDSSJ0939+4136',
            #'SDSSJ0939+5019',
            #'SDSSJ1040+2407',
            #'SDSSJ1043+3516',
            #'SDSSJ1055+3725',
            #'SDSSJ1229+0743',
            #'SDSSJ1234+5208',
            #'SDSSJ1321-0237',
            #'SDSSJ1336+3547',
            #'SDSSJ1345+1153',
            #'SDSSJ1411+3410',
            #'SDSSJ1524+4049',
            #'SDSSJ2230+1905',
            #'PG0843+516Gaensicke',
            #'PG0843+516Xu',
            #'PG1015+161',
            #'SDSS1228+1040',
            #'GALEX1931+0117',
            #'GD362',
            #'GD40',
            #'G241-6',
            'GD61Corr',
            'NLTT43806Corr',
            #'HS2253+8023',
            #'PG1225-079',
            #'SDSSJ0738+1835',
            #'SDSSJ1228+1040',
            #'G29-38',
            #'SDSSJ1242+5226',
            #'SDSSJ0845+2257',
            #'WD1536+520',
            #'SDSSJ1043+0855',
            #'WD1425+540',
            'WD0446-255',
            #'WD2216-657',
            #'WD2157-574',
            #'WD2115-560',
            #'WD1232+563',
            #'WD2207+121',
            #'WD1145+017',
            #'WD0122-227',
            #'WD0449-259',
            'WD1350-162NoNaCorr',
            'WD2105-820',
            #'SDSSJ1043+0855Corr',
            #'GD424'
        ],
        'nifemgfe': [
            #'PG0843+516',
            #'SDSSJ2047-1259'
            #'PG0843+516Gaensicke',
            #'PG0843+516Xu',
            #'PG1015+161',
            #'SDSS1228+1040',
            #'GALEX1931+0117',
            #'GD362',
            #'GD40',
            #'G241-6',
            'GD61Corr',
            'NLTT43806Corr',
            #'HS2253+8023',
            #'PG1225-079',
            #'SDSSJ0738+1835',
            #'SDSSJ1228+1040',
            #'G29-38',
            #'SDSSJ1242+5226',
            #'SDSSJ0845+2257',
            #'WD1536+520',
            #'SDSSJ1043+0855',
            #'WD1425+540',
            'WD0446-255',
            #'WD2216-657',
            #'WD2157-574',
            #'WD2115-560',
            #'WD1232+563',
            #'WD2207+121',
            #'WD1145+017',
            #'WD0122-227',
            #'WD0449-259',
            #'WD1350-162',
            'WD2105-820',
            'WD0449-259NoNaCorr',
            'WD1350-162NoNaCorr',
            #'SDSSJ1043+0855Corr',
            #'GD424'
        ],
        'crfenife': [
            'PG0843+516'
            'PG0843+516Gaensicke',
            'PG0843+516Xu',
            'PG1015+161',
            'SDSS1228+1040',
            'GALEX1931+0117',
            'GD362',
            'GD40',
            'G241-6',
            'GD61',
            'NLTT43806',
            'HS2253+8023',
            'PG1225-079',
            'SDSSJ0738+1835',
            'SDSSJ1228+1040',
            'G29-38',
            'SDSSJ1242+5226',
            'SDSSJ0845+2257',
            'WD1536+520',
            'SDSSJ1043+0855',
            'WD1425+540',
            'WD0446-255',
            'WD2216-657',
            'WD2157-574',
            'WD2115-560',
            'WD1232+563',
            'WD2207+121',
            'WD1145+017',
            'WD0122-227',
            'WD0449-259',
            'WD1350-162',
            'WD2105-820'
        ],
        'crmgmgfe': [
            'PG0843+516',
            #'WD1551+175',
            #'SDSSJ1535+1247',
            #'SDSSJ1430-0151',
            #'SDSSJ0150+1354',
            #'SDSSJ0512-0505',
            #'SDSSJ0736+4118',
            #'SDSSJ0956+5912',
            #'SDSSJ1144+1218',
            #'SDSSJ0252-0401',
            #'SDSSJ1055+3725',
            #'SDSSJ1627+4646',
            #'SDSSJ0046+2717',
            #'SDSSJ2047-1259',
            #'WD1551+175',
            #'WD1929+011',
            #'SDSSJ0116+2050',
            #'SDSSJ0741+3146',
            #'SDSSJ0806+3055',
            #'SDSSJ0901+0752',
            #'SDSSJ0916+2540',
            #'SDSSJ0929+4247',
            #'SDSSJ0939+4136',
            #'SDSSJ0939+5019',
            #'SDSSJ1040+2407',
            #'SDSSJ1043+3516',
            #'SDSSJ1055+3725',
            #'SDSSJ1229+0743',
            #'SDSSJ1234+5208',
            #'SDSSJ1321-0237',
            'SDSSJ1336+3547',
            #'SDSSJ1345+1153',
            #'SDSSJ1411+3410',
            #'SDSSJ1524+4049',
            #'SDSSJ2230+1905'
            'PG0843+516Gaensicke',
            'PG0843+516Xu',
            'PG1015+161',
            'SDSS1228+1040',
            'GALEX1931+0117',
            'GD362',
            'GD40',
            'G241-6',
            'GD61',
            'NLTT43806',
            'HS2253+8023',
            'PG1225-079',
            'SDSSJ0738+1835',
            'SDSSJ1228+1040',
            'G29-38',
            'SDSSJ1242+5226',
            'SDSSJ0845+2257',
            'WD1536+520',
            'SDSSJ1043+0855',
            'WD1425+540',
            'WD0446-255',
            'WD2216-657',
            'WD2157-574',
            'WD2115-560',
            'WD1232+563',
            'WD2207+121',
            'WD1145+017',
            'WD0122-227',
            'WD0449-259',
            'WD1350-162',
            'WD2105-820'
        ],
        'sifemgfe': [
            #'PG0843+516',
            #'WD1551+175',
            #'SDSSJ1535+1247',
            #'SDSSJ1430-0151',
            #'SDSSJ0150+1354',
            #'SDSSJ0512-0505',
            #'SDSSJ0736+4118',
            #'SDSSJ0956+5912',
            #'SDSSJ1144+1218',
            #'SDSSJ0252-0401',
            #'SDSSJ1055+3725',
            #'SDSSJ1627+4646',
            #'SDSSJ0046+2717',
            #'SDSSJ2047-1259',
            #'WD1551+175',
            #'WD1929+011',
            #'SDSSJ0116+2050',
            #'SDSSJ0741+3146',
            #'SDSSJ0806+3055',
            #'SDSSJ0901+0752',
            #'SDSSJ0916+2540',
            #'SDSSJ0929+4247',
            #'SDSSJ0939+4136',
            #'SDSSJ0939+5019',
            #'SDSSJ1040+2407',
            'SDSSJ1043+3516Corr',
            #'SDSSJ1055+3725',
            #'SDSSJ1229+0743',
            #'SDSSJ1234+5208',
            #'SDSSJ1321-0237',
            #'SDSSJ1336+3547',
            #'SDSSJ1345+1153',
            #'SDSSJ1411+3410',
            #'SDSSJ1524+4049',
            #'SDSSJ2230+1905',
            #'PG0843+516Gaensicke',
            'PG0843+516XCorr',
            #'PG1015+161',
            #'SDSS1228+1040',
            #'GALEX1931+0117',
            #'GD362',
            #'GD40',
            #'G241-6',
            'GD61Corr',
            'NLTT43806Corr', #NEW
            #'HS2253+8023',
            #'PG1225-079',
            #'SDSSJ0738+1835',
            #'SDSSJ1228+1040',
            #'G29-38',
            #'SDSSJ1242+5226',
            #'SDSSJ0845+2257',
            #'WD1536+520',
            #'SDSSJ1043+0855',
            #'WD1425+540',
            'WD0446-255',
            #'WD2216-657',
            #'WD2157-574',
            #'WD2115-560',
            #'WD1232+563',
            #'WD2207+121',
            #'WD1145+017',
            #'WD0122-227',
            'WD0449-259NoNaCorr',
            'WD1350-162NoNaCorr',
            'WD2105-820', #NEW
            #'SDSSJ1043+0855Corr',
            #'GD424'
        ]
    }
    return systems, 7, system_obs, system_obs_err, elel_by_system, system_ps, system_fcfs, text_offset_dict, stellar_dict, synth_stellar_dict, synth_wd_dict, systems_to_show_ellipse_dict, system_obs_lb, system_obs_ub

def calculate_sinking_vectors(y_el_1, y_el_2, x_el_1, x_el_2):
    # where the y axis is log(y_el_1/y_el_2)
    # and the x axis is log(x_el_1/x_el_2)
    # eg for CrFe_MgFe we have y_el_1 = Cr, y_el_2 = Fe, x_el_1 = Mg, x_el_2 = Fe
    el_indices = {
        ci.Element.Al: 0,
        ci.Element.Ti: 1,
        ci.Element.Ca: 2,
        ci.Element.Ni: 3,
        ci.Element.Fe: 4,
        ci.Element.Cr: 5,
        ci.Element.Mg: 6,
        ci.Element.Si: 7,
        ci.Element.Na: 8,
        ci.Element.O: 9,
        ci.Element.C: 10,
        ci.Element.N: 11
    }
    wd_timescales = np.array([19194795.09,10663368.1,14083479.59,9241341.812,9430702.751,9965597.748,22230674.45,19801266.37,21802996.49,33976746.42,44367224.76,38432208.65])  # In order of ci.usual_elements. Just an average of all the real timescale
    dummy_planetesimal_abundance = np.array([0.077449061,0.002751825,0.0536328,0.03313212,0.489834545,0.006169951,1,1.07237162,0.02754622,45.67899096,9.310358117,1.000270853])  # Just the first star in the sample. This shouldn't matter
    t_disc = 50000000  # 50 Myr, should ensure we reach a steady state
    
    prev_coords = None
    time = 0
    timestep = 0.05
    max_time = 55  # Myr
    arrow_list = list()
    while time < max_time:
        abundance = wdm.process_abundances(time, t_disc, dummy_planetesimal_abundance, dummy_planetesimal_abundance, wd_timescales, wd_timescales, 0)
        coords = [abundance[el_indices[x_el_1]] - abundance[el_indices[x_el_2]], abundance[el_indices[y_el_1]] - abundance[el_indices[y_el_2]]]
        if prev_coords is None:
            prev_coords = coords
        else:
            arrow = prev_coords + [coords[0] - prev_coords[0]] + [coords[1] - prev_coords[1]]
            prev_coords = coords
            arrow_list.append(arrow)
        time += timestep
    el_key = str(y_el_1) + str(y_el_2) + str(x_el_1) + str(x_el_2)
    offset_dict = {
        'CrFeMgFe': [0.55, -0.85],
        'CrFeNiFe': [1.5, -0.8],
        'SiFeMgFe': [0.3, -1.5],
        'NiFeMgFe': [0.5, 1.1],
        'CrMgMgFe': [0.6, 1.5],
    }
    offset = offset_dict.get(el_key, [0, 0])
    for arrow in arrow_list:
        arrow[0] += offset[0]
        arrow[1] += offset[1]
    return arrow_list

def calculate_heating_vectors(y_el_1, y_el_2, x_el_1, x_el_2):
    # where the y axis is log(y_el_1/y_el_2)
    # and the x axis is log(x_el_1/x_el_2)
    # eg for CrFe_MgFe we have y_el_1 = Cr, y_el_2 = Fe, x_el_1 = Mg, x_el_2 = Fe
    
    fe_star = 428  # Pick the median value (only matters for oxygen anyway)
    z_formation = 0  # This might be a bad choice? Also 0.05 is the default! But then we do want to show the effect of heating, not the feeding zone
    t_formation = 1.5        
        
    prev_coords = None
    dist = 1
    diststep = 0.001
    min_dist = 0.238  # AU
    arrow_list = list()
    while dist > min_dist:
        abundances = am.get_all_abundances(ci.usual_elements, dist, z_formation, t_formation, fe_star)
        coords = [np.log10(abundances[x_el_1]/abundances[x_el_2]), np.log10(abundances[y_el_1]/abundances[y_el_2])]
        if prev_coords is None:
            prev_coords = coords
        else:
            arrow = prev_coords + [coords[0] - prev_coords[0]] + [coords[1] - prev_coords[1]]
            prev_coords = coords
            arrow_list.append(arrow)
        dist -= diststep
    el_key = str(y_el_1) + str(y_el_2) + str(x_el_1) + str(x_el_2)
    offset_dict = {
        'CrFeMgFe': [-1, -0.5],
        'CrFeNiFe': [0.1, -1.7],
        'SiFeMgFe': [-1.1, 0.7],
        'NiFeMgFe': [-0.85, 0.2],
        'CrMgMgFe': [-1.55, -1.7],
    }
    offset = offset_dict.get(el_key, [0, 0])
    for arrow in arrow_list:
        arrow[0] += offset[0]
        arrow[1] += offset[1]
    return arrow_list


def main():
    list_of_video_ps = list()
    video_p = 0
    while video_p < 60.01:
        list_of_video_ps.append(video_p)
        video_p += 1
    # Next line is to avoid rounding errors
    # list_of_video_ps = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19, 19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20, 20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21, 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9, 22, 22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23, 23.1, 23.2, 23.3, 23.4, 23.5, 23.6, 23.7, 23.8, 23.9, 24, 24.1, 24.2, 24.3, 24.4, 24.5, 24.6, 24.7, 24.8, 24.9, 25, 25.1, 25.2, 25.3, 25.4, 25.5, 25.6, 25.7, 25.8, 25.9, 26, 26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9, 27, 27.1, 27.2, 27.3, 27.4, 27.5, 27.6, 27.7, 27.8, 27.9, 28, 28.1, 28.2, 28.3, 28.4, 28.5, 28.6, 28.7, 28.8, 28.9, 29, 29.1, 29.2, 29.3, 29.4, 29.5, 29.6, 29.7, 29.8, 29.9, 30]
    trial_fits = {
        'PG0843+516': [(0.616654946, 0, -2.744307251), (0.616654946, 30, -2.744307251), (0.616654946, 60, -2.744307251), (0.616654946, 100, -2.744307251)],#, (0.9, 400, -4)],
        'WD1551+175': [(0.050019429804295315, p, -1.6937455473588918) for p in [1,45,60]],
        'SDSSJ1535+1247': [(0.072592648, p, -2.351406914) for p in [5,25,50]],
        #'SDSSJ1430-0151': [],
        #'SDSSJ0150+1354': [],#[(0.090547456, 0, -1), (0.090547456, 0, -3), (0.02, 0, -3), (0.06, 0, -3)],
        #'SDSSJ0512-0505': [],
        #'SDSSJ0736+4118': [],
        #'SDSSJ1144+1218': [(0.000067514, p, -1.642344232) for p in list_of_video_ps],
        #'SDSSJ1055+3725': [],#[(0.422469808, 400, -3), (0.422469808, 400, -1), (0.5, 400, -3), (0.6, 400, -3)],
        #'SDSSJ1627+4646': [],
        #'SDSSJ0046+2717': []
        'Earthfcf0': [(0, p, -1.5) for p in [1, 5, 10, 60, 100]], #0.020520230211611815   #-2.3096788622895144
        'GD424': [(0.114364826808444, p, -2.36726679401981) for p in [1, 60]] + [(0.114364826808444, 25, -2.8)] + [(0.05, 1, -2.36726679401981)] + [(0.114364826808444, 54, -2)], # list(range(1, 61))
        'GD61': [(0.026231629101023, p, -2.47677571461788) for p in [0, 60]] + [(0.026231629101023, 42.9965877524006, f) for f in [-1, -3]],
        'PG0843+516Gaensicke': [(0.639718848123253, p, -2.56976334306055) for p in [0, 60]] + [(0.639718848123253, 27.0819223781586, f) for f in [-1, -3]],
        'PG0843+516Xu': [(0.728745771883592, p, -1.76765362383884) for p in [0, 60]] + [(0.728745771883592, 23.0904826689188, f) for f in [-1, -3]],
        'NLTT43806': [(0.027875631165411, p, -2.67509568637878) for p in [0, 60]] + [(0.027875631165411, 46.2001486427106, f) for f in [-1, -3]],
        'SDSSJ0738+1835': [(0.43722029937676, p, -2.54048938334322) for p in [0, 60]] + [(0.43722029937676, 31.2042441423048, f) for f in [-1, -3]],
        'SDSSJ0845+2257': [(0.331040611889366, p, -2.01081246790329) for p in [0, 60]] + [(0.331040611889366, 29.5502718185378, f) for f in [-1, -3]],
        'GD 378': [(0.054500320061958, p, -2.2182777080587) for p in [0, 60]] + [(0.054500320061958, 37.0105226526938, f) for f in [-1, -3]],
        
        'WD0449-259': [(0.66877586362797, p, -1.36618477730749) for p in [0, 60]] + [(0.66877586362797, 19.157755678069, f) for f in [-1, -3]],
        'WD1350-162': [(0.39945856964945, p, -1.95588390243148) for p in [0, 60]] + [(0.39945856964945, 20.0942317714276, f) for f in [-1, -3]],
        'WD0122-227': [(0.567305307981525, p, -1.99238401733305) for p in [0, 60]] + [(0.567305307981525, 29.4379454699773, f) for f in [-1, -3]],
        'SDSSJ0512-0505': [(0.545002949374465, p, -1.77609214975635) for p in [0, 60]] + [(0.545002949374465, 25.2549663467436, f) for f in [-1, -3]],
        'SDSSJ0823+0546': [(0.956795923975906, p, -2.10081082736568) for p in [0, 60]] + [(0.956795923975906, 29.0070252532306, f) for f in [-1, -3]],
        'WD0611-6931Full': [(0.768404825227812, p, -1.91571596482581) for p in [0, 60]] + [(0.768404825227812, 29.661637997752, f) for f in [-1, -3]],
        'LHS2534': [(0.518208753763491, p, -1.28532198095239) for p in [0, 60]] + [(0.518208753763491, 35.415573421677, f) for f in [-1, -3]],
        #'WD2105-820': [(0.838348492316552, p, -1.6750110526307) for p in [0, 60]] + [(0.838348492316552, 18.7277178945314, f) for f in [-1, -3]],
        'GALEX1931+0117Melis': [(0.302430342279287, p, -2.63547009336713) for p in [0, 60]] + [(0.302430342279287, 15.9961525686632, f) for f in [-1, -3]],
        'IdealHPMantle1': [(0.002505105355787, p, -2.88191070564759) for p in [0, 60]] + [(0.002505105355787, 56.4533456889211, f) for f in [-1, -3]],
        'IdealHPCore1': [(0.993936868451227, p, -1.34753838938756) for p in [0, 60]] + [(0.993936868451227, 50.2971106091221, f) for f in [-1, -3]],
        'IdealLPMantle1': [(1.19943771354248E-06, p, -2.98008745573218) for p in [0, 60]] + [(1.19943771354248E-06, 0.024131980675524, f) for f in [-1, -3]],
        'GALEX1931+0117': [],
        'GaiaJ0510+2315': [],
        'GaiaJ0644-0352': [],
        'SynthEarthfcf0': [(0.013857052365627, p, -2.35147228174733) for p in [0, 60]] + [(0.013857052365627, 23.4507397312912, f) for f in [-1, -3]],
        #'SynthEarthfcf0': [(0, p, -2.35147228174733) for p in [0, 60]] + [(0, 23.4507397312912, f) for f in [-1, -3]] + [(0, 54, -2)], # The slight offset in fcf is enough to throw off the P/fO2 fit - rerun fixing fcf=0? And similarly for GD61?!?
        #'WD0449-259NoNa': [(0.637580010112239, p, -1.44238021207705) for p in [0, 60]] + [(0.637580010112239, 20.1788399789155, f) for f in [-1, -3]] + [(0.637580010112239, 0, -1)],
        'WD1350-162NoNa': [(0.414381955635034, p, -2.08319750488267) for p in [0, 60]] + [(0.414381955635034, 23.0031345846313, f) for f in [-1, -3]]
    }
    extended_fits = { # For changing any and all parameters. Usually this is overkill...
                 #Stellar metallicity indices	Time since Accretion/Myrs	log(Formation Distance/AU)	Feeding Zone Size/AU	Fragment Core Fraction	log(Pollution Fraction)	log(Accretion Event Timescale/Yrs)	Pressure /GPa	Oxygen Fugacity /ΔIW
        'WD0449-259NoNa': {
            'High fcf + heating': [800,6.97618833434471,-1.5,0.1,None,None,0.78,None,-7.5,2,10,-1.1],
            'Declining': [800,6.97618833434471,-0.4,0.1,None,None,0.42,None,-7.35,7.9,10,-1.1]
        },
        'IdealHPMantle4': {
            'low fo2': [616.187306876256,0.101288647261848,-0.509276468517034,0.05,None,None,0.012473279169699,None,-4.0068215773501,3.92043001595668,46.6028000969223,-3],
            'Ideal': [616.187306876256,0.101288647261848,-0.509276468517034,0.05,None,None,0.005,None,-4.0068215773501,3.92043001595668,60,-3]
        },
        'IdealLPMantle2': {
            'mod': [320.542433461403,0.647239361906416,-0.750639925657234,0.141993949232361,None,None,0.003,None,-3.98617291312008,2.7204957285566,0.087532511419774,-2.95460193289795],
            'mod2': [320.542433461403,0.647239361906416,-0.750639925657234,0.141993949232361,None,None,0.003,None,-3.98617291312008,2.7204957285566,60,-2.95460193289795]
        },
        'IdealHPCfcf75err5': {
            #'original': [356.91492883401,11.8351315556828,-0.584288132998062,0.097038753008786,None,None,0.824414901894526,None,-4.02152647155219,7.4333379374253,11.7852127666372,-2.40027437365142]
            'originaltweak': [356.91492883401,10.8351315556828,-0.584288132998062,0.097038753008786,None,None,0.824414901894526,None,-4.02152647155219,7.4333379374253,11.7852127666372,-2.40027437365142],
            'intended': [356.91492883401,0.01,-0.584288132998062,0.097038753008786,None,None,0.75,None,-4.02152647155219,7.4333379374253,60,-1],
            'intendedtweak': [356.91492883401,1,-0.584288132998062,0.097038753008786,None,None,0.75,None,-4.02152647155219,7.4333379374253,60,-1]
        },
        'SynthEarthfcf0.8': {
            #'original': [383.01136690618,4.99677630667471,-0.584687934079343,0.099466147558644,None,None,0.855818734724399,None,-4.00660297203982,7.14425290494588,24.2043203306747,-2.77307189857887]
            'originaltweak': [383.01136690618,3.99677630667471,-0.584687934079343,0.099466147558644,None,None,0.855818734724399,None,-4.00660297203982,7.14425290494588,24.2043203306747,-2.77307189857887],
            'intended': [383.01136690618,0.01,-0.584687934079343,0.099466147558644,None,None,0.8,None,-4.00660297203982,7.14425290494588,54,-2],
            'intendedtweak': [383.01136690618,1,-0.584687934079343,0.099466147558644,None,None,0.8,None,-4.00660297203982,7.14425290494588,54,-2]
        },
        'SynthEarthfcf0.7': {
            #'original': [596.938597489607,0.150314806568319,-0.510798429826005,0.05,None,None,0.66802004674132,None,-4.01984252881437,4.25676799308911,46.1640507207147,-2.75832533857166],
            'sinking': [596.938597489607,10.150314806568319,-0.560798429826005,0.1,None,None,0.75,None,-4.01984252881437,7.25676799308911,20,-2.75832533857166],
        },
        'WD1350-162NoNaCorr': {
            'original': [162.019231367967,6.83503111540107,0.01033097266606,0.05,None,None,0.414644299207587,None,-6.03802276078376,3.12952116775951,22.4471729043176,-2.05890354283997],
            'soln1': [162.019231367967,6.83503111540107,0.01033097266606,0.05,None,None,0.414644299207587,None,-6.03802276078376,3.12952116775951,5,-2.65],
            'soln2': [162.019231367967,6.83503111540107,-0.25,0.05,None,None,0.47,None,-6.03802276078376,3.12952116775951,50,-1.2],
            #'soln1': [162.019231367967,6.83503111540107,0.15,0.05,None,None,0.414644299207587,None,-6.03802276078376,1.8,5,-2.65],
            #'soln2': [162.019231367967,6.83503111540107,0,0.05,None,None,0.414644299207587,None,-6.03802276078376,2.5,50,-1.2]
        },
        'GALEX1931+0117GCorr': {
            #'original': [459.436688813984,0.002795356013718,2,0.05,None,None,None,None,-3.74182233851499,3.85652188295597,54,-2],
            'highfecr': [459,0.002795356013718,0.25,0.05,None,None,0.5,None,-4.03802276078376,3.85652188295597,60,-3],
        },
        'GALEX1931+0117VCorr': {
            #'original': [137.895272219274,0.003483767472115,0.552328519115545,0.05,None,None,0.158660318708795,None,-3.46325777513603,3.93522350145396,45.0079465185255,-2.62218607860118],
            'sinking': [137.895272219274,0.00861438,0.552328519115545,0.05,None,None,0.158660318708795,None,-3.46325777513603,3.93522350145396,45.0079465185255,-2.62218607860118],
        },
        'GALEX1931+0117MCorr': {
            #'original': [561.534557730359,0.000702881452596,0.433069991893536,0.070534707971139,None,None,0.314143312681015,None,-3.58125092627202,3.06579389404442,51.5925458254742,-2.73855772669102],
            'highfecrmg': [561.534557730359,0.000702881452596,0.5,0.070534707971139,None,None,0.4,None,-3.58125092627202,3.06579389404442,60,-1],
        },
        'WD2105-820': {
            'original': [209.665823846459,9.58421401303261E-05,0.070881042833234,0.05,None,None,0.837580130963948,None,-5.95833726628792,1.53608665508434,18.8654573745211,-1.67716051960656],
            'high P': [209.665823846459,9.58421401303261E-05,0.070881042833234,0.05,None,None,0.837580130963948,None,-5.95833726628792,1.53608665508434,60,-1.67716051960656],
        },
        'WD1232+563Corr': {
            #'original': [559.52104594338,13.2447022443142,2,0.05,None,None,None,None,-5.10106592894453,7.48555150602173,54,-2],
            'highfecrmg': [559.52104594338,11,2,0.05,None,None,None,None,-4.9,7,54,-2]
        },
        'SDSSJ0738+1835Corr': {
            #'original': [375.900885063282, 0.49247702176291, 0.099327565963983, 0.05, None, None, 0.440526840524824, None, -3.78919557507339, 2.38503940644603, 32.3298145531553, -2.53861020691926],
            '15, -2.8': [375.900885063282, 0.49247702176291, 0.099327565963983, 0.05, None, None, 0.440526840524824, None, -3.78919557507339, 2.38503940644603, 15, -2.8],
            '47, -2.8': [375.900885063282, 0.49247702176291, 0.099327565963983, 0.05, None, None, 0.440526840524824, None, -3.78919557507339, 2.38503940644603, 47, -2.8],
            '55, -1.2': [375.900885063282, 0.49247702176291, 0.099327565963983, 0.05, None, None, 0.440526840524824, None, -3.78919557507339, 2.38503940644603, 55, -1.2],
            '15, -1.2': [375.900885063282, 0.49247702176291, 0.099327565963983, 0.05, None, None, 0.440526840524824, None, -3.78919557507339, 2.38503940644603, 15, -1.2],
            '33, -2.8': [375.900885063282, 0.49247702176291, 0.099327565963983, 0.05, None, None, 0.440526840524824, None, -3.78919557507339, 2.38503940644603, 33, -2.8],
        }
    }
    print(trial_fits)
    data, abundances_dict, errors_dict, stellar_data, abundance_lower_bounds_dict, abundance_upper_bounds_dict = run_complete_model(trial_fits, extended_fits)
    synthetic_stellar_data = generate_synthetic_stellar_data(stellar_data)
    synthetic_wd_data = generate_synthetic_wd_data(abundances_dict, errors_dict)
    #pf_keys = list(data['SDSSJ1144+1218'].keys())  # This is a hack
    
    #For these, assume fO2 = IW - 2
    ni_over_mg_vals = dict()
    cr_over_mg_vals = dict()
    
    graph_fac = gf.GraphFactory()
    
    #graph_fac.make_composition_plot('PG0843+516', all_observations['PG0843+516'], all_errors['PG0843+516'], collect_model_output('PG0843+516', all_observations, data, trial_fits))
    #graph_fac.make_composition_plot('SDSSJ1535+1247', all_observations['SDSSJ1535+1247'], all_errors['SDSSJ1535+1247'], collect_model_output('SDSSJ1535+1247', all_observations, data, trial_fits))
    #graph_fac.make_composition_plot('SDSSJ0046+2717', all_observations['SDSSJ0046+2717'], all_errors['SDSSJ0046+2717'], collect_model_output('SDSSJ0046+2717', all_observations, data, trial_fits))
    
    manager = mn.Manager(
        Namespace(
            wd_data_filename='BlouinConglomNewTimescales.csv',
            stellar_compositions_filename='StellarCompositionsSortFE.csv',
            n_live_points = 0, # This argument shouldn't matter,
            enhancement_model = 'NonEarthlike',
            base_dir = '.',
            pollution_model_names=['Model_24']
        )
    )
    
    system_names, num_ref_systems, system_obs, system_obs_err, elel_by_system, system_ps, system_fcfs, text_offset_dict, stellar_dict, synth_stellar_dict, synth_wd_dict, systems_to_show_ellipse_dict, system_obs_lb, system_obs_ub = collect_best_fit_ratios(data, abundances_dict, errors_dict, stellar_data, synthetic_stellar_data, synthetic_wd_data, abundance_lower_bounds_dict, abundance_upper_bounds_dict) # Sorry
    made_up_crfe_contour_vals = list()
    made_up_mgfe_contour_vals = list()
    made_up_nife_contour_vals = list()
    made_up_crmg_contour_vals = list()
    made_up_sife_contour_vals = list()
    made_up_p_contour_vals = list()

    geo_model = gi.GeologyModel()
    example_fcf_lines = dict()
    fcf_list = list()
    fcf = 0
    while fcf < 0.95:
        fcf_list.append(fcf)
        if fcf > 0.9:
            fcf += 0.0002
        else:
            fcf += 0.001
        #fcf += 0.01
    do_big_loop = True
    if do_big_loop:
        for fragment_core_number_fraction in fcf_list:
            example_fcf_lines[fragment_core_number_fraction] = {'CrFe': list(), 'NiFe': list(), 'MgFe': list(), 'CrMg': list(), 'SiFe': list(), 'p': list()}
            fragment_mantle_number_fraction = 1 - fragment_core_number_fraction
            for p in list(range(0, 65, 5)):
                for fO2 in [-2]:
                    number_abundances, parent_core_number_fraction, Ds, all_Ds = geo_model.form_a_planet_iteratively(p, fO2)  # Could speed up if don't repeat this unnecessarily
                    ni_abundance = (fragment_mantle_number_fraction*number_abundances[ci.Element.Ni][gi.Layer.mantle]) + (fragment_core_number_fraction*number_abundances[ci.Element.Ni][gi.Layer.core])
                    fe_abundance = (fragment_mantle_number_fraction*number_abundances[ci.Element.Fe][gi.Layer.mantle]) + (fragment_core_number_fraction*number_abundances[ci.Element.Fe][gi.Layer.core])
                    cr_abundance = (fragment_mantle_number_fraction*number_abundances[ci.Element.Cr][gi.Layer.mantle]) + (fragment_core_number_fraction*number_abundances[ci.Element.Cr][gi.Layer.core])
                    mg_abundance = (fragment_mantle_number_fraction*number_abundances[ci.Element.Mg][gi.Layer.mantle]) + (fragment_core_number_fraction*number_abundances[ci.Element.Mg][gi.Layer.core])
                    si_abundance = (fragment_mantle_number_fraction*number_abundances[ci.Element.Si][gi.Layer.mantle]) + (fragment_core_number_fraction*number_abundances[ci.Element.Si][gi.Layer.core])
                    made_up_p_contour_vals.append(p)
                    made_up_crfe_contour_vals.append(np.log10(cr_abundance/fe_abundance))
                    made_up_mgfe_contour_vals.append(np.log10(mg_abundance/fe_abundance))
                    made_up_nife_contour_vals.append(np.log10(ni_abundance/fe_abundance))
                    made_up_crmg_contour_vals.append(np.log10(cr_abundance/mg_abundance))
                    made_up_sife_contour_vals.append(np.log10(si_abundance/fe_abundance))
                    example_fcf_lines[fragment_core_number_fraction]['CrFe'].append(np.log10(cr_abundance/fe_abundance))  # This will get messy if we ever do more than 1 fO2 value!
                    example_fcf_lines[fragment_core_number_fraction]['NiFe'].append(np.log10(ni_abundance/fe_abundance))
                    example_fcf_lines[fragment_core_number_fraction]['MgFe'].append(np.log10(mg_abundance/fe_abundance))
                    example_fcf_lines[fragment_core_number_fraction]['CrMg'].append(np.log10(cr_abundance/mg_abundance))
                    example_fcf_lines[fragment_core_number_fraction]['SiFe'].append(np.log10(si_abundance/fe_abundance))
                    example_fcf_lines[fragment_core_number_fraction]['p'].append(p)

    cr_plot_dict = graph_fac.plot_all_system_ratios_poster_version(system_names, num_ref_systems, system_obs['crfe'], system_obs['mgfe'], system_obs_err['crfe'], system_obs_err['mgfe'], elel_by_system['crfe'], elel_by_system['mgfe'], ci.Element.Cr, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe, system_ps, system_fcfs, made_up_crfe_contour_vals, made_up_mgfe_contour_vals, made_up_p_contour_vals, text_offset_dict['crfemgfe'], example_fcf_lines, stellar_dict['crfe'], stellar_dict['mgfe'], synth_stellar_dict['crfe'], synth_stellar_dict['mgfe'], calculate_sinking_vectors(ci.Element.Cr, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe), calculate_heating_vectors(ci.Element.Cr, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe), synth_wd_dict['crfe'], synth_wd_dict['mgfe'], systems_to_show_ellipse_dict['crfemgfe'], system_obs_lb['crfe'], system_obs_lb['mgfe'], system_obs_ub['crfe'], system_obs_ub['mgfe'], system_categories)
    ni_plot_dict = graph_fac.plot_all_system_ratios_poster_version(system_names, num_ref_systems, system_obs['nife'], system_obs['mgfe'], system_obs_err['nife'], system_obs_err['mgfe'], elel_by_system['nife'], elel_by_system['mgfe'], ci.Element.Ni, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe, system_ps, system_fcfs, made_up_nife_contour_vals, made_up_mgfe_contour_vals, made_up_p_contour_vals, text_offset_dict['nifemgfe'], example_fcf_lines, stellar_dict['nife'], stellar_dict['mgfe'], synth_stellar_dict['nife'], synth_stellar_dict['mgfe'], calculate_sinking_vectors(ci.Element.Ni, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe), calculate_heating_vectors(ci.Element.Ni, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe), synth_wd_dict['nife'], synth_wd_dict['mgfe'], systems_to_show_ellipse_dict['nifemgfe'], system_obs_lb['nife'], system_obs_lb['mgfe'], system_obs_ub['nife'], system_obs_ub['mgfe'], system_categories)    
    si_plot_dict = graph_fac.plot_all_system_ratios_poster_version(system_names, num_ref_systems, system_obs['sife'], system_obs['mgfe'], system_obs_err['sife'], system_obs_err['mgfe'], elel_by_system['sife'], elel_by_system['mgfe'], ci.Element.Si, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe, system_ps, system_fcfs, made_up_sife_contour_vals, made_up_mgfe_contour_vals, made_up_p_contour_vals, text_offset_dict['sifemgfe'], example_fcf_lines, stellar_dict['sife'], stellar_dict['mgfe'], synth_stellar_dict['sife'], synth_stellar_dict['mgfe'], calculate_sinking_vectors(ci.Element.Si, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe), calculate_heating_vectors(ci.Element.Si, ci.Element.Fe, ci.Element.Mg, ci.Element.Fe), synth_wd_dict['sife'], synth_wd_dict['mgfe'], systems_to_show_ellipse_dict['sifemgfe'], system_obs_lb['sife'], system_obs_lb['mgfe'], system_obs_ub['sife'], system_obs_ub['mgfe'], system_categories)
    graph_fac.multipanelise([cr_plot_dict, ni_plot_dict, si_plot_dict], 3, 1, ['bowtie_multipanel.png', 'bowtie_multipanel.pdf'], 15, 10)

    make_video = False
    if make_video:
        for system in ['GD61']:
            manager.publish_live_data(fise_run_args[system][0])
            collected_output = collect_model_output_mk2(system, data, trial_fits, make_video)
            for co_key, co_val in collected_output.items():
                if co_key == 'LP run':
                    continue
                fit_dict = {co_key: co_val}
                graph_fac.make_composition_plot_mk2(system, system + '_forvid_p' + str(int(100*co_key[1])).zfill(4) + '_', abundances_dict[system], errors_dict[system], fit_dict, None, None, abundance_upper_bounds_dict[system], abundance_lower_bounds_dict[system], make_video)
        return
    
    make_comp_and_elel_plots = False
    make_elel_plots = False
    if make_comp_and_elel_plots:
        for system in ['WD2105-820']:
            manager.publish_live_data(fise_run_args[system][0])
            only_pressure_in_legend = False
            collected_output = collect_model_output_mk2(system, manager.wd_types[system], data, trial_fits, extended_fits, only_pressure_in_legend)
            raw_output = collect_model_output_mk2(system, manager.wd_types[system], data, trial_fits, extended_fits, only_pressure_in_legend, True)
            collected_xfe_data = collect_xfe_data(system, data, trial_fits)
            print(abundances_dict[system])
            print(collected_output)
            print(raw_output)
            graph_fac.make_composition_plot_mk2(system, system + '_sisi_', abundances_dict[system], errors_dict[system], collected_output, None, None, abundance_upper_bounds_dict[system], abundance_lower_bounds_dict[system], False, only_pressure_in_legend)
            graph_fac.make_composition_plot_raw(system, manager.wd_types[system], system + '_sisi_', abundances_dict[system], errors_dict[system], raw_output, None, None, abundance_upper_bounds_dict[system], abundance_lower_bounds_dict[system], False, only_pressure_in_legend)
            if make_elel_plots:
                for element in [ci.Element.Cr, ci.Element.Ni]:
                    graph_fac.plot_log_elel_ratio(
                        collected_xfe_data[1][element],  # xfe vals
                        collected_xfe_data[0],  # fcf vals
                        'fcf',
                        system,
                        element,
                        ci.Element.Fe,
                        None,
                        abundances_dict[system][element] - abundances_dict[system][ci.Element.Fe],  # observed x/fe for this system  index 5 = Cr, 4 = Fe
                        np.sqrt(((errors_dict[system][element])**2) + ((errors_dict[system][ci.Element.Fe])**2)), # observed error
                        'nt',
                        None,
                        colour_dict,
                        False,
                        collected_xfe_data[2],  # fe / mg vals
                        abundances_dict[system][ci.Element.Fe] - abundances_dict[system][ci.Element.Mg],  # observed fe/mg for this system  index 6 = Mg, 4 = Fe
                        np.sqrt(((errors_dict[system][ci.Element.Fe])**2) + ((errors_dict[system][ci.Element.Mg])**2)) # observed error
                    )
                    graph_fac.plot_log_elel_ratio(
                        collected_xfe_data[4][element],  # xfe vals
                        collected_xfe_data[3],  # p vals
                        'Pressure',
                        system,
                        element,
                        ci.Element.Fe,
                        None,
                        abundances_dict[system][element] - abundances_dict[system][ci.Element.Fe],  # observed x/fe for this system  index 5 = Cr, 4 = Fe
                        np.sqrt(((errors_dict[system][element])**2) + ((errors_dict[system][ci.Element.Fe])**2)), # observed error
                        'nt',
                        None,
                        colour_dict,
                        False,
                        collected_xfe_data[5],  # fe / mg vals
                        abundances_dict[system][ci.Element.Fe] - abundances_dict[system][ci.Element.Mg],  # observed fe/mg for this system  index 6 = Mg, 4 = Fe
                        np.sqrt(((errors_dict[system][ci.Element.Fe])**2) + ((errors_dict[system][ci.Element.Mg])**2)) # observed error
                    )
                    if not np.isnan(abundances_dict[system][element]):
                        graph_fac.plot_3d_log_elel_ratio(
                            collected_xfe_data[6][element],
                            generate_pressure_vals(),
                            generate_fcf_vals(),
                            system,
                            element,
                            ci.Element.Fe,
                            fise_run_args[system][12],
                            abundances_dict[system][element] - abundances_dict[system][ci.Element.Fe],  # observed x/fe for this system  index 5 = Cr, 4 = Fe
                            np.sqrt(((errors_dict[system][element])**2) + ((errors_dict[system][ci.Element.Fe])**2)), # observed error
                            'nt',
                            'Pressure',
                            'fcf',
                            'fO2'
                        )
        
        #graph_fac.plot_3d_log_elel_ratio(wd_props['log_crfe_array'], generate_pressure_vals(), generate_fO2_vals(), wd, ci.Element.Cr, ci.Element.Fe, wd_props['cnf'], wd_props['observed_log_elfe_ratios'][ci.Element.Cr], wd_props['observed_log_elfe_ratio_errors'][ci.Element.Cr], tag)
    #graph_fac.make_composition_plot('SDSSJ1535+1247_r2000', all_observations['SDSSJ1535+1247'], all_errors['SDSSJ1535+1247'], collect_model_output('SDSSJ1535+1247_r2000', all_observations, data))
    #graph_fac.make_composition_plot('SDSSJ0046+2717_r2000', all_observations['SDSSJ0046+2717'], all_errors['SDSSJ0046+2717'], collect_model_output('SDSSJ0046+2717_r2000', all_observations, data))
    #graph_fac.make_composition_plot('PG0843+516_r2000', all_observations['PG0843+516'], all_errors['PG0843+516'], collect_model_output('PG0843+516_r2000', all_observations, data))
    
if __name__ == '__main__':
    main()
