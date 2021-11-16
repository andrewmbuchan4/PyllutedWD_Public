#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import numpy as np
import os

import model_parameters as mp

whitelist = [
    'SDSSJ0002+3209',
    'SDSSJ0004+0819',
    'SDSSJ0006+0520',
    'SDSSJ0010-0430',
    'SDSSJ0013+1109',
    'SDSSJ0019+2209',
    'SDSSJ0044+0418',
    'SDSSJ0046+2717',
    'SDSSJ0047+1628',
    'SDSSJ0052+1846',
    'SDSSJ0053+3115',
    'SDSSJ0056+2453',
    'SDSSJ0108-0537',
    'SDSSJ0114+3505',
    'SDSSJ0116+2050',
    'SDSSJ0117+0021',
    'SDSSJ0126+2534',
    'SDSSJ0135+1302',
    'SDSSJ0143+0113',
    'SDSSJ0144+1920',
    'SDSSJ0144+0305',
    'SDSSJ0148-0112',
    'SDSSJ0150+1354',
    'SDSSJ0158-0942',
    'SDSSJ0201+2015',
    'SDSSJ0208-0542',
    'SDSSJ0234-0510',
    'SDSSJ0252-0401',
    'SDSSJ0252+0054',
    'SDSSJ0447+1124',
    'SDSSJ0721+3928',
    'SDSSJ0736+4118',
    'SDSSJ0739+3112',
    'SDSSJ0741+3146',
    'SDSSJ0744+4649',
    'SDSSJ0744+4408',
    'SDSSJ0744+2701',
    'SDSSJ0744+1640',
    'SDSSJ0758+1013',
    'SDSSJ0800+2242',
    'SDSSJ0806+3055',
    'SDSSJ0807+4930',
    'SDSSJ0816+2330',
    'SDSSJ0818+1247',
    'SDSSJ0830-0319',
    'SDSSJ0838+2322',
    'SDSSJ0842+1406',
    'SDSSJ0842+1536',
    'SDSSJ0843+5614',
    'SDSSJ0851+1543',
    'SDSSJ0852+3402',
    'SDSSJ0901+0752',
    'SDSSJ0906+1141',
    'SDSSJ0908+5136',
    'SDSSJ0908+4119',
    'SDSSJ0913+2627',
    'SDSSJ0913+4127',
    'SDSSJ0916+2540',
    'SDSSJ0924+4301',
    'SDSSJ0925+3130',
    'SDSSJ0929+4247',
    'SDSSJ0933+6334',
    'SDSSJ0937+5228',
    'SDSSJ0939+4136',
    'SDSSJ0939+5019',
    'SDSSJ0946+2024',
    'SDSSJ0948+3008',
    'SDSSJ0956+5912',
    'SDSSJ1005+2244',
    'SDSSJ1006+1752',
    'SDSSJ1014+2827',
    'SDSSJ1017+3447',
    'SDSSJ1017+2419',
    'SDSSJ1019+3535',
    'SDSSJ1019+2045',
    'SDSSJ1024+4531',
    'SDSSJ1024+1014',
    'SDSSJ1032+1338',
    'SDSSJ1033+1809',
    'SDSSJ1038-0036',
    'SDSSJ1038+0432',
    'SDSSJ1040+2407',
    'SDSSJ1041+3432',
    'SDSSJ1043+3516',
    'SDSSJ1046+1329',
    'SDSSJ1055+3725',
    'SDSSJ1058+3143',
    'SDSSJ1102+2827',
    'SDSSJ1102+0214',
    'SDSSJ1103+4144',
    'SDSSJ1105+0228',
    'SDSSJ1112+0700',
    'SDSSJ1132+3323',
    'SDSSJ1134+1542',
    'SDSSJ1144+3720',
    'SDSSJ1144+1218',
    'SDSSJ1147+5429',
    'SDSSJ1149+0519',
    'SDSSJ1150+4928',
    'SDSSJ1152+5101',
    'SDSSJ1157+6138',
    'SDSSJ1158+0454',
    'SDSSJ1158+1845',
    'SDSSJ1158+4712',
    'SDSSJ1158+5448',
    'SDSSJ1158+5942',
    'SDSSJ1205+3536',
    'SDSSJ1211+2326',
    'SDSSJ1217+1157',
    'SDSSJ1218+0023',
    'SDSSJ1220+0929',
    'SDSSJ1224+2838',
    'SDSSJ1229+0743',
    'SDSSJ1230+3143',
    'SDSSJ1234+5208',
    'SDSSJ1238+2149',
    'SDSSJ1245+0822',
    'SDSSJ1254+3551',
    'SDSSJ1257+3238',
    'SDSSJ1257-0310',
    'SDSSJ1259+3112',
    'SDSSJ1259+4729',
    'SDSSJ1303+4055',
    'SDSSJ1308+0957',
    'SDSSJ1308+0258',
    'SDSSJ1314+3748',
    'SDSSJ1316+1918',
    'SDSSJ1319+3641',
    'SDSSJ1320+0204',
    'SDSSJ1321-0237',
    'SDSSJ1329+1301',
    'SDSSJ1336+3547',
    'SDSSJ1339+2643',
    'SDSSJ1340+2702',
    'SDSSJ1342+1813',
    'SDSSJ1345+1153',
    'SDSSJ1347+1415',
    'SDSSJ1350+1058',
    'SDSSJ1351+2645',
    'SDSSJ1356+2416',
    'SDSSJ1356+0236',
    'SDSSJ1401+3659',
    'SDSSJ1404+3620',
    'SDSSJ1405+2542',
    'SDSSJ1405+1549',
    'SDSSJ1411+3410',
    'SDSSJ1421+1843',
    'SDSSJ1428+4403',
    'SDSSJ1429+3841',
    'SDSSJ1430-0151',
    'SDSSJ1443+5833',
    'SDSSJ1443+3014',
    'SDSSJ1445+0913',
    'SDSSJ1448+1047',
    'SDSSJ1500+2315',
    'SDSSJ1502+3744',
    'SDSSJ1507+4034',
    'SDSSJ1518+0506',
    'SDSSJ1524+4049',
    'SDSSJ1534+1242',
    'SDSSJ1535+1247',
    'SDSSJ1537+3608',
    'SDSSJ1540+5352',
    'SDSSJ1542+4650',
    'SDSSJ1543+2024',
    'SDSSJ1545+5236',
    'SDSSJ1549+2633',
    'SDSSJ1549+1906',
    'SDSSJ1554+1735',
    'SDSSJ1604+1830',
    'SDSSJ1610+4006',
    'SDSSJ1612+3534',
    'SDSSJ1616+3303',
    'SDSSJ1624+3310',
    'SDSSJ1626+3303',
    'SDSSJ1627+4646',
    'SDSSJ1636+1619',
    'SDSSJ1641+1856',
    'SDSSJ1649+2238',
    'SDSSJ1706+2541',
    'SDSSJ2109-0039',
    'SDSSJ2110+0512',
    'SDSSJ2123+0016',
    'SDSSJ2157+1206',
    'SDSSJ2225+2338',
    'SDSSJ2230+1905',
    'SDSSJ2231+0906',
    'SDSSJ2235-0056',
    'SDSSJ2238+0213',
    'SDSSJ2238-0113',
    'SDSSJ2304+2415',
    'SDSSJ2319+3018',
    'SDSSJ2328+0830',
    'SDSSJ2330+2805',
    'SDSSJ2333+1058',
    'SDSSJ2340+0124',
    'SDSSJ2340+0817',
    'SDSSJ2343-0010',
    'SDSSJ2352+1922',
    'SDSSJ2352+3344',
    'SDSSJ2357+2348',
    #'GaiaJ1814-7355',
    #'SDSSJ2047-1259',
    #'PG0843+516',
    #'WD1551+175',
    #'WD1929+011',
    #'Earthfcf0',
    #'Earthfcf0.1',
    #'Earthfcf0.2',
    #'Earthfcf0.3',
    #'Earthfcf0.4',
    #'Earthfcf0.5',
    #'Earthfcf0.6',
    #'Earthfcf0.7',
    #'Earthfcf0.8',
    #'Earthfcf0.9',
    #'Earthfcf0.99',
    #'SynthMantle',
    #'SynthCore',
    #'PG0843+516Gaensicke',
    #'PG0843+516Xu',
    #'PG1015+161',
    #'GALEX1931+0117',
    #'GD362',
    #'GD40',
    #'G241-6',
    #'GD61',
    #'NLTT43806',
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
    #'GD 378',
    #'WD2216-657',
    'WD2157-574',
    #'WD2115-560',
    #'WD1232+563',
    #'WD2207+121',
    #'WD1145+017',
    #'WD0122-227',
    #'WD0449-259',
    #'WD1350-162',
    'WD2105-820',
    'SDSSJ0512-0505',
    'SDSSJ0823+0546',
    #'G29-38_CaFeMgSiO_1D',
    #'G29-38_CaFeMgSiO_3D',
    'GD424',
    #'WDJ1644-0449',
    #'WDJ2356–209',
    #'SDSSJ1330+6435',
    'LHS2534',
    #'WDJ2317+1830',
    #'WDJ1824+1213',
    'GALEXJ2339',
    'GD378',
    #'GD61noSi',
    #'GD61noO',
    #'GD61noSiO',
    #'GD424noSi',
    #'GD424noO',
    #'GD424noSiO',
    #'Marsfcf0',
    #'Marsfcf0.1',
    #'Marsfcf0.2',
    #'Marsfcf0.3',
    #'Marsfcf0.4',
    #'Marsfcf0.5',
    #'Marsfcf0.6',
    #'Marsfcf0.7',
    #'Marsfcf0.8',
    #'Marsfcf0.9',
    #'Marsfcf0.99',
    #'EarthFcf0',
    #'EarthFcf0.1',
    #'EarthFcf0.2',
    #'EarthFcf0.3',
    #'EarthFcf0.4',
    #'EarthFcf0.5',
    #'EarthFcf0.6',
    #'EarthFcf0.7',
    #'EarthFcf0.8',
    #'EarthFcf0.9',
    #'EarthFcf0.99',
    #'MarsFcf0',
    #'MarsFcf0.1',
    #'MarsFcf0.2',
    #'MarsFcf0.3',
    #'MarsFcf0.4',
    #'MarsFcf0.5',
    #'MarsFcf0.6',
    #'MarsFcf0.7',
    #'MarsFcf0.8',
    #'MarsFcf0.9',
    #'MarsFcf0.99',
    #'GaiaJ0510+2315',  # Can be included once published
    #'GaiaJ0644-0352',  # Can be included once published
    #'SDSSJ0006+2858',  # Can be included once published
    #'WD0611-6931',
    #'WD0611-6931Full',  # Can be included once published
    #'GALEX1931+0117Melis',
    #'SynthEarthfcf0',
    #'SynthEarthfcf0.1',
    #'SynthEarthfcf0.2',
    #'SynthEarthfcf0.3',
    #'SynthEarthfcf0.4',
    #'SynthEarthfcf0.5',
    #'SynthEarthfcf0.6',
    #'SynthEarthfcf0.7',
    #'SynthEarthfcf0.8',
    #'SynthEarthfcf0.9',
    #'SynthEarthfcf0.99',
    #'SynthMarsfcf0',
    #'SynthMarsfcf0.1',
    #'SynthMarsfcf0.2',
    #'SynthMarsfcf0.3',
    #'SynthMarsfcf0.4',
    #'SynthMarsfcf0.5',
    #'SynthMarsfcf0.6',
    #'SynthMarsfcf0.7',
    #'SynthMarsfcf0.8',
    #'SynthMarsfcf0.9',
    #'SynthMarsfcf0.99',
    #'IdealHPMantle1',
    #'IdealHPCore1',
    #'IdealLPMantle1',
    #'IdealLPCore1',
    #'IdealHPCore2',
    #'IdealHPCore3',
    #'IdealHPCore4',
    #'IdealHPMantle2',
    #'IdealLPMantle2',
    #'IdealLPCore2',
    #'IdealHPMantle3',
    #'IdealHPMantle4',
    #'IdealHPMantle5',
    #'WD0449-259NoNa',
    #'WD1350-162NoNa',
    #'WD2216-657Si',
    'WD0122-227SiO',
    #'WD2230-125',
    'WD0446-255',
    'WDJ1814-7354',
    #'PG0843+516GCorr',
    'PG0843+516XCorr',
    'GALEX1931+0117GCorr',
    #'GALEX1931+0117VCorr',
    #'GALEX1931+0117MCorr',
    'SDSSJ2047-1259Corr',
    'WD1551+175Corr',
    #'PG1015+161Corr', exclude  because this was based on Gaensicke data
    'GD362Corr',
    'GD40Corr',
    'G241-6Corr',
    'GD61Corr',
    'NLTT43806Corr',
    'HS2253+8023Corr',
    'SDSSJ0738+1835Corr',
    #'SDSSJ1228+1040OptCorr',
    'SDSSJ1228+1040Corr',
    'G29-38Corr',
    'SDSSJ1242+5226Corr',
    'SDSSJ0845+2257Corr',
    'WD1536+520Corr',
    'SDSSJ1043+0855Corr',
    'WD1425+540Corr',
    'WD2115-560Corr',
    'WD1232+563Corr',
    'WD2207+121Corr',
    'WD1145+017Corr',
    'WD0449-259NoNaCorr',
    'WD1350-162NoNaCorr',
    'WD2216-657SiCorr',
    #'IdealLPMantle2',
    'PG1225-079Corr',
    #'IdealHPMfcf0err5',
    #'IdealHPMfcf2err5',
    #'IdealHPMfcf0err10',
    #'IdealHPMfcf2err10',
    #'IdealHPMfcf0err20',
    #'IdealHPMfcf2err20',
    #'IdealHPMfcf0err40',
    #'IdealHPMfcf2err40',
    #'IdealLPMfcf0err5',
    #'IdealLPMfcf2err5',
    #'IdealLPMfcf0err10',
    #'IdealLPMfcf2err10',
    #'IdealLPMfcf0err20',
    #'IdealLPMfcf2err20',
    #'IdealLPMfcf0err40',
    #'IdealLPMfcf2err40',
    #'IdealHPCfcf75err5',
    #'IdealHPCfcf99err5',
    #'IdealHPCfcf75err10',
    #'IdealHPCfcf99err10',
    #'IdealHPCfcf75err20',
    #'IdealHPCfcf99err20',
    #'IdealHPCfcf75err40',
    #'IdealHPCfcf99err40',
    #'IdealLPCfcf75err5',
    #'IdealLPCfcf99err5',
    #'IdealLPCfcf75err10',
    #'IdealLPCfcf99err10',
    #'IdealLPCfcf75err20',
    #'IdealLPCfcf99err20',
    #'IdealLPCfcf75err40',
    #'IdealLPCfcf99err40',
    'G166-58',
    #'WD2221-165',
    #'WD0307+077',
    #'PG1541+651',
    'WD1145+288',
    #'WD1150-153',
    'GD56',
    #'WD0107-192',
    'HE0106-3253',
    'PG1015+161Xu',
    #'PG1457-086',
    #'WD1226+110',
    #'PG1018+411',
    #'WD1341+036',
    #'PG0010+280'
]

systems_in_sample_including_dupes = [
    'G166-58', # No evidence of differentiation
    'G241-6Corr',
    'G29-38Corr',
    'GALEX1931+0117MCorr', # Poor fit
    'GALEX1931+0117GCorr', # Poor fit
    'GALEX1931+0117VCorr', # Poor fit
    'GALEXJ2339',
    'GD362Corr',
    'GD378',
    'GD40Corr',
    'GD424', # High pressure mantle (not any more! NED now)
    'GD56',
    'GD61Corr',
    'HE0106-3253', # Pressure degenerate
    'HS2253+8023Corr',
    'LHS2534',
    'NLTT43806Corr',
    'PG0843+516GCorr',
    'PG0843+516XCorr',
    'PG1015+161Corr',
    'PG1015+161Xu', # Pressure unconstrained
    'PG1225-079Corr',
    'SDSSJ0512-0505',
    'SDSSJ0738+1835Corr',
    'SDSSJ0823+0546',
    'SDSSJ0845+2257Corr',
    'SDSSJ1043+0855Corr',
    'SDSSJ1228+1040OptCorr',
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
    'GD424': 'NED', # High pressure mantle
    'GD56': 'NED',
    'GD61Corr': 'HPM',
    'HE0106-3253': 'Pdegen', # Pressure degenerate
    'HS2253+8023Corr': 'NED',
    'LHS2534': 'Unphysical',
    'NLTT43806Corr': 'HPM*',
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
    'WD1350-162NoNaCorr': 'LPC$^\dagger$',
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

system_sources = {
    'G166-58': '\citet{Xu2019}',
    'G241-6Corr': '\citet{Jura2012}',
    'G29-38Corr': '\citet{Xu2014}',
    'GALEX1931+0117MCorr': '\citet{Melis2011}',
    'GALEX1931+0117GCorr': '\citet{Gaensicke2012}',
    'GALEX1931+0117VCorr': '\citet{Vennes2011b}',
    'GALEXJ2339': '\citet{Klein2021}',
    'GD362Corr': '\citet{Zuckerman2007}',
    'GD378': '\citet{Klein2021}',
    'GD40Corr': '\citet{Jura2012}',
    'GD424': '\citet{Izquierdo2020}',
    'GD56': '\citet{Xu2019}',
    'GD61Corr': '\citet{Farihi2013}',
    'HE0106-3253': '\citet{Xu2019}',
    'HS2253+8023Corr': '\citet{Klein2011}',
    'LHS2534': '\citet{Hollands2021}',
    'NLTT43806Corr': '\citet{Zuckerman2011}',
    'PG0843+516XCorr': '\citet{Xu2019}',
    'PG0843+516GCorr': '\citet{Gaensicke2012}',
    'PG1015+161Xu': '\citet{Xu2019}', 
    'PG1015+161Corr': '\citet{Gaensicke2012}', 
    'PG1225-079Corr': '\citet{Klein2011} and \citet{Xu2013}',
    'SDSSJ0512-0505': '\citet{Harrison2021}, although Si not previously included',
    'SDSSJ0738+1835Corr': '\citet{Dufour2012}',
    'SDSSJ0823+0546': '\citet{Harrison2021}, although Ti not previously included',
    'SDSSJ0845+2257Corr': '\citet{Wilson2015}',
    'SDSSJ1043+0855Corr': '\citet{Melis2016}',
    'SDSSJ1228+1040Corr': '\citet{Gaensicke2012}',
    'SDSSJ1228+1040OptCorr': '\citet{Gaensicke2012}',
    'SDSSJ1242+5226Corr': '\citet{Raddi2015}',
    'SDSSJ2047-1259Corr': '\citet{Hoskin2020}',
    'WD0122-227SiO': '\citet{Swan2019}',
    'WD0446-255': '\citet{Swan2019}',
    'WD0449-259NoNaCorr': '\citet{Swan2019}',
    'WD1145+017Corr': '\citet{Fortin-Archambault2020}',
    'WD1145+288': '\citet{Xu2019}',
    'WD1232+563Corr': '\citet{Xu2019}',
    'WD1350-162NoNaCorr': '\citet{Swan2019}',
    'WD1425+540Corr': '\citet{Xu2017}',
    'WD1536+520Corr': '\citet{Farihi2016}',
    'WD1551+175Corr': '\citet{Xu2019}',
    'WD2105-820': '\citet{Swan2019}',
    'WD2115-560Corr': '\citet{Swan2019}',
    'WD2157-574': '\citet{Swan2019}',
    'WD2207+121Corr': '\citet{Xu2019}',
    'WD2216-657SiCorr': '\citet{Swan2019}',
    'WDJ1814-7354': '\citet{GonzalezEgea2020}'
}

def extract_obs_number(stats_file):
    obs_number = stats_file.split('obs')[1].split('_')[0]
    return int(obs_number)

def find_stats_files(path=None):
    stats_files = list()
    for filename in os.listdir(path):
        if filename.startswith('stats') and filename.endswith('.csv'):
            stats_files.append(filename)
    stats_files.sort(key=extract_obs_number)
    return stats_files

def create_empty_variable_dict():
    toret = dict()
    for k, v in mp.model_parameter_strings.items():
        toret[v] = None
    return toret

def collect_system_stats(path=None):
    system_stats_dict = dict()
    na_string = 'N/A'
    volatile_rich_temp = 1000  # Below this temp, we're volatile rich
    volatile_poor_temp = 1400 # Above this temp, we're volatile poor (including moderate volatiles) (loosely basing this on Lodders 2003 table 8)
    for stats_file in find_stats_files(path):
        path_to_use = stats_file
        if path is not None:
            path_to_use = path + path_to_use
        system_name = None
        diff_sigma = None
        pcnf = None
        temperature = None
        best_model_index = None
        
        variable_indices = create_empty_variable_dict()
        medians = create_empty_variable_dict()
        percentile_16s = create_empty_variable_dict()
        percentile_84s = create_empty_variable_dict()
        upper_errors = create_empty_variable_dict()
        lower_errors = create_empty_variable_dict()

        good_fit_index = None
        good_fit = None
        with open(path_to_use, encoding='utf-8') as stats_csv:
            for row in csv.reader(stats_csv):
                if len(row) > 1:
                    if row[0] == 'System Name:':
                        system_name = row[1]
                        print(system_name)
                    if row[0] == 'Diff Sigma':
                        diff_sigma = row[1]
                    if row[0] == 'Parent Core Number Fraction:':
                        pcnf = row[1]
                    if row[0] == 'Temperature /K, +error, -error:':
                        temperature = row[1]
                    if row[0] == 'Parameter:':
                        for variable in variable_indices.keys():
                            try:
                                variable_indices[variable] = row.index(variable)
                            except ValueError:
                                # Then variable was not invoked for this system, ignore
                                pass
                    if row[0] == 'Median:':
                        for variable in medians.keys():
                            if variable_indices[variable] is not None:
                                medians[variable] = float(row[variable_indices[variable]])
                    if row[0] == '16th percentile:':
                        for variable in percentile_16s.keys():
                            if variable_indices[variable] is not None:
                                percentile_16s[variable] = float(row[variable_indices[variable]])
                    if row[0] == '84th percentile:':
                        for variable in percentile_84s.keys():
                            if variable_indices[variable] is not None:
                                percentile_84s[variable] = float(row[variable_indices[variable]])
                    if row[0] == 'Model':
                        best_model_index = row.index('Best model?')
                        good_fit_index = row.index('Good fit?')
                    if best_model_index is not None and len(row) > best_model_index + 1:
                        if row[best_model_index] == 'True':
                            good_fit = row[good_fit_index]
        
        for variable, median in medians.items():
            if median is not None:
                upper_error = percentile_84s[variable] - median
                lower_error = median - percentile_16s[variable]
                upper_errors[variable] = upper_error
                lower_errors[variable] = lower_error
        
        if temperature == '' or temperature is None:
            temperature = 0
        else:
            temperature = float(temperature)
        
        if diff_sigma is not None and diff_sigma != na_string:
            if diff_sigma == 'inf':
                diff_sigma = np.inf
            else:
                diff_sigma = float(diff_sigma)
        if good_fit is not None:
            good_fit = good_fit == 'True'
        
        if pcnf is not None:
            pcnf = float(pcnf)
        fcf = medians['Fragment Core Fraction']
            
        core_rich = fcf > pcnf if fcf is not None else False
        mantle_rich = fcf < pcnf if fcf is not None else False
        primitive = not core_rich and not mantle_rich
        
        volatile_rich = temperature < volatile_rich_temp
        volatile_depleted = temperature >= volatile_rich_temp and temperature <= volatile_poor_temp
        moderate_volatile_depleted = temperature > volatile_poor_temp

        system_stats_dict[system_name] = {
            'GoodFit': good_fit,
            'Primitive': primitive,
            'CoreRich': core_rich,
            'MantleRich': mantle_rich,
            'VolatileRich': volatile_rich,
            'VolatileDepleted': volatile_depleted,
            'ModerateVolatileDepleted': moderate_volatile_depleted,
            'DifferentiationSigma': diff_sigma,
            'Temperature': temperature,
            'ParentCoreNumberFraction': pcnf,
            'Medians': medians,
            'Upper1Sigmas': percentile_84s,
            'Lower1Sigmas': percentile_16s,
            'UpperErrors': upper_errors,
            'LowerErrors': lower_errors
        }
    return system_stats_dict

def latex_error_format(value, lower_error, upper_error, dp=2):
    if value is None or lower_error is None or upper_error is None:
        return 'N/A'
    if dp == 0:
        toret = str(int(value)) + '_{-' + str(int(lower_error)) + '}^{+' + str(int(upper_error)) + '}'
    else:
        toret = '{:.{dp}f}'.format(round(value, dp), dp=dp) + '_{-' + '{:.{dp}f}'.format(round(lower_error, dp), dp=dp) + '}^{+' + '{:.{dp}f}'.format(round(upper_error, dp), dp=dp) + '}'
    return '$'+toret+'$'

def strip_system_suffix(system_name):
    obs_data_series_name = system_name
    strippable_suffixes = ['Corr', 'NoNa', 'X', 'Xu', 'G', 'V', 'M', 'SiO', 'Si', 'Bounds', 'Opt']
    for suffix in strippable_suffixes:
        if obs_data_series_name.endswith(suffix):
            len_to_strip = len(suffix)
            obs_data_series_name = obs_data_series_name[:-len_to_strip]
    return obs_data_series_name

def compile_stats_table(system_stats_dict, path=None):
    path_to_use = 'system_stats.csv'
    latex_path_to_use = 'results.tex'
    if path is not None:
        path_to_use = path + path_to_use
        latex_path_to_use = path + latex_path_to_use
    with open(path_to_use, 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            'System',
            'Stellar metallicity',
            'Stellar metallicity Lower Error',
            'Stellar metallicity Upper Error',
            'Time since Accretion/Myrs',
            'Time since Accretion/Myrs Lower Error',
            'Time since Accretion/Myrs Upper Error',
            'log(Formation Distance/AU)',
            'log(Formation Distance/AU) Lower Error',
            'log(Formation Distance/AU) Upper Error',
            'Feeding Zone Size/AU',
            'Feeding Zone Size/AU Lower Error',
            'Feeding Zone Size/AU Upper Error',
            'Parent Core Number Fraction (fitted)',
            'Parent Core Number Fraction (fitted) Lower Error',
            'Parent Core Number Fraction (fitted) Upper Error',
            'Parent Crust Number Fraction',
            'Parent Crust Number Fraction Lower Error',
            'Parent Crust Number Fraction Upper Error',
            'Fragment Core Number Fraction',
            'Fragment Core Number Fraction Lower Error',
            'Fragment Core Number Fraction Upper Error',
            'Fragment Crust Number Fraction',
            'Fragment Crust Number Fraction Lower Error',
            'Fragment Crust Number Fraction Upper Error',
            'log(Pollution Fraction)',
            'log(Pollution Fraction) Lower Error',
            'log(Pollution Fraction) Upper Error',
            'log(Accretion Event Timescale/Yrs)',
            'log(Accretion Event Timescale/Yrs) Lower Error',
            'log(Accretion Event Timescale/Yrs) Upper Error',
            'Pressure /GPa',
            'Pressure /GPa Lower Error',
            'Pressure /GPa Upper Error',
            'Oxygen Fugacity /ΔIW',
            'Oxygen Fugacity /ΔIW Lower Error',
            'Oxygen Fugacity /ΔIW Upper Error',
            'Good Fit?',
            'Primitive?',
            'Core Rich?',
            'Mantle Rich?',
            'Volatile Rich?',
            'Volatile Depleted?',
            'Moderate Volatile Depleted?',
            'Temperature',
            'Inferred Parent Core Number Fraction',
            'Differentiation Sigma',
            'Include In Sample?',
        ])
        for sys, vals in system_stats_dict.items():
            to_write.writerow([
                sys,
                vals['Medians']['Stellar metallicity indices'] if vals['Medians']['Stellar metallicity indices'] is not None else 'N/A',
                vals['LowerErrors']['Stellar metallicity indices'] if vals['LowerErrors']['Stellar metallicity indices'] is not None else 'N/A',
                vals['UpperErrors']['Stellar metallicity indices'] if vals['UpperErrors']['Stellar metallicity indices'] is not None else 'N/A',
                vals['Medians']['Time since Accretion/Myrs'] if vals['Medians']['Time since Accretion/Myrs'] is not None else 'N/A',
                vals['LowerErrors']['Time since Accretion/Myrs'] if vals['LowerErrors']['Time since Accretion/Myrs'] is not None else 'N/A',
                vals['UpperErrors']['Time since Accretion/Myrs'] if vals['UpperErrors']['Time since Accretion/Myrs'] is not None else 'N/A',
                vals['Medians']['log(Formation Distance/AU)'] if vals['Medians']['log(Formation Distance/AU)'] is not None else 'N/A',
                vals['LowerErrors']['log(Formation Distance/AU)'] if vals['LowerErrors']['log(Formation Distance/AU)'] is not None else 'N/A',
                vals['UpperErrors']['log(Formation Distance/AU)'] if vals['UpperErrors']['log(Formation Distance/AU)'] is not None else 'N/A',
                vals['Medians']['Feeding Zone Size/AU'] if vals['Medians']['Feeding Zone Size/AU'] is not None else 'N/A',
                vals['LowerErrors']['Feeding Zone Size/AU'] if vals['LowerErrors']['Feeding Zone Size/AU'] is not None else 'N/A',
                vals['UpperErrors']['Feeding Zone Size/AU'] if vals['UpperErrors']['Feeding Zone Size/AU'] is not None else 'N/A',
                vals['Medians']['Parent Core Fraction'] if vals['Medians']['Parent Core Fraction'] is not None else 'N/A',
                vals['LowerErrors']['Parent Core Fraction'] if vals['LowerErrors']['Parent Core Fraction'] is not None else 'N/A',
                vals['UpperErrors']['Parent Core Fraction'] if vals['UpperErrors']['Parent Core Fraction'] is not None else 'N/A',
                vals['Medians']['Parent Crust Fraction'] if vals['Medians']['Parent Crust Fraction'] is not None else 'N/A',
                vals['LowerErrors']['Parent Crust Fraction'] if vals['LowerErrors']['Parent Crust Fraction'] is not None else 'N/A',
                vals['UpperErrors']['Parent Crust Fraction'] if vals['UpperErrors']['Parent Crust Fraction'] is not None else 'N/A',
                vals['Medians']['Fragment Core Fraction'] if vals['Medians']['Fragment Core Fraction'] is not None else 'N/A',
                vals['LowerErrors']['Fragment Core Fraction'] if vals['LowerErrors']['Fragment Core Fraction'] is not None else 'N/A',
                vals['UpperErrors']['Fragment Core Fraction'] if vals['UpperErrors']['Fragment Core Fraction'] is not None else 'N/A',
                vals['Medians']['Fragment Crust Fraction'] if vals['Medians']['Fragment Crust Fraction'] is not None else 'N/A',
                vals['LowerErrors']['Fragment Crust Fraction'] if vals['LowerErrors']['Fragment Crust Fraction'] is not None else 'N/A',
                vals['UpperErrors']['Fragment Crust Fraction'] if vals['UpperErrors']['Fragment Crust Fraction'] is not None else 'N/A',
                vals['Medians']['log(Pollution Fraction)'] if vals['Medians']['log(Pollution Fraction)'] is not None else 'N/A',
                vals['LowerErrors']['log(Pollution Fraction)'] if vals['LowerErrors']['log(Pollution Fraction)'] is not None else 'N/A',
                vals['UpperErrors']['log(Pollution Fraction)'] if vals['UpperErrors']['log(Pollution Fraction)'] is not None else 'N/A',
                vals['Medians']['log(Accretion Event Timescale/Yrs)'] if vals['Medians']['log(Accretion Event Timescale/Yrs)'] is not None else 'N/A',
                vals['LowerErrors']['log(Accretion Event Timescale/Yrs)'] if vals['LowerErrors']['log(Accretion Event Timescale/Yrs)'] is not None else 'N/A',
                vals['UpperErrors']['log(Accretion Event Timescale/Yrs)'] if vals['UpperErrors']['log(Accretion Event Timescale/Yrs)'] is not None else 'N/A',
                vals['Medians']['Pressure /GPa'] if vals['Medians']['Pressure /GPa'] is not None else 'N/A',
                vals['LowerErrors']['Pressure /GPa'] if vals['LowerErrors']['Pressure /GPa'] is not None else 'N/A',
                vals['UpperErrors']['Pressure /GPa'] if vals['UpperErrors']['Pressure /GPa'] is not None else 'N/A',
                vals['Medians']['Oxygen Fugacity /ΔIW'] if vals['Medians']['Oxygen Fugacity /ΔIW'] is not None else 'N/A',
                vals['LowerErrors']['Oxygen Fugacity /ΔIW'] if vals['LowerErrors']['Oxygen Fugacity /ΔIW'] is not None else 'N/A',
                vals['UpperErrors']['Oxygen Fugacity /ΔIW'] if vals['UpperErrors']['Oxygen Fugacity /ΔIW'] is not None else 'N/A',
                vals['GoodFit'],
                vals['Primitive'],
                vals['CoreRich'],
                vals['MantleRich'],
                vals['VolatileRich'],
                vals['VolatileDepleted'],
                vals['ModerateVolatileDepleted'],
                vals['Temperature'] if vals['Temperature'] != 0 else 'N/A',
                vals['ParentCoreNumberFraction'] if vals['ParentCoreNumberFraction'] is not None else 'N/A',
                vals['DifferentiationSigma'],
                sys in whitelist
            ])
    with open(latex_path_to_use, 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            '%System',
            'Source (abundances)',
            'Stellar metallicity',
            'Time since Accretion/Myrs',
            'log(Formation Distance/AU)',
            'Feeding Zone Size/AU',
            'Fragment Core Number Fraction',
            'log(Pollution Fraction)',
            'log(Accretion Event Timescale/Yrs)',
            'Pressure /GPa',
            'Oxygen Fugacity /ΔIW',
            'Good Fit?',
            'Primitive?',
            'Core Rich?',
            'Mantle Rich?',
            'Volatile Rich?',
            'Volatile Depleted?',
            'Moderate Volatile Depleted?',
            'Temperature',
            'Inferred Parent Core Number Fraction',
            'Differentiation Sigma',
            'Category'
        ])
        list_of_all_rows = list()
        for sys, vals in system_stats_dict.items():
            print(sys)
            print(vals['DifferentiationSigma'])
            if sys in systems_in_sample_including_dupes:
                list_of_vals = [
                    strip_system_suffix(sys),
                    system_sources[sys],
                    latex_error_format(
                        vals['Medians']['Stellar metallicity indices'],
                        vals['LowerErrors']['Stellar metallicity indices'],
                        vals['UpperErrors']['Stellar metallicity indices'],
                        0
                    ),
                    latex_error_format(
                        vals['Medians']['Time since Accretion/Myrs'],
                        vals['LowerErrors']['Time since Accretion/Myrs'],
                        vals['UpperErrors']['Time since Accretion/Myrs']
                    ),
                    latex_error_format(
                        vals['Medians']['log(Formation Distance/AU)'],
                        vals['LowerErrors']['log(Formation Distance/AU)'],
                        vals['UpperErrors']['log(Formation Distance/AU)']
                    ),
                    latex_error_format(
                        vals['Medians']['Feeding Zone Size/AU'],
                        vals['LowerErrors']['Feeding Zone Size/AU'],
                        vals['UpperErrors']['Feeding Zone Size/AU']
                    ),
                    latex_error_format(
                        vals['Medians']['Fragment Core Fraction'],
                        vals['LowerErrors']['Fragment Core Fraction'],
                        vals['UpperErrors']['Fragment Core Fraction'],
                        2
                    ),
                    latex_error_format(
                        vals['Medians']['log(Pollution Fraction)'],
                        vals['LowerErrors']['log(Pollution Fraction)'],
                        vals['UpperErrors']['log(Pollution Fraction)']
                    ),
                    latex_error_format(
                        vals['Medians']['log(Accretion Event Timescale/Yrs)'],
                        vals['LowerErrors']['log(Accretion Event Timescale/Yrs)'],
                        vals['UpperErrors']['log(Accretion Event Timescale/Yrs)']
                    ),
                    latex_error_format(
                        vals['Medians']['Pressure /GPa'],
                        vals['LowerErrors']['Pressure /GPa'],
                        vals['UpperErrors']['Pressure /GPa'],
                        1
                    ),
                    latex_error_format(
                        vals['Medians']['Oxygen Fugacity /ΔIW'],
                        vals['LowerErrors']['Oxygen Fugacity /ΔIW'],
                        vals['UpperErrors']['Oxygen Fugacity /ΔIW'],
                        1
                    ),
                    'Y' if vals['GoodFit'] else 'N',
                    'Y' if vals['Primitive'] else 'N',
                    'Y' if vals['CoreRich'] else 'N',
                    'Y' if vals['MantleRich'] else 'N',
                    'Y' if vals['VolatileRich'] else 'N',
                    'Y' if vals['VolatileDepleted'] else 'N',
                    'Y' if vals['ModerateVolatileDepleted'] else 'N',
                    str(int(vals['Temperature'])) if vals['Temperature'] != 0 else 'N/A',
                    str(round(vals['ParentCoreNumberFraction'],2)) if vals['ParentCoreNumberFraction'] is not None else 'N/A',
                    str(round(vals['DifferentiationSigma'], 1)) if vals['DifferentiationSigma'] is not None and vals['DifferentiationSigma'] != 'N/A' else 'N/A',
                    system_categories.get(sys, 'N/A')
                ]
                list_of_all_rows.append([' & '.join(list_of_vals) + ' \\\\'])
        for row in sorted(list_of_all_rows):
            to_write.writerow(row)
    with open(latex_path_to_use + '.mod1', 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            '%System',
            'Stellar metallicity',
            'Time since Accretion/Myrs',
            'log(Formation Distance/AU)',
            'Feeding Zone Size/AU',
            'Fragment Core Number Fraction',
            'log(Pollution Fraction)',
            'log(Accretion Event Timescale/Yrs)',
            'Pressure /GPa',
            'Oxygen Fugacity /ΔIW'
        ])
        list_of_all_rows = list()
        for sys, vals in system_stats_dict.items():
            print(sys)
            print(vals['DifferentiationSigma'])
            if sys in systems_in_sample_including_dupes:
                list_of_vals = [
                    strip_system_suffix(sys),
                    latex_error_format(
                        vals['Medians']['Stellar metallicity indices'],
                        vals['LowerErrors']['Stellar metallicity indices'],
                        vals['UpperErrors']['Stellar metallicity indices'],
                        0
                    ),
                    latex_error_format(
                        vals['Medians']['Time since Accretion/Myrs'],
                        vals['LowerErrors']['Time since Accretion/Myrs'],
                        vals['UpperErrors']['Time since Accretion/Myrs']
                    ),
                    latex_error_format(
                        vals['Medians']['log(Formation Distance/AU)'],
                        vals['LowerErrors']['log(Formation Distance/AU)'],
                        vals['UpperErrors']['log(Formation Distance/AU)']
                    ),
                    latex_error_format(
                        vals['Medians']['Feeding Zone Size/AU'],
                        vals['LowerErrors']['Feeding Zone Size/AU'],
                        vals['UpperErrors']['Feeding Zone Size/AU']
                    ),
                    latex_error_format(
                        vals['Medians']['Fragment Core Fraction'],
                        vals['LowerErrors']['Fragment Core Fraction'],
                        vals['UpperErrors']['Fragment Core Fraction'],
                        2
                    ),
                    latex_error_format(
                        vals['Medians']['log(Pollution Fraction)'],
                        vals['LowerErrors']['log(Pollution Fraction)'],
                        vals['UpperErrors']['log(Pollution Fraction)']
                    ),
                    latex_error_format(
                        vals['Medians']['log(Accretion Event Timescale/Yrs)'],
                        vals['LowerErrors']['log(Accretion Event Timescale/Yrs)'],
                        vals['UpperErrors']['log(Accretion Event Timescale/Yrs)']
                    ),
                    latex_error_format(
                        vals['Medians']['Pressure /GPa'],
                        vals['LowerErrors']['Pressure /GPa'],
                        vals['UpperErrors']['Pressure /GPa'],
                        1
                    ),
                    latex_error_format(
                        vals['Medians']['Oxygen Fugacity /ΔIW'],
                        vals['LowerErrors']['Oxygen Fugacity /ΔIW'],
                        vals['UpperErrors']['Oxygen Fugacity /ΔIW'],
                        1
                    )
                ]
                list_of_all_rows.append([' & '.join(list_of_vals) + ' \\\\'])
        for row in sorted(list_of_all_rows):
            to_write.writerow(row)
    with open(latex_path_to_use + '.mod2', 'w', newline='', encoding='utf-8') as f:
        to_write = csv.writer(f)
        to_write.writerow([
            '%System',
            'Good Fit?',
            'Primitive?',
            'Core Rich?',
            'Mantle Rich?',
            'Volatile Rich?',
            'Volatile Depleted?',
            'Moderate Volatile Depleted?',
            'Temperature',
            'Inferred Parent Core Number Fraction',
            'Differentiation Sigma',
            'Category'
        ])
        list_of_all_rows = list()
        for sys, vals in system_stats_dict.items():
            print(sys)
            print(vals['DifferentiationSigma'])
            if sys in systems_in_sample_including_dupes:
                list_of_vals = [
                    strip_system_suffix(sys),
                    'Y' if vals['GoodFit'] else 'N',
                    'Y' if vals['Primitive'] else 'N',
                    'Y' if vals['CoreRich'] else 'N',
                    'Y' if vals['MantleRich'] else 'N',
                    'Y' if vals['VolatileRich'] else 'N',
                    'Y' if vals['VolatileDepleted'] else 'N',
                    'Y' if vals['ModerateVolatileDepleted'] else 'N',
                    str(int(vals['Temperature'])) if vals['Temperature'] != 0 else 'N/A',
                    str(round(vals['ParentCoreNumberFraction'],2)) if vals['ParentCoreNumberFraction'] is not None else 'N/A',
                    str(round(vals['DifferentiationSigma'], 1)) if vals['DifferentiationSigma'] is not None and vals['DifferentiationSigma'] != 'N/A' else 'N/A',
                    system_categories.get(sys, 'N/A')
                ]
                list_of_all_rows.append([' & '.join(list_of_vals) + ' \\\\'])
        for row in sorted(list_of_all_rows):
            to_write.writerow(row)

def compile_stats_summary(system_stats_dict, path=None):
    unclassified = 0
    differentiated = 0
    differentiated_sigma_requirement = 3
    classifications = {
        'VolatileRich': {
            'Primitive': 0,
            'CoreRich': 0,
            'MantleRich': 0
        },
        'VolatileDepleted': {
            'Primitive': 0,
            'CoreRich': 0,
            'MantleRich': 0
        },
        'ModerateVolatileDepleted': {
            'Primitive': 0,
            'CoreRich': 0,
            'MantleRich': 0
        }
    }
    for sys, vals in system_stats_dict.items():
        if sys in whitelist:
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
                    
                classifications[t_classification][c_classification] += 1
                try:
                    if vals['DifferentiationSigma'] > differentiated_sigma_requirement:
                        differentiated += 1
                except TypeError:
                    # Then diff sigma was None or N/A
                    pass

    path_to_use = 'system_stats_summary.csv'
    if path is not None:
        path_to_use = path + path_to_use
    with open(path_to_use, 'w', newline='', encoding='utf-8') as f:
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
            'Mantle Rich'
        ])
        for t_classification, vals in classifications.items():
            to_write.writerow([
                t_classification,
                vals['Primitive'],
                vals['CoreRich'],
                vals['MantleRich']
            ])
        to_write.writerow([])
        to_write.writerow([
            'Differentiated to >' + str(differentiated_sigma_requirement) + ' sigma:',
            differentiated
        ])

def main():
    path = '<your_filepath_here>'
    system_stats_dict = collect_system_stats(path)
    compile_stats_table(system_stats_dict, path)
    compile_stats_summary(system_stats_dict, path)
    
if __name__ == '__main__':
    main()
