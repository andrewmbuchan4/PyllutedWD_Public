#!/usr/bin/env python

import math
import numpy as np

import chemistry_info as ci
import graph_factory as gf
import pewdd_sample_selector as pss
import pewdd_values as pv
import timescale_interpolator as ti
import thermohaline_interpolator as thi

all_plausible_he_atm_wds = {
    'LHS2534_Hollands2021': (4780.0, 7.97),
    'SDSSJ122943.92+074311.8_Blouin2020': (6014.0, 7.133),
    'SDSSJ094813.74+300851.2_Blouin2020': (6284.0, 7.21),
    'SDSSJ010825.79-053755.6_Blouin2020': (6250.0, 7.256),
    'WDJ1824+1213_Hollands2021': (3350.0, 7.41),
    'SDSSJ121731.31+115715.9_Blouin2020': (6012.0, 7.401),
    'SDSSJ072144.23+392843.3_Blouin2020': (6022.0, 7.433),
    'SDSSJ234307.67-001016.3_Blouin2020': (5778.0, 7.475),
    'SDSSJ110304.15+414434.9_Blouin2020': (5728.0, 7.521),
    'SDSSJ105533.73+372542.7_Blouin2020': (5614.0, 7.605),
    'SDSSJ001949.26+220926.5_Blouin2020': (5797.0, 7.608),
    'SDSSJ131612.87+191806.5_Blouin2020': (5160.0, 7.639),
    'SDSSJ142939.38+384113.2_Blouin2020': (5636.0, 7.634),
    'SDSSJ1330+6435_Hollands2021': (3660.0, 7.65),
    'SDSSJ132941.79+130131.9_Blouin2020': (6706.0, 7.636),
    'SDSSJ123024.04+314339.7_Blouin2020': (6310.0, 7.676),
    'PG1225-079Model1_Klein2011': (10500.0, 7.7),
    'SDSSJ135123.86+264546.5_Blouin2020': (5640.0, 7.73),
    'SDSSJ085217.60+340211.3_Blouin2020': (4806.0, 7.743),
    'SDSSJ083033.66-031911.0_Blouin2020': (6424.0, 7.736),
    'SDSSJ114408.05+372007.8_Blouin2020': (7280.0, 7.768),
    'SDSSJ090803.35+513633.1_Blouin2020': (5779.0, 7.78),
    'SDSSJ162703.34+464658.2_Blouin2020': (6000.0, 7.785),
    'WDJ1644-0449_Kaiser2021': (3830.0, 7.77),
    'SDSSJ125945.32+472953.6_Blouin2020': (5682.0, 7.795),
    'SDSSJ000614.53+052039.0_Blouin2020': (5783.0, 7.839),
    'SDSSJ092932.50+424757.9_Blouin2020': (6676.0, 7.853),
    'SDSSJ044751.21+112403.8_Blouin2020': (6966.0, 7.858),
    'SDSSJ111215.06+070052.4_Blouin2020': (6891.0, 7.867),
    'SDSSJ014415.12+192021.4_Blouin2020': (6024.0, 7.877),
    'WD1622+587(Spec;Opt)_Rogers2023': (23430.0, 7.9),
    'WD1622+587(Spec;UV)_Rogers2023': (23430.0, 7.9),
    'SDSSJ104658.12+132911.3_Blouin2020': (5177.0, 7.895),
    'SDSSJ074456.21+164041.8_Blouin2020': (4703.0, 7.91),
    'J0848_Swan2023': (11700.0, 7.89),
    'SDSSJ114944.95+051947.7_Blouin2020': (7173.0, 7.914),
    'SDSSJ025206.06-040130.3_Blouin2020': (6773.0, 7.924),
    'SDSSJ230414.49+241554.0_Blouin2020': (5102.0, 7.935),
    'SDSSJ215752.30+120603.1_Blouin2020': (6042.0, 7.93),
    'SDSSJ222503.71+233855.1_Blouin2020': (6029.0, 7.937),
    'SDSSJ093944.58+501917.6_Blouin2020': (6030.0, 7.939),
    'SDSSJ164104.95+185602.1_Blouin2020': (5470.0, 7.941),
    'J1618_Swan2023': (9870.0, 7.94),
    'GALEXJ2339_Klein2021': (13735.0, 7.93),
    'GD378_Klein2021': (15620.0, 7.93),
    'SDSSJ140137.31+365909.9_Blouin2020': (5931.0, 7.973),
    'WDJ2356-209_Kaiser2021': (4040.0, 7.98),
    'SDSSJ142120.11+184351.6_Blouin2020': (6517.0, 7.978),
    'SDSSJ115818.76+045447.5_Blouin2020': (5344.0, 7.991),
    'SDSSJ015849.02-094225.3_Blouin2020': (6115.0, 7.99),
    'SDSSJ084223.14+140615.9_Blouin2020': (7075.0, 7.985),
    'SDSSJ223101.12+090635.1_Blouin2020': (5679.0, 7.992),
    'SDSSJ161603.00+330301.2_Blouin2020': (6491.0, 7.991),
    'SDSSJ090146.87+075206.8_Blouin2020': (7263.0, 7.989),
    'WDJ1922+0233_Elms2022': (3343.0, 8.0),
    'J1313_Swan2023': (8800.0, 7.99),
    'SDSSJ210916.51-003921.6_Blouin2020': (6132.0, 7.997),
    'SDSSJ080740.69+493059.7_Blouin2020': (5172.0, 8.005),
    'SDSSJ234034.61+012416.6_Blouin2020': (6072.0, 8.0),
    'SDSSJ101711.54+344710.5_Blouin2020': (6089.0, 8.004),
    'SDSSJ122437.07+283853.0_Blouin2020': (4991.0, 8.01),
    'SDSSJ114441.92+121829.2_Blouin2020': (5320.0, 8.009),
    'GD17_GentileFusillo2017': (8300.0, 8.0),
    'SDSSJ135632.63+241606.0_Blouin2020': (6173.0, 8.01),
    'WD1536+520_Farihi2016': (20800.0, 7.96),
    'PG1225-079Model2_Klein2011': (10800.0, 8.0),
    'GaiaJ0644-0352(Phot;Opt)_Rogers2023': (17000.0, 7.98),
    'WD0446-255_Swan2019': (10120.0, 8.0),
    'PG1225-079_Xu2013': (10800.0, 8.0),
    'PG1225-079Updated_Xu2013': (10800.0, 8.0),
    'GD16_GentileFusillo2017': (11000.0, 8.0),
    'SDSSJ093719.14+522802.4_Blouin2020': (6660.0, 8.015),
    'SDSSJ115809.88+594210.3_Blouin2020': (6046.0, 8.022),
    'SDSSJ125922.03+311215.2_Blouin2020': (5664.0, 8.027),
    'SDSSJ0956+5912(NoOvershoot)_Hollands2022': (8100.0, 8.02),
    'SDSSJ0956+5912(Overshoot)_Hollands2022': (8100.0, 8.02),
    'GD40_Klein2010': (15300.0, 8.0),
    'SDSSJ124231.07+522626.6_Raddi2015': (13000.0, 8.0),
    'WD1622+587(Phot;Opt)_Rogers2023': (21530.0, 7.98),
    'WD1622+587(Phot;UV)_Rogers2023': (21530.0, 7.98),
    'SDSSJ135654.01+023641.5_Blouin2020': (7662.0, 8.028),
    'SDSSJ025253.20+005439.4_Blouin2020': (7478.0, 8.031),
    'SDSSJ125454.93+355145.7_Blouin2020': (6417.0, 8.036),
    'SDSSJ115809.38+184557.3_Blouin2020': (6696.0, 8.04),
    'WD1350-162_Swan2019': (11640.0, 8.02),
    'SDSSJ074414.66+464912.5_Blouin2020': (4861.0, 8.052),
    'SDSSJ011759.81+002138.1_Blouin2020': (6994.0, 8.052),
    'SDSSJ133624.26+354751.2_Blouin2020': (6172.0, 8.057),
    'SDSSJ1038-0036_Hollands2022': (7560.0, 8.06),
    'SDSSJ1535+1247_Hollands2022': (5680.0, 8.06),
    'WD0122-227_Swan2019': (8380.0, 8.06),
    'WD2216-657_Swan2019': (9190.0, 8.05),
    'WD0449-259_Swan2019': (9850.0, 8.04),
    'SDSSJ091621.36+254028.4_Blouin2020': (5497.0, 8.066),
    'SDSSJ142833.78+440346.2_Blouin2020': (6574.0, 8.064),
    'SDSSJ154913.46+263301.1_Blouin2020': (6794.0, 8.063),
    'WDJ2047-1259_Hoskin2020': (17970.0, 8.04),
    'SDSSJ130328.07+405545.5_Blouin2020': (6481.0, 8.071),
    'J0939_Swan2023': (9020.0, 8.07),
    'J1227_Swan2023': (7420.0, 8.07),
    'SDSSJ140557.10+154940.5_Blouin2020': (7055.0, 8.073),
    'SDSSJ143007.15-015129.5_Blouin2020': (6344.0, 8.076),
    'SDSSJ133905.98+264322.9_Blouin2020': (6452.0, 8.077),
    'SDSSJ150228.71+374452.9_Blouin2020': (5525.0, 8.084),
    'SDSSJ020128.66+201521.8_Blouin2020': (6250.0, 8.091),
    'SDSSJ123826.93+214937.7_Blouin2020': (5437.0, 8.097),
    'SDSSJ013504.11+130240.4_Blouin2020': (6013.0, 8.097),
    'SDSSJ164939.23+223807.2_Blouin2020': (5261.0, 8.101),
    'SDSSJ160429.80+183035.4_Blouin2020': (6421.0, 8.1),
    'SDSSJ123415.21+520808.1_Blouin2020': (7627.0, 8.098),
    'SDSSJ154510.31+523618.4_Blouin2020': (6068.0, 8.107),
    'SDSSJ085100.23+154301.6_Blouin2020': (6284.0, 8.12),
    'HS2253+8023Model1_Klein2011': (14000.0, 8.1),
    'SDSSJ001052.56-043014.3_Blouin2020': (6903.0, 8.122),
    'WD1145+017_Fortin-Archambault2020': (14500.0, 8.11),
    'SDSSJ104319.84+351641.6_Blouin2020': (6069.0, 8.135),
    'SDSSJ110234.21+021459.2_Blouin2020': (5699.0, 8.137),
    'J0956_Swan2023': (8720.0, 8.13),
    'SDSSJ234048.74+081753.3_Blouin2020': (6151.0, 8.141),
    'SDSSJ151835.63+050627.5_Blouin2020': (5187.0, 8.161),
    'HE1349-2305_Melis2012': (18173.0, 8.13),
    'SDSSJ105826.60+314358.3_Blouin2020': (7173.0, 8.159),
    'SDSSJ212312.20+001653.5_Blouin2020': (5463.0, 8.166),
    'SDSSJ074450.24+270132.6_Blouin2020': (7829.0, 8.16),
    'Ton345_Wilson2015': (19780.0, 8.18),
    'SDSSJ144804.49+104709.1_Blouin2020': (6331.0, 8.173),
    'SDSSJ153505.75+124744.2_Blouin2020': (5950.0, 8.175),
    'SDSSJ095645.15+591240.6_Blouin2020': (8843.0, 8.168),
    'SDSSJ120548.97+353642.4_Blouin2020': (6000.0, 8.178),
    'SDSSJ104046.48+240759.5_Blouin2020': (6023.0, 8.181),
    'SDSSJ014300.52+011356.8_Blouin2020': (7229.0, 8.178),
    'SDSSJ134520.99+115357.6_Blouin2020': (6000.0, 8.19),
    'SDSSJ115015.74+492843.7_Blouin2020': (7417.0, 8.189),
    'SDSSJ134711.47+141528.0_Blouin2020': (6520.0, 8.194),
    'GaiaJ0644-0352(Spec;Opt)_Rogers2023': (18350.0, 8.18),
    'SDSSJ074153.45+314620.4_Blouin2020': (5974.0, 8.205),
    'SDSSJ152449.58+404938.1_Blouin2020': (6203.0, 8.203),
    'SDSSJ140410.72+362056.8_Blouin2020': (6284.0, 8.216),
    'GD61_Farihi2011': (17820.0, 8.2),
    'SDSSJ004451.69+041819.2_Blouin2020': (6104.0, 8.22),
    'SDSSJ073635.22+411828.2_Blouin2020': (5010.0, 8.222),
    'GD362Updated_Xu2013': (10540.0, 8.24),
    'GD362_Xu2013': (10540.0, 8.24),
    'SDSSJ155429.01+173545.9_Blouin2020': (6847.0, 8.231),
    'SDSSJ083858.56+232252.9_Blouin2020': (6048.0, 8.244),
    'SDSSJ100609.16+175221.3_Blouin2020': (5605.0, 8.251),
    'SDSSJ000215.65+320914.1_Blouin2020': (6466.0, 8.257),
    'SDSSJ1330+6435_Kaiser2021': (4310.0, 8.26),
    'SDSSJ015008.55+135433.9_Blouin2020': (6953.0, 8.262),
    'SDSSJ103352.89+180935.3_Blouin2020': (6147.0, 8.267),
    'SDSSJ154201.75+465020.2_Blouin2020': (6130.0, 8.269),
    'SDSSJ223507.65-005607.7_Blouin2020': (6514.0, 8.288),
    'SDSSJ084300.23+561452.8_Blouin2020': (6624.0, 8.292),
    'SDSSJ131900.19+364149.8_Blouin2020': (7464.0, 8.294),
    'SDSSJ104130.65+343240.7_Blouin2020': (7728.0, 8.301),
    'GD424(WHT)_Izquierdo2021': (16560.0, 8.25),
    'GD424(Keck)_Izquierdo2021': (16560.0, 8.25),
    'PG1225-079Model3_Klein2011': (11100.0, 8.3),
    'SDSSJ121837.12+002304.0_Blouin2020': (6500.0, 8.314),
    'SDSSJ102414.83+453109.9_Blouin2020': (6339.0, 8.356),
    'SDSSJ141140.27+341039.4_Blouin2020': (5500.0, 8.363),
    'SDSSJ131420.49+374806.5_Blouin2020': (6201.0, 8.379),
    'SDSSJ084239.85+153628.8_Blouin2020': (5966.0, 8.385),
    'HS2253+8023Model2_Klein2011': (14400.0, 8.4),
    'SDSSJ101750.24+241911.6_Blouin2020': (6772.0, 8.409),
    'SDSSJ073842.56+183509.06_Dufour2012': (13950.0, 8.4),
    'SDSSJ093916.04+413612.9_Blouin2020': (6321.0, 8.422),
    'J2050_Swan2023': (10200.0, 8.44),
    'SDSSJ012620.48+253433.6_Blouin2020': (5588.0, 8.456),
    'SDSSJ004634.22+271737.6_Blouin2020': (8053.0, 8.465),
    'SDSSJ135054.02+105808.0_Blouin2020': (5176.0, 8.505),
    'SDSSJ0738+1835_Dufour2010': (13600.0, 8.5),
    'WDJ2317+1830_Hollands2021': (4210.0, 8.64),
    'SDSSJ121106.43+232623.0_Blouin2020': (6609.0, 8.663),
    'HS2253+8023Model3_Klein2011': (14800.0, 8.7),
    'WD1425+540(Model1)_Xu2017': (14490.0, 7.95),
    'WD1425+540(Model2)_Xu2017': (14490.0, 7.95),
    'HS2253+8023_Friedrich1999': (14700.0, 8.0),
    'L745-46A_Koester+Wolff2000': (7500.0, 8.0),
    'WDJ1333-6751_OBrien2023': (5640.0, 8.17),
    'WD1232+563_Badenas-Agusti2024': (11777.8, 8.24),
    'LP658-2_Blouin2018': (4430.0, 7.967),
    'SDSSJ080440.63+223948.6_Blouin2018': (4970.0, 7.98),
    'Ross640_Blouin2018': (8070.0, 7.923),
    'SDSSJ1330+6435_Blouin2019': (4310.0, 8.26),
    'WDJ2356-209_Blouin2019': (4040.0, 7.98),
    'SDSSJ000418.68+081929.9_Blouin2020': (5843.0, 8.0),
    'SDSSJ001309.43+110949.0_Blouin2020': (6090.0, 8.0),
    'SDSSJ004757.07+162836.5_Blouin2020': (6300.0, 8.0),
    'SDSSJ005247.16+184649.5_Blouin2020': (5305.0, 8.0),
    'SDSSJ005304.15+311555.8_Blouin2020': (6548.0, 8.0),
    'SDSSJ005649.27+245335.4_Blouin2020': (5061.0, 8.0),
    'SDSSJ011421.17+350547.1_Blouin2020': (6209.0, 8.0),
    'SDSSJ011646.10+205001.9_Blouin2020': (6245.0, 8.0),
    'SDSSJ014441.64+030536.2_Blouin2020': (6753.0, 8.0),
    'SDSSJ014834.00-011235.9_Blouin2020': (6760.0, 8.0),
    'SDSSJ020809.93-054258.3_Blouin2020': (6085.0, 8.0),
    'SDSSJ023407.46-051028.1_Blouin2020': (6601.0, 8.0),
    'SDSSJ073953.46+311204.2_Blouin2020': (5221.0, 8.0),
    'SDSSJ074444.03+440845.6_Blouin2020': (6612.0, 8.0),
    'SDSSJ075853.46+101347.4_Blouin2020': (5585.0, 8.0),
    'SDSSJ080003.90+224210.1_Blouin2020': (6049.0, 8.0),
    'SDSSJ080626.70+305555.7_Blouin2020': (7017.0, 8.0),
    'SDSSJ081606.19+233030.1_Blouin2020': (7642.0, 8.0),
    'SDSSJ081828.12+124717.2_Blouin2020': (6895.0, 8.0),
    'SDSSJ090652.19+114149.9_Blouin2020': (6556.0, 8.0),
    'SDSSJ090814.52+411918.3_Blouin2020': (6746.0, 8.0),
    'SDSSJ091322.38+262752.0_Blouin2020': (5252.0, 8.161),
    'SDSSJ091356.05+412728.6_Blouin2020': (6010.0, 8.0),
    'SDSSJ092450.03+430136.4_Blouin2020': (5500.0, 8.0),
    'SDSSJ092523.10+313019.0_Blouin2020': (6050.0, 8.0),
    'SDSSJ093320.01+633441.2_Blouin2020': (6337.0, 8.0),
    'SDSSJ094648.94+202423.2_Blouin2020': (7540.0, 8.0),
    'SDSSJ100537.43+224403.1_Blouin2020': (6165.0, 8.0),
    'SDSSJ101451.15+282701.6_Blouin2020': (6269.0, 8.0),
    'SDSSJ101924.73+353527.6_Blouin2020': (6224.0, 8.0),
    'SDSSJ101959.51+204553.4_Blouin2020': (5515.0, 8.0),
    'SDSSJ102438.05+101410.5_Blouin2020': (6105.0, 8.0),
    'SDSSJ103205.15+133833.4_Blouin2020': (5479.0, 8.0),
    'SDSSJ103809.19-003622.5_Blouin2020': (7996.0, 8.0),
    'SDSSJ103839.01+043223.8_Blouin2020': (6363.0, 8.0),
    'SDSSJ110216.09+282730.7_Blouin2020': (6320.0, 8.0),
    'SDSSJ110556.17+022849.0_Blouin2020': (5842.0, 8.0),
    'SDSSJ113209.59+332353.0_Blouin2020': (6062.0, 8.0),
    'SDSSJ113410.85+154245.9_Blouin2020': (6806.0, 8.0),
    'SDSSJ114709.09+542940.8_Blouin2020': (5000.0, 8.0),
    'SDSSJ115207.15+510126.2_Blouin2020': (4790.0, 8.0),
    'SDSSJ115748.35+613845.9_Blouin2020': (6607.0, 8.0),
    'SDSSJ115822.32+471214.9_Blouin2020': (7840.0, 8.0),
    'SDSSJ115844.35+544837.5_Blouin2020': (6213.0, 8.0),
    'SDSSJ122035.76+092948.1_Blouin2020': (6677.0, 8.0),
    'SDSSJ124547.11+082231.4_Blouin2020': (6074.0, 8.0),
    'SDSSJ125710.13+323848.5_Blouin2020': (5376.0, 8.0),
    'SDSSJ125720.87-031025.1_Blouin2020': (6269.0, 8.0),
    'SDSSJ130826.36+095724.0_Blouin2020': (7692.0, 8.0),
    'SDSSJ130830.03+025844.5_Blouin2020': (6003.0, 8.0),
    'SDSSJ132005.53+020419.0_Blouin2020': (7356.0, 8.0),
    'SDSSJ132144.04-023751.4_Blouin2020': (5592.0, 8.0),
    'SDSSJ134050.31+270219.0_Blouin2020': (8413.0, 8.0),
    'SDSSJ134203.60+181332.8_Blouin2020': (5524.0, 8.0),
    'SDSSJ140525.20+254212.4_Blouin2020': (5890.0, 8.0),
    'SDSSJ144301.55+583301.6_Blouin2020': (7061.0, 8.0),
    'SDSSJ144354.13+301413.3_Blouin2020': (6955.0, 8.0),
    'SDSSJ144535.03+091340.4_Blouin2020': (7035.0, 8.0),
    'SDSSJ150028.02+231554.0_Blouin2020': (6630.0, 8.0),
    'SDSSJ150739.03+403408.9_Blouin2020': (7304.0, 8.0),
    'SDSSJ153407.58+124254.4_Blouin2020': (6197.0, 8.0),
    'SDSSJ153745.52+360818.6_Blouin2020': (5519.0, 8.0),
    'SDSSJ154022.80+535239.8_Blouin2020': (6500.0, 8.0),
    'SDSSJ154349.79+202442.9_Blouin2020': (6206.0, 8.0),
    'SDSSJ154933.23+190646.7_Blouin2020': (6246.0, 8.0),
    'SDSSJ161026.10+400619.7_Blouin2020': (6552.0, 8.0),
    'SDSSJ161248.17+353434.8_Blouin2020': (7181.0, 8.0),
    'SDSSJ162408.57+331019.0_Blouin2020': (6654.0, 8.0),
    'SDSSJ162612.73+330308.2_Blouin2020': (6715.0, 8.0),
    'SDSSJ170638.11+254111.7_Blouin2020': (5813.0, 8.0),
    'SDSSJ211045.34+051214.8_Blouin2020': (5828.0, 8.0),
    'SDSSJ223014.70+190514.4_Blouin2020': (5631.0, 8.0),
    'SDSSJ223811.10+021352.9_Blouin2020': (6986.0, 8.0),
    'SDSSJ223815.97-011336.9_Blouin2020': (6228.0, 8.0),
    'SDSSJ231937.39+301848.4_Blouin2020': (7478.0, 8.0),
    'SDSSJ232833.31+083028.4_Blouin2020': (5566.0, 8.0),
    'SDSSJ233054.31+280517.4_Blouin2020': (6344.0, 8.0),
    'SDSSJ233320.38+105830.2_Blouin2020': (6515.0, 8.0),
    'SDSSJ235224.26+192247.3_Blouin2020': (6067.0, 8.0),
    'SDSSJ235249.13+334439.2_Blouin2020': (7200.0, 8.0),
    'SDSSJ235715.03+234848.9_Blouin2020': (6117.0, 8.0),
    'G270-124_Desharnais2008': (20820.0, 7.94),
    'GD408_Desharnais2008': (14660.0, 8.31),
    'GD61_Desharnais2008': (17280.0, 8.2),
    'GD378_Desharnais2008': (16600.0, 8.03),
    'SDSSJ1734+6052_Doyle2023': (16340.0, 8.04),
    'SDSSJ2248+2632_Doyle2023': (17370.0, 8.02),
    'EC22211-?2525_Doyle2023': (14740.0, 7.89),
    'GaiaJ0218+3625_Doyle2023': (14700.0, 7.86),
    'WD1415+234_Doyle2023': (17300.0, 8.17),
    'WD1244+498_Doyle2023': (15150.0, 7.97),
    'GaiaJ1922+4709_Doyle2023': (15500.0, 7.95),
    'SDSSJ1248+1005_Doyle2023': (15180.0, 8.11),
    'GD61_Farihi2013': (17280.0, 8.2),
    'GD40_Friedrich1999': (15150.0, 8.0),
    'HE0446-2531_Friedrich2000': (12600.0, 7.77),
    'G241-6_Jura2012': (15300.0, 8.0),
    'GD40_Jura2012': (15300.0, 8.0),
    'Ton345_Jura2015': (18700.0, 8.0),
    'Ross640_Koester+Wolff2000': (8500.0, 8.0),
    'HS0146+1847_Koester2005': (11500.0, 8.0),
    'PG1115+158_Koester2014b': (25000.0, 7.91),
    'WDJ1241-2434_OBrien2023': (6310.0, 8.13),
    'WDJ1057-0413_OBrien2023': (6500.0, 8.03),
    'WDJ2236-5548_OBrien2023': (5350.0, 8.17),
    'WDJ1927-0355_OBrien2024': (6540.0, 7.99),
    'WDJ2141-3300_OBrien2024': (6870.0, 7.96),
    'WD2138-332_Subasavage2007': (7188.0, 8.0),
    'WDJ191246.12+024239.11_Tremblay2020': (6050.0, 8.15),
    'WDJ183352.68+321757.25_Tremblay2020': (7650.0, 8.05),
    'WDJ235750.73+194905.90_Tremblay2020': (5700.0, 7.95),
    'GD378_Wolff2002': (17000.0, 7.9),
    'GD408_Wolff2002': (13750.0, 8.0),
    'vMa2_Wolff2002': (5700.0, 7.9),
    'GD303_Wolff2002': (18000.0, 7.8),
    'K789-37_Wolff2002': (10500.0, 8.0),
    'L119-34_Wolff2002': (9200.0, 8.0),
    'WD1551+175_Xu2019': (14756.0, 8.02),
    'WD1232+563_Xu2019': (11787.0, 8.3),
    'WD2207+121_Xu2019': (14752.0, 7.97),
    'G241-6_Zuckerman2010': (15300.0, 8.0),
    'WD1145+017(Opt)_LeBourdais2024': (15500.0, 8.19),
    'WD1145+017(UV)_LeBourdais2024': (15500.0, 8.19),
    'SDSSJ082019.49+253035.3_Aguilera-Gomez2024': (11900.0, 8.06),
    'WDJ1644-0449_Kaiser2024': (3830.0, 7.85),
    'SDSSJ1330+6435_Kaiser2024': (4310.0, 8.26),
    'WDJ1824+1213_Kaiser2024': (3540.0, 7.53),
    'WDJ2317+1830_Kaiser2024': (4430.0, 8.74),
    'LHS2534_Kaiser2024': (5020.0, 8.1)
}

def get_logg_values(grid_steps=100, grid_range=(None, None), timescale_types=list()): # Bedard goes from 7.5 to 9, Koester from 7.5 to 8.5
    # grid_range = min logg, max logg
    min_logg = grid_range[0]
    max_logg = grid_range[1]
    if min_logg is None:
        min_logg = max([logg_ranges[timescale_type][0] for timescale_type in timescale_types])
    if max_logg is None:
        max_logg = min([logg_ranges[timescale_type][1] for timescale_type in timescale_types])
    print(min_logg)
    print(max_logg)
    return np.linspace(min_logg, max_logg, grid_steps + 1)

def get_Teff_values(grid_steps=100, grid_range=(None, None), timescale_types=list()): # Bedard goes from 11000 to 20000 (actually it varies! 5000 to 30000 by default, or 11000 to 20000 for 3O), Koester from 5000 to 20000
    min_teff = grid_range[0]
    max_teff = grid_range[1]
    if min_teff is None:
        min_teff = max([teff_ranges[timescale_type][0] for timescale_type in timescale_types])
    if max_teff is None:
        max_teff = min([teff_ranges[timescale_type][1] for timescale_type in timescale_types])
    print(min_teff)
    print(max_teff)
    return np.linspace(min_teff, max_teff, grid_steps + 1)

def get_CaHe_values():
    return [-6.5, -7, -7.5, -8, -8.5, -9, -9.5, -10, -10.5, -11, -11.5, -12, -12.5, -13, -13.5, -14, -14.5, -15, -15.5, -16]

def get_variable_vals(variable):
    if variable == 'logg':
        return get_logg_values()
    if variable == 'Teff':
        return get_Teff_values()
    if variable == 'CaHe':
        return get_CaHe_values()
    return None

def get_elements_for_2D_plot():
    #return [ci.Element.Ca, ci.Element.Fe, ci.Element.Mg, ci.Element.O]
    #return [ci.Element.C, ci.Element.Si]
    return [ci.Element.Ca, ci.Element.Fe]

def get_elements_for_model_comparison(): # needs to be a subset of get_elements_for_2D_plot().
    return get_elements_for_2D_plot()
    #return [ci.Element.Ca, ci.Element.Mg]
    #return [ci.Element.C, ci.Element.Si]
    #return [ci.Element.Fe, ci.Element.Ca]

def generate_timescales():
    timescale_interpolator = ti.TimescaleInterpolator()
    logg_vals = get_logg_values()
    Teff_vals = get_Teff_values()
    CaHe_vals = get_CaHe_values()
    timescale_vals = {
        'H': dict(),
        'He': dict()
    }
    elements_we_care_about = get_elements_for_2D_plot()
    for element in elements_we_care_about:
        timescale_vals['H'][element] = np.zeros((len(logg_vals), len(Teff_vals)))
        timescale_vals['He'][element] = np.zeros((len(logg_vals), len(Teff_vals), len(CaHe_vals)))
    for i in range(len(logg_vals)):
        logg = logg_vals[i]
        for j in range(len(Teff_vals)):
            Teff = Teff_vals[j]
            H_timescales = timescale_interpolator.get_wd_timescales('H', logg, Teff, None)
            for element in elements_we_care_about:
                timescale_vals['H'][element][i,j] = H_timescales[element]
            for k in range(len(CaHe_vals)):
                CaHe = CaHe_vals[k]
                He_timescales = timescale_interpolator.get_wd_timescales('He', logg, Teff, CaHe)
                for element in elements_we_care_about:
                    timescale_vals['He'][element][i,j,k] = He_timescales[element]
    return timescale_vals

def extract_2D_vals_to_plot(timescale_vals, HorHe, logg, Teff, CaHe, elements_to_extract):
    # Whichever of these variables is None is the one we need to plot against
    num_variables = 0
    variable = None
    if logg is None:
        num_variables += 1
        variable = 'logg'
    if Teff is None:
        num_variables += 1
        variable = 'Teff'
    if CaHe is None:
        num_variables += 1
        variable = 'CaHe'
    if num_variables != 1:
        print('Error: Exactly 1 of logg, Teff and CaHe must be None. (This indicates the variable to plot against)')
    assert num_variables == 1
    if HorHe == 'H' and variable == 'CaHe':
        print('Error: H does not vary with CaHe')
        return None
    assert logg is None or logg in get_logg_values()
    assert Teff is None or Teff in get_Teff_values()
    assert CaHe is None or CaHe in get_CaHe_values()
    try:
        logg_index = get_logg_values().index(logg)
    except ValueError:
        logg_index = None
    try:
        Teff_index = get_Teff_values().index(Teff)
    except ValueError:
        Teff_index = None
    try:
        CaHe_index = get_CaHe_values().index(CaHe)
    except ValueError:
        CaHe_index = None
    to_plot = dict()
    for element in elements_to_extract:
        to_plot[element] = list()
        to_slice_into = timescale_vals[HorHe][element]
        if HorHe == 'H':
            if variable == 'logg':
                to_plot[element] = to_slice_into[:,Teff_index]
            if variable == 'Teff':
                to_plot[element] = to_slice_into[logg_index,:]
        elif HorHe == 'He':
            if variable == 'logg':
                to_plot[element] = to_slice_into[:,Teff_index,CaHe_index]
            if variable == 'Teff':
                to_plot[element] = to_slice_into[logg_index,:,CaHe_index]
            if variable == 'CaHe':
                to_plot[element] = to_slice_into[logg_index,Teff_index,:]
        else:
            print('Error: HorHe must be H or He')
            return None
    return to_plot, variable

def plot_2D_timescale(timescale_vals, HorHe, logg, Teff, CaHe):
    to_plot, variable = extract_2D_vals_to_plot(timescale_vals, HorHe, logg, Teff, CaHe, get_elements_for_2D_plot())
    graph_factory = gf.GraphFactory()
    graph_factory.plot_2D_timescales(to_plot, variable, get_variable_vals(variable), HorHe, logg, Teff, CaHe)

def plot_model_comparison(timescale_vals, HorHe, logg, Teff, CaHe):
    to_plot_new, variable = extract_2D_vals_to_plot(timescale_vals, HorHe, logg, Teff, CaHe, get_elements_for_model_comparison())
    if variable == 'logg':
        # These are at (or near) Teff = 6250 according to Simon (but we take the timescales themselves from Hollands) (also excluding those where logg was just assumed to be 8 by default)
        # System      logg     Ca       Fe
        #J0108-0537 7.256 1230268.771 1230268.771
        #J0201+2015 8.091 1348962.883 1348962.883
        #J0851+1543 8.12 1202264.435 1230268.771
        #J0939+4136 8.422 1122018.454 1148153.621
        #J0948+3008 7.21 1513561.248 1380384.265
        #J1024+4531 8.356 1258925.412 1258925.412
        #J1230+3143 7.676 1584893.192 1548816.619
        #J1314+3748 8.379 1071519.305 1071519.305
        #J1336+3547 8.057 1318256.739 1348962.883
        #J1356+2416 8.01 1584893.192 1412537.545
        #J1404+3620 8.216 1584893.192 1348962.883
        #J1430-0151 8.076 831763.7711 870963.59
        #J1448+1047 8.173 1445439.771 1445439.771
        #J1524+4049 8.203 1230268.771 1230268.771
        #J2340+0817 8.141 1288249.552 1174897.555
        old_vals_dict = {
            7.256:( 1230268.771, 1230268.771),
            8.091:( 1348962.883, 1348962.883),
            8.12: (1202264.435, 1230268.771) ,
            8.422:( 1122018.454, 1148153.621),
            7.21: (1513561.248, 1380384.265) ,
            8.356: (1258925.412, 1258925.412),
            7.676: (1584893.192, 1548816.619),
            8.379: (1071519.305, 1071519.305),
            8.057: (1318256.739, 1348962.883),
            8.01: (1584893.192 ,1412537.545 ),
            8.216: (1584893.192, 1348962.883),
            8.076: (831763.7711, 870963.59  ),
            8.173: (1445439.771, 1445439.771),
            8.203: (1230268.771, 1230268.771),
            8.141: (1288249.552, 1174897.555)
        }
        x_vals_old = list(old_vals_dict.keys())
        x_vals_old.sort()
        ca_vals = list()
        fe_vals = list()
        for x_val in x_vals_old:
            ca_vals.append(old_vals_dict[x_val][0])
            fe_vals.append(old_vals_dict[x_val][1])
        to_plot_old = {
            ci.Element.Ca: np.array(ca_vals),
            ci.Element.Fe: np.array(fe_vals)
        }
        print(x_vals_old)
        print(to_plot_old)
    elif variable == 'Teff':
        # These are at (or near) logg = 8 according to Simon (but we take the timescales themselves from Hollands) (also excluding those where logg was just assumed to be 8 by default)
        # System   Teff     Ca       Fe
        #J0158-0942 6115 2089296.131 1513561.248
        #J0252+0054 7478 1479108.388 1513561.248
        #J0807+4930 5172 851138.0382 851138.0382
        #J0842+1406 7075 1380384.265 1445439.771
        #J0901+0752 7263 954992.586 1000000
        #J0937+5228 6660 1288249.552 1318256.739
        #J1017+3447 6089 1819700.859 1513561.248
        #J1144+1218 5320 1737800.829 1288249.552
        #J1158+0454 5344 1023292.992 977237.221
        #J1158+1845 6696 1230268.771 1318256.739
        #J1158+5942 6046 1318256.739 1288249.552
        #J1224+2838 4991 2884031.503 1548816.619
        #J1254+3551 6417 1513561.248 1548816.619
        #J1259+3112 5664 2187761.624 1548816.619
        #J1356+2416 6173 1584893.192 1412537.545
        #J1356+0236 7662 1380384.265 1479108.388
        #J1401+3659 5931 2398832.919 1659586.907
        #J1421+1843 6517 1096478.196 1148153.621
        #J1616+3303 6491 1148153.621 1174897.555
        #J2109-0039 6132 1230268.771 1230268.771
        #J2231+0906 5679 2454708.916 1659586.907
        old_vals_dict = {
            6115: ( 2089296.131, 1513561.248),
            7478: ( 1479108.388, 1513561.248),
            5172: ( 851138.0382, 851138.0382),
            7075: ( 1380384.265, 1445439.771),
            7263: ( 954992.586, 1000000     ),
            6660: ( 1288249.552, 1318256.739),
            6089: ( 1819700.859, 1513561.248),
            5320: ( 1737800.829, 1288249.552),
            5344: ( 1023292.992, 977237.221 ),
            6696: ( 1230268.771, 1318256.739),
            6046: ( 1318256.739, 1288249.552),
            4991: ( 2884031.503, 1548816.619),
            6417: ( 1513561.248, 1548816.619),
            5664: ( 2187761.624, 1548816.619),
            6173: ( 1584893.192, 1412537.545),
            7662: ( 1380384.265, 1479108.388),
            5931: ( 2398832.919, 1659586.907),
            6517: ( 1096478.196, 1148153.621),
            6491: ( 1148153.621, 1174897.555),
            6132: ( 1230268.771, 1230268.771),
            5679: ( 2454708.916, 1659586.907)
        }
        x_vals_old = list(old_vals_dict.keys())
        x_vals_old.sort()
        ca_vals = list()
        fe_vals = list()
        for x_val in x_vals_old:
            ca_vals.append(old_vals_dict[x_val][0])
            fe_vals.append(old_vals_dict[x_val][1])
        to_plot_old = {
            ci.Element.Ca: np.array(ca_vals),
            ci.Element.Fe: np.array(fe_vals)
        }
        print(x_vals_old)
        print(to_plot_old)
    else:
        print('Old timescales only varied with Teff and logg')
        x_vals_old = None
        to_plot_old = None
    if HorHe != 'He':
        print('Old timescales only present for He')
        x_vals_old = None
        to_plot_old = None
    graph_factory = gf.GraphFactory()
    graph_factory.plot_timescale_model_comparison(to_plot_new, to_plot_old, variable, get_variable_vals(variable), x_vals_old, HorHe, logg, Teff, CaHe, get_elements_for_model_comparison())

def extract_timescale_ratios(element1, element2, Hx, logg_vals, Teff_vals):
    timescale_interpolator = ti.TimescaleInterpolator()
    #CaHe_vals = get_CaHe_values()
    timescale_vals = {tt: {element1: np.zeros((len(logg_vals), len(Teff_vals))), element2: np.zeros((len(logg_vals), len(Teff_vals))), 'SS': np.zeros((len(logg_vals), len(Teff_vals))), 'Dec': np.zeros((len(logg_vals), len(Teff_vals)))} for tt in ti.TimescaleType}
    for i in range(len(logg_vals)):
        logg = logg_vals[i]
        for j in range(len(Teff_vals)):
            Teff = Teff_vals[j]
            print('Getting timescales for logg = ' + str(logg) + ', Teff = ' + str(Teff))
            timescales = timescale_interpolator.get_wd_timescales(Hx, logg, Teff, -15)
            for tt in timescales:
                if timescales[tt] is not None:
                    t1 = timescales[tt][element1]
                    t2 = timescales[tt][element2]
                    timescale_vals[tt][element1][i,j] = t1
                    timescale_vals[tt][element2][i,j] = t2
                    timescale_vals[tt]['SS'][i,j] = t1/t2
                    t_e = 10*max(t1, t2) # Let's say accretion lasts for 10 sinking timescales of elements 1 and 2, so that the system has reached steady state
                    t_obs = 15*max(t1, t2) # ...And let's say you observe it a few sinking timescales after that
                    timescale_vals[tt]['Dec'][i,j] = t1/t2 * np.exp(t_obs/t2 - t_obs/t1) * ((np.exp(t_e/t1) - 1)/(np.exp(t_e/t2) - 1))
    return timescale_vals

def plot_timescale_ratios(element1, element2, timescale_type_pairs, Hx, correction_types=['SS'], grid_steps=100, grid_range=((None, None), (None, None))):
    #grid_range is a 2-tuple of 2-tuples: ( (min Teff, max Teff) , (min logg, max logg) )

    graph_fac = gf.GraphFactory()
    include_wd_markers = True
    #reference_systems = {
    #    #'GD56': (15270, 8.09),
    #    #'WD0107-192': (15440, 7.95),
    #    #'GD133': (12600, 8.1),
    #    #'WD2105-820': (10890, 8.41),
    #    #'WD1145+288': (12140, 8.14),
    #    #'GaiaJ0611': (16530, 7.81),
    #    #'GaiaJ0611-uv-P': (16530, 7.81),
    #    #'Synthetic System': (14500, 7.9),
    #    #'WD2058+181': (17308, 7.92)
    #}
    reference_systems = dict()
    full_sample = dict()
    if include_wd_markers:
        if Hx == ci.Element.H:
            #reference_systems = {
            #    #DA Overshoot targets:
            #    'Gaia J0611-6931 (Phot; Opt) Rogers 2023': (16530, 7.81),
            #    'Gaia J0611-6931 (Phot; UV) Rogers 2023': (16530, 7.81),
            #    'GD133 Xu 2014': (12600, 8.1),
            #    'WD0145+234 Melis 2020': (12720, 8.1),
            #    'G29-38 Xu 2014': (11820, 8.4),
            #    'G 29-38 Koester 1997':	(11600, 8.05),
            #    'GD 56 Xu 2019': (15270, 8.09),
            #    'WD 1150-153 Xu 2019': (12640, 8.22),
            #    'WD 1145+288 Xu 2019': (12140, 8.14),
            #    'WD 0310-688 Limbach 2024':	(15865, 8.076)
            #}
            #reference_systems = {
            #    #DA BVK targets:
            #    'NLTT 1675 Kawka 2012': (6020, 8.04),
            #    'NLTT 43806 Zuckerman 2011': (5900, 8),
            #    'WD 1124-293 Steele 2021': (9367, 7.99),
            #    'WD2115-560 Swan 2019': (9600, 7.97),
            #    'NLTT 25792 Vennes + Kawka 2013': (7903, 8.04),
            #    'GD133 Xu 2014': (12600, 8.1),
            #    'WD0354+463 Vennes + Kawka 2013': (8240, 7.96),
            #    'G 29-38 Koester 1997': (11600, 8.05),
            #    'WDJ113444.64+610826.68 Tremblay 2020': (7590, 7.96),
            #    'WD1455+298 Vennes + Kawka 2013': (7383, 7.97),
            #    'G74-7 Vennes + Kawka 2013': (7306, 8.06),
            #    'WD1257+278 Vennes + Kawka 2013': (8609, 8.24),
            #    'G 166-58 Xu 2019': (7390, 7.99),
            #    'WD 1145+288 Xu 2019': (12140, 8.14),
            #    'G149-28 Zuckerman 2011': (8600, 8.1)
            #}
            name_of_sample = 'BVK_DA'
            #name_of_sample = 'OVERSHOOT_DA'
            reference_systems = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pss.pick_out_sample(name_of_sample)}
            full_sample = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pss.pick_out_all_das()}
        elif Hx == ci.Element.He:
            #reference_systems = {
            #    #DB Overshoot candidates:
            #    'WD 1622+587 (Spec; Opt) Rogers 2023': (23430, 7.9),
            #    'J0848 Swan 2023': (11700, 7.89),
            #    'J1618 Swan 2023': (9870, 7.94),
            #    'GALEXJ2339 Klein 2021': (13735, 7.93),
            #    'WD0446-255 Swan 2019': (10120, 8),
            #    'GD 16 Gentile Fusillo 2017': (11000, 8),
            #    'SDSS J0956+5912 (No Overshoot) Hollands 2022': (8100, 8.02),
            #    'WD2216-657 Swan 2019': (9190, 8.05),
            #    'WD0449-259 Swan 2019': (9850, 8.04),
            #    'J0956 Swan 2023': (8720, 8.13),
            #    'L 745-46A Koester and Wolff 2000': (7500, 8),
            #    'K 789-37 Wolff 2002': (10500, 8),
            #    'L 119-34 Wolff 2002': (9200, 8)
            #}
            #name_of_sample = 'BVK_DB'
            name_of_sample = 'OVERSHOOT_DB'
            reference_systems = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pss.pick_out_sample(name_of_sample)}
            full_sample = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pss.pick_out_all_dbs()}
        else:
            print(Hx)
            raise
    all_timescale_types = list()
    for timescale_type1, timescale_type2 in timescale_type_pairs:
        if timescale_type1 not in all_timescale_types:
            all_timescale_types.append(timescale_type1)
        if timescale_type2 not in all_timescale_types:
            all_timescale_types.append(timescale_type2)
    print(all_timescale_types)
    print([ti.teff_ranges[Hx][timescale_type][0] for timescale_type in all_timescale_types])
    print(min([ti.teff_ranges[Hx][timescale_type][0] for timescale_type in all_timescale_types]))
    min_teff = min([ti.teff_ranges[Hx][timescale_type][0] for timescale_type in all_timescale_types]) if grid_range[0][0] is None else grid_range[0][0]
    max_teff = max([ti.teff_ranges[Hx][timescale_type][1] for timescale_type in all_timescale_types]) if grid_range[0][1] is None else grid_range[0][1]
    min_logg = min([ti.logg_ranges[Hx][timescale_type][0] for timescale_type in all_timescale_types]) if grid_range[1][0] is None else grid_range[1][0]
    max_logg = max([ti.logg_ranges[Hx][timescale_type][1] for timescale_type in all_timescale_types]) if grid_range[1][1] is None else grid_range[1][1]
    full_logg_range = np.linspace(min_logg, max_logg, grid_steps + 1)
    full_teff_range = np.linspace(min_teff, max_teff, grid_steps + 1)
    print(full_logg_range)
    print(full_teff_range)
    data_pool = extract_timescale_ratios(element1, element2, Hx, full_logg_range, full_teff_range)
    overshoot_test_timescale_types = [ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot, ti.TimescaleType.BedardVariableOvershoot, ti.TimescaleType.Bedard3DOvershoot]
    for correction_type in correction_types:
        plot_dicts = list()
        overshoot_plot_dicts = list()
        for timescale_type1, timescale_type2 in timescale_type_pairs:
            #Teff_vals = get_Teff_values(grid_steps, grid_range[0], [timescale_type_1, timescale_type2])
            #logg_vals = get_logg_values(grid_steps, grid_range[1], [timescale_type_1, timescale_type2])
            grid_range_to_plot = (
                (   #Find the teff/logg range that is within BOTH of the interpolation limits (ie for both timescale types)
                    max([ti.teff_ranges[Hx][timescale_type][0] for timescale_type in [timescale_type1, timescale_type2]]) if grid_range[0][0] is None else grid_range[0][0],
                    min([ti.teff_ranges[Hx][timescale_type][1] for timescale_type in [timescale_type1, timescale_type2]]) if grid_range[0][1] is None else grid_range[0][1],
                ),
                (
                    max([ti.logg_ranges[Hx][timescale_type][0] for timescale_type in [timescale_type1, timescale_type2]]) if grid_range[1][0] is None else grid_range[1][0],
                    min([ti.logg_ranges[Hx][timescale_type][1] for timescale_type in [timescale_type1, timescale_type2]]) if grid_range[1][1] is None else grid_range[1][1],
                ),
            )
            data_to_plot = np.log10(data_pool[timescale_type1][correction_type]/data_pool[timescale_type2][correction_type])
            plot_dict = graph_fac.plot_timescale_type_ratios(element1, element2, timescale_type1, timescale_type2, Hx, correction_type, full_teff_range, full_logg_range, data_to_plot, reference_systems, full_sample, grid_range_to_plot)
            plot_dicts.append(plot_dict)
            if timescale_type1 in overshoot_test_timescale_types and timescale_type2 in overshoot_test_timescale_types:
                overshoot_plot_dicts.append(plot_dict)
        graph_fac.multipanelise(plot_dicts, math.ceil(len(plot_dicts)/2), 2, 'timescale_type_comparison_' + str(element1) + '_' + str(element2) + '_' + str(Hx) + '_' + correction_type + '.pdf', 16, 24, 0.07, 0, True, True)
        if len(overshoot_plot_dicts) >= 1:
            graph_fac.multipanelise(overshoot_plot_dicts, 1, len(overshoot_plot_dicts), 'timescale_type_overshoot_comparison_' + str(element1) + '_' + str(element2) + '_' + str(Hx) + '_' + correction_type + '.pdf', 16, 24, 0.07, 0, True, True)

def plot_timescale_comparison_benchmarked(element1, element2, timescale_type_pairs, timescale_benchmark, Hx, correction_types=['SS'], grid_steps=10, grid_range=((None, None), (None, None))):
    graph_fac = gf.GraphFactory()
    include_wd_markers = False
    reference_systems = dict()
    full_sample = dict()
    if include_wd_markers:
        if Hx == ci.Element.H:
            #name_of_sample = 'BVK_DA'
            name_of_sample = 'OVERSHOOT_DA'
            reference_systems = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pss.pick_out_sample(name_of_sample)}
            full_sample = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pss.pick_out_all_das()}
        elif Hx == ci.Element.He:
            #name_of_sample = 'BVK_DB'
            name_of_sample = 'OVERSHOOT_DB'
            reference_systems = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pss.pick_out_sample(name_of_sample)}
            full_sample = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pss.pick_out_all_dbs()}
        else:
            print(Hx)
            raise
    all_timescale_types = list()
    for timescale_type1, timescale_type2 in timescale_type_pairs:
        if timescale_type1 not in all_timescale_types:
            all_timescale_types.append(timescale_type1)
        if timescale_type2 not in all_timescale_types:
            all_timescale_types.append(timescale_type2)
    min_teff = min([ti.teff_ranges[Hx][timescale_type][0] for timescale_type in all_timescale_types]) if grid_range[0][0] is None else grid_range[0][0]
    max_teff = max([ti.teff_ranges[Hx][timescale_type][1] for timescale_type in all_timescale_types]) if grid_range[0][1] is None else grid_range[0][1]
    min_logg = min([ti.logg_ranges[Hx][timescale_type][0] for timescale_type in all_timescale_types]) if grid_range[1][0] is None else grid_range[1][0]
    max_logg = max([ti.logg_ranges[Hx][timescale_type][1] for timescale_type in all_timescale_types]) if grid_range[1][1] is None else grid_range[1][1]
    full_logg_range = np.linspace(min_logg, max_logg, grid_steps + 1)
    full_teff_range = np.linspace(min_teff, max_teff, grid_steps + 1)
    print(full_logg_range)
    print(full_teff_range)
    data_pool = extract_timescale_ratios(element1, element2, Hx, full_logg_range, full_teff_range)
    for correction_type in correction_types:
        plot_dicts = list()
        overshoot_plot_dicts = list()
        for timescale_type1, timescale_type2 in timescale_type_pairs:
            #Teff_vals = get_Teff_values(grid_steps, grid_range[0], [timescale_type_1, timescale_type2])
            #logg_vals = get_logg_values(grid_steps, grid_range[1], [timescale_type_1, timescale_type2])
            grid_range_to_plot = (
                (   #Find the teff/logg range that is within BOTH of the interpolation limits (ie for both timescale types)
                    max([ti.teff_ranges[Hx][timescale_type][0] for timescale_type in [timescale_type1, timescale_type2]]) if grid_range[0][0] is None else grid_range[0][0],
                    min([ti.teff_ranges[Hx][timescale_type][1] for timescale_type in [timescale_type1, timescale_type2]]) if grid_range[0][1] is None else grid_range[0][1],
                ),
                (
                    max([ti.logg_ranges[Hx][timescale_type][0] for timescale_type in [timescale_type1, timescale_type2]]) if grid_range[1][0] is None else grid_range[1][0],
                    min([ti.logg_ranges[Hx][timescale_type][1] for timescale_type in [timescale_type1, timescale_type2]]) if grid_range[1][1] is None else grid_range[1][1],
                ),
            )

            log_t1 = np.log10(data_pool[timescale_type1][correction_type])
            log_t2 = np.log10(data_pool[timescale_type2][correction_type])
            log_tb = np.log10(data_pool[timescale_benchmark][correction_type])
            diff_t1b = abs(log_t1 - log_tb)
            diff_t2b = abs(log_t2 - log_tb)
            #data_to_plot = np.where(diff_t1b < diff_t2b, diff_t2b-diff_t1b, diff_t2b-diff_t1b) #Introduce sign convention that t2 is negative, t1 is positive
            # I realised this just simplifies to this:
            data_to_plot = diff_t2b-diff_t1b
            plot_dict = graph_fac.plot_timescale_type_benchmark(element1, element2, timescale_type1, timescale_type2, timescale_benchmark, Hx, correction_type, full_teff_range, full_logg_range, data_to_plot, reference_systems, full_sample, grid_range_to_plot)
            plot_dicts.append(plot_dict)
        graph_fac.multipanelise(plot_dicts, math.ceil(len(plot_dicts)/2), 2, 'timescale_type_comparison_benchmark_' + str(element1) + '_' + str(element2) + '_' + str(Hx) + '_' + correction_type + '.pdf', 16, 24, 0.07, 0, True, True)


def extract_therm_factors(therm_factor_interpolator, logg_vals, Teff_vals, logMdot):
    #CaHe_vals = get_CaHe_values()
    therm_factors = np.zeros((len(logg_vals), len(Teff_vals)))
    for i in range(len(logg_vals)):
        logg = logg_vals[i]
        for j in range(len(Teff_vals)):
            Teff = Teff_vals[j]
            print('Getting therm factor for logg = ' + str(logg) + ', Teff = ' + str(Teff))
            therm_factors[i,j] = therm_factor_interpolator((logMdot, Teff, logg))
    return therm_factors

def plot_therm_factors():
    #import sys
    #import pwd_utils as pu
    #sys.path.append(pu.get_path_to_da_pollution_tables_dir())
    #
    #import InterpThermohaline as intt
    therm_factor_interpolator = thi.ThermohalineInterpolator()
    logg_vals = np.linspace(7.5, 8.5, 101)
    Teff_vals = np.linspace(6000, 30000, 101)
    logg_vals = [7.5, 7.6, 7.7, 7.8, 7.9, 8, 8.1, 8.2, 8.3, 8.4, 8.5]
    Teff_vals = [6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000]
    logMdots = [6, 7, 8, 9, 10, 11]
    reference_systems = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pv.pick_out_thermohaline_das()}
    full_sample = {wd.full_name(): (wd.get_teff().value, wd.get_logg().value) for wd in pv.pick_out_all_das()}
    graph_fac = gf.GraphFactory()
    for logMdot in logMdots:
        therm_factors = extract_therm_factors(therm_factor_interpolator, logg_vals, Teff_vals, logMdot)
        #reference_systems = {
        #    #'WD2058+181': (17308, 7.92),
        #    'HE0106-3253': (17350, 8.12),
        #    'SDSSJ1043+0855': (18330, 8.05),
        #    #'HS2229+2335': (18538, 7.92),
        #    'PG1015+161G': (19200, 8.22),
        #    'PG1015+161X': (20420, 8.11),
        #    #'WD1943+163': (19250, 7.87),
        #    #'WD1953-715': (19254, 8.12)
        #    #'WD2105-820': (10890, 8.41),
        #    #'WD1145+288': (12140, 8.14),
        #    #'WD2221-165': (10130, 8.15),
        #    #'WDJ1814-7354': (10090, 8),
        #    #'WD0307+077': (10230, 7.96),
        #    #'GD362': (10057, 7.95)
        #    #'G29-38': (11800, 8.4)
        #}
        plot_dict = graph_fac.plot_therm_factors(Teff_vals, logg_vals, logMdot, therm_factors, reference_systems, full_sample)

def main():
    plot_thermohaline_factors = False
    plot_timescale_pairs = False

    plot_timescale_comparison_benchmarked(ci.Element.Ca, ci.Element.Fe, [(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot)], ti.TimescaleType.Bedard3DOvershootPatched, ci.Element.H)
    plot_timescale_comparison_benchmarked(ci.Element.Mg, ci.Element.Fe, [(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot)], ti.TimescaleType.Bedard3DOvershootPatched, ci.Element.H)
    plot_timescale_comparison_benchmarked(ci.Element.Ca, ci.Element.Mg, [(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot)], ti.TimescaleType.Bedard3DOvershootPatched, ci.Element.H)
    plot_timescale_comparison_benchmarked(ci.Element.Ca, ci.Element.O, [(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot)], ti.TimescaleType.Bedard3DOvershootPatched, ci.Element.H)

    #plot_timescale_comparison_benchmarked(ci.Element.Ca, ci.Element.Fe, [(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot)], ti.TimescaleType.BedardVariableOvershoot, ci.Element.He)
    #plot_timescale_comparison_benchmarked(ci.Element.Mg, ci.Element.Fe, [(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot)], ti.TimescaleType.BedardVariableOvershoot, ci.Element.He)
    #plot_timescale_comparison_benchmarked(ci.Element.Ca, ci.Element.Mg, [(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot)], ti.TimescaleType.BedardVariableOvershoot, ci.Element.He)
    #plot_timescale_comparison_benchmarked(ci.Element.Ca, ci.Element.O, [(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.BedardOvershoot)], ti.TimescaleType.BedardVariableOvershoot, ci.Element.He)
    if plot_thermohaline_factors:
        plot_therm_factors()

    if plot_timescale_pairs:
        # Possible timescale types:
        #     KoesterNoOvershoot
        #     KoesterOvershoot
        #     BedardNoOvershoot
        #     BedardOvershoot
        #     BedardVariableOvershoot
        #     Bedard3DOvershoot
        #     MWDD
        timescale_type_pairs = [
            #(ti.TimescaleType.KoesterOvershoot, ti.TimescaleType.KoesterNoOvershoot),
            #(ti.TimescaleType.BedardOvershoot, ti.TimescaleType.BedardNoOvershoot),

            #(ti.TimescaleType.BedardVariableOvershoot, ti.TimescaleType.BedardNoOvershoot),
            #(ti.TimescaleType.Bedard3DOvershoot, ti.TimescaleType.BedardNoOvershoot)
            #(ti.TimescaleType.Bedard3DOvershootPatched, ti.TimescaleType.BedardNoOvershoot)

            (ti.TimescaleType.KoesterNoOvershoot, ti.TimescaleType.BedardNoOvershoot),
            #(ti.TimescaleType.KoesterOvershoot, ti.TimescaleType.BedardOvershoot),

            #(ti.TimescaleType.BedardNoOvershoot, ti.TimescaleType.MWDD)
        ]

        #for i, reference_element in enumerate(ci.usual_elements):
        #    j = i + 1
        #    while j < len(ci.usual_elements):
        #        el = ci.usual_elements[j]
        #        #plot_timescale_ratios(el, reference_element, timescale_type_pairs, ci.Element.H)
        #        plot_timescale_ratios(el, reference_element, timescale_type_pairs, ci.Element.He, ['SS'])
        #        j += 1

        #plot_timescale_ratios(ci.Element.C, ci.Element.Ca, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.O, ci.Element.Ca, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.Fe, ci.Element.Ca, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.Mg, ci.Element.Ca, timescale_type_pairs, ci.Element.H)

        plot_timescale_ratios(ci.Element.Ca, ci.Element.Mg, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.Fe, ci.Element.Mg, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.Si, ci.Element.Mg, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.O, ci.Element.Mg, timescale_type_pairs, ci.Element.H)

        #plot_timescale_ratios(ci.Element.Ca, ci.Element.Mg, timescale_type_pairs, ci.Element.He, ['SS'])
        #plot_timescale_ratios(ci.Element.Fe, ci.Element.Mg, timescale_type_pairs, ci.Element.He, ['SS'])
        #plot_timescale_ratios(ci.Element.Si, ci.Element.Mg, timescale_type_pairs, ci.Element.He, ['SS'])
        #plot_timescale_ratios(ci.Element.O, ci.Element.Mg, timescale_type_pairs,  ci.Element.He, ['SS'])

        #plot_timescale_ratios(ci.Element.Fe, ci.Element.Ca, timescale_type_pairs, ci.Element.He, ['SS'])

        #plot_timescale_ratios(ci.Element.Ca, ci.Element.Fe, timescale_type_pairs, ci.Element.He, ['SS'])

        #plot_timescale_ratios(ci.Element.Ca, ci.Element.Fe, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.Ca, ci.Element.Fe, timescale_type_pairs, ci.Element.He, ['SS', 'Dec'])
        #plot_timescale_ratios(ci.Element.Ca, ci.Element.Fe, timescale_type_pairs, ci.Element.He, ['SS'])

        #plot_timescale_ratios(ci.Element.Mg, ci.Element.Fe, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.Fe, ci.Element.Mg, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.Mg, ci.Element.Fe, timescale_type_pairs, ci.Element.He, ['SS'])

        #plot_timescale_ratios(ci.Element.Ca, ci.Element.Mg, timescale_type_pairs, ci.Element.H)
        #plot_timescale_ratios(ci.Element.Ca, ci.Element.Mg, timescale_type_pairs, ci.Element.He, ['SS'])

        #plot_timescale_ratios(ci.Element.Si, ci.Element.Fe, timescale_type_pairs, ci.Element.H)

        #    timescale_vals = generate_timescales()
        #    print(timescale_vals)
        #    #plot_2D_timescale(timescale_vals, 'He', None, 8000, -9.5)
        #    #plot_2D_timescale(timescale_vals, 'He', 8, None, -9.5)
        #    #plot_2D_timescale(timescale_vals, 'He', 8, 8000, None)
        #    #plot_2D_timescale(timescale_vals, 'H', None, 8000, -9.5)
        #    #plot_2D_timescale(timescale_vals, 'H', 8, None, -9.5)
        #
        #    #plot_model_comparison(timescale_vals, 'H', None, 6250, -9.5) # At the moment, we're limited to Teff == 6250 only  (the only val. I got old model values for)
        #    #plot_model_comparison(timescale_vals, 'He', None, 6250, -9.5) # At the moment, we're limited to Teff == 6250 only  (the only val. I got old model values for)
        #    plot_model_comparison(timescale_vals, 'He', 8, 6500, None)
        #    plot_model_comparison(timescale_vals, 'He', 8, None, -9.5)
        #    plot_model_comparison(timescale_vals, 'He', None, 6500, -9.5)
        #    #plot_model_comparison(timescale_vals, 'H', 8, None, -9.5)  # At the moment, we're limited to logg == 8 only

if __name__ == '__main__':
    main()
