#!/usr/bin/env python

import csv
import os
import sys

import pwd_utils as pu

sys.path.append(pu.get_path_to_utils())

import dict_plotter as dp
import timescale_interpolator as ti

test_dict = {
    "Test 1": {
        "Grid1": ti.TimescaleType.Bedard3DOvershoot,
        "Grid2": ti.TimescaleType.BedardNoOvershoot,
        "Entries": [
            "GaiaJ061-693(P;O)_R23",
            "GaiaJ061-693(P;U)_R23",
            "GD133_Xu2014",
            "GaiaJ061-693(S;O)_R23",
            "GD56_Xu2019",
            "HE0106-3253_Xu2019",
            "WD0310-688_Limbach24",
        ],
    },
    "Test 2": {
        "Grid1": ti.TimescaleType.BedardVariableOvershoot,
        "Grid2": ti.TimescaleType.BedardNoOvershoot,
        "Entries": [
            "PG1225-079Model1_Kl11",
            "GaiaJ064-035(P;O)_R23",
            "SDSSJ095+591(N)_H22",
            "SDSSJ095+591(O)_H22",
            "WDJ2047-1259_Hoskin20",
            "GD61_Farihi2011",
            "GD424(WHT)_Izquierd21",
            "GD424(Keck)_Izquier21",
            "HS2253+8023Model2_K11",
            "HS2253+8023Model3_K11",
            "Ross640_Blouin2018",
            "WD1232+563_Xu2019",
        ],
    },
    "Test 3": {
        "Grid1": ti.TimescaleType.KoesterNoOvershoot,
        "Grid2": ti.TimescaleType.BedardNoOvershoot,
        "Entries": [
            "NLTT1675_Kawka2012",
            "NLTT43806_Zuckerman11",
            "WDJ1935-3252_O'Brie23",
            "NLTT7547Spectral_Ka19",
            "NLTT7547(Balmer)_Ka19",
            "NLTT7547Balmer_Kawk19",
            "WD1124-293_Steele2021",
            "WD2115-560_Swan2019",
            "NLTT25792_Vennes&Ka13",
            "WD2157-574_Swan2019",
            "GD133_Xu2014",
            "WD0354+463_Vennes&K13",
            "WDJ113444+610826_T20",
            "WD1455+298_Vennes&K13",
            "G74-7_Vennes&Kawka13",
            "WD1257+278_Vennes&K13",
            "G166-58_Xu2019",
        ],
    },
    "Test 4": {
        "Grid1": ti.TimescaleType.KoesterNoOvershoot,
        "Grid2": ti.TimescaleType.BedardNoOvershoot,
        "Entries": [
            "PG1225-079Model1_Kl11",
            "GALEXJ2339_Klein2021",
            "GD17_GentileFusillo17",
            "PG1225-079Model2_Kl11",
            "GD16_GentileFusillo17",
            "SDSSJ095+591(N)_H22",
            "SDSSJ095+591(O)_H22",
            "SDSSJ12423+52262_R15",
            "WD1350-162_Swan2019",
            "SDSSJ1038-0036_Holl22",
            "WD2216-657_Swan2019",
            "J0939_Swan2023",
            "J1227_Swan2023",
            "HS2253+8023Model1_K11",
            "J0956_Swan2023",
            "PG1225-079Model3_Kl11",
            "HS2253+8023Model2_K11",
            "WD142+54(Model1)_X17",
            "WD142+54(Model2)_X17",
            "L745-46A_Koester&Wo00",
            "WD1232+563_Badenas-24",
            "Ross640_Blouin2018",
            "GaiaJ0218+3625_Doyl23",
            "WDJ183352+321757_T20",
            "WD1232+563_Xu2019",
        ],
    },
    "Test 5": {
        "Grid1": ti.TimescaleType.Bedard3DOvershootPatched,
        "Grid2": "DUMMY",
        "Entries": [
            "GaiaJ061-693(P;O)_R23",
            "GaiaJ061-693(P;U)_R23",
            "SDSSJ1043+0855_Meli17",
            "GD133_Xu2014",
            "GaiaJ051+231(PO)_R23",
            "GaiaJ051+231(PU)_R23",
            "GaiaJ061-693(S;O)_R23",
            "GaiaJ061-693(S;U)_R23",
            "G29-38_Xu2014",
            "PG1015+161_Gansicke12",
            "GD56_Xu2019",
            "HE0106-3253_Xu2019",
            "PG1015+161_Xu2019",
        ],
    },
}

output_dir = "/specify/path/to/results/here/"


def get_results_dir_and_file_name(entry, entry_dir, grid, thermohaline):
    print(entry)
    print(entry_dir)
    print(grid)
    thermohaline_suffix = "/t" if thermohaline else "/n"
    results_dir = entry_dir + grid.short_str() + thermohaline_suffix
    stats_file = None
    for file_name in os.listdir(results_dir):
        if file_name.endswith("stats.csv"):
            stats_file = file_name
    # file_name = entry + '_' + grid.short_str() + '_p2000_HD_NEL_D_stats.csv'
    print(results_dir)
    print(stats_file)
    print()
    return results_dir, stats_file


def retrieve_stats_from_file(results_dir, file_name):
    file_path = f"{results_dir}/{file_name}"
    print(file_path)
    with open(file_path, encoding="utf-8") as input_csv:
        models_outputted = 0
        good_fit = False
        model_section = False
        unphysical_fit = False
        for row in csv.reader(input_csv):
            if len(row) > 1:
                if row[0] == "Differentiation Sigma:":
                    diff_sigma = row[1]
                if row[0] == "Heating Sigma:":
                    heat_sigma = row[1]
                if row[0] == "Mg":
                    tau_Mg = float(row[1])
                if row[0] == "delta time/Myrs, +error, -error:":
                    delta_time = float(row[1])
                if row[0] == "Results from model:":
                    models_outputted += 1
                if row[0] == "Best model name:":
                    model_section = False
                if model_section:
                    if row[2] == "True":
                        good_fit = True
                if row[0] == "Model":
                    model_section = True
                if row[0] == "Parent Core Number Fraction:":
                    if row[1] != "":
                        if float(row[1]) < 0.01:
                            unphysical_fit = True
        assert models_outputted == 1
    if not good_fit or unphysical_fit:
        return None
    return diff_sigma, heat_sigma, tau_Mg, delta_time


def combine_stats(stats_1, stats_2):
    print()
    print()
    print()
    print()
    print()
    print(stats_1)
    print(stats_2)
    if stats_1 is None or stats_2 is None:
        toret = ["?", "?"]
    else:
        # Overturned sigmas are ones that only appear once
        toret = list()
        overturned = list()
        unoverturned = list()
        for stat_index in [0, 1]:
            if (stats_1[stat_index] == "N/A") ^ (stats_2[stat_index] == "N/A"):
                # Exclusive OR ---------------^
                # Then one of these is overturned
                overturned.append(stats_1[stat_index])
                overturned.append(stats_2[stat_index])
            else:
                unoverturned.append(stats_1[stat_index])
                unoverturned.append(stats_2[stat_index])

        valid_over_str = [str(sig) for sig in overturned if sig != "N/A"]
        overturned_str = ";".join(valid_over_str)
        toret.append(overturned_str)
        valid_unover_str = [str(sig) for sig in unoverturned if sig != "N/A"]
        unoverturned_str = ";".join(valid_unover_str)
        toret.append(unoverturned_str)
    for stats_obj in [stats_1, stats_2]:
        if stats_obj is None:
            toret.append(None)
            toret.append(None)
            toret.append(None)
        else:
            toret.append(stats_obj[3])
            toret.append(stats_obj[2])
            toret.append(1000000 * (stats_obj[3] / stats_obj[2]))
    return toret


def process_entry(entry, grid1, grid2, thermohaline):
    entry_dir = f"{output_dir}{entry}/"
    results_dir_1, file_1 = get_results_dir_and_file_name(
        entry, entry_dir, grid1, thermohaline
    )
    if not thermohaline:
        results_dir_2, file_2 = get_results_dir_and_file_name(
            entry, entry_dir, grid2, thermohaline
        )
    else:
        results_dir_2, file_2 = get_results_dir_and_file_name(
            entry, entry_dir, grid1, False
        )
    stats_1 = retrieve_stats_from_file(results_dir_1, file_1)
    print()
    print(stats_1)
    stats_2 = retrieve_stats_from_file(results_dir_2, file_2)
    print(stats_2)
    row_to_write = combine_stats(stats_1, stats_2)
    return row_to_write


def output_rows(list_of_rows):
    with open(
        "timescale_project_postprocess.csv", "w", newline="", encoding="utf-8"
    ) as f:
        to_write = csv.writer(f)
        to_write.writerow(
            [
                "Entry",
                "Sigma Overturned",
                "Sigma Unoverturned",
                "Delta Time 1 (Myr)",
                "tau Mg 1 (yr)",
                "D 1",
                "Delta Time 2 (Myr)",
                "tau Mg 2 (yr)",
                "D 2",
            ]
        )
        for row in list_of_rows:
            to_write.writerow(row)


def plot_metric_prediction():
    metric_overturned = [
        1.8,
        1.8,
        1.8,
        1.3,
        1,
        1.8,
        1.3,
        2.6,
        7.3,
        9.1,
        9.1,
        9.1,
        9.1,
        3.3,
        3.4,
        3.6,
        3.1,
        3.5,
        7.9,
        7.9,
        3.9,
        1.2,
        1,
        2,
        1.5,
        1.1,
    ]
    sigma_overturned = [
        1.35966916,
        1.35966916,
        1.49955924,
        2.562796205,
        1.407278135,
        1.401395077,
        2.140207376,
        1.26684019,
        2.192559182,
        4.408530133,
        4.408530133,
        4.274483102,
        4.274483102,
        4.685945573,
        3.636892944,
        2.399284553,
        1.245786305,
        2.717599212,
        1.407278135,
        2.680665767,
        1.401395077,
        1.661239502,
        2.439914937,
        1.921626685,
        2.009237626,
        1.270904128,
    ]
    metric_unoverturned = [
        1.3,
        1.3,
        1.9,
        1.9,
        1.9,
        1.9,
        1.1,
        1.1,
        1.3,
        1.3,
        1.3,
        1.3,
        1.5,
        1.5,
        1.5,
        1.5,
        2,
        2,
        1,
        1,
        2.8,
        2.8,
        1.3,
        1.3,
        1,
        1,
        2.3,
        2.3,
        2.3,
        2.3,
        6.8,
        6.8,
        6.8,
        6.8,
        1.1,
        1.1,
        1.8,
        1.8,
        1,
        1,
        1.5,
        1.5,
        2.8,
        2.8,
        2.7,
        2.7,
        3,
        3,
        3.3,
        3.3,
        3.4,
        3.4,
        4.2,
        4.2,
        3.6,
        3.6,
        3.2,
        3.2,
        2.6,
        2.6,
        2.6,
        2.6,
        2.1,
        2.1,
        11.3,
        11.3,
        11.3,
        11.3,
        1,
        1,
        3.6,
        3.6,
        3.6,
        3.6,
        1.3,
        1.3,
        1.3,
        1.3,
        1.4,
        1.4,
        1.6,
        1.6,
        1.6,
        1.6,
    ]
    sigma_unoverturned = [
        5.93830883,
        6.070770433,
        8.852872431,
        9.254754982,
        8.065222436,
        7.5112284,
        16.03050551,
        16.73824838,
        2.904947043,
        2.614783563,
        6.026608289,
        5.286614302,
        3.794235581,
        4.408530133,
        3.794235581,
        4.408530133,
        3.490238854,
        2.723185747,
        3.361971651,
        3.159852496,
        10.97029122,
        11.06997159,
        6.332186616,
        6.493544509,
        2.770069303,
        2.680665767,
        3.121450785,
        2.017806324,
        3.437825175,
        2.017806324,
        5.944501108,
        6.373208649,
        5.944501108,
        6.373208649,
        1.439628473,
        1.851318457,
        2.057168595,
        2.487228715,
        5.382709731,
        6.070770433,
        2.366261061,
        1.752728128,
        16.6477865,
        16.73824838,
        16.35001836,
        16.98119396,
        1.337698271,
        2.007598521,
        3.346083396,
        4.685945573,
        8.781783749,
        9.580226671,
        5.226961048,
        4.696653494,
        4.587500593,
        2.615399381,
        16.57337666,
        16.99112206,
        3.361032456,
        4.951153208,
        3.498240605,
        4.951153208,
        1.649503024,
        1.128966,
        31.28482972,
        31.24975509,
        11.07853148,
        10.78170478,
        3.861457164,
        2.773061662,
        7.825539286,
        8.977025944,
        7.825539286,
        8.977025944,
        2.128206891,
        2.995845198,
        2.1675833,
        2.995845198,
        6.0862972948038,
        5.50673590356409,
        1.67688789460071,
        5.91218518640607,
        2.22030121464052,
        5.88465037756019,
    ]
    series_dict = {
        "Unambiguous": {
            "type": dp.SeriesType.scatter_2d,
            "x_data": metric_unoverturned,
            "y_data": sigma_unoverturned,
            "line_color": "#3a2aba",
            "line_marker": ".",
            "line_style": "None",
            "line_markersize": 10,
            "legend": True,
        },
        "Ambiguous": {
            "type": dp.SeriesType.scatter_2d,
            "x_data": metric_overturned,
            "y_data": sigma_overturned,
            "line_color": "#bb2a2a",
            "line_marker": "x",
            "line_style": "None",
            "line_markersize": 10,
            "legend": True,
        },
        "y=x": {
            "type": dp.SeriesType.scatter_2d,
            "x_data": [0, 30],
            "y_data": [0, 30],
            "line_color": "silver",
            "line_style": "--",
            "legend": False,
        },
    }
    plot_dict1 = {
        "metricplot": {
            "show": False,
            "filenames": ["metric_plot.pdf", "metric_plot.png"],
            "fig_height": 5,
            "fig_width": 7,
            "subplots": {
                "subplot1": {
                    "subplot_region": 111,
                    "legend": True,
                    "legend_loc": "best",
                    "legend_text_size": 18,
                    "title_fontsize": 12,
                    "title_fontweight": "bold",
                    "xlabel_text": "Discrepancy Metric",
                    "xlabel_fontsize": 20,
                    "xlabel_fontweight": "bold",
                    "ylabel_text": "Sigma significance of result",
                    "ylabel_fontsize": 20,
                    "ylabel_fontweight": "bold",
                    "x_min": 0,
                    "x_max": 12,
                    "y_min": 0,
                    # 'y_max': 1,
                    "font": "STIXGeneral",
                    "series": series_dict,
                }
            },
        }
    }
    plotter = dp.DictPlotter(plot_dict1)
    plotter.draw()
    plotter.yield_output()


def make_stats_table():
    list_of_rows = list()
    for test_name, test_description in test_dict.items():
        grid1 = test_description["Grid1"]
        grid2 = test_description["Grid2"]
        entries = test_description["Entries"]
        if test_name == "Test 5":
            thermohaline = True
        else:
            thermohaline = False
        for entry in entries:
            results = process_entry(entry, grid1, grid2, thermohaline)
            print()
            print(entry)
            print(results)
            list_of_rows.append([entry] + results)
        list_of_rows.append([])
    output_rows(list_of_rows)


def main():
    make_stats_table()
    plot_metric_prediction()


if __name__ == "__main__":
    main()
