#!/usr/bin/env python

import csv
import numpy as np

import graph_factory as gf
import pwd_utils as pu

def load_raw_water_content(file_including_path):
    toret = list()
    with open(file_including_path, encoding='utf-8') as input_csv:
        for row in csv.reader(input_csv):
            toret.append(float(row[0]))
    return toret

def process_wc(raw_wc):
    linear = [max(0, wc) for wc in raw_wc]
    linear_positive = [wc for wc in raw_wc if wc > 0]
    log = [0.0-np.log10(wc) for wc in raw_wc if wc > 0]
    return linear, linear_positive, log


def load_water_content(file_including_path):
    raw_wc = load_raw_water_content(file_including_path)
    toret_linear, toret_linear_positive, toret_log = process_wc(raw_wc)
    return toret_linear, toret_linear_positive, toret_log

def main():
    directory = pu.get_path_to_pipeline_base_dir() + 'volatiles/'
    print(directory)
    wet_linear, wet_linear_positive, wet_log = load_water_content(directory + 'raw_water_content_wet.csv')
    dry_linear, dry_linear_positive, dry_log = load_water_content(directory + 'raw_water_content_dry.csv')
    wc_to_plot_dict_linear = {
        'Wet_Linear': wet_linear,
        'Dry_Linear': dry_linear,
        'Wet_Linear_Nonzero': wet_linear_positive,
        'Dry_Linear_Nonzero': dry_linear_positive,
    }
    wc_to_plot_dict_log = {
        'Wet_Log': wet_log,
        'Dry_Log': dry_log
    }

    linear_bin_size = 0.05
    linear_bins = np.arange(0.0, 1.0 + linear_bin_size, linear_bin_size)
    linear_bin_centres = np.arange(0.0 + (linear_bin_size/2), 1.0, linear_bin_size)

    max_log_bin = max(max(wet_log), max(dry_log))
    print(max_log_bin)
    log_bin_size = 0.05
    log_bins = np.arange(0.0, max_log_bin + log_bin_size, log_bin_size)
    log_bin_centres = np.arange(0.0 + (log_bin_size/2), max_log_bin, log_bin_size)

    print(len(linear_bins))
    print(len(linear_bin_centres))

    print(log_bins)
    print(log_bin_centres)

    graph_fac = gf.GraphFactory(directory)
    for pop_name, vals in wc_to_plot_dict_linear.items():
        heights, bins2 = np.histogram(
            vals,
            linear_bins,
            density=True
        )

        graph_fac.make_histogram(
            linear_bin_centres,
            [heights],
            [pop_name],
            pop_name,
            linear_bin_size,
            1.1,
            'Water Mass Fraction',
            '_water',
            None,
            None
        )

    for pop_name, vals in wc_to_plot_dict_log.items():
        heights, bins2 = np.histogram(
            vals,
            log_bins,
            density=True
        )

        print(len(log_bins))
        print(len(log_bin_centres))

        graph_fac.make_histogram(
            log_bin_centres,
            [heights],
            [pop_name],
            pop_name,
            log_bin_size,
            1.1,
            'log_10(1 - Water Mass Fraction)',
            '_water',
            None,
            None
        )

if __name__ == '__main__':
    main()
