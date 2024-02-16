#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum

import csv
import geology_info as gi
import model_parameters as mp
import numpy as np

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
                    teff = int(row[3])
                except ValueError:
                    teff = None
                if teff is not None:
                    teffs.append(teff)
                try:
                    logg = float(row[4])
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

def plot_demo_distributions():
    import graph_factory as gf
    graph_fac = gf.GraphFactory()
    xbar = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
    all_heights1 = [
        [1, 1, 0, 0, 0, 0, 0, 0, 0, 1]
    ]
    all_heights2 = [
        [10, 9, 8, 7, 7, 6, 3, 3, 3, 5]
    ]
    all_heights3 = [
        [10, 11, 12, 20, 15, 4, 3, 3, 3, 3]
    ]

    graph_fac.make_histogram(
        xbar,
        all_heights1,
        ['Example 1'],
        'DemoDist',
        0.1,
        1,
        'Fragment Core Number Fraction',
        '1',
        dict(),
        None
    )

    graph_fac.make_histogram(
        xbar,
        all_heights2,
        ['Example 2'],
        'DemoDist',
        0.1,
        1,
        'Fragment Core Number Fraction',
        '2',
        dict(),
        None
    )

    graph_fac.make_histogram(
        xbar,
        all_heights3,
        ['Example 3'],
        'DemoDist',
        0.1,
        1,
        'Fragment Core Number Fraction',
        '3',
        dict(),
        None
    )

def main():
    plot_demo_distributions()

if __name__ == '__main__':
    main()

