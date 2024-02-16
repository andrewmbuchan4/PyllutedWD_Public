#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import scipy.stats as st
import synthetic_pipeline as spi

def ks_test(data1, data2):
    pipeline = spi.Pipeline('TestPipeline', dict(), dict(), dict())
    min_bin_edge = 0.54
    max_bin_edge = 0.64
    half_bin_size = 0.01
    ks_statistic, p_value = pipeline.run_ks_test_on_samples(data1, data2, min_bin_edge, max_bin_edge, half_bin_size)
    scipy_result = st.ks_2samp(data1, data2, "two-sided", "asymp")
    print(ks_statistic)
    print(p_value)
    print(scipy_result)

def read_data_from_txt(filename):
    # This function assumes that the input file has no header, and two columns (of equal length) separated by a space
    data1 = list()
    data2 = list()
    with open(filename) as csvfile:
        read = csv.reader(csvfile, delimiter=' ')
        for row in read:
            data1.append(float(row[0]))
            data2.append(float(row[1]))
    return data1, data2

def main():
    #data1, data2 = read_data_from_txt('KS_test_data.txt')
    data1 = [0.6, 0.61, 0.63, 0.59, 0.55]
    data2 = [0.62, 0.61, 0.63, 0.59, 0.55]
    ks_test(data1, data2)

if __name__ == '__main__':
    main()
