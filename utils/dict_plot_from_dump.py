#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections as cn
import dict_plotter as dp
import numpy as np
import os
import sys

# Usage: python dict_plot_from_dump.py name-of-dump-file
# The dump file should be a .txt file produced by dict_plotter.py
# This will read in the .txt file and plot the graph it describes
# This is untested on "function" SeriesTypes - they probably don't work (update: it is, and they don't)

def load_dict_dump(dump_file):
    container_dict = dict()
    dict_string = open(dump_file, 'r').read()
    dict_to_plot = eval(dict_string)
    container_dict['reconstructed_plot'] = dict_to_plot
    return container_dict

def main():
    text_dump_inc_path = sys.argv[1]
    output_dir = os.path.dirname(text_dump_inc_path)
    dict_to_plot = load_dict_dump(text_dump_inc_path)
    plotter = dp.DictPlotter(dict_to_plot)
    plotter.draw()
    plotter.yield_output(output_dir, False)

if __name__ == '__main__':
    main()
