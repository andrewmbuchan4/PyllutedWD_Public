#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pprint
import matplotlib as mpl
import matplotlib.colors as clr
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import gridspec
import matplotlib.transforms as transforms
import numpy as np
from enum import Enum
from matplotlib import ticker
from ternary_diagram import TernaryDiagram
plt.rcParams['text.usetex'] = 'True'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
#plt.rcParams["text.latex.preamble"].append(r"\usepackage{mathabx}")
# TODO: Add comments. Remove repeated strings. Think of an alternative to the big if/elif/else blocks. Decide whether it's color or colour!

class SeriesType(Enum):
    function_2d = 0
    scatter_2d = 1
    function_3d = 2
    scatter_3d = 3
    vline = 4
    text = 5
    shade = 6
    contour_scatter = 7
    scatter_2d_error = 8
    tricontour_scatter = 9
    ellipse = 10
    hline = 11
    arrow = 12
    hist2d = 13
    bar = 14
    pie = 15
    ternary_scatter = 16 # Could remove this and make scatter_3d function as a ternary_scatter if the subplottype is ternary?!? Already done something similar for text. TODO

class SubplotType(Enum):
    cartesian = 0
    ternary = 1
    polar = 2 #TODO!

class DictPlotter:

    def __init__(self, input_dict):
        self.unpack_input_dict(input_dict)

    def unpack_input_dict(self, input_dict):
        self.plots = dict()
        for plot_name, plot_parameters in input_dict.items():
            self.plots[plot_name] = Plot(plot_parameters)

    def draw(self):
        for plot_name, plot in self.plots.items():
            print('Drawing ' + plot_name)
            plot.draw()

    def yield_output(self, output_dir=None, dump=True):
        for plot_name, plot in self.plots.items():
            print('Yielding output for ' + plot_name)
            plot.yield_output(output_dir, dump)

class Plot:

    def __init__(self, input_dict):
        self.raw_input_dict = input_dict
        self.unpack_input_dict(input_dict)

    def unpack_input_dict(self, input_dict):
        self.show = input_dict.get('show')
        self.filenames = input_dict.get('filenames')
        self.fig_height = input_dict.get('fig_height')
        self.fig_width = input_dict.get('fig_width')
        self.dpi = input_dict.get('dpi')
        self.gridspec = None
        if input_dict.get('gridspec_y_x') is not None:
            self.gridspec = gridspec.GridSpec(
                input_dict.get('gridspec_y_x')[0],
                input_dict.get('gridspec_y_x')[1],
                width_ratios=input_dict.get('gridspec_width_ratios'),
                height_ratios=input_dict.get('gridspec_height_ratios'),
                wspace=input_dict.get('gridspec_wspace'),
                hspace=input_dict.get('gridspec_hspace')
            )
        self.subplots = dict()
        for subplot_name, subplot_parameters in input_dict.get('subplots').items():
            self.subplots[subplot_name] = Subplot(subplot_parameters)

    def draw(self):
        if self.fig_width is not None and self.fig_height is not None:
            self.fig = plt.figure(figsize=(self.fig_width, self.fig_height))
        else:
            self.fig = plt.figure()
        for subplot_name, subplot in self.subplots.items():
            print('Plotting subplot ' + subplot_name)
            if subplot.twin_subplot is None:
                try:
                    sharex_sub = self.subplots[subplot.sharex_subplot]
                except KeyError:
                    sharex_sub = None
                try:
                    sharey_sub = self.subplots[subplot.sharey_subplot]
                except KeyError:
                    sharey_sub = None
                subplot.draw(self.fig, self.gridspec, None, sharex_sub, sharey_sub)
                if sharex_sub is not None:
                    plt.setp(subplot.ax.get_xticklabels(), visible=False)
                    subplot.ax.xaxis.set_tick_params(direction='inout')
                    subplot.ax.xaxis.set_label_text(None)
                if sharey_sub is not None:
                    plt.setp(subplot.ax.get_yticklabels(), visible=False)
                    subplot.ax.yaxis.set_tick_params(direction='inout')
                    subplot.ax.yaxis.set_label_text(None)
            else:
                subplot.draw(self.fig, self.gridspec, self.subplots[subplot.twin_subplot])  # TODO: This does assume we've already plotted the twin, which is not guaranteed
        plt.tight_layout()
        #self.fig.subplots_adjust(wspace=0)
        #self.fig.subplots_adjust(hspace=0)

    def yield_output(self, output_dir=None, dump=True, close=True):
        if output_dir is not None:
            if not output_dir.endswith('/'):
                output_dir += '/'
            try:
                os.makedirs(output_dir)
            except FileExistsError:
                pass
        if self.show:
            plt.show()
        if self.filenames is not None:
            for filename in self.filenames:
                if output_dir is None:
                    print('Attempting to save to ' + filename)
                    self.fig.savefig(filename, dpi=self.dpi, bbox_inches='tight', pad_inches=0.1)
                else:
                    print('Attempting to save to: ' + output_dir + filename)
                    self.fig.savefig(output_dir + filename, dpi=self.dpi, bbox_inches='tight', pad_inches=0.1)
        if dump:
            s = pprint.pformat(self.raw_input_dict)
            s = s.replace('array', 'np.array')
            s = s.replace('nan', 'np.nan')
            s = s.replace('inf', 'np.inf')
            # TODO: Do this logic by cycling through all possible enum values? Might be hard to read
            s = s.replace('<SeriesType.function_2d: 0>', 'dp.SeriesType.function_2d')
            s = s.replace('<SeriesType.scatter_2d: 1>', 'dp.SeriesType.scatter_2d')
            s = s.replace('<SeriesType.function_3d: 2>', 'dp.SeriesType.function_3d')
            s = s.replace('<SeriesType.scatter_3d: 3>', 'dp.SeriesType.scatter_3d')
            s = s.replace('<SeriesType.vline: 4>', 'dp.SeriesType.vline')
            s = s.replace('<SeriesType.text: 5>', 'dp.SeriesType.text')
            s = s.replace('<SeriesType.shade: 6>', 'dp.SeriesType.shade')
            s = s.replace('<SeriesType.contour_scatter: 7>', 'dp.SeriesType.contour_scatter')
            s = s.replace('<SeriesType.scatter_2d_error: 8>', 'dp.SeriesType.scatter_2d_error')
            s = s.replace('<SeriesType.tricontour_scatter: 9>', 'dp.SeriesType.tricontour_scatter')
            s = s.replace('<SeriesType.ellipse: 10>', 'dp.SeriesType.ellipse')
            s = s.replace('<SeriesType.hline: 11>', 'dp.SeriesType.hline')
            s = s.replace('<SeriesType.arrow: 12>', 'dp.SeriesType.arrow')
            s = s.replace('<SeriesType.hist2d: 13>', 'dp.SeriesType.hist2d')
            s = s.replace('<SeriesType.bar: 14>', 'dp.SeriesType.bar')
            s = s.replace('<SeriesType.pie: 15>', 'dp.SeriesType.pie')
            s = s.replace('<SeriesType.ternary_scatter: 16>', 'dp.SeriesType.ternary_scatter')
            s = s.replace('<SubplotType.cartesian: 0>', 'dp.SubplotType.cartesian')
            s = s.replace('<SubplotType.ternary: 1>', 'dp.SubplotType.ternary')
            s = s.replace('<SubplotType.polar: 2>', 'dp.SubplotType.polar')
            s = s.replace('OrderedDict', 'cn.OrderedDict')
            destination = self.filenames[0] + '.txt'
            if output_dir is not None:
                destination = output_dir + destination
            with open(destination, 'w') as dump_file:
                print(s, file=dump_file)
            print('Dumped text version to ' + destination)
        if close:
            plt.clf()
            plt.close(self.fig)

class Subplot:

    def __init__(self, input_dict):
        self.unpack_input_dict(input_dict)

    def unpack_input_dict(self, input_dict):
        self.subplot_region = input_dict.get('subplot_region')
        self.subplot_type = input_dict.get('type', SubplotType.cartesian) # cartesian by default
        self.legend = input_dict.get('legend')
        self.legend_loc = input_dict.get('legend_loc')
        self.legend_coord_x = input_dict.get('legend_coord_x')
        self.legend_coord_y = input_dict.get('legend_coord_y')
        self.legend_text_size = input_dict.get('legend_text_size')
        self.title_text = input_dict.get('title_text')
        #try:
        #    self.title_text = self.title_text.replace('_', '\_')
        #except AttributeError:
        #    pass
        self.title_fontsize = input_dict.get('title_fontsize')
        self.title_fontweight = input_dict.get('title_fontweight')
        self.xlabel_text = input_dict.get('xlabel_text')
        #try:
        #    self.xlabel_text = input_dict.get('xlabel_text').replace('_', '\_') # Is this broken??
        #except AttributeError:
        #    self.xlabel_text = None
        self.ylabel_text = input_dict.get('ylabel_text')
        #try:
        #    self.ylabel_text = input_dict.get('ylabel_text').replace('_', '\_')
        #except AttributeError:
        #    self.ylabel_text = None
        self.xlabel_fontsize = input_dict.get('xlabel_fontsize')
        self.xlabel_fontweight = input_dict.get('xlabel_fontweight')
        self.ylabel_fontsize = input_dict.get('ylabel_fontsize')
        self.zlabel_fontsize = input_dict.get('zlabel_fontsize')
        self.ylabel_fontweight = input_dict.get('ylabel_fontweight')
        self.xlabel_pad = input_dict.get('xlabel_pad')
        self.ylabel_pad = input_dict.get('ylabel_pad')
        self.x_min = input_dict.get('x_min')
        self.x_max = input_dict.get('x_max')
        self.y_min = input_dict.get('y_min')
        self.y_max = input_dict.get('y_max')
        self.x_scale = input_dict.get('x_scale')
        self.y_scale = input_dict.get('y_scale')
        self.font = input_dict.get('font')
        self.x_tick_labels = input_dict.get('x_tick_labels')
        #try:
        #    self.x_tick_labels = self.x_tick_labels.replace('_', '\_')
        #except AttributeError:
        #    pass
        self.y_tick_labels = input_dict.get('y_tick_labels')
        #try:
        #    self.y_tick_labels = self.y_tick_labels.replace('_', '\_')
        #except AttributeError:
        #    pass
        self.x_tick_locations = input_dict.get('x_tick_locations')
        self.y_tick_locations = input_dict.get('y_tick_locations')
        self.x_tick_fontsize = input_dict.get('x_tick_fontsize')
        self.y_tick_fontsize = input_dict.get('y_tick_fontsize')
        self.x_tick_colour = input_dict.get('x_tick_colour')
        self.y_tick_colour = input_dict.get('y_tick_colour')
        self.x_label_colour = input_dict.get('x_label_colour')
        self.y_label_colour = input_dict.get('y_label_colour')
        self.z_label_colour = input_dict.get('z_label_colour')
        self.twin_subplot = input_dict.get('twin_subplot')
        self.twin_on_x = input_dict.get('twin_on_x')  # Sorry that this setup is a bit silly. TODO: Fix
        self.sharex_subplot = input_dict.get('sharex_subplot')
        self.sharey_subplot = input_dict.get('sharey_subplot')
        self.gridspec_index = input_dict.get('gridspec_index')
        self.x_hide_ticks = input_dict.get('x_hide_ticks')  # indices of ticks/labels to hide
        self.y_hide_ticks = input_dict.get('y_hide_ticks')  # indices of ticks/labels to hide
        self.x_label_shade_colour = input_dict.get('x_label_shade_colour')
        self.x_label_alpha = input_dict.get('x_label_alpha')
        self.y_label_shade_colour = input_dict.get('y_label_shade_colour')
        self.y_label_alpha = input_dict.get('y_label_alpha')
        self.z_label_shade_colour = input_dict.get('z_label_shade_colour')
        self.z_label_alpha = input_dict.get('z_label_alpha')
        self.series = dict()
        if input_dict.get('series') is not None:
            for series_name, series_parameters in input_dict.get('series').items():
                self.series[series_name] = Series(series_parameters)
        if self.subplot_type == SubplotType.ternary:
            # I don't like the implementation of labels in ternary diagram (ie putting them at the corners). Doing this instead:
            if input_dict.get('x_label') is not None:
                self.series['x_label_text'] = Series({
                    'type': SeriesType.text,
                    'x_pos': 0.85,
                    'y_pos': 0.5,
                    'text_string': input_dict.get('x_label'),
                    'fontsize': self.xlabel_fontsize,
                    'rotation': 0,
                    'text_colour': self.x_label_colour,
                    'shade_colour': self.x_label_shade_colour,
                    'text_alpha': self.x_label_alpha,
                    'horizontalalignment': 'left',
                    'verticalalignment': 'center',
                })
            if input_dict.get('y_label') is not None:
                self.series['y_label_text'] = Series({
                    'type': SeriesType.text,
                    'x_pos': 0.15,
                    'y_pos': 0.5,
                    'text_string': input_dict.get('y_label'),
                    'fontsize': self.ylabel_fontsize,
                    'rotation': -60,
                    'text_colour': self.y_label_colour,
                    'shade_colour': self.y_label_shade_colour,
                    'text_alpha': self.y_label_alpha,
                    'horizontalalignment': 'right',
                    'verticalalignment': 'center',
                })
            if input_dict.get('z_label') is not None:
                self.series['z_label_text'] = Series({
                    'type': SeriesType.text,
                    'x_pos': 0.5,
                    'y_pos': -0.1,
                    'text_string': input_dict.get('z_label'),
                    'fontsize': self.zlabel_fontsize,
                    'rotation': 60,
                    'text_colour': self.z_label_colour,
                    'shade_colour': self.z_label_shade_colour,
                    'text_alpha': self.z_label_alpha,
                    'horizontalalignment': 'center',
                    'verticalalignment': 'center',
                })

    def draw(self, figure, gridspec=None, twin_subplot=None, sharex_subplot=None, sharey_subplot=None):
        if twin_subplot is None:
            if gridspec is not None and self.gridspec_index is not None:
                try:
                    sharex_sub_ax = sharex_subplot.ax
                except AttributeError:
                    sharex_sub_ax = None
                try:
                    sharey_sub_ax = sharey_subplot.ax
                except AttributeError:
                    sharey_sub_ax = None
                self.ax = figure.add_subplot(gridspec[self.gridspec_index], sharex=sharex_sub_ax, sharey=sharey_sub_ax)
            else:
                try:
                    # Assume the subplot region is a 3-tuple or list of len >= 3...
                    #TODO: This is a bit buggy. 5, 2, 9 doesn't appear!
                    self.ax = figure.add_subplot(self.subplot_region[0], self.subplot_region[1], self.subplot_region[2], sharex=sharex_subplot.ax, sharey=sharey_subplot.ax)
                except (TypeError, IndexError) as e:
                    # ...But if not, assume we just want to use the int version of this function
                    try:
                        sharex_sub_ax = sharex_subplot.ax
                    except AttributeError:
                        sharex_sub_ax = None
                    try:
                        sharey_sub_ax = sharey_subplot.ax
                    except AttributeError:
                        sharey_sub_ax = None
                    self.ax = figure.add_subplot(self.subplot_region, sharex=sharex_sub_ax, sharey=sharey_sub_ax)
        else:
            if self.twin_on_x:
                self.ax = twin_subplot.ax.twiny()
                self.ax.xaxis.tick_top()
                self.ax.xaxis.set_label_position('top')
            else:
                self.ax = twin_subplot.ax.twinx()
                self.ax.yaxis.tick_right()
                self.ax.yaxis.set_label_position('right')
        self.legend_lines = list()
        self.legend_labels = list()
        if self.subplot_type == SubplotType.ternary:
            self.td = TernaryDiagram(['', '', ''], ax=self.ax) # We will trick the Series object into thinking that this is the actual axis object
        for series_name, series in self.series.items():
            print('Plotting series ' + series_name)
            if self.subplot_type == SubplotType.ternary:
                line = series.draw(self.td, figure)
            else:
                line = series.draw(self.ax, figure)
            if series.legend:
                self.legend_labels.append(series_name.replace('_', '\_'))
                # Possibly would be better to just check whether line is wrapped in a list rather than doing this:
                if series.plot_type in [SeriesType.shade, SeriesType.hline, SeriesType.vline, SeriesType.arrow]:
                    self.legend_lines.append(line)
                elif series.plot_type == SeriesType.ternary_scatter and series.line_style is None:
                    self.legend_lines.append(line)
                else:
                    self.legend_lines.append(line[0])  # Take 0th element because the line gets wrapped in a list
        if self.legend:
            if twin_subplot is None:
                if self.legend_coord_x is None or self.legend_coord_y is None:
                    self.ax.legend(self.legend_lines, self.legend_labels, loc=self.legend_loc, prop={'size': self.legend_text_size, 'family': self.font})
                else:
                    self.ax.legend(self.legend_lines, self.legend_labels, loc=self.legend_loc, prop={'size': self.legend_text_size, 'family': self.font}, bbox_to_anchor=(self.legend_coord_x, self.legend_coord_y))
            else:
                if twin_subplot.legend_coord_x is None or twin_subplot.legend_coord_y is None:
                    twin_subplot.ax.legend(twin_subplot.legend_lines + self.legend_lines, twin_subplot.legend_labels + self.legend_labels, loc=twin_subplot.legend_loc, prop={'size': twin_subplot.legend_text_size, 'family': twin_subplot.font})
                else:
                    twin_subplot.ax.legend(twin_subplot.legend_lines + self.legend_lines, twin_subplot.legend_labels + self.legend_labels, loc=twin_subplot.legend_loc, prop={'size': twin_subplot.legend_text_size, 'family': twin_subplot.font}, bbox_to_anchor=(twin_subplot.legend_coord_x, twin_subplot.legend_coord_y))
        if self.title_text:
            self.ax.set_title(self.title_text, fontsize=self.title_fontsize, fontweight=self.title_fontweight)#, fontfamily=self.font)
        if self.xlabel_text:
            self.ax.set_xlabel(self.xlabel_text, fontsize=self.xlabel_fontsize, fontweight=self.xlabel_fontweight, labelpad=self.xlabel_pad)#, fontfamily=self.font)
        if self.ylabel_text:
            self.ax.set_ylabel(self.ylabel_text, fontsize=self.ylabel_fontsize, fontweight=self.ylabel_fontweight, labelpad=self.ylabel_pad)#, fontfamily=self.font)
        if self.x_min is not None:
            self.ax.set_xlim(left=self.x_min)
        if self.x_max is not None:
            self.ax.set_xlim(right=self.x_max)
        if self.y_min is not None:
            self.ax.set_ylim(bottom=self.y_min)
        if self.y_max is not None:
            self.ax.set_ylim(top=self.y_max)
        if self.x_scale is not None:
            self.ax.set_xscale(self.x_scale)
        if self.y_scale is not None:
            self.ax.set_yscale(self.y_scale)
        if self.x_tick_fontsize is not None:
            self.ax.xaxis.set_tick_params(labelsize=self.x_tick_fontsize)
        if self.x_tick_labels is not None:
            #self.ax.set_xticks(list(range(0, len(self.x_tick_labels)))) # Not sure why this is here? This silently messes with your tick locations
            self.ax.set_xticklabels(self.x_tick_labels)
        if self.y_tick_fontsize is not None:
            self.ax.yaxis.set_tick_params(labelsize=self.y_tick_fontsize)
        if self.y_tick_labels is not None:
            #self.ax.set_yticks(list(range(0, len(self.y_tick_labels))))
            self.ax.set_yticklabels(self.y_tick_labels)
        if self.x_tick_locations is not None:
            self.ax.set_xticks(self.x_tick_locations)
        if self.y_tick_locations is not None:
            self.ax.set_yticks(self.y_tick_locations)
        if self.x_tick_colour is not None:
            self.ax.xaxis.set_tick_params(colors=self.x_tick_colour)
        if self.y_tick_colour is not None:
            self.ax.yaxis.set_tick_params(colors=self.y_tick_colour)
        if self.x_label_colour is not None:
            self.ax.xaxis.label.set_color(self.x_label_colour)
        if self.y_label_colour is not None:
            self.ax.yaxis.label.set_color(self.y_label_colour)
        if self.x_hide_ticks is not None:
            for x_tick_index in self.x_hide_ticks:
                #self.ax.xaxis.ticks[x_tick_index].set_visible(False)
                plt.setp(self.ax.get_xticklabels()[x_tick_index], visible=False)
        if self.y_hide_ticks is not None:
            for y_tick_index in self.y_hide_ticks:
                plt.setp(self.ax.get_yticklabels()[y_tick_index], visible=False)

class Series:

    def __init__(self, input_dict):
        self.mandatory_args = {
            SeriesType.function_2d: [
                'x_start',
                'x_end',
                'x_points',
                'function',
                'line_type'
            ],
            SeriesType.scatter_2d: [
                'x_data',
                'y_data'
            ],
            SeriesType.function_3d: [
                'x_start',
                'x_end',
                'x_points',
                'y_start',
                'y_end',
                'y_points',
                'function',
                'fill'
            ],
            SeriesType.contour_scatter: [
                'x_data',
                'y_data',
                'z_data',
                'fill'
            ],
            SeriesType.tricontour_scatter: [
                'x_data',
                'y_data',
                'z_data',
                'fill'
            ],
            SeriesType.scatter_3d: [
                'x_data',
                'y_data',
                'z_data'
            ],
            SeriesType.vline: [
                'x_start',
                'y_min',
                'y_max'
            ],
            SeriesType.hline: [
                'y_start',
                'x_min',
                'x_max'
            ],
            SeriesType.text: [
                'x_pos',
                'y_pos',
                'text_string'
            ],
            SeriesType.shade: [
                'x_data',
                'y_data',
                'y_shade_data'
            ],
            SeriesType.scatter_2d_error: [
                'x_data',
                'y_data'
            ],
            SeriesType.ellipse: [
                'x_data',
                'y_data',
                'n_std'
            ],
            SeriesType.arrow: [
                'x_start',
                'y_start',
                'x_change',
                'y_change'
            ],
            SeriesType.hist2d: [
                'x_data',
                'y_data'
            ],
            SeriesType.bar: [
                'x_data',
                'y_data'
            ],
            SeriesType.pie: [
                'x_data'
            ],
            SeriesType.ternary_scatter: [
                'x_data',
                'y_data',
                'z_data'
            ]
        }
        self.unpack_input_dict(input_dict)

    def unpack_input_dict(self, input_dict):
        self.plot_type = input_dict['type']  #TODO: Error handling
        self.function = input_dict.get('function')
        self.function_args = input_dict.get('function_args')
        self.line_color = input_dict.get('line_color')
        self.line_style = input_dict.get('line_style')
        self.line_marker = input_dict.get('line_marker')
        self.line_markersize = input_dict.get('line_markersize')
        self.line_linewidth = input_dict.get('line_linewidth')
        self.line_type = input_dict.get('line_type', '')
        self.x_start = input_dict.get('x_start')
        self.x_end = input_dict.get('x_end')
        self.x_points = input_dict.get('x_points')
        self.y_start = input_dict.get('y_start')
        self.y_end = input_dict.get('y_end')
        self.y_points = input_dict.get('y_points')
        self.x_data = input_dict.get('x_data')
        self.y_data = input_dict.get('y_data')
        self.z_data = input_dict.get('z_data')
        self.c_data = input_dict.get('c_data') # Intended as a 4th dimension (eg the colour on a ternary plot)
        self.cbar_label = input_dict.get('cbar_label')
        try:
            self.cbar_label = self.cbar_label.replace('_', '\_')
        except AttributeError:
            pass
        self.cbar_ticks = input_dict.get('cbar_ticks')
        self.cbar_ticklabels = input_dict.get('cbar_ticklabels')
        self.cbar_labelfontsize = input_dict.get('cbar_labelfontsize')
        self.cbar_tickfontsize = input_dict.get('cbar_tickfontsize')
        self.cbar_labelrotation = input_dict.get('cbar_labelrotation')
        self.cbar_labelpad = input_dict.get('cbar_labelpad')
        self.cbar_shrink = input_dict.get('cbar_shrink')
        self.colour_map_colours = input_dict.get('colour_map_colours')
        self.colour_map_min = input_dict.get('colour_map_min')
        self.colour_map_max = input_dict.get('colour_map_max')
        self.legend = input_dict.get('legend')
        self.vline_label = input_dict.get('vline_label')
        self.y_min = input_dict.get('y_min')
        self.y_max = input_dict.get('y_max')
        self.hline_label = input_dict.get('hline_label')
        self.x_min = input_dict.get('x_min')
        self.x_max = input_dict.get('x_max')
        self.x_pos = input_dict.get('x_pos')
        self.y_pos = input_dict.get('y_pos')
        self.text_string = input_dict.get('text_string')
        try:
            self.text_string = self.text_string.replace('_', '\_')
        except AttributeError:
            pass
        self.y_shade_data = input_dict.get('y_shade_data')
        self.y_error_data = input_dict.get('y_error_data')
        self.x_error_data = input_dict.get('x_error_data')
        self.shade_colour = input_dict.get('shade_colour')
        self.shade_alpha = input_dict.get('shade_alpha')
        self.shade_hatch = input_dict.get('shade_hatch')
        self.shade_where = input_dict.get('shade_where')
        self.fill = input_dict.get('fill')
        self.levels = input_dict.get('levels')
        self.colours = input_dict.get('colours')
        self.colour_below_min = input_dict.get('colour_below_min')
        self.colour_above_max = input_dict.get('colour_above_max')
        self.zorder = input_dict.get('zorder')
        self.capsize = input_dict.get('capsize')
        self.capthick = input_dict.get('capthick')
        self.fontsize = input_dict.get('fontsize')
        self.horizontalalignment = input_dict.get('horizontalalignment')
        self.verticalalignment = input_dict.get('verticalalignment')
        self.rotation = input_dict.get('rotation')
        self.text_colour = input_dict.get('text_colour', 'k')
        self.lower_limits = input_dict.get('lower_limits')
        self.upper_limits = input_dict.get('upper_limits')
        self.error_linewidth = input_dict.get('error_linewidth')
        self.n_std = input_dict.get('n_std')
        self.line_styles = input_dict.get('line_styles', 'solid')
        self.face_colour = input_dict.get('face_colour')
        self.edge_colour = input_dict.get('edge_colour')
        self.face_alpha = input_dict.get('face_alpha')
        self.head_width = input_dict.get('head_width')
        self.x_change = input_dict.get('x_change')
        self.y_change = input_dict.get('y_change')
        self.fill = input_dict.get('fill')
        self.length_includes_head = input_dict.get('length_includes_head')
        self.bins = input_dict.get('bins')
        self.normed = input_dict.get('normed')
        self.cmap = input_dict.get('cmap')
        self.bar_width = input_dict.get('bar_width')
        self.bar_bottom = input_dict.get('bar_bottom')
        self.bar_align = input_dict.get('bar_align')
        self.edge_colours = input_dict.get('edge_colours')
        self.step = input_dict.get('step', False)
        self.labels = input_dict.get('labels')
        self.text_alpha = input_dict.get('text_alpha')
        self.position_text_relative_to_plot = input_dict.get('position_text_relative_to_plot', False)
        self.tick_locator = self.map_z_scale_to_tick_locator(input_dict.get('z_scale'))

        missing_mandatory_args = list()
        for ma in self.mandatory_args[self.plot_type]:
            if getattr(self, ma) is None:
                missing_mandatory_args.append(ma)
        if len(missing_mandatory_args) > 0:
            mma_string = ', '.join(missing_mandatory_args)
            raise TypeError('Missing the following mandatory args: ' + mma_string)  # TODO: this should be passed upwards to create an overall report

    def call_2d_function(self, x, func_args=None):
        if func_args is None:
            return self.function(x)
        return self.function(x, func_args)

    def map_z_scale_to_tick_locator(self, z_scale=None):
        if z_scale is None:
            return None
        elif z_scale == 'log':
            return ticker.LogLocator()
        else:
            return None

    def call_3d_function(self, x, y, func_args=None):
        if func_args is None:
            return self.function(x, y)
        return self.function(x, y, func_args)

    # This function is awkwardly placed, but it needs axis as an input...
    def generate_confidence_ellipse(self, ax):
        # Adapted from https://matplotlib.org/3.1.1/gallery/statistics/confidence_ellipse.html#sphx-glr-gallery-statistics-confidence-ellipse-py
        if len(self.x_data) != len(self.y_data):
            raise ValueError('Inputs mismatched: x and y have different lengths')

        covariance = np.cov(self.x_data, self.y_data)
        pearson = covariance[0, 1]/np.sqrt(covariance[0, 0] * covariance[1, 1])
        ell_radius_x = np.sqrt(1 + pearson)
        ell_radius_y = np.sqrt(1 - pearson)
        ellipse = Ellipse(
            (0, 0),
            width=2*ell_radius_x,
            height=2*ell_radius_y,
            color=self.face_colour,
            alpha=self.face_alpha,
            fill=self.fill
        )
        scale_x = np.sqrt(covariance[0, 0]) * self.n_std
        mean_x = np.mean(self.x_data)

        scale_y = np.sqrt(covariance[1, 1]) * self.n_std
        mean_y = np.mean(self.y_data)

        transformation = transforms.Affine2D().rotate_deg(45).scale(scale_x, scale_y).translate(mean_x, mean_y)

        ellipse.set_transform(transformation + ax.transData)   # ... for this bit
        return ellipse

    def draw(self, ax, figure):
        if self.plot_type == SeriesType.function_2d:
            x = np.linspace(self.x_start, self.x_end, self.x_points)
            y = np.zeros(self.x_points)
            for index, x_point in enumerate(x):
                y[index] = self.call_2d_function(x_point, self.function_args)
            return ax.plot(x, y, color=self.line_color, ls=self.line_style, marker=self.line_marker, markersize=self.line_markersize, linewidth=self.line_linewidth)
        elif self.plot_type == SeriesType.scatter_2d:
            return ax.plot(self.x_data, self.y_data, color=self.line_color, marker=self.line_marker, ls=self.line_style, markersize=self.line_markersize, linewidth=self.line_linewidth, zorder=self.zorder)
        elif self.plot_type == SeriesType.function_3d:
            x = np.linspace(self.x_start, self.x_end, self.x_points)
            y = np.linspace(self.y_start, self.y_end, self.y_points)
            X, Y = np.meshgrid(x, y)
            Z = self.call_3d_function(X, Y, self.function_args)
            if self.fill:
                contours = ax.contourf(X, Y, Z, locator=self.tick_locator)
            else:
                contours = ax.contour(X, Y, Z, locator=self.tick_locator)
            if self.cbar_label is not None:
                figure.colorbar(contours, label=self.cbar_label)
            return contours
        elif self.plot_type == SeriesType.contour_scatter:
            if self.colour_map_colours is not None:
                cmap = clr.LinearSegmentedColormap.from_list('dummy', self.colour_map_colours, N=256)
            else:
                cmap = None
            #norm = clr.Normalize(vmin=-2, vmax = 1)
            if self.fill:
                contours = ax.contourf(self.x_data, self.y_data, self.z_data, locator=self.tick_locator, levels=self.levels, colors=self.colours, cmap=cmap, vmin=self.colour_map_min, vmax=self.colour_map_max)
            else:
                contours = ax.contour(self.x_data, self.y_data, self.z_data, locator=self.tick_locator, levels=self.levels, colors=self.colours, cmap=cmap, vmin=self.colour_map_min, vmax=self.colour_map_max)
            # TODO This bit doesn't work:
            if self.colour_below_min is not None:
                contours.get_cmap().set_under(self.colour_below_min)
            if self.colour_above_max is not None:
                contours.get_cmap().set_over(self.colour_above_max)
            if self.cbar_label is not None:
                if self.cbar_shrink is not None:
                    cbar = figure.colorbar(contours, spacing='proportional', shrink=self.cbar_shrink)
                else:
                    cbar = figure.colorbar(contours, spacing='proportional')
            if self.cbar_label is not None:
                cbar.set_label(self.cbar_label, fontsize=self.cbar_labelfontsize, rotation=self.cbar_labelrotation, labelpad=self.cbar_labelpad)
                cbar.ax.tick_params(labelsize=self.cbar_tickfontsize)
            if self.cbar_ticks is not None:
                cbar.set_ticks(self.cbar_ticks)
            if self.cbar_ticklabels is not None:
                cbar.set_ticklabels(self.cbar_ticklabels)
            return contours
        elif self.plot_type == SeriesType.tricontour_scatter:
            if self.fill:
                contours = ax.tricontourf(self.x_data, self.y_data, self.z_data, locator=self.tick_locator, levels=self.levels, colors=self.colours)
            else:
                contours = ax.tricontour(self.x_data, self.y_data, self.z_data, locator=self.tick_locator, levels=self.levels, colors=self.colours)
            # TODO This bit doesn't work:
            if self.colour_below_min is not None:
                contours.get_cmap().set_under(self.colour_below_min)
            if self.colour_above_max is not None:
                contours.get_cmap().set_over(self.colour_above_max)
            if self.cbar_label is not None:
                if self.cbar_shrink is not None:
                    cbar = figure.colorbar(contours, spacing='proportional', shrink=self.cbar_shrink)
                else:
                    cbar = figure.colorbar(contours, spacing='proportional')
            if self.cbar_label is not None:
                cbar.set_label(self.cbar_label, fontsize=self.cbar_labelfontsize, rotation=self.cbar_labelrotation, labelpad=self.cbar_labelpad)
                cbar.ax.tick_params(labelsize=self.cbar_tickfontsize)
            if self.cbar_ticks is not None:
                cbar.set_ticks(self.cbar_ticks)
            if self.cbar_ticklabels is not None:
                cbar.set_ticklabels(self.cbar_ticklabels)
            #if self.cbar_label is not None:
            #    figure.colorbar(contours, label=self.cbar_label)
            return contours
        elif self.plot_type == SeriesType.scatter_3d:
            scatter = ax.scatter(self.x_data, self.y_data, c=self.z_data)
            if self.cbar_label is not None:
                figure.colorbar(scatter, label=self.cbar_label)
            return scatter
        elif self.plot_type == SeriesType.vline:
            vline = ax.vlines(self.x_start, self.y_min, self.y_max, label=self.vline_label, linewidth=self.line_linewidth, colors=self.colours, linestyles=self.line_styles)
            return vline
        elif self.plot_type == SeriesType.hline:
            hline = ax.hlines(self.y_start, self.x_min, self.x_max, label=self.hline_label, linewidth=self.line_linewidth, colors=self.colours, linestyles=self.line_styles)
            return hline
        elif self.plot_type == SeriesType.text:
            text_args = {
                'fontsize': self.fontsize,
                'rotation': self.rotation,
                'color': self.text_colour,
                'backgroundcolor': self.shade_colour,
                'alpha': self.text_alpha,
                'horizontalalignment': self.horizontalalignment,
                'verticalalignment': self.verticalalignment,
                'transform': ax.transAxes if self.position_text_relative_to_plot else None
            }
            filtered_text_args = {key: val for key, val in text_args.items() if val is not None}
            try:
                text = ax.text(self.x_pos, self.y_pos, self.text_string, **filtered_text_args)
            except AttributeError:
                # happens if this is actually a TernaryDiagram. luckily it has an internal ax
                text = ax.ax.text(self.x_pos, self.y_pos, self.text_string, **filtered_text_args)
            if self.shade_alpha is not None:
                text.set_bbox(dict(facecolor=self.shade_colour, alpha=self.shade_alpha, edgecolor=self.shade_colour))
            return text
        elif self.plot_type == SeriesType.shade:
            shade = ax.fill_between(self.x_data, self.y_data, self.y_shade_data, where=self.shade_where, facecolor=self.shade_colour, alpha=self.shade_alpha, hatch=self.shade_hatch, zorder=self.zorder)
            return shade
        elif self.plot_type == SeriesType.scatter_2d_error:
            errors = ax.errorbar(self.x_data, self.y_data, yerr=self.y_error_data, xerr=self.x_error_data, fmt=self.line_type, capsize=self.capsize, capthick=self.capthick, color=self.line_color, marker=self.line_marker, ls=self.line_style, markersize=self.line_markersize, linewidth=self.line_linewidth, zorder=self.zorder, uplims=self.upper_limits, lolims=self.lower_limits, elinewidth=self.error_linewidth)
            return errors
        elif self.plot_type == SeriesType.ellipse:
            patch = ax.add_patch(self.generate_confidence_ellipse(ax))
            return patch
        elif self.plot_type == SeriesType.arrow:
            arrow = ax.arrow(self.x_start, self.y_start, self.x_change, self.y_change, facecolor=self.face_colour, alpha=self.face_alpha, head_width=self.head_width, length_includes_head=self.length_includes_head, linewidth=self.line_linewidth, edgecolor=self.edge_colour)
            return arrow
        elif self.plot_type == SeriesType.hist2d:
            mpl_version = mpl.__version__
            mpl_version_parse = mpl_version.split('.')
            mpl_major_version = int(mpl_version_parse[0])
            mpl_minor_version = int(mpl_version_parse[1])
            newversion = mpl_major_version >= 3 and mpl_minor_version >= 1
            # From version 3.1.0 onwards, normed was renamed to density
            if newversion:
                hist = ax.hist2d(self.x_data, self.y_data, bins=self.bins, density=self.normed, cmap=self.cmap)
            else:
                hist = ax.hist2d(self.x_data, self.y_data, bins=self.bins, normed=self.normed, cmap=self.cmap)
            if (self.cbar_label is not None) or (self.cbar_ticks is not None) or (self.cbar_ticklabels is not None):
                cbar = figure.colorbar(hist[3], ax=ax, format='%.2f')
                if self.cbar_label is not None:
                    cbar.set_label(self.cbar_label, fontsize=self.cbar_labelfontsize, rotation=self.cbar_labelrotation, labelpad=self.cbar_labelpad)
                    cbar.ax.tick_params(labelsize=self.cbar_tickfontsize)
                if self.cbar_ticks is not None:
                    cbar.set_ticks(self.cbar_ticks)
                if self.cbar_ticklabels is not None:
                    cbar.set_ticklabels(self.cbar_ticklabels)
            return hist
        elif self.plot_type == SeriesType.bar:
            if self.step:
                if self.line_style is None:
                    bar = ax.step(self.x_data, self.y_data, where='mid', color=self.colours)
                else:
                    bar = ax.step(self.x_data, self.y_data, where='mid', color=self.colours, linestyle=self.line_style)
            else:
                if self.bar_align is None:
                    bar = ax.bar(self.x_data, self.y_data, width=self.bar_width, bottom=self.bar_bottom, color=self.colours, edgecolor=self.edge_colours, alpha=self.face_alpha, hatch=self.shade_hatch, fill=self.fill, linewidth=self.line_linewidth)
                else:
                    bar = ax.bar(self.x_data, self.y_data, width=self.bar_width, bottom=self.bar_bottom, align=self.bar_align, color=self.colours, edgecolor=self.edge_colours, hatch=self.shade_hatch, fill=self.fill, linewidth=self.line_linewidth)
            return bar
        elif self.plot_type == SeriesType.pie:
            pie = ax.pie(self.x_data, labels=self.labels, colors=self.colours, startangle=90, textprops={'fontsize': self.fontsize})
            return pie
        elif self.plot_type == SeriesType.ternary_scatter:
            data_vector = list(zip(self.x_data, self.y_data, self.z_data))
            if self.line_style is None:
                if self.line_color is None:
                    scatter = ax.scatter(vector=data_vector, z=self.c_data, marker=self.line_marker, s=self.line_markersize)
                else:
                    scatter = ax.scatter(vector=data_vector, z=self.c_data, marker=self.line_marker, c=self.line_color, s=self.line_markersize)
            else:
                scatter = ax.plot(vector=data_vector, color=self.line_color)
            return scatter
        else:
            raise TypeError('Unknown plot type')
