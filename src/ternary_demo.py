import matplotlib.pyplot as plt
import sys
from ternary_diagram import TernaryDiagram

import pwd_utils as pu

sys.path.append(pu.get_path_to_utils())

import dict_plotter as dp

def golden_example():
    # You can set `ax` to select which axes to draw. If not, the current axes will be used.
    td = TernaryDiagram(["Li2O", "La2O3", "TiO2"])

    # scatter
    td.scatter(vector=[[1, 1, 1], [1, 2, 3]], z=[0, 1])
    # You can set some options in `plt.scatter` like `marker`, `c` etc.
    td.scatter(vector=[[2, 1, 3], [3, 2, 1]], marker="s", c="#022c5e", s=30)

    # line plot
    # You can set some options in `plt.plot` like `lw`, `c`, and so on.
    td.plot([[1, 1, 1], [1, 2, 3]], color="black")

    # save figure
    td.fig.savefig(pu.get_path_to_default_graphs() + "testternarygolden.pdf", dpi=144)

def mimic_golden_example():
    x_data = [1, 1]
    y_data = [1, 2]
    z_data = [1, 3]
    c_data = [0, 1]
    series_dict = {
        'test_ternary1': {
            'type': dp.SeriesType.ternary_scatter,
            'x_data': [1, 1],
            'y_data': [1, 2],
            'z_data': [1, 3],
            'c_data': [0, 1],
            'legend': True
        },
        'test_ternary2': {
            'type': dp.SeriesType.ternary_scatter,
            'x_data': [2, 3],
            'y_data': [1, 2],
            'z_data': [3, 1],
            'line_color': "#022c5e",
            'line_markersize': 30,
            'line_marker': 's',
            'legend': True
        },
        'test_ternary3': {
            'type': dp.SeriesType.ternary_scatter,
            'x_data': [1, 1],
            'y_data': [1, 2],
            'z_data': [1, 3],
            'line_style': '-',
            'line_color': 'black',
            'legend': True
        }
    }
    plot_dict = {
        'test_ternary_plot': {
            'show': False,
            'filenames': ['testternary.pdf'],
            'subplots': {
                'subplot1': {
                    'subplot_region': 121,
                    'legend': True,
                    'type': dp.SubplotType.ternary,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'font': 'STIXGeneral',
                    'series': series_dict,
                    'x_label': 'Li2O',
                    'y_label': 'La2O3',
                    'z_label': 'TiO2'
                },
                'subplot2': {
                    'subplot_region': 122,
                    'legend': True,
                    'type': dp.SubplotType.cartesian,
                    'legend_loc': 'best',
                    'legend_text_size': 8,
                    'font': 'STIXGeneral',
                    'series': {
                        'dummy': {
                            'type': dp.SeriesType.scatter_2d,
                            'x_data': [0,1],
                            'y_data': [4,5],
                            'legend': True
                        }
                    }
                }
            }
        }
    }
    plotter = dp.DictPlotter(plot_dict)
    plotter.draw()
    plotter.yield_output(pu.get_path_to_default_graphs())

def main():
    golden_example()
    mimic_golden_example()

if __name__ == '__main__':
    main()
