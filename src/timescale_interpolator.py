#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum

import collections
import csv
import numpy as np
import scipy.interpolate as si

import chemistry_info as ci
import pwd_utils as pu

class TimescaleType(Enum):
    KoesterNoOvershoot = 0
    KoesterOvershoot  = 1
    BedardNoOvershoot = 2
    BedardOvershoot = 3
    BedardVariableOvershoot = 4
    Bedard3DOvershoot = 5
    MWDD = 6
    Bedard3DOvershootPatched = 7

    def __str__(self):
        return self.name

    def short_str(self):
        abbreviations = {
            'KoesterNoOvershoot': 'kn',
            'KoesterOvershoot': 'ko',
            'BedardNoOvershoot': 'bn',
            'BedardOvershoot': 'bo',
            'BedardVariableOvershoot': 'vo',
            'Bedard3DOvershoot': '3o',
            'MWDD': 'm',
            'Bedard3DOvershootPatched': '3p',
        }
        return abbreviations[self.name]

    def __len__(self):
        return 1

    def __lt__(self, other):
        return (self.value < other.value)

    def __le__(self, other):
        return(self.value <= other.value)

    def __gt__(self, other):
        return(self.value > other.value)

    def __ge__(self, other):
        return(self.value >= other.value)

    def __eq__(self, other):
        if not isinstance(other, TimescaleType):
            return False
        return (self.value == other.value)

    def __ne__(self, other):
        return not(self.__eq__(other))

    def __hash__(self):
        return self.value

logg_ranges = {
    ci.Element.H: {
        TimescaleType.BedardNoOvershoot: (7.5, 9.0),
        TimescaleType.BedardOvershoot: (7.5, 9.0),
        TimescaleType.BedardVariableOvershoot: (7.5, 9.0),
        TimescaleType.Bedard3DOvershoot: (7.5, 9.0),
        TimescaleType.KoesterNoOvershoot: (7.5, 8.5),
        TimescaleType.KoesterOvershoot: (7.5, 8.5),
        TimescaleType.MWDD: (7.5, 9.0),
        TimescaleType.Bedard3DOvershootPatched: (7.5, 9.0)
    },
    ci.Element.He: {
        TimescaleType.BedardNoOvershoot: (7.5, 9.0),
        TimescaleType.BedardOvershoot: (7.5, 9.0),
        TimescaleType.BedardVariableOvershoot: (7.5, 9.0),
        TimescaleType.Bedard3DOvershoot: (7.5, 9.0),
        TimescaleType.KoesterNoOvershoot: (7.5, 8.5),
        TimescaleType.KoesterOvershoot: (7.5, 8.5),
        TimescaleType.MWDD: (7.5, 9.0),
        TimescaleType.Bedard3DOvershootPatched: (7.5, 9.0)
    }
}

teff_ranges = {
    ci.Element.H: {
        TimescaleType.BedardNoOvershoot: (5000, 30000),
        TimescaleType.BedardOvershoot: (5000, 30000),
        TimescaleType.BedardVariableOvershoot: (5000, 30000),
        TimescaleType.Bedard3DOvershoot: (11000, 20000),
        TimescaleType.KoesterNoOvershoot: (5000, 20000),
        TimescaleType.KoesterOvershoot: (5000, 20000),
        TimescaleType.MWDD: (5000, 30000),
        TimescaleType.Bedard3DOvershootPatched: (5000, 30000)
    },
    ci.Element.He: {
        TimescaleType.BedardNoOvershoot: (7000, 30000),
        TimescaleType.BedardOvershoot: (7000, 30000),
        TimescaleType.BedardVariableOvershoot: (7000, 30000),
        TimescaleType.Bedard3DOvershoot: (np.nan, np.nan),
        TimescaleType.KoesterNoOvershoot: (3000, 15000),
        TimescaleType.KoesterOvershoot: (3000, 15000),
        TimescaleType.MWDD: (7000, 30000),
        TimescaleType.Bedard3DOvershootPatched: (7000, 30000)
    }
}

timescale_types_with_CaHe = [
    TimescaleType.KoesterNoOvershoot,
    TimescaleType.KoesterOvershoot
]

class TimescaleInterpolator():

    def __init__(self):
        # For each key in this, we will now calculate a different set of timescales
        self.file_dict = {
            TimescaleType.KoesterOvershoot: {
                ci.Element.H: {
                    7.5: 'timescales_H_g750_ov1.csv',
                    8.0: 'timescales_H_g800_ov1.csv',
                    8.5: 'timescales_H_g850_ov1.csv'
                },
                ci.Element.He: {
                    7.5: 'timescales_He_g750_ov1.csv',
                    7.75: 'timescales_He_g775_ov1.csv',
                    8: 'timescales_He_g800_ov1.csv',
                    8.25: 'timescales_He_g825_ov1.csv',
                    8.5: 'timescales_He_g850_ov1.csv'
                }
            },
            TimescaleType.KoesterNoOvershoot: {
                ci.Element.H: {
                    7.5: 'timescales_H_g750_ov0.csv',
                    8.0: 'timescales_H_g800_ov0.csv',
                    8.5: 'timescales_H_g850_ov0.csv'
                },
                ci.Element.He: {
                    7.5: 'timescales_He_g750_ov0.csv',
                    7.75: 'timescales_He_g775_ov0.csv',
                    8: 'timescales_He_g800_ov0.csv',
                    8.25: 'timescales_He_g825_ov0.csv',
                    8.5: 'timescales_He_g850_ov0.csv'
                }
            },
            TimescaleType.BedardNoOvershoot: {
                ci.Element.H: {
                    7.5: 'tdif_H_no_overshoot_abedop1_750.csv',
                    8.0: 'tdif_H_no_overshoot_abedop1_800.csv',
                    8.5: 'tdif_H_no_overshoot_abedop1_850.csv',
                    9.0: 'tdif_H_no_overshoot_abedop1_900.csv'
                },
                ci.Element.He: {
                    7.5: 'tdif_He_no_overshoot_abedop1_750.csv',
                    8.0: 'tdif_He_no_overshoot_abedop1_800.csv',
                    8.5: 'tdif_He_no_overshoot_abedop1_850.csv',
                    9.0: 'tdif_He_no_overshoot_abedop1_900.csv'
                }
            },
            TimescaleType.BedardOvershoot: {
                ci.Element.H: {
                    7.5: 'tdif_H_1hp_overshoot_abedop2_750.csv',
                    8.0: 'tdif_H_1hp_overshoot_abedop2_800.csv',
                    8.5: 'tdif_H_1hp_overshoot_abedop2_850.csv',
                    9.0: 'tdif_H_1hp_overshoot_abedop2_900.csv'
                },
                ci.Element.He: {
                    7.5: 'tdif_He_1hp_overshoot_abedop2_750.csv',
                    8.0: 'tdif_He_1hp_overshoot_abedop2_800.csv',
                    8.5: 'tdif_He_1hp_overshoot_abedop2_850.csv',
                    9.0: 'tdif_He_1hp_overshoot_abedop2_900.csv'
                }
            },
            TimescaleType.BedardVariableOvershoot: {
                ci.Element.H: {
                    7.5: 'tdif_H_variable_overshoot_abedop3_750.csv',
                    8.0: 'tdif_H_variable_overshoot_abedop3_800.csv',
                    8.5: 'tdif_H_variable_overshoot_abedop3_850.csv',
                    9.0: 'tdif_H_variable_overshoot_abedop3_900.csv'
                },
                ci.Element.He: {
                    7.5: 'tdif_He_variable_overshoot_abedop3_750.csv',
                    8.0: 'tdif_He_variable_overshoot_abedop3_800.csv',
                    8.5: 'tdif_He_variable_overshoot_abedop3_850.csv',
                    9.0: 'tdif_He_variable_overshoot_abedop3_900.csv'
                }
            },
            TimescaleType.Bedard3DOvershoot: {
                ci.Element.H: {
                    7.5: 'tdif_H_3D_overshoot_abedop4_750.csv',
                    8.0: 'tdif_H_3D_overshoot_abedop4_800.csv',
                    8.5: 'tdif_H_3D_overshoot_abedop4_850.csv',
                    9.0: 'tdif_H_3D_overshoot_abedop4_900.csv'
                }
                #ci.Element.He: dict() # 3D data available for H-dominated only
            },
            TimescaleType.MWDD: {
                ci.Element.H: {
                    7.5: 'MWDD_compiled_H_750.csv',
                    8.0: 'MWDD_compiled_H_800.csv',
                    8.5: 'MWDD_compiled_H_850.csv',
                    9.0: 'MWDD_compiled_H_900.csv'
                },
                ci.Element.He: {
                    7.5: 'MWDD_compiled_He_750.csv',
                    8.0: 'MWDD_compiled_He_800.csv',
                    8.5: 'MWDD_compiled_He_850.csv',
                    9.0: 'MWDD_compiled_He_900.csv'
                }
            }
        }
        print('Reminder! For the 3D option, should we default to the no overshoot case if we land outside the 3D grid?')
        self.all_timescale_data = self.load_data()
        self.wd_data = None
        self.expected_vals = dict()
        for timescale_name in self.all_timescale_data.keys():
            self.expected_vals[timescale_name] = {ci.Element.H: dict(), ci.Element.He: dict()}
        self.set_up_interpolators()

    @staticmethod
    def perform_format(input_file_name, output_file_name):
        header_written = False
        file_logg = None
        file_overshoot = None
        print('Opening input file ' + input_file_name + ' and output file ' + output_file_name)
        with open(input_file_name, 'r') as input_file, open(output_file_name, 'w', newline='', encoding='utf-8') as output_file:
            read = csv.reader(input_file, delimiter='|')  # Intentionally using a delimiter that isn't present
            to_write = csv.writer(output_file)
            for row in read:
                try:
                    line = row[0].lstrip()  # This should contain the whole line as a str, minus leading whitespace
                except IndexError:
                    # Row was empty
                    continue
                if line.startswith('EL'):
                    # Read Teff, log(g) and overshoot from this line
                    Teff = line.split('Teff[K] =')[1].split('log')[0].strip()
                    logg = line.split('log g =')[1].split('log')[0].strip()
                    overshoot = line.split('including')[1].split('Hp')[0].strip()
                    if not header_written:
                        # Write a header saying what values of log(g) and overshoot this corresponds to
                        to_write.writerow(['logg', logg, 'overshoot', overshoot])
                        file_logg = logg
                        file_overshoot = overshoot
                        header_written = True
                    else:
                        # Check that log(g) and overshoot are consistent
                        assert logg == file_logg
                        assert overshoot == file_overshoot
                    to_write.writerow(['Teff', Teff])
                elif line.startswith('Z'): # This line is actually showing the [Ca/He] values
                    cahe_values = line.replace('  ', ' ').split(' ')
                    cahe_values[0] = 'CaHe'
                    to_write.writerow(cahe_values)
                elif line.startswith('qcvz'): # This line is showing the log(q) values
                    logq_values = line.split('  ')
                    to_write.writerow(logq_values)
                else:
                    # Test to see if first element is an int
                    values = line.split('   ')
                    try:
                        test_val = int(values[0])
                        to_write.writerow(values)
                    except ValueError:
                        # Then this wasn't a row containing element data
                        pass
        print('Done!')

    @staticmethod
    def perform_format_abed(input_file_name):
        # NB These files contain all log(g) values in one, whereas we need a separate file for each
        output_dict = dict() # keys will be log(g) values
        numerical_entries = list()
        string_entries = list()
        print('Opening input file ' + input_file_name)
        with open(input_file_name, 'r') as input_file:
            read = csv.reader(input_file, delimiter=' ')
            #to_write = csv.writer(output_file)
            for row in read:
                # Strategy is to read each line up to the next blank line, saving strings and numbers in 2 separate lists, then pair them afterwards
                if len(row) == 0:
                    pairs_dict = dict(zip(string_entries, numerical_entries))
                    TimescaleInterpolator.process_pairs_dict(pairs_dict, output_dict)
                    numerical_entries = list()
                    string_entries = list()
                for entry in row:
                    if entry == '':
                        continue
                    else:
                        is_numerical = True
                        try:
                            numerical_entry = float(entry)
                        except ValueError:
                            is_numerical = False
                        if is_numerical:
                            numerical_entries.append(numerical_entry)
                        else:
                            string_entries.append(entry)
            pairs_dict = dict(zip(string_entries, numerical_entries))
            TimescaleInterpolator.process_pairs_dict(pairs_dict, output_dict) # Call this one last time because there's no extra blank line at EOF
            for logg, teff_dict in output_dict.items():
                output_file_name = input_file_name.split('raw')[0] + str(int(100*logg)) + '.csv'
                list_of_Teffs = list()
                qcvz_list = list()
                element_value_lists = dict()
                for teff, el_dict in teff_dict.items():
                    list_of_Teffs.append(int(teff))
                    qcvz_list.append(el_dict['qcvz'])
                    for element, value in el_dict.items():
                        if isinstance(element, ci.Element):
                            if element.value not in element_value_lists:
                                element_value_lists[element.value] = list()
                            element_value_lists[element.value].append(value)
                print('Writing to ' + output_file_name)
                with open(output_file_name, 'w', newline='', encoding='utf-8') as output_file:
                    to_write = csv.writer(output_file)
                    to_write.writerow(['T:'] + list_of_Teffs)
                    to_write.writerow(['qcvz:'] + qcvz_list)
                    for element_no, list_of_values in element_value_lists.items():
                        to_write.writerow([str(element_no)] + list_of_values)
        print('Done!')

    @staticmethod
    def process_pairs_dict(pairs_dict, output_dict): # Alters output_dict in place
        if len(pairs_dict) == 0:
            return
        current_logg = pairs_dict['lg']
        current_teff = pairs_dict['Teff']
        if current_logg not in output_dict:
            output_dict[current_logg] = dict()
        if current_teff not in output_dict[current_logg]:
            output_dict[current_logg][current_teff] = dict()
        for param, value in pairs_dict.items():
            if param not in ['lg', 'Teff']:
                if param == 'lqc':
                    output_dict[current_logg][current_teff]['qcvz'] = value
                else:
                    #param should be of the form lt(<Element Symbol>)
                    element_symbol = param.split('lt(')[1].split(')')[0]
                    element = ci.Element[element_symbol]
                    output_dict[current_logg][current_teff][element] = value

    def get_csv_format_type(self, timescale_type, HorHe):
        # Different csv files are formatted differently - this just returns a flag to let us know how the particular request csv is formatted
        # 0 is basically default
        # 1 has CaHe as an additional variable
        if timescale_type in timescale_types_with_CaHe and HorHe == ci.Element.He:
            return 1
        return 0

    def get_arbitrary_val(self, timescale_name, HorHe, variable):
        return self.expected_vals[timescale_name][HorHe][variable][0]

    def get_values_for_variable(self, timescale_name, HorHe, variable):
        if variable == 'g':
            vals = np.array(list(self.all_timescale_data[timescale_name][HorHe].keys()))
        elif variable == 't':
            vals = np.array(list(self.all_timescale_data[timescale_name][HorHe][self.get_arbitrary_val(timescale_name, HorHe, 'g')].keys()))
        elif variable == 'c':
            vals_list = list(self.all_timescale_data[timescale_name][HorHe][self.get_arbitrary_val(timescale_name, HorHe, 'g')][self.get_arbitrary_val(timescale_name, HorHe, 't')].keys())
            if 'logq' in vals_list:
                # Then this was a file set which didn't have CaHe as a variable! We went straight to the elements.
                print('Warning! Variable CaHe was not present. Returning None')
                return None
            vals = np.array(vals_list)
        else:
            print('Warning! Unrecognised variable ' + str(variable))
            return None
        vals.sort()
        if len(vals) < 2:
            # Can't really do interpolation
            print('Warning! received invalid interpolation values for variable ' + str(variable) + ':')
            print(vals)
            return None
        if self.expected_vals[timescale_name][HorHe].get(variable) is None:
            self.expected_vals[timescale_name][HorHe][variable] = vals
        else:
            if self.expected_vals[timescale_name][HorHe][variable] != vals:
                print('Warning! received invalid interpolation values for variable ' + str(variable) + ':')
                print(vals)
                return None
        return vals

    def set_up_interpolators(self):
        # This assumes that self.timescale_data has no missing keys anywhere or anything like that
        self.interpolators = dict()
        for timescale_name, timescale_data in self.all_timescale_data.items():
            self.interpolators[timescale_name] = dict()
            for HorHe in timescale_data.keys():
                self.interpolators[timescale_name][HorHe] = dict()
                # You have to call them in this order!
                g_vals = self.get_values_for_variable(timescale_name, HorHe, 'g')
                if g_vals is None:
                    continue  # This needs to be present
                t_vals = self.get_values_for_variable(timescale_name, HorHe, 't')
                if t_vals is None:
                    continue  # This needs to be present
                c_vals = self.get_values_for_variable(timescale_name, HorHe, 'c')
                if c_vals is None:
                    # Then CaHe was not a variable. Set up a 2D interpolator
                    for element in timescale_data[HorHe][self.get_arbitrary_val(timescale_name, HorHe, 'g')][self.get_arbitrary_val(timescale_name, HorHe, 't')].keys():
                        grid_vals = np.zeros((len(g_vals), len(t_vals)))
                        for i in range(len(g_vals)):
                            for j in range(len(t_vals)):
                                grid_vals[i,j] = timescale_data[HorHe][g_vals[i]][t_vals[j]][element]
                        # 'linear', False, None means linear interpolation, no error if we go out of bounds, extrapolate in that case
                        self.interpolators[timescale_name][HorHe][element] = si.RegularGridInterpolator((g_vals, t_vals), grid_vals, 'linear', False, None)
                else:
                    for element in timescale_data[HorHe][self.get_arbitrary_val(timescale_name, HorHe, 'g')][self.get_arbitrary_val(timescale_name, HorHe, 't')][self.get_arbitrary_val(timescale_name, HorHe, 'c')].keys():
                        grid_vals = np.zeros((len(g_vals), len(t_vals), len(c_vals)))
                        for i in range(len(g_vals)):
                            for j in range(len(t_vals)):
                                for k in range(len(c_vals)):
                                    grid_vals[i,j,k] = timescale_data[HorHe][g_vals[i]][t_vals[j]][c_vals[k]][element]
                        # 'linear', False, None means linear interpolation, no error if we go out of bounds, extrapolate in that case
                        self.interpolators[timescale_name][HorHe][element] = si.RegularGridInterpolator((g_vals, t_vals, c_vals), grid_vals, 'linear', False, None)

    def load_data(self):
        toret = dict()
        # Structure to build: data[HorHe][g][Teff][CaHe] = {'logq': <q>, Element.Li: <t_Li>, ... Element.Zn: <t_Zn> }
        for timescale_set_name, HorHe_dict in self.file_dict.items():
            toret[timescale_set_name] = dict()
            for HorHe, HorHe_files in HorHe_dict.items():
                toret[timescale_set_name][HorHe] = dict()
                for g, input_file in HorHe_files.items():
                    toret[timescale_set_name][HorHe][g] = dict()
                    current_Teff = None
                    current_CaHe_vals = None
                    try:
                        with open('../data/' + input_file) as csvfile:
                            read = csv.reader(csvfile, delimiter=',')
                            if self.get_csv_format_type(timescale_set_name, HorHe) == 0:
                                i = 0
                                for row in read:
                                    if i == 0:
                                        T_vals = row[1:]
                                        for j in range(1, len(row)):
                                            toret[timescale_set_name][HorHe][g][int(T_vals[j-1])] = dict()
                                    elif i == 1:
                                        q_vals = row[1:]
                                        for j in range(1, len(row)):
                                            toret[timescale_set_name][HorHe][g][int(T_vals[j-1])]['logq'] = float(q_vals[j-1])
                                    else:
                                        element = ci.Element(int(row[0]))
                                        for j in range(1, len(row)):
                                            toret[timescale_set_name][HorHe][g][int(T_vals[j-1])][element] = float(row[j])
                                    i += 1
                            else:
                                for row in read:
                                    if row[0] == 'Teff':
                                        current_Teff = int(row[1])
                                        toret[timescale_set_name][HorHe][g][current_Teff] = dict()
                                    elif row[0] == 'CaHe':
                                        current_CaHe_vals = row[1:]
                                        for j in range(1, len(row)):
                                            toret[timescale_set_name][HorHe][g][current_Teff][float(current_CaHe_vals[j-1])] = dict()
                                    elif row[0] == 'qcvz':
                                        q_vals = row[1:]
                                        for j in range(1, len(row)):
                                            toret[timescale_set_name][HorHe][g][current_Teff][float(current_CaHe_vals[j-1])]['logq'] = float(q_vals[j-1])
                                    elif row[0] == 'logg':
                                        claimed_logg = float(row[1])
                                        assert g == claimed_logg
                                    else:
                                        element = ci.Element(int(row[0]))
                                        for j in range(1, len(row)):
                                            toret[timescale_set_name][HorHe][g][current_Teff][float(current_CaHe_vals[j-1])][element] = float(row[j])
                    except (TypeError, FileNotFoundError):
                        print('Warning! Could not open an input file')
        toret = self.load_3d_patch(toret)
        return toret

    def load_3d_patch(self, timescale_data):
        to_add = dict()
        for g_val in timescale_data[TimescaleType.Bedard3DOvershoot][ci.Element.H].keys():
            to_add[g_val] = dict()
            for temp, data_vals in timescale_data[TimescaleType.Bedard3DOvershoot][ci.Element.H][g_val].items():
                to_add[g_val][temp] = data_vals
            for temp, data_vals in timescale_data[TimescaleType.BedardVariableOvershoot][ci.Element.H][g_val].items():
                if temp not in to_add[g_val].keys():
                    to_add[g_val][temp] = data_vals
        timescale_data[TimescaleType.Bedard3DOvershootPatched] = {ci.Element.H: to_add}
        return timescale_data

    def load_wd_data(self, wd_data_file='WDInputData.csv'):
        toret = collections.OrderedDict()
        with open(pu.get_path_to_data() + pu.get_wd_input_file(), encoding='utf-8') as wdcsv:
            for row in csv.DictReader(wdcsv):
                wd_name = row['Star']
                wd_name_tag = row['Variant'] # TODO: Add this to the white dwarf class, so that we can either use the name or the full name on demand
                if wd_name_tag is not None and wd_name_tag != '':
                    wd_name += '_' + wd_name_tag
                wd_type_raw = row['atmosphere']
                try:
                    wd_type = ci.Element[wd_type_raw]
                except KeyError:
                    raise KeyError('Invalid atmospheric type: must be H or He, received ' + str(wd_type_raw))
                Teff = int(row['T_eff'])
                try:
                    logg = float(row['logg'])
                except ValueError:
                    logg = pu.get_default_logg()
                try:
                    CaHe = float(row['log(' + str(ci.Element.Ca) + '/H(e))'])
                except ValueError:
                    CaHe = None
                toret[wd_name] = {'Type': wd_type, 'Teff': Teff, 'logg': logg, 'CaHe': CaHe}
        return toret

    def process_wd_data(self, wd_data_file='WDInputData.csv'):
        self.wd_data = self.load_wd_data(wd_data_file)
        toret = collections.OrderedDict()
        for wd_name, wd_entry in self.wd_data.items():
            HorHe = wd_entry['Type']
            Teff = wd_entry['Teff']
            logg = wd_entry['logg']
            CaHe = wd_entry['CaHe']
            timescales = self.extract_timescales(HorHe, logg, Teff, CaHe)
            toret[wd_name] = timescales
        return toret

    def extract_timescales(self, HorHe, logg, Teff, CaHe=None):
        toret = dict()
        for timescale_name, timescale_grids in self.all_timescale_data.items():
            try:
            # Firstly, we check to see if we happen to land on an exact grid point (i.e., no need to interpolate)
                if self.get_csv_format_type(timescale_name, HorHe) == 1:
                    toret[timescale_name] = timescale_grids[timescale_name][HorHe][logg][Teff][CaHe]
                else:
                    toret[timescale_name] = timescale_grids[timescale_name][HorHe][logg][Teff]
            except KeyError:
                pass
            # Need to interpolate!
            # Firstly, check we actually have an appropriate interpolator!
            if HorHe not in self.interpolators[timescale_name].keys() or self.interpolators[timescale_name][HorHe] == dict():
                #print('Warning! Could not interpolate. Did not find an interpolator for ' + str(HorHe))
                toret[timescale_name] = None
            toret[timescale_name] = dict()
            if self.get_csv_format_type(timescale_name, HorHe) == 1:
                if CaHe is None:
                    print('Warning! Could not find timescales for ' + str(timescale_name) + ', ' + str(HorHe) + ' because CaHe was None')
                    toret[timescale_name] = None
                else:
                    try:
                        timescale_grid = timescale_grids[HorHe]
                    except KeyError:
                        toret[timescale_name] = None
                        continue
                    for element in timescale_grid[self.get_arbitrary_val(timescale_name, HorHe, 'g')][self.get_arbitrary_val(timescale_name, HorHe, 't')][self.get_arbitrary_val(timescale_name, HorHe, 'c')].keys():
                        interpolator = self.interpolators[timescale_name][HorHe][element]
                        point_to_sample = np.array([[logg, Teff, CaHe]])
                        toret[timescale_name][element] = interpolator(point_to_sample)[0]  # The float value is wrapped in a numpy array so indexing [0] extracts what we actually want
            else:
                try:
                    timescale_grid = timescale_grids[HorHe]
                except KeyError:
                    toret[timescale_name] = None
                    continue
                for element in timescale_grid[self.get_arbitrary_val(timescale_name, HorHe, 'g')][self.get_arbitrary_val(timescale_name, HorHe, 't')].keys():
                    interpolator = self.interpolators[timescale_name][HorHe][element]
                    point_to_sample = np.array([[logg, Teff]])
                    toret[timescale_name][element] = interpolator(point_to_sample)[0]  # The float value is wrapped in a numpy array so indexing [0] extracts what we actually want
        return toret

    def get_wd_timescales(self, HorHe, logg, Teff, CaHe=None, all_timescales=False): # TODO: could HorHe be a ci.Element rather than a str?
        timescales = self.extract_timescales(HorHe, logg, Teff, CaHe)
        return self.return_wd_timescales_as_dict(timescales, all_timescales)

    def return_wd_timescales_as_dict(self, wd_entry, all_timescales=False):
        toret = dict()
        for timescale_name in self.all_timescale_data.keys():
            if wd_entry[timescale_name] is None:
                toret[timescale_name] = None
            else:
                toret[timescale_name] = {'logq': wd_entry[timescale_name]['logq']} # TODO: replace with mp.WDParameter.logq?
                for el in ci.all_elements:
                    if all_timescales or el in ci.usual_elements:
                        try:
                            toret[timescale_name][el] = 10**wd_entry[timescale_name][el]
                        except KeyError:
                            pass
        return(toret)

    def return_wd_timescales_as_dict_of_lists(self, wd_entry):
        toret = dict()
        for timescale_name in self.timescale_data.keys():
            toret[timescale_name] = list()
            toret[timescale_name].append(wd_entry[timescale_name]['logq'])
            for el in ci.usual_elements:
                try:
                    entry_to_use = 10**wd_entry[timescale_name][el]
                except KeyError:
                    entry_to_use = 0.0
                toret[timescale_name].append(entry_to_use)
        return(toret)

    def return_wd_timescales_as_list(self, wd_entry):
        toret = list()
        toret.append(wd_entry['logq'])
        for el in ci.usual_elements:
            try:
                entry_to_use = 10**wd_entry[el]
            except KeyError:
                entry_to_use = 0.0
            toret.append(entry_to_use)
        return(toret)

    def dump_wd_timescales(self, wd_timescales, outfile='wd_timescales'):
        if outfile is not None:
            print('Writing to ' + outfile + '.csv')
            with open(outfile + '.csv', 'w', newline='', encoding='utf-8') as f:
                to_write = csv.writer(f)
                to_write.writerow(['WD Name', 'Timescale Name', 'log(q)', 't_Al', 't_Ti', 't_Ca', 't_Ni', 't_Fe', 't_Cr', 't_Mg', 't_Si', 't_Na', 't_O', 't_C', 't_N'])
                for wd_name, timescale_data in wd_timescales.items():
                    for timescale_name, wd_entry in timescale_data.items():
                        if wd_entry is None:
                            to_write.writerow([
                                wd_name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                            ])
                        else:
                            output_row = [wd_name] + [timescale_name] + self.return_wd_timescales_as_list(wd_entry)
                            to_write.writerow(output_row)

def create_grid():
    timescale_interpolator = TimescaleInterpolator()
    timescale_type = TimescaleType.KoesterNoOvershoot
    elements = [ci.Element.Ca, ci.Element.Fe, ci.Element.Mg]
    Hx_values = [ci.Element.H, ci.Element.He]
    logg_values = [7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5]
    teff_values = [4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000]
    CaHe = -15
    timescales = dict()
    for Hx in Hx_values:
        timescales[Hx] = dict()
        for logg in logg_values:
            timescales[Hx][logg] = dict()
            for teff in teff_values:
                timescales[Hx][logg][teff] = dict()
                timescales_e = timescale_interpolator.get_wd_timescales(Hx, logg, teff, CaHe)
                for test_el in elements:
                    timescales[Hx][logg][teff][test_el] = timescales_e[timescale_type][test_el]
    for Hx in Hx_values:
        print()
        print(Hx)
        for element in elements:
            print(','.join([str(element)] + [str(teff) for teff in teff_values]))
            for logg in logg_values:
                row_to_write = [str(logg)] + [str(timescales[Hx][logg][teff][element]) for teff in teff_values]
                print(','.join(row_to_write))
            print()
            print()

def main():
    timescale_interpolator = TimescaleInterpolator()
    HorHe = ci.Element.H
    logg = 8
    Teff = 6000
    CaHe = -15
    print('Example timescales for a ' + str(HorHe) + ' WD with log(g) = ' + str(logg) + ', Teff = ' + str(Teff) + ' and Ca/He = ' + str(CaHe))
    result = timescale_interpolator.get_wd_timescales(HorHe, logg, Teff, CaHe, True)
    print(result)
    for timescale_type in TimescaleType:
        print()
        print(timescale_type)
        for element in ci.usual_elements:
            try:
                print(str(element) + ': ' + str(np.log10(result[timescale_type][element])))
            except TypeError:
                print(str(element) + ': None')
    #for timescale_type in TimescaleType:
    #    print()
    #    print(timescale_type)
    #    print(result[timescale_type][ci.Element.Mg]/result[timescale_type][ci.Element.Si])
    #    print(result[timescale_type][ci.Element.Fe]/result[timescale_type][ci.Element.Si])
    #    print(result[timescale_type][ci.Element.Mg]/result[timescale_type][ci.Element.Fe])
    #interpolated_timescales = timescale_interpolator.process_wd_data()
    #print(interpolated_timescales)
    #print('Call timescale_interpolator.dump_wd_timescales(interpolated_timescales) to dump these timescales into a file')

if __name__ == '__main__':
    main()
