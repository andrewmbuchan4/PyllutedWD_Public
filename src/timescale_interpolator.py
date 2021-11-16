#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import csv
import numpy as np
import scipy.interpolate as si

import chemistry_info as ci

class TimescaleInterpolator():
    
    def __init__(self):
        self.file_dict = {
            'H': {
                7.5: 'timescales_H_g750_ov1.csv',
                8.0: 'timescales_H_g800_ov1.csv',
                8.5: 'timescales_H_g850_ov1.csv'
            },
            'He': {
                7.5: 'timescales_He_g750_ov1.csv',
                7.75: 'timescales_He_g775_ov1.csv',
                8: 'timescales_He_g800_ov1.csv',
                8.25: 'timescales_He_g825_ov1.csv',
                8.5: 'timescales_He_g850_ov1.csv'
            }
        }
        self.timescale_data = self.load_data()
        self.wd_data = None
        self.expected_vals = {'H': dict(), 'He': dict()}
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

    def get_arbitrary_val(self, HorHe, variable):
        return self.expected_vals[HorHe][variable][0]
    
    def get_values_for_variable(self, HorHe, variable):
        if variable == 'g':
            vals = np.array(list(self.timescale_data[HorHe].keys()))
        elif variable == 't':
            vals = np.array(list(self.timescale_data[HorHe][self.get_arbitrary_val(HorHe, 'g')].keys()))
        elif variable == 'c':
            vals_list = list(self.timescale_data[HorHe][self.get_arbitrary_val(HorHe, 'g')][self.get_arbitrary_val(HorHe, 't')].keys())
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
            print('Warning! received invalid vals for variable ' + str(variable) + ':')
            print(vals)
            return None
        if self.expected_vals[HorHe].get(variable, None) is None:
            self.expected_vals[HorHe][variable] = vals
        else:
            if self.expected_vals[HorHe][variable] != vals:
                print('Warning! received invalid vals for variable ' + str(variable) + ':')
                print(vals)
                return None
        return vals
    
    def set_up_interpolators(self):
        # This assumes that self.timescale_data has no missing keys anywhere or anything like that
        self.interpolators = dict()
        for HorHe in self.timescale_data.keys():
            self.interpolators[HorHe] = dict()
            # You have to call them in this order!
            g_vals = self.get_values_for_variable(HorHe, 'g')
            if g_vals is None:
                continue  # This needs to be present
            t_vals = self.get_values_for_variable(HorHe, 't')
            if t_vals is None:
                continue  # This needs to be present
            c_vals = self.get_values_for_variable(HorHe, 'c')
            if c_vals is None:
                # Then CaHe was not a variable. Set up a 2D interpolator
                for element in self.timescale_data[HorHe][self.get_arbitrary_val(HorHe, 'g')][self.get_arbitrary_val(HorHe, 't')].keys():
                    grid_vals = np.zeros((len(g_vals), len(t_vals)))
                    for i in range(len(g_vals)):
                        for j in range(len(t_vals)):
                            grid_vals[i,j] = self.timescale_data[HorHe][g_vals[i]][t_vals[j]][element]
                    # 'linear', False, None means linear interpolation, no error if we go out of bounds, extrapolate in that case
                    self.interpolators[HorHe][element] = si.RegularGridInterpolator((g_vals, t_vals), grid_vals, 'linear', False, None)
            else:
                for element in self.timescale_data[HorHe][self.get_arbitrary_val(HorHe, 'g')][self.get_arbitrary_val(HorHe, 't')][self.get_arbitrary_val(HorHe, 'c')].keys():
                    grid_vals = np.zeros((len(g_vals), len(t_vals), len(c_vals)))
                    for i in range(len(g_vals)):
                        for j in range(len(t_vals)):
                            for k in range(len(c_vals)):
                                grid_vals[i,j,k] = self.timescale_data[HorHe][g_vals[i]][t_vals[j]][c_vals[k]][element]
                    # 'linear', False, None means linear interpolation, no error if we go out of bounds, extrapolate in that case
                    self.interpolators[HorHe][element] = si.RegularGridInterpolator((g_vals, t_vals, c_vals), grid_vals, 'linear', False, None)

    def load_data(self):
        toret = dict()
        # Structure to build: data[HorHe][g][Teff][CaHe] = {'logq': <q>, Element.Li: <t_Li>, ... Element.Zn: <t_Zn> }
        for HorHe, HorHe_files in self.file_dict.items():
            toret[HorHe] = dict()
            for g, input_file in HorHe_files.items():
                toret[HorHe][g] = dict()
                current_Teff = None
                current_CaHe_vals = None
                try:
                    with open('../data/' + input_file) as csvfile:
                        read = csv.reader(csvfile, delimiter=',')
                        if HorHe == 'H':
                            i = 0
                            for row in read:
                                if i == 0:
                                    T_vals = row[1:]
                                    for j in range(1, len(row)):
                                        toret[HorHe][g][int(T_vals[j-1])] = dict()
                                elif i == 1:
                                    q_vals = row[1:]
                                    for j in range(1, len(row)):
                                        toret[HorHe][g][int(T_vals[j-1])]['logq'] = float(q_vals[j-1])
                                else:
                                    element = ci.Element(int(row[0]))
                                    for j in range(1, len(row)):
                                        toret[HorHe][g][int(T_vals[j-1])][element] = float(row[j]) 
                                i += 1
                        elif HorHe == 'He':
                            for row in read:
                                if row[0] == 'Teff':
                                    current_Teff = int(row[1])
                                    toret[HorHe][g][current_Teff] = dict()
                                elif row[0] == 'CaHe':
                                    current_CaHe_vals = row[1:]
                                    for j in range(1, len(row)):
                                        toret[HorHe][g][current_Teff][float(current_CaHe_vals[j-1])] = dict()
                                elif row[0] == 'qcvz':
                                    q_vals = row[1:]
                                    for j in range(1, len(row)):
                                        toret[HorHe][g][current_Teff][float(current_CaHe_vals[j-1])]['logq'] = float(q_vals[j-1])
                                elif row[0] == 'logg':
                                    claimed_logg = float(row[1])
                                    assert g == claimed_logg
                                else:
                                    element = ci.Element(int(row[0]))
                                    for j in range(1, len(row)):
                                        toret[HorHe][g][current_Teff][float(current_CaHe_vals[j-1])][element] = float(row[j])
                        else:
                            print('Unrecognised HorHe value: ' + str(HorHe) + ', could not read input')
                except (TypeError, FileNotFoundError):
                    print('Warning! Could not open an input file')
        return toret
        
    def load_wd_data(self, wd_data_file='BlouinConglomNewTimescales.csv'):
        toret = collections.OrderedDict()
        with open('../data/' + wd_data_file) as csvfile:
            read = csv.reader(csvfile, delimiter=',')
            i = 0
            for row in read:
                if i == 0:
                    pass  # Header row
                else:
                    wd_name = row[0]
                    wd_type = row[1]
                    Teff = int(row[3])
                    try:
                        CaHe = float(row[10])
                    except ValueError:
                        CaHe = None # Likely means Ca was an upper bound
                    if CaHe == 0:
                        CaHe = None
                    try:
                        logg = float(row[4])
                    except ValueError:
                        # Assume 8 by default
                        logg = 8.0
                        print('Warning: ' + wd_name + ' has no log(g) value. Assuming ' + str(logg))
                    toret[wd_name] = {'Type': wd_type, 'Teff': Teff, 'logg': logg, 'CaHe': CaHe}
                i += 1
        return toret
        
    def process_wd_data(self, wd_data_file='BlouinConglomNewTimescales.csv'):
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
    
    def extract_timescales(self, HorHe, logg, Teff, CaHe=None):  # CaHe can be None for HorHe == H
        if HorHe == 'He' and CaHe == 0:
            # We need to careful! CaHe == 0 means that Ca is actually not present and we shouldn't attempt to find timescales
            # This only applies to He though because H doesn't include CaHe as a dimension
            print('This is a He WD but no Ca is present. Therefore interpolation is not possible. Returning None')
            return None
        try:
            if HorHe == 'He':
                return self.timescale_data[HorHe][logg][Teff][CaHe]
            elif HorHe == 'H':
                return self.timescale_data[HorHe][logg][Teff]
            else:
                print('Unrecognised HorHe value: ' + str(HorHe) + ', could not find timescales')
                return None
        except KeyError:
            pass
        # Need to interpolate!
        # Firstly, check we actually have an appropriate interpolator!
        if HorHe not in self.interpolators.keys() or self.interpolators[HorHe] == dict():
            print('Warning! Could not interpolate. Did not find an interpolator for ' + str(HorHe))
            return None
        toret = dict()
        if HorHe == 'H':
            for element in self.timescale_data[HorHe][self.get_arbitrary_val(HorHe, 'g')][self.get_arbitrary_val(HorHe, 't')].keys():
                interpolator = self.interpolators[HorHe][element]
                point_to_sample = np.array([[logg, Teff]])
                toret[element] = interpolator(point_to_sample)[0]  # The float value is wrapped in a numpy array so indexing [0] extracts what we actually want
        elif HorHe == 'He':
            for element in self.timescale_data[HorHe][self.get_arbitrary_val(HorHe, 'g')][self.get_arbitrary_val(HorHe, 't')][self.get_arbitrary_val(HorHe, 'c')].keys():
                interpolator = self.interpolators[HorHe][element]
                point_to_sample = np.array([[logg, Teff, CaHe]])
                toret[element] = interpolator(point_to_sample)[0]  # The float value is wrapped in a numpy array so indexing [0] extracts what we actually want
        else:
            print('Unrecognised HorHe value: ' + str(HorHe) + ', could not find timescales')
            return None
        return toret

    def get_wd_timescales(self, HorHe, logg, Teff, CaHe): # TODO: Consider making CaHe=None by default? Also, could HorHe be a ci.Element rather than a str? (Or better yet, a duck type str)
        timescales = self.extract_timescales(HorHe, logg, Teff, CaHe)
        if timescales is None:
            return None
        return self.return_wd_timescales_as_dict(timescales)

    def return_wd_timescales_as_dict(self, wd_entry):
        toret = dict()
        toret['logq'] = wd_entry['logq']
        for el in ci.usual_elements:
            try:
                entry_to_use = 10**wd_entry[el]
            except KeyError:
                entry_to_use = 0.0
            toret[el] = entry_to_use
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
                to_write.writerow(['WD Name', 'log(q)', 't_Al', 't_Ti', 't_Ca', 't_Ni', 't_Fe', 't_Cr', 't_Mg', 't_Si', 't_Na', 't_O', 't_C', 't_N'])
                for wd_name, wd_entry in wd_timescales.items():
                    if wd_entry is None:
                        to_write.writerow([
                            wd_name, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                        ])
                    else:
                        output_row = [wd_name] + self.return_wd_timescales_as_list(wd_entry)
                        to_write.writerow(output_row)

def main():
    timescale_interpolator = TimescaleInterpolator()
    HorHe = 'He'
    logg = 7.93
    Teff = 15620
    CaHe = -8.7
    print('Example timescales for a ' + HorHe + ' WD with log(g) = ' + str(logg) + ', Teff = ' + str(Teff) + ' and Ca/He = ' + str(CaHe))
    print(timescale_interpolator.get_wd_timescales(HorHe, logg, Teff, CaHe))
    #interpolated_timescales = timescale_interpolator.process_wd_data()
    #print(interpolated_timescales)
    #print('Call timescale_interpolator.dump_wd_timescales(interpolated_timescales) to dump these timescales into a file')
    
    print('Call TimescaleInterpolator.perform_format(detlev_file, output_file) to write files into a machine readable format')

if __name__ == '__main__':
    main()
