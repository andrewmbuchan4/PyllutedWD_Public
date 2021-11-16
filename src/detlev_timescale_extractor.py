#!/usr/bin/env python
# -*- coding: utf-8 -*-

import collections
import csv
import chemistry_info as ci

def load_detlev_data():
    toret = dict()
    g_values = [7.5, 8, 8.5]
    file_dict = {
        'H': {
            7.5: 'Detlev_timescales_g7.5_H_overshoot.csv',
            8: 'Detlev_timescales_g8_H_overshoot.csv',
            8.5: 'Detlev_timescales_g8.5_H_overshoot.csv'
        },
        'He': {
            7.5: 'Detlev_timescales_g7.5_He_overshoot.csv',
            8: 'Detlev_timescales_g8_He_overshoot.csv',
            8.5: 'Detlev_timescales_g8.5_He_overshoot.csv'
        }
    }
    for HorHe in ['H', 'He']:
        toret[HorHe] = dict()
        for g in g_values:
            toret[HorHe][g] = dict()
            T_vals = list()
            q_vals = list()
            with open('../data/' + file_dict[HorHe][g]) as csvfile:
                read = csv.reader(csvfile, delimiter=',')
                i = 0
                for row in read:
                    if i == 0:
                        T_vals = row[1:]
                        for j in range(1, len(row)):
                            toret[HorHe][g][int(T_vals[j-1])]= dict()
                    elif i == 1:
                        q_vals = row[1:]
                        for j in range(1, len(row)):
                            toret[HorHe][g][int(T_vals[j-1])]['logq'] = float(q_vals[j-1])
                    else:
                        element = ci.Element(int(row[0]))
                        for j in range(1, len(row)):
                            toret[HorHe][g][int(T_vals[j-1])][element] = float(row[j]) 
                    i += 1
    return toret

def perform_interpolation(data, x1, x2, y1, y2, x, y):
    toret = dict()
    for key in data[x1][y1].keys():  # This is just a way of looping through all the elements (and q)
        if x1 == x2:
            # Then just do 1D interpolation in y
            z1 = data[x1][y1][key]
            z2 = data[x1][y2][key]
            numerator = z1*(y2 - y) + z2*(y - y1)
            denominator = y2 - y1
            result = numerator / denominator
        elif y1 == y2:
            # Then just do 1D interpolation in x
            z1 = data[x1][y1][key]
            z2 = data[x2][y1][key]
            numerator = z1*(x2 - x) + z2*(x - x1)
            denominator = x2 - x1
            result = numerator / denominator
        else:
            # Bilinear interpolation in x and y, following notation of https://en.wikipedia.org/wiki/Bilinear_interpolation
            fQ11 = data[x1][y1][key]
            fQ12 = data[x1][y2][key]
            fQ21 = data[x2][y1][key]
            fQ22 = data[x2][y2][key]
            numerator = fQ11*(x2 - x)*(y2 - y) + fQ21*(x - x1)*(y2 - y) + fQ12*(x2 - x)*(y - y1) + fQ22*(x - x1)*(y - y1)
            denominator = (x2 - x1)*(y2 - y1)
            result = numerator / denominator
        toret[key] = result
    return toret

def extract_timescales(data, HorHe, Teff, logg):
    try:
        return data[HorHe][logg][Teff]
    except KeyError:
        pass
    # Need to interpolate!
    g_lower_bound = None
    current_best_lower_diff = None
    g_upper_bound = None
    current_best_upper_diff = None
    for g_test in data[HorHe].keys():
        if g_test <= logg:
            lower_diff = logg - g_test
            if (current_best_lower_diff is None) or (lower_diff < current_best_lower_diff):
                current_best_lower_diff = lower_diff
                g_lower_bound = g_test
        if g_test >= logg:
            upper_diff = g_test - logg
            if (current_best_upper_diff is None) or (upper_diff < current_best_upper_diff):
                current_best_upper_diff = upper_diff
                g_upper_bound = g_test
    # At this point, if either of the g bounds is None, then we are outside the range of the data, so just use the min/max value of g as appropriate
    if g_upper_bound is None:
        g_upper_bound = g_lower_bound
    if g_lower_bound is None:
        g_lower_bound = g_upper_bound
    
    # The next section is a copy/paste. TODO: make this into a function
    g_value = g_upper_bound if g_lower_bound is None else g_lower_bound  # This doesn't matter, it's just to read the T values which should be the same for all g
    t_lower_bound = None
    current_best_lower_diff = None
    t_upper_bound = None
    current_best_upper_diff = None
    for t_test in data[HorHe][g_value].keys():
        if t_test <= Teff:
            lower_diff = Teff - t_test
            if (current_best_lower_diff is None) or (lower_diff < current_best_lower_diff):
                current_best_lower_diff = lower_diff
                t_lower_bound = t_test
        if t_test >= Teff:
            upper_diff = t_test - Teff
            if (current_best_upper_diff is None) or (upper_diff < current_best_upper_diff):
                current_best_upper_diff = upper_diff
                t_upper_bound = t_test
    
    if t_upper_bound is None:
        t_upper_bound = t_lower_bound
    if t_lower_bound is None:
        t_lower_bound = t_upper_bound
    if (g_lower_bound == g_upper_bound) and (t_lower_bound == t_upper_bound):
        # This means that one of the parameters was out of bounds and the other was in the keys already (or both out of bounds?)
        # Can now just look up the appropriate value
        return data[HorHe][g_lower_bound][t_lower_bound]
    return perform_interpolation(data[HorHe], g_lower_bound, g_upper_bound, t_lower_bound, t_upper_bound, logg, Teff)

def load_wd_data():
    toret = collections.OrderedDict()
    with open('../data/BlouinConglomNewTimescales.csv') as csvfile:
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
                    logg = float(row[4])
                except ValueError:
                    # One of them has logg missing! Assume 8 by default
                    logg = 8.0
                    print('Warning: ' + wd_name + ' has no log(g) value. Assuming ' + str(logg))
                toret[wd_name] = {'Type': wd_type, 'Teff': Teff, 'logg': logg}
            i += 1
    return toret

def process_wd_data(wd_data, timescale_data):
    toret = collections.OrderedDict()
    for wd_name, wd_entry in wd_data.items():
        HorHe = wd_entry['Type']
        Teff = wd_entry['Teff']
        logg = wd_entry['logg']
        timescales = extract_timescales(timescale_data, HorHe, Teff, logg)
        toret[wd_name] = timescales
    return toret

def dump_wd_timescales(wd_timescales, outfile='wd_timescales'):
    if outfile is not None:
        with open(outfile + '.csv', 'w', newline='', encoding='utf-8') as f:
            to_write = csv.writer(f)
            to_write.writerow(['WD Name', 'log(q)', 't_Al', 't_Ti', 't_Ca', 't_Ni', 't_Fe', 't_Cr', 't_Mg', 't_Si', 't_Na', 't_O', 't_C', 't_N'])
            for wd_name, wd_entry in wd_timescales.items():
                to_write.writerow([
                    wd_name,
                    wd_entry['logq'],
                    10**wd_entry[ci.Element.Al],  # We only care about these ones
                    10**wd_entry[ci.Element.Ti],
                    10**wd_entry[ci.Element.Ca],
                    10**wd_entry[ci.Element.Ni],
                    10**wd_entry[ci.Element.Fe],
                    10**wd_entry[ci.Element.Cr],
                    10**wd_entry[ci.Element.Mg],
                    10**wd_entry[ci.Element.Si],
                    10**wd_entry[ci.Element.Na],
                    10**wd_entry[ci.Element.O],
                    10**wd_entry[ci.Element.C],
                    10**wd_entry[ci.Element.N]
                ])

def main():
    timescale_data = load_detlev_data()  # Structure: data[g][Teff] = {'logq': <q>, Element.Li: <t_Li>, ... Element.Zn: <t_Zn> }
    wd_data = load_wd_data()
    wd_timescales = process_wd_data(wd_data, timescale_data)
    dump_wd_timescales(wd_timescales)
    

if __name__ == '__main__':
    main()
