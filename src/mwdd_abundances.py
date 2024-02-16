#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv

import chemistry_info as ci
import pwd_utils as pu

def get_mwdd_abundances(element, Hx, include_teff=False, include_upper_bounds=False):
    element_indices = {
        (ci.Element.H, ci.Element.He): 5,
        (ci.Element.C, ci.Element.He): 6,
        (ci.Element.C, ci.Element.H): 7,
        (ci.Element.N, ci.Element.H): 8,
        (ci.Element.O, ci.Element.He): 9,
        (ci.Element.O, ci.Element.H): 10,
        (ci.Element.Na, ci.Element.He): 11,
        (ci.Element.Na, ci.Element.H): 12,
        (ci.Element.Mg, ci.Element.He): 13,
        (ci.Element.Mg, ci.Element.H): 14,
        (ci.Element.Al, ci.Element.He): 15,
        (ci.Element.Al, ci.Element.H): 16,
        (ci.Element.Si, ci.Element.He): 17,
        (ci.Element.Si, ci.Element.H): 18,
        (ci.Element.Ca, ci.Element.He): 19,
        (ci.Element.Ca, ci.Element.H): 20,
        (ci.Element.Ti, ci.Element.He): 21,
        (ci.Element.Ti, ci.Element.H): 22,
        (ci.Element.Cr, ci.Element.He): 23,
        (ci.Element.Cr, ci.Element.H): 24,
        (ci.Element.Fe, ci.Element.He): 25,
        (ci.Element.Fe, ci.Element.H): 26,
        (ci.Element.Ni, ci.Element.He): 27,
        (ci.Element.Ni, ci.Element.H): 28
    }
    teff_index = 4
    try:
        element_index = element_indices[(element, Hx)]
    except KeyError:
        if include_upper_bounds:
            return None, None
        else:
            return None
    toret = list()
    toret_ub = list()
    with open(pu.get_path_to_data() + 'MWDD-export-allpollution.csv', encoding='utf-8') as mwddcsv:
        row_count = 0
        for row in csv.reader(mwddcsv):
            if row_count == 0: # Heading row
                pass
            else:
                teff = None
                upper_bound = None
                abundance = None
                try:
                    teff = float(row[teff_index])
                except ValueError:
                    pass
                if row[element_index].startswith('<'):
                    try:
                        upper_bound = float(row[element_index][1:])
                    except ValueError:
                        pass
                else:
                    try:
                        abundance = float(row[element_index])
                    except ValueError:
                        pass
                if include_teff:
                    if teff is not None:
                        if abundance is not None:
                            toret.append((teff, abundance))
                        if upper_bound is not None:
                            toret_ub.append((teff, upper_bound))
                else:
                    if abundance is not None:
                        toret.append(abundance)
                    if upper_bound is not None:
                        toret_ub.append(upper_bound)
            row_count += 1
    if include_upper_bounds:
        return toret, toret_ub
    else:
        return toret

def main():
    el1 = ci.Element.Mg
    el2 = ci.Element.H
    values = get_mwdd_abundances(el1, el2)
    print(str(len(values)) + ' values for ' + str(el1) + '/' + str(el2) + ' found in MWDD:')
    print(values)

if __name__ == '__main__':
    main()
