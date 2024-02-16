#!/usr/bin/env python
# -*- coding: utf-8 -*-
import csv

import chemistry_info as ci
import pwd_utils as pu

def get_40pc_abundance_values(element, fp_hx):
    desired_header = 'log' + str(element).lower() + str(fp_hx).lower()
    header_index = None
    toret = list()
    with open(pu.get_path_to_data() + 'MWDD_40pc_complete_clean.csv', encoding='utf-8') as fpcsv:
        header_next = True
        for row in csv.reader(fpcsv, delimiter='|'):
            if header_next: # Heading row
                header_index = row.index(desired_header)
                header_next = False
            else:
                raw_abundance = row[header_index]
                try:
                    abundance = float(raw_abundance)
                except ValueError:
                    # This excludes blanks, as well as any upper bounds
                    continue
                toret.append(abundance)
    return toret

def get_40pc_el_el_values_hx(element1, element2, fp_hx):
    desired_header1 = 'log' + str(element1).lower() + str(fp_hx).lower()
    desired_header2 = 'log' + str(element2).lower() + str(fp_hx).lower()
    header_index1 = None
    header_index2 = None
    toret = list()
    with open(pu.get_path_to_data() + 'MWDD_40pc_complete_clean.csv', encoding='utf-8') as fpcsv:
        header_next = True
        for row in csv.reader(fpcsv, delimiter='|'):
            if header_next: # Heading row
                header_index1 = row.index(desired_header1)
                header_index2 = row.index(desired_header2)
                header_next = False
            else:
                raw_abundance1 = row[header_index1]
                raw_abundance2 = row[header_index2]
                try:
                    ratio = float(raw_abundance1) - float(raw_abundance2)
                except ValueError:
                    # This excludes blanks, as well as any upper bounds
                    continue
                toret.append(ratio)
    return toret
    
def get_40pc_el_el_values(element1, element2):
    return get_40pc_el_el_values_hx(element1, element2, ci.Element.H) + get_40pc_el_el_values_hx(element1, element2, ci.Element.He)
    #return get_40pc_el_el_values_hx(element1, element2, ci.Element.H)
        
def main():
    print(get_40pc_el_el_values(ci.Element.Ca, ci.Element.Fe))

if __name__ == '__main__':
    main()
