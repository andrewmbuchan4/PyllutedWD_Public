#!/usr/bin/env python
# -*- coding: utf-8 -*-

import chemistry_info as ci
import geology_info as gi
import graph_factory as gf
import model_parameters as mp
import timescale_interpolator as ti

import csv
import numpy as np

wd_dict = {
    'SDSSJ1228+1040Corr': { #UV, ID 365
        'Properties': {
            mp.WDParameter.spectral_type: 'H',
            mp.WDParameter.temperature: 20900,
            mp.WDParameter.logg: 8.15
        },
        'Abundances': {
            ci.Element.Al: -5.75,
            #ci.Element.Ti: 0,
            ci.Element.Ca: -5.94,
            ci.Element.Ni: -6.5,
            ci.Element.Fe: -5.2,
            ci.Element.Cr: -6.00,
            ci.Element.Mg: -5.2,
            ci.Element.Si: -5.2,
            #ci.Element.Na: 0,
            ci.Element.O: -4.55,
            ci.Element.C: -7.5,
            #ci.Element.N: 0
        },
        'Errors': {
            ci.Element.Al: 0.2,
            #ci.Element.Ti: 0,
            ci.Element.Ca: 0.2,
            #ci.Element.Ni: 0,
            ci.Element.Fe: 0.3,
            #ci.Element.Cr: 0,
            ci.Element.Mg: 0.2,
            ci.Element.Si: 0.2,
            #ci.Element.Na: 0,
            ci.Element.O: 0.2,
            ci.Element.C: 0.2,
            #ci.Element.N: 0
        },
        'UpperBounds': {
            ci.Element.Ni: True,
            ci.Element.Cr: True
        }
    },
    'SDSSJ1228+1040OptCorr': { #Optical, ID 364
        'Properties': {
            mp.WDParameter.spectral_type: 'H',
            mp.WDParameter.temperature: 20900,
            mp.WDParameter.logg: 8.15
        },
        'Abundances': {
            ci.Element.Al: -5.75,
            #ci.Element.Ti: 0,
            ci.Element.Ca: -5.94,
            ci.Element.Ni: -6.5,
            ci.Element.Fe: -5.2,
            ci.Element.Cr: -6.00,
            ci.Element.Mg: -5.2,
            ci.Element.Si: -4.7,
            #ci.Element.Na: 0,
            ci.Element.O: -4.55,
            ci.Element.C: -7.5,
            #ci.Element.N: 0
        },
        'Errors': {
            ci.Element.Al: 0.2,
            #ci.Element.Ti: 0,
            ci.Element.Ca: 0.2,
            #ci.Element.Ni: 0,
            ci.Element.Fe: 0.3,
            #ci.Element.Cr: 0,
            ci.Element.Mg: 0.2,
            ci.Element.Si: 0.2,
            #ci.Element.Na: 0,
            ci.Element.O: 0.2,
            ci.Element.C: 0.2,
            #ci.Element.N: 0
        },
        'UpperBounds': {
            ci.Element.Ni: True,
            ci.Element.Cr: True
        }
    },
    'PG1225-079Corr': { # ID 380
        'Properties': {
            mp.WDParameter.spectral_type: 'He',
            mp.WDParameter.temperature: 10800,
            mp.WDParameter.logg: 8
        },
        'Abundances': {
            ci.Element.Al: -7.84,
            ci.Element.Ti: -9.45,
            ci.Element.Ca: -8.06,
            ci.Element.Ni: -8.76,
            ci.Element.Fe: -7.42,
            ci.Element.Cr: -9.27,
            ci.Element.Mg: -7.5,
            ci.Element.Si: -7.45,
            ci.Element.Na: -8.26,
            ci.Element.O: -5.54,
            ci.Element.C: -7.8,
            #ci.Element.N: 0
        },
        'Errors': {
            #ci.Element.Al: 0,
            ci.Element.Ti: 0.231,
            ci.Element.Ca: 0.192,
            ci.Element.Ni: 0.184,
            ci.Element.Fe: 0.231,
            ci.Element.Cr: 0.228,
            ci.Element.Mg: 0.2,
            ci.Element.Si: 0.1,
            #ci.Element.Na: 0,
            #ci.Element.O: 0,
            ci.Element.C: 0.1,
            #ci.Element.N: 0
        },
        'UpperBounds': {
            ci.Element.Al: True,
            ci.Element.Na: True,
            ci.Element.O: True
        }
    },
    'WD1145+017Corr': { #Fortin-Archambault, ID 375
        'Properties': {
            mp.WDParameter.spectral_type: 'He',
            mp.WDParameter.temperature: 14500,
            mp.WDParameter.logg: 8.11
        },
        'Abundances': {
            ci.Element.Al: -6.89,
            ci.Element.Ti: -8.57,
            ci.Element.Ca: -7,
            ci.Element.Ni: -7.02,
            ci.Element.Fe: -5.61,
            ci.Element.Cr: -7.92,
            ci.Element.Mg: -5.91,
            ci.Element.Si: -5.89,
            #ci.Element.Na: 0,
            ci.Element.O: -5.12,
            ci.Element.C: -7.5,
            ci.Element.N: -7.00
        },
        'Errors': {
            ci.Element.Al: 0.2,
            ci.Element.Ti: 0.2,
            ci.Element.Ca: 0.2,
            ci.Element.Ni: 0.3,
            ci.Element.Fe: 0.2,
            ci.Element.Cr: 0.4,
            ci.Element.Mg: 0.2,
            ci.Element.Si: 0.2,
            #ci.Element.Na: 0,
            ci.Element.O: 0.35,
            ci.Element.C: 0.4,
            #ci.Element.N: 0
        },
        'UpperBounds': {
            ci.Element.N: True
        }
    },
    'WD1145+017Budaj': { #Budaj, ID 472
        'Properties': {
            mp.WDParameter.spectral_type: 'He',
            mp.WDParameter.temperature: 15000,
            mp.WDParameter.logg: 8
        },
        'Abundances': {
            ci.Element.Al: -6.89,
            ci.Element.Ti: -8.66,
            ci.Element.Ca: -7.32,
            ci.Element.Ni: -7.62,
            ci.Element.Fe: -5.92,
            ci.Element.Cr: -8.1,
            ci.Element.Mg: -6.05,
            ci.Element.Si: -5.94,
            #ci.Element.Na: 0,
            ci.Element.O: -5.19,
            ci.Element.C: -8.03,
            ci.Element.N: -8.49
        },
        'Errors': { # No errors were given
            ci.Element.Al: 0.2,
            ci.Element.Ti: 0.2,
            ci.Element.Ca: 0.2,
            ci.Element.Ni: 0.2,
            ci.Element.Fe: 0.2,
            ci.Element.Cr: 0.2,
            ci.Element.Mg: 0.2,
            ci.Element.Si: 0.2,
            #ci.Element.Na: 0,
            ci.Element.O: 0.2,
            ci.Element.C: 0.2,
            #ci.Element.N: 0
        },
        'UpperBounds': {
            ci.Element.N: True
        }
    },
    'PG1015+161Xu': { #Xu, ID 422
        'Properties': {
            mp.WDParameter.spectral_type: 'H',
            mp.WDParameter.temperature: 19226,
            mp.WDParameter.logg: 8.04
        },
        'Abundances': {
            #ci.Element.Al: 0,
            #ci.Element.Ti: 0,
            ci.Element.Ca: -6.4,
            #ci.Element.Ni: 0,
            ci.Element.Fe: -4.92,
            #ci.Element.Cr: 0,
            ci.Element.Mg: -5.6,
            ci.Element.Si: -5.42,
            #ci.Element.Na: 0,
            #ci.Element.O: 0,
            #ci.Element.C: 0,
            #ci.Element.N: 0
        },
        'Errors': {
            #ci.Element.Al: 0,
            #ci.Element.Ti: 0,
            ci.Element.Ca: 0.2,
            #ci.Element.Ni: 0,
            ci.Element.Fe: 0.2,
            #ci.Element.Cr: 0,
            ci.Element.Mg: 0.2,
            ci.Element.Si: 0.21,
            #ci.Element.Na: 0,
            #ci.Element.O: 0,
            #ci.Element.C: 0,
            #ci.Element.N: 0
        },
        'UpperBounds': dict()
    },
    'PG1015+161Corr': { #Gaensicke, ID 356
        'Properties': {
            mp.WDParameter.spectral_type: 'H',
            mp.WDParameter.temperature: 19226,
            mp.WDParameter.logg: 8.04
        },
        'Abundances': {
            #ci.Element.Al: 0,
            #ci.Element.Ti: 0,
            ci.Element.Ca: -6.45,
            #ci.Element.Ni: 0,
            ci.Element.Fe: -5.5,
            ci.Element.Cr: -5.8,
            ci.Element.Mg: -5.3,
            ci.Element.Si: -6.4,
            #ci.Element.Na: 0,
            ci.Element.O: -5.5,
            ci.Element.C: -8,
            #ci.Element.N: 0
        },
        'Errors': {
            #ci.Element.Al: 0,
            #ci.Element.Ti: 0,
            ci.Element.Ca: 0.2,
            #ci.Element.Ni: 0,
            ci.Element.Fe: 0.3,
            #ci.Element.Cr: 0,
            ci.Element.Mg: 0.2,
            ci.Element.Si: 0.2,
            #ci.Element.Na: 0,
            ci.Element.O: 0.2,
            #ci.Element.C: 0,
            #ci.Element.N: 0
        },
        'UpperBounds': {
            ci.Element.Cr: True,
            ci.Element.C: True
        }
    }
}

# Use forward = True if you know the intrinsic (ie planetesimal) composition and want to find out
# what you would observe if you put that material into a WD (in steady state)
# Use forward = False if you know the WD observations and want to work backwards to get
# the intrinsic planetesimal composition, assuming the WD accretion is in steady state
def find_steady_state_adjusted_abundances(planetesimal_abundance_dict, HorHe, logg, Teff, CaHe=None, forward=True):
    timescale_interpolator = ti.TimescaleInterpolator()
    all_timescales = False
    for el in planetesimal_abundance_dict:
        if el not in ci.usual_elements:
            all_timescales = True
    sinking_timescales = timescale_interpolator.get_wd_timescales(HorHe, logg, Teff, CaHe, all_timescales)
    total = 0
    toret = dict()
    print(sinking_timescales)
    for el, sub_dict in planetesimal_abundance_dict.items():
        unadjusted = sub_dict[gi.Layer.bulk]
        if forward:
            adjusted = unadjusted*sinking_timescales[el]
        else:
            adjusted = unadjusted/sinking_timescales[el]
        toret[el] = dict()
        toret[el][gi.Layer.bulk] = adjusted
        total += adjusted
    for el in toret:
        toret[el][gi.Layer.bulk] /= total
    return toret

def read_wd_csv_into_dict(wd_data_filename):
    output_dict = dict()
    input_abundance_indices = {
        ci.Element.Al: 9,
        ci.Element.Ti: 11,
        ci.Element.Ca: 13,
        ci.Element.Ni: 15,
        ci.Element.Fe: 17,
        ci.Element.Cr: 19,
        ci.Element.Mg: 21,
        ci.Element.Si: 23,
        ci.Element.Na: 25,
        ci.Element.O: 27,
        ci.Element.C: 29,
        ci.Element.N: 31
    }
    with open(wd_data_filename, encoding='utf-8') as wdcsv:
        row_count = 0
        for row in csv.reader(wdcsv):
            print()
            print(row)
            if row_count == 0: # Heading row
                pass  # TODO: Find column indices dynamically using the heading row
            else:
                wd_name = row[0]
                output_dict[wd_name] = {
                    'Properties': dict(),
                    'Abundances': dict(),
                    'Errors': dict(),
                    'UpperBounds': dict()
                }
                output_dict[wd_name]['Properties'][mp.WDParameter.spectral_type] = row[1]
                output_dict[wd_name]['Properties'][mp.WDParameter.temperature] = int(row[4])
                output_dict[wd_name]['Properties'][mp.WDParameter.logg] = float(row[6])
                for el in ci.usual_elements:
                    input_abundance_index = input_abundance_indices[el]
                    input_error_index = input_abundance_index + 1
                    print(el)
                    print(row[input_abundance_index])
                    try:
                        abundance = float(row[input_abundance_index])
                        output_dict[wd_name]['Abundances'][el] = abundance
                        output_dict[wd_name]['Errors'][el] = float(row[input_error_index])
                    except ValueError:
                        # Occurs if this is an upper bound so the value starts with <
                        if row[input_abundance_index].startswith('<'):
                            # Then this is an upper bound
                            output_dict[wd_name]['Abundances'][el] = float(row[input_abundance_index][1:])
                            output_dict[wd_name]['UpperBounds'][el] = True
                        elif row[input_abundance_index] == '':
                            pass
                        else:
                            # Then we don't know what it is
                            raise ValueError('Unable to parse abundances for ' + wd_name)
            row_count += 1
    return output_dict

def batch_convert(wd_dict):
    geo_model = gi.GeologyModel()
    default_Ca = -9
    forward = False
    all_wd_ss_number_abundances_dict = dict()
    all_errors_dict = dict()
    ss_upper_limits_dict = dict()
    all_wd_ss_wt_linear_abundances_dict = dict()
    for wd_name, wd_info in wd_dict.items():
        ss_wd_abundances_dict = dict()
        all_errors_dict[wd_name] = wd_info['Errors']
        wd_abundances = wd_info['Abundances']
        HorHe = wd_info['Properties'][mp.WDParameter.spectral_type]
        Teff = wd_info['Properties'][mp.WDParameter.temperature]
        logg = wd_info['Properties'][mp.WDParameter.logg]
        pol_frac = np.log10(sum([10**a for a in wd_abundances.values()]))
        planetesimal_abundance_dict = geo_model.convert_log_abundances_to_linear(wd_abundances)
        ca_abundance = wd_abundances.get(ci.Element.Ca, default_Ca)
        adjusted_planetesimal_dict = find_steady_state_adjusted_abundances(planetesimal_abundance_dict, HorHe, logg, Teff, ca_abundance, forward)
        ss_wd_abundances = geo_model.convert_linear_abundances_to_log(adjusted_planetesimal_dict, pol_frac)
        all_wd_ss_number_abundances_dict[wd_name] = dict()
        ss_upper_limits_dict[wd_name] = dict()
        for el, abundance in ss_wd_abundances.items():
            if wd_info['UpperBounds'].get(el, False):
                # Then this was actually an upper bound
                ss_upper_limits_dict[wd_name][el] = abundance
            else:
                all_wd_ss_number_abundances_dict[wd_name][el] = abundance
        present_ss_abundances = {el: adjusted_planetesimal_dict[el] for el in all_wd_ss_number_abundances_dict[wd_name]}
        all_wd_ss_wt_linear_abundances_dict[wd_name] = geo_model.convert_number_abundances_to_mass(present_ss_abundances)
    systems_to_print = [
        'G29-38Corr',
        'GaiaJ0644-0352',
        'GaiaJ2100+2122SpecMarch',
        'GALEX1931+0117VCorr',
        'GALEX1931+0117MCorr',
        'GD61Corr',
        'HE0106-3253',
        'PG0843+516XCorr',
        'PG0843+516GCorr',
        'PG1015+161Xu',
        'PG1015+161Corr',
        'PG1225-079Corr',
        'SDSSJ0738+1835Corr',
        'SDSSJ0845+2257Corr',
        'SDSSJ1043+0855Corr',
        'SDSSJ1228+1040OptCorr',
        'SDSSJ1228+1040Corr',
        'WD1145+017Corr',
        'WD1232+563Corr',
        'WD1536+520Corr',
        'WD1551+175Corr',
        'WD2115-560Corr',
        'WD2207+121Corr'
    ]
    for system in systems_to_print:
        print()
        print(system)
        print(all_wd_ss_number_abundances_dict[system])
        print(all_wd_ss_wt_linear_abundances_dict[system])
        print(all_errors_dict[system])
        print(ss_upper_limits_dict[system])
        print(100*all_wd_ss_wt_linear_abundances_dict[system][ci.Element.Fe][gi.Layer.bulk])
        print(100*all_wd_ss_wt_linear_abundances_dict[system][ci.Element.Mg][gi.Layer.bulk])
        print(100*all_wd_ss_wt_linear_abundances_dict[system][ci.Element.Si][gi.Layer.bulk])
    make_plot = False
    if make_plot:
        graph_fac = gf.GraphFactory()
        graph_fac.make_composition_plot_multi_system(
            'testmulti',
            all_wd_abundances_dict,
            all_errors_dict,
            upper_limits_dict
        )

def example():
    wd_abundances = {
        ci.Element.Al: -7.216657942,
        ci.Element.Ti: -8.880391333,
        ci.Element.Ca: -7.284162087,
        ci.Element.Ni: -7.474560357,
        ci.Element.Fe: -6.192683939,
        ci.Element.Cr: -8.050817538,
        ci.Element.Mg: -6,
        ci.Element.Si: -6.025018216,
        ci.Element.Na: -7.379933142,
        ci.Element.O: -5.484776702,
        ci.Element.C: -100,
        ci.Element.N: -100
    }
    HorHe = 'H'
    Teff = 20000
    logg = 8

    override_with_synthetic = False
    override_with_real = False
    override_with_ideal = False
    planet = 'LPC'
    fcf = 0.99
    error = 0.05
    geo_model = gi.GeologyModel()
    planetesimal_abundance_dict = geo_model.convert_log_abundances_to_linear(wd_abundances)
    if (not override_with_synthetic) and (not override_with_real) and (not override_with_ideal):
        planetesimal_abundance_dict = geo_model.convert_log_abundances_to_linear(wd_abundances)
    else:
        if override_with_synthetic:
            if planet == 'Earth':
                raw_abundances = geo_model.form_a_planet_iteratively(geo_model.get_earth_differentiation_pressure(), geo_model.get_earth_oxygen_fugacity())[0]
            if planet == 'Mars':
                raw_abundances = geo_model.form_a_planet_iteratively(geo_model.get_mars_differentiation_pressure(), geo_model.get_mars_oxygen_fugacity())[0]
        if override_with_real:
            if planet == 'Earth':
                raw_abundances = geo_model.element_info
            if planet == 'Mars':
                raw_abundances = geo_model.mars_abundances
        if override_with_ideal:
            if planet == 'HPM':
                raw_abundances = geo_model.form_a_planet_iteratively(60, -3)[0]
            if planet == 'HPC':
                raw_abundances = geo_model.form_a_planet_iteratively(60, -1)[0]
            if planet == 'LPM':
                raw_abundances = geo_model.form_a_planet_iteratively(0, -3)[0]
            if planet == 'LPC':
                raw_abundances = geo_model.form_a_planet_iteratively(0, -1)[0]
        mixed_dict = geo_model.get_earth_mix(fcf, raw_abundances)
        for ue in ci.usual_elements:
            planetesimal_abundance_dict[ue] = dict()
            planetesimal_abundance_dict[ue][gi.Layer.bulk] = mixed_dict[ue]
    print()
    print('Initial abundances')
    print(planetesimal_abundance_dict)
    adjusted_planetesimal_dict = find_steady_state_adjusted_abundances(planetesimal_abundance_dict, HorHe, logg, Teff, wd_abundances[ci.Element.Ca])
    print()
    print('Steady state adjusted abundances (by number):')
    print(adjusted_planetesimal_dict)
    # Now also output the abundances by weight, not number
    wt_adjusted_planetesimal_dict = geo_model.convert_number_abundances_to_mass(adjusted_planetesimal_dict)
    print()
    print('Steady state adjusted abundances (by mass):')
    print(wt_adjusted_planetesimal_dict)

    #Finally, the main thing, steady state adjusted abundances:
    pol_frac = np.log10(sum([10**a for a in wd_abundances.values()]))
    #print(pol_frac)
    ss_wd_abundances = geo_model.convert_linear_abundances_to_log(adjusted_planetesimal_dict, pol_frac)
    print()
    print('Steady state adjusted abundances (log scale, by number):')
    print(ss_wd_abundances)
    print()
    print('Steady state adjusted abundances (log scale, by number, comma-separated with errors):')
    to_print = ','.join(str(ss_wd_abundances[ue]) + ',' + str(error) for ue in ci.usual_elements)
    print(to_print)

def main():
    example()
    wd_dict = read_wd_csv_into_dict('../data/WDInputData.csv')
    batch_convert(wd_dict)

if __name__ == '__main__':
    main()
