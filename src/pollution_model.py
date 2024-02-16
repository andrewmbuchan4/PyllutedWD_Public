#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path

import json

import loglike_functions as lf
import model_parameters as mp
import prior_functions as pf
import pwd_utils as pu
import pymultinest as pn

class PollutionModel:
    
    def __init__(self, basename, enhancement_model_name, prior_name, live_points=2000, seed=-1, verbose=True):
        self.basename = basename
        print(self.basename)
        print(mp.model_definitions_dict.keys())
        try:
            model_definition = mp.model_definitions_dict[self.basename]
        except KeyError:
            raise KeyError('Unrecognised model name: ' + self.basename)
        self.set_model_params(model_definition)
        self.prior = pf.universal_prior
        self.loglike = lf.universal_loglike
        self.verbose = verbose
        self.live_points = live_points
        self.seed = seed
        self.result = None
        self.comparison = dict()
        self.best_model = None
        self.best_model_without_parameter = dict()
        self.best_heated_model = False
        self.enhancement_model_name = enhancement_model_name
        self.prior_name = prior_name
        self.minimum_likelihood = mp.minimum_likelihood

    def __repr__(self):
        return 'Model ' + self.basename
    
    def get_n_dims(self):
        return len(self.params)
        
    def get_output_path(self, output_dir, observation_number):
        return output_dir + self.get_identifier1(observation_number) + self.basename
    
    def get_identifier1(self, observation_number):
        identifier =  'p' + str(self.live_points) + '_'
        return identifier
    
    def get_identifier2(self):
        try:
            identifier =  '_' + pu.abbreviations[self.enhancement_model_name] + '_' + pu.abbreviations[self.prior_name] + '_'
        except KeyError:
            raise KeyError('Could not identify abbreviation for at least one of ' + self.prior_name + ', or ' + self.enhancement_model_name + ' (may need to add it to abbreviations in pwd_utils.py)')
        return identifier
    
    def get_prefix(self, observation_number):
        return self.get_identifier1(observation_number) + self.basename + self.get_identifier2()
    
    def get_full_prefix(self, output_dir, observation_number):
        return self.get_output_path(output_dir, observation_number) + self.get_identifier2()
    
    def set_model_params(self, model_definition):
        #TODO: Certain params require other params to be present (see universal_prior) --> raise error if not the case
        to_set = list()
        for parameter in mp.get_model_params_in_order():
            parameter_present = model_definition.get(parameter, False)
            if parameter_present:
                to_set.append(mp.model_parameter_strings[parameter])
        self.params = to_set
    
    def get_model_params(self):
        return self.params
    
    def execute(self, observation_number, output_dir):
        # What other parameters should (could) this have?
        # - Core/mantle/crust compositions
        # - Core/crust fractions of fragment/parent
        # - Distance of formation + width of feeding zone
        # - Accretion: timescale and time since
        # - Phase: Build-up, steady state or declining
        # - Ice
        # - Heating
        print('Executing ' + str(self))
        print('pn args:')
        print(self.enhancement_model_name)
        print(self.prior_name)
        print(self.get_n_dims())
        print(self.get_full_prefix(output_dir, observation_number))
        print(self.verbose)
        print(self.live_points)
        print(self.seed)
        print(self.minimum_likelihood)
        
        disable_filename_error = False
        max_file_length = 78 # 100 minus 22 characters that PyMultiNest will add
        
        if not disable_filename_error and len(self.get_full_prefix(output_dir, observation_number)) > max_file_length:
            raise IOError('\nOutput path ' + self.get_full_prefix(output_dir, observation_number) + ' is too long!\nLength was ' + str(len(self.get_full_prefix(output_dir, observation_number))) + ', but max length is ' + str(max_file_length) + ' so that the full names fit in 100 characters\nThe 100 character limit is hardcoded in MultiNest v3.10.\nIf you are using MultiNest 3.11 or later, you should be able to just disable this error (line 97 in ' + str(Path(__file__).resolve()) + ')\n(This is untested though!)')
        
        #progress_plotter = pn.ProgressPlotter(n_params = self.get_n_dims(), outputfiles_basename = self.get_full_prefix(observation_number))
        #progress_plotter.start()
        
        self.result = pn.solve(
            LogLikelihood=self.loglike,
            Prior=self.prior,
            n_dims=self.get_n_dims(),
            outputfiles_basename = self.get_full_prefix(output_dir, observation_number),
            verbose=self.verbose,
            n_live_points=self.live_points,
            seed=self.seed,
            #evidence_tolerance = 0.5,
            #sampling_efficiency = 0.8,
            #multimodal = False,
            log_zero = self.minimum_likelihood  # if likelihood < minimum_likelihood, point gets ignored (use this for errors)
        )
        #progress_plotter.stop()
        print('evidence: %(logZ).1f +- %(logZerr).1f' % self.result)
        print("MultiNest Model " + self.get_full_prefix(output_dir, observation_number) + " Completed")
        with open('%sparams.json' % self.get_full_prefix(output_dir, observation_number), 'w') as f:
            json.dump(self.params, f, indent=2)
        #self.results[str(stellar_composition)] = 'Result of executing ' + str(self) + ' on composition ' + str(stellar_composition)
