#!/usr/bin/env python

from pathlib import Path

import json
import numpy as np
import pymultinest as pn
# import sys

import loglike_functions as lf
import model_parameters as mp
import prior_functions as pf
import pwd_utils as pu

# sys.path.append(pu.get_path_to_pocomc_dir())
# import pocomc as pmc


class PollutionModel:

    def __init__(
        self,
        basename,
        timescale_type,
        enhancement_model_name,
        prior_name,
        consider_thermohaline=False,
        live_points=2000,
        seed=-1,
        verbose=True,
        resume=True,
    ):
        self.basename = basename
        try:
            model_definition = mp.model_definitions_dict[self.basename]
        except KeyError:
            raise KeyError(f"Unrecognised model name: {self.basename}")
        self.set_model_params(model_definition)
        self.prior = pf.universal_prior
        self.loglike = lf.universal_loglike
        self.verbose = verbose
        self.resume = resume
        self.live_points = live_points
        self.seed = seed
        self.result = None
        self.comparison = dict()
        self.best_model = None
        self.best_model_without_parameter = dict()
        self.best_heated_model = False
        self.enhancement_model_name = enhancement_model_name
        self.prior_name = prior_name
        self.consider_thermohaline = consider_thermohaline
        self.minimum_likelihood = mp.minimum_likelihood
        self.timescale_type = timescale_type
        self.actual_output_dir = None
        self.use_pymultinest = True
        if self.use_pymultinest:
            self.execute = self.execute_pymultinest
        else:
            self.execute = self.execute_pocomc

    def __repr__(self):
        return f"Model {self.basename}"

    def get_n_dims(self):
        return len(self.params)

    # def get_output_path(self, output_dir, observation_number):
    #     return output_dir + self.get_identifier1(observation_number) + self.basename

    def get_identifier1(self):
        if self.consider_thermohaline:
            identifier = "t_"
        else:
            identifier = "n_"
        identifier += f"p{self.live_points}_"
        return identifier

    def get_identifier2(self):
        try:
            identifier = (
                f"_{pu.abbreviations[self.enhancement_model_name]}_"
                f"{pu.abbreviations[self.prior_name]}_"
            )
        except KeyError:
            raise KeyError(
                "Could not identify abbreviation for at least one of "
                + f"{self.prior_name}, or {self.enhancement_model_name}"
                + " (may need to add it to abbreviations in pwd_utils.py)"
            )
        return identifier

    def get_prefix(self):
        return (
            self.timescale_type.short_str()
            + f"_{self.get_identifier1()}"
            + self.get_identifier1()
            + self.basename
            + self.get_identifier2()
        )

    def get_full_prefix(self, output_dir):
        return output_dir + self.get_prefix()

    def set_model_params(self, model_definition):
        to_set = list()
        for parameter in mp.get_model_params_in_order():
            parameter_present = model_definition.get(parameter, False)
            if parameter_present:
                to_set.append(mp.model_parameter_strings[parameter])
        self.params = to_set

    def get_model_params(self):
        return self.params

    def get_max_file_length(self):
        if self.use_pymultinest:
            return 78  # 100 minus 22 characters that PyMultiNest will add
        else:
            return np.inf

    def get_file_name_overflow(self, output_dir):
        return len(self.get_full_prefix(output_dir)) - self.get_max_file_length()

    def execute_pymultinest(self, output_dir):
        # What other parameters should (could) this have?
        # - Core/mantle/crust compositions
        # - Core/crust fractions of fragment/parent
        # - Distance of formation + width of feeding zone
        # - Accretion: timescale and time since
        # - Phase: Build-up, steady state or declining
        # - Ice
        # - Heating
        print(f"Executing {self}")
        print("pn args:")
        print(self.enhancement_model_name)
        print(self.prior_name)
        print(self.get_n_dims())
        print(self.get_full_prefix(output_dir))
        print(self.verbose)
        print(self.live_points)
        print(self.seed)
        print(self.minimum_likelihood)
        print(self.consider_thermohaline)

        disable_filename_error = False
        self.actual_output_dir = output_dir

        if (
            not disable_filename_error
            and len(self.get_full_prefix(output_dir)) > self.get_max_file_length()
        ):
            # TODO: In this case, we can now make use of the new functions in the
            # WhiteDwarf class to abbreviate the system_name and variant!
            raise IOError(
                f"\nOutput path {self.get_full_prefix(output_dir)} is too long!"
                + f"\nLength was {len(self.get_full_prefix(output_dir))}, but max"
                + f" length is {max_file_length} so that the full names fit in 100"
                + " characters"
                + "\nThe 100 character limit is hardcoded in MultiNest v3.10."
                + "\nIf you are using MultiNest 3.11 or later, you should be able to"
                + "just disable this error (go to execute_pymultinest in"
                + f"{str(Path(__file__).resolve())})\n(This is untested though!)"
            )

        # progress_plotter = pn.ProgressPlotter(
        #     n_params = self.get_n_dims(),
        #     outputfiles_basename = self.get_full_prefix(observation_number)
        # )
        # progress_plotter.start()

        self.result = pn.solve(
            LogLikelihood=self.loglike,
            Prior=self.prior,
            n_dims=self.get_n_dims(),
            outputfiles_basename=self.get_full_prefix(output_dir),
            verbose=self.verbose,
            n_live_points=self.live_points,
            seed=self.seed,
            resume=self.resume,
            # evidence_tolerance = 0.5,
            # sampling_efficiency = 0.8,
            # multimodal = False,
            log_zero=self.minimum_likelihood,
            # if likelihood < minimum_likelihood, point gets ignored
            # use this for errors)
        )
        # progress_plotter.stop()
        print(f"evidence: {self.result['logZ']:.1f} +- {self.result['logZerr']:.1f}")
        print(f"MultiNest Model {self.get_full_prefix(output_dir)} Completed")
        with open(f"{self.get_full_prefix(output_dir)}params.json", "w") as f:
            json.dump(self.params, f, indent=2)
        # self.results[str(stellar_composition)] = (
        #     'Result of executing '
        #     + str(self)
        #     + ' on composition '
        #     + str(stellar_composition)
        # )

    def execute_pocomc(self, output_dir):
        import pocomc as pmc
        # What other parameters should (could) this have?
        # - Core/mantle/crust compositions
        # - Core/crust fractions of fragment/parent
        # - Distance of formation + width of feeding zone
        # - Accretion: timescale and time since
        # - Phase: Build-up, steady state or declining
        # - Ice
        # - Heating
        print(f"Executing {self}")
        print("pocomc args:")
        print(self.enhancement_model_name)
        print(self.prior_name)
        print(self.get_n_dims())
        print(self.get_full_prefix(output_dir))
        print(self.verbose)
        print(self.live_points)
        print(self.seed)
        print(self.minimum_likelihood)
        print(self.consider_thermohaline)
        self.actual_output_dir = output_dir

        from scipy.stats import uniform

        hacky_prior = pmc.Prior(
            [
                uniform(loc=0.0, scale=958.0),  # metallicity
                uniform(loc=0.0, scale=10.0),  # t, in Myr
                uniform(loc=-1.0, scale=3.0),  # distance
                uniform(loc=0.0, scale=0.15),  # feeding zone
                # uniform(loc=0.0, scale=1.0), # fc
                uniform(loc=10.0, scale=15.0),  # mass
                uniform(loc=5.0, scale=3.0),  # t_event, log(yr)
                # uniform(loc=0.0, scale=60.0), # pressure
                # uniform(loc=-3.0, scale=2.0) # fO2
            ]
        )

        sampler = pmc.Sampler(
            prior=hacky_prior,
            likelihood=self.loglike,
            # vectorize=False
            # ^ Switching this to True could be good for speed-up purposes, but requires
            # some rewiring...
            random_state=0,
            # ^ This should be self.seed - pocomc doesn't accept the default value of -1
            # as a valid input unlike pymultinest though
            dynamic=True,
            # ^ Seems good? Dynamically adjusts particle count as necessary
            output_dir=output_dir,
            output_label=self.get_prefix(),
        )
        # To do checkpointing, need to call resume_state_path = "states/pmc_final.state"
        # or pmc_i.state where final should be preferred, otherwise higher values of i
        sampler.run(save_every=3)
        samples, weights, logl, logp = sampler.posterior()
        import matplotlib.pyplot as plt
        import corner

        fig = corner.corner(samples, weights=weights, color="C0")
        plt.show()
        print(logl)
        print(logp)
        logz, logz_err = sampler.evidence()
        print(f"logZ: {logz} +- {logz_err}")
        raise
