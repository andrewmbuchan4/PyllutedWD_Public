#!/usr/bin/env python

from scipy import interpolate

import csv
import numpy as np
import os
import pwd_utils as pu


class ThermohalineInterpolator:
    # Partially based on Keith Williams' solution:
    # https://stackoverflow.com/questions/20516762/extrapolate-with-linearndinterpolator
    # But with a different solution for extrapolation since the nearest neighbour
    # solution had weird behaviour

    def __init__(self):
        self.accretion_rate_limits_dict = dict()
        self.temperature_limits_dict = dict()
        self.setup()

    def setup(self):
        # Setting up a 3D interpolator that finds logMCa_therm, logMCa_diff, and
        # therm_factor as a function of logMdot, Teff and logg

        base_path = pu.get_path_to_da_pollution_tables_dir()

        input_file_name = self.compile_MESA_outputs(base_path)

        all_logmdots = list()
        # all_logXCa_therms = list()
        # all_logXCa_diffs = list()
        all_therm_factors = list()
        all_teffs = list()
        all_logg = list()

        with open(input_file_name, "r") as input_file:
            read = csv.reader(input_file, delimiter=",")
            next(read)  # skip header
            for row in read:
                all_logmdots.append(float(row[0]))
                all_teffs.append(float(row[1]))
                all_logg.append(float(row[2]))
                # all_logXCa_therms.append(float(row[3]))
                # all_logXCa_diffs.append(float(row[4]))
                all_therm_factors.append(float(row[5]))

        points = list(zip(all_logmdots, all_teffs, all_logg))

        self.interpolate = interpolate.LinearNDInterpolator(points, all_therm_factors)

    def compile_MESA_outputs(self, base_path):
        dirname = base_path + "thermohaline_surfaceX/"
        subdirnames = list(os.walk(dirname))[0][1]
        all_logmdots = np.array([])
        all_logXCa_therms = np.array([])
        all_logXCa_diffs = np.array([])
        all_teffs = np.array([])
        all_logg = np.array([])
        all_therm_factors = np.array([])
        # ^ We will use this array to save the thermohaline correction factor, on a
        # log scale. Generally negative - thermohaline mixing decreases surface mass
        for subdirname in subdirnames:
            mass = float(subdirname.split("_")[2])
            subdir = dirname + subdirname
            Teff_file_names = list(os.walk(subdir))[0][2]
            for Teff_file_name in Teff_file_names:
                Teff_file = subdir + "/" + Teff_file_name
                Teff = float(Teff_file.split("K.data")[0].split("/")[-1])
                d = np.genfromtxt(Teff_file)
                logmdot = np.log10(d[:, 0])
                logXCa_therm = np.log10(d[:, 1])
                # ^ This is the Ca surface abundance with thermohaline mixing
                logXCa_diff = np.log10(d[:, 2])
                # ^ This is the Ca surface abundance without thermohaline mixing, just
                # diffusion
                teff = d[:, 3]
                logg = d[:, 4]
                therm_factors = logXCa_therm - logXCa_diff
                if not np.isnan(therm_factors).any():
                    all_logmdots = np.append(all_logmdots, logmdot)
                    all_logXCa_therms = np.append(all_logXCa_therms, logXCa_therm)
                    all_logXCa_diffs = np.append(all_logXCa_diffs, logXCa_diff)
                    all_teffs = np.append(all_teffs, teff)
                    all_logg = np.append(all_logg, logg)
                    all_therm_factors = np.append(all_therm_factors, therm_factors)
                else:
                    print(subdirname)
                    print(Teff_file_name)
                    print(
                        "This logic is hopefully irrelevant now - there used to be some invalid values in the raw tables"
                    )
        output_file_name = dirname + "all_mesa_outputs.csv"
        print("Writing to " + output_file_name)
        with open(output_file_name, "w", newline="", encoding="utf-8") as output_file:
            to_write = csv.writer(output_file)
            to_write.writerow(
                [
                    "log(Mdot)",
                    "Teff",
                    "log(g)",
                    "log(XCa)_therm",
                    "log(XCa)_diff",
                    "Therm_factor",
                ]
            )
            i = 0
            while i < len(all_logmdots):
                if np.isnan(all_therm_factors[i]):
                    print(all_logmdots[i])
                    print(all_teffs[i])
                    print(all_logg[i])
                    print(all_logXCa_therms[i])
                    print(all_logXCa_diffs[i])
                    print(all_therm_factors[i])
                    raise
                to_write.writerow(
                    [
                        all_logmdots[i],
                        all_teffs[i],
                        all_logg[i],
                        all_logXCa_therms[i],
                        all_logXCa_diffs[i],
                        all_therm_factors[i],
                    ]
                )
                i += 1
        return output_file_name

    def find_accretion_rate_limits(self, Teff, logg):
        if (Teff, logg) in self.accretion_rate_limits_dict:
            pass
        else:
            # At a given Teff/logg, there will be a maximum and a mininum value of
            # logMdot for which interpolation (rather than extrapolation) is possible -
            # we want to find these limits (and cache them)
            # Then if we come across a logMdot outside these limits, snap to the closest
            # limit.
            # Returns a tuple: (minimum accretion rate in log Mdot, maximum accretion
            # rate in log Mdot)
            logmdot_bound_high = 20  # Assume there's no chance we'll ever be above this
            logmdot_bound_low = 0  # Assume there's no chance we'll ever be below this

            iterations_to_try = 100
            # ^ The parameter space will halve each iteration, so this effectively sets
            # the tolerance level. 100 seems like ample

            current_upper_maximum_logmdot_limit = logmdot_bound_high
            current_lower_maximum_logmdot_limit = logmdot_bound_low

            current_upper_minimum_logmdot_limit = logmdot_bound_high
            current_lower_minimum_logmdot_limit = logmdot_bound_low

            iteration_count = 0

            valid_max_limit = False
            valid_min_limit = False

            while iteration_count < iterations_to_try:
                trial_upper_logmdot = (
                    current_upper_maximum_logmdot_limit
                    + current_lower_maximum_logmdot_limit
                ) / 2
                trial_lower_logmdot = (
                    current_upper_minimum_logmdot_limit
                    + current_lower_minimum_logmdot_limit
                ) / 2
                interpolated_upper_result = self.interpolate(
                    (trial_upper_logmdot, Teff, logg)
                )
                interpolated_lower_result = self.interpolate(
                    (trial_lower_logmdot, Teff, logg)
                )
                if np.isnan(interpolated_upper_result):
                    # Then we were out of bounds, and the trial value becomes the new upper
                    # limit:
                    current_upper_maximum_logmdot_limit = trial_upper_logmdot
                else:
                    valid_max_limit = True
                    current_lower_maximum_logmdot_limit = trial_upper_logmdot
                if np.isnan(interpolated_lower_result):
                    # Then we were out of bounds, and the trial value becomes the new
                    # lower limit:
                    current_lower_minimum_logmdot_limit = trial_lower_logmdot
                else:
                    valid_min_limit = True
                    current_upper_minimum_logmdot_limit = trial_lower_logmdot
                iteration_count += 1
            if valid_min_limit and valid_max_limit:
                self.accretion_rate_limits_dict[(Teff, logg)] = (
                    current_upper_minimum_logmdot_limit,
                    current_lower_maximum_logmdot_limit,
                )
                # ^ Conservatively these as the min/max values (to make sure we can only
                # ever interpolate within bounds). There's some inefficiency here (i.e.
                # some info from final iteration can potentially be unused) but I'm not
                # too concerned
            else:
                self.accretion_rate_limits_dict[(Teff, logg)] = np.nan, np.nan
        return self.accretion_rate_limits_dict[(Teff, logg)]

    def find_max_temperature(self, logg):
        if logg in self.temperature_limits_dict:
            pass
        else:
            # I think it should be safe to assume that the max temperature will always
            # be between 15000 and 25000K
            lower_bound = 15000
            upper_bound = 25000

            iterations_to_try = 100
            # ^ The parameter space will halve each iteration, so this effectively sets
            # the tolerance level. 100 seems like ample

            current_lower_bound = lower_bound
            current_upper_bound = upper_bound

            iteration_count = 0

            valid_max_limit = False

            low_logmdot_edge = 4
            high_logmdot_edge = 11

            while iteration_count < iterations_to_try:
                trial_upper_bound = (current_upper_bound + current_lower_bound) / 2
                interpolated_upper_result_1 = self.interpolate(
                    (low_logmdot_edge, trial_upper_bound, logg)
                )
                interpolated_upper_result_2 = self.interpolate(
                    (high_logmdot_edge, trial_upper_bound, logg)
                )
                # ^^ Need BOTH of these to be valid
                if np.isnan(interpolated_upper_result_1) or np.isnan(
                    interpolated_upper_result_2
                ):
                    # TOO HOT! we were out of bounds, and the trial value becomes the
                    # new upper limit:
                    current_upper_bound = trial_upper_bound
                else:
                    # Too cold - but this is a valid endpoint
                    valid_max_limit = True
                    current_lower_bound = trial_upper_bound
                iteration_count += 1
            if valid_max_limit:
                # Conservatively use the lower limit (to make sure we can only ever
                # interpolate within bounds)
                self.temperature_limits_dict[logg] = current_lower_bound
            else:
                self.temperature_limits_dict[logg] = np.nan
        return self.temperature_limits_dict[logg]

    def __call__(self, *args):
        interpolated_result = self.interpolate(*args)
        if not np.isnan(interpolated_result):
            return interpolated_result
        logMdot = args[0][0]
        Teff = args[0][1]
        logg = args[0][2]

        # We firstly need to find the max usable teff, for this log(g) - define this as
        # the highest teff which gives valid therm factors at both logmdot = 4 and
        # logmdot = 11 (basically the full range of the grid). Then Step 1 is to snap to
        # this teff if we're above it.
        # After that, proceed as before
        max_teff = self.find_max_temperature(logg)
        if Teff > max_teff:
            Teff = max_teff

        min_accretion_rate, max_accretion_rate = self.find_accretion_rate_limits(
            Teff, logg
        )
        if logMdot > max_accretion_rate:
            # Snap to edge of grid
            return self.interpolate((max_accretion_rate, Teff, logg))
        elif logMdot < min_accretion_rate:
            # We'll do something v simple! Take two points and interpolate linearly
            # between them (BUT only for low Mdot!!!)
            point1_mdot = min_accretion_rate
            point1_therm_factor = self.interpolate((min_accretion_rate, Teff, logg))

            point2_mdot = 0
            # By assumption (i.e., we'll force thermohaline mixing to be completely
            # inactive for logMdot = 0)
            point2_therm_factor = 0.0

            if logMdot < point2_mdot:
                return 0.0

            fraction_to_interpolate = (min_accretion_rate - logMdot) / (
                min_accretion_rate - point2_mdot
            )

            interpolated_therm_factor = (
                point1_therm_factor
                + fraction_to_interpolate * (point2_therm_factor - point1_therm_factor)
            )

            return interpolated_therm_factor
        else:
            # No special action needed, hopefully
            return self.interpolate((logMdot, Teff, logg))


def example():
    therm_interpolator = ThermohalineInterpolator()
    therm_factor = therm_interpolator((0, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((1, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((2, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((3, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((4, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((5, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((6, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((7, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((8, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((9, 15270, 8.09))
    print(therm_factor)
    therm_factor = therm_interpolator((10, 15270, 8.09))
    print(therm_factor)


def main():
    example()


if __name__ == "__main__":
    main()
