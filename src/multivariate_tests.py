#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import Enum

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.packages import PackageNotInstalledError

imported_correctly = True
print('Importing R modules...')
try:
    _cramer = importr('cramer')  # Making this global - they take a while to import
    #_fasano_franceschini = importr('fasano.franceschini.test') # Feel free to add this back in, it seemed insensitive to the differences I saw in the data though
except PackageNotInstalledError:
    imported_correctly = False
print('Done')

mv_test_display_strings = {
    0: 'Cramér'
#    1: 'Fasano Franceschini'
}

class MultiVariateTest(Enum):
    cramer = 0
    ff = 1

    def __str__(self):
        return self.name

    def get_display_string(self):
        return mv_test_display_strings[self.value]

def convert_data_to_r_format(list_of_data_lists):
    # Assuming the input is in this format (using a 3d example with 5 data points)
    #list_of_data_lists = [
    #    [1, 2, 3, 4, 5], # <-- all the x values
    #    [10, 20, 30, 40, 50], # <-- all the y values
    #    [11, 21, 31, 41, 51] # <-- all the z values
    #] # so the first data point was (1, 10, 11) etc
    len_of_first_list = len(list_of_data_lists[0])
    concatenated_list = list()
    for data_list in list_of_data_lists:
        assert len(data_list) == len_of_first_list
        concatenated_list += data_list
    v = robjects.FloatVector(concatenated_list)
    m = robjects.r['matrix'](v, nrow=len_of_first_list)
    return m

def convert_rlistvector_to_dict(rlistvector):
    return {key : rlistvector.rx2(key)[0] for key in rlistvector.names}

def cramer_test_p_value(data1, data2, confidence_level=0.95, replicate_count=1000, bootstrap_type="ordinary"):
    # bootstrap_type can be "ordinary", "permutation" or "eigenvalue"
    just_statistic = False
    result = _cramer.cramer_test(
        data1,
        data2,
        **{'conf.level': confidence_level, 'replicates': replicate_count, 'sim': bootstrap_type, 'just.statistic': just_statistic}
    )
    result_dict = convert_rlistvector_to_dict(result)
    return result_dict['p.value']

def fasano_franceschini_test_p_value(data1, data2, n_permutations=100):
    verbose = False
    result = _fasano_franceschini.fasano_franceschini_test(
        data1,
        data2,
        **{'nPermute': n_permutations, 'verbose': verbose}
    )
    result_dict = convert_rlistvector_to_dict(result)
    return result_dict['p.value']

def demo():
    data_set1 = [
        [1, 2, 3, 4, 5],
        [10, 20, 30, 40, 50],
        [11, 21, 31, 41, 51]
    ]
    data_set2 = [
        [2, 1, 6, 1.7],
        [11, 28, 20, 50],
        [11, 22, 31, 42]
    ]
    data_set1_r = convert_data_to_r_format(data_set1)
    data_set2_r = convert_data_to_r_format(data_set2)
    p_value_cramer = cramer_test_p_value(
        data_set1_r,
        data_set2_r
    )
    print(p_value_cramer)
    #p_value_fasano_franceschini = fasano_franceschini_test_p_value(
    #    data_set1_r,
    #    data_set2_r
    #)
    #print(p_value_fasano_franceschini)

def get_mv_test_method_dict():
    #return {MultiVariateTest.cramer: cramer_test_p_value, MultiVariateTest.ff: fasano_franceschini_test_p_value}
    return {MultiVariateTest.cramer: cramer_test_p_value}

def main():
    demo()

if __name__ == '__main__':
    main()
