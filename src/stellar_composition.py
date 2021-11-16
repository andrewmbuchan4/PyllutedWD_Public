#!/usr/bin/env python
# -*- coding: utf-8 -*-

class StellarComposition:
    
    def __init__(self, input_list_of_floats):
        self.composition = input_list_of_floats
    
    def __repr__(self):
        return 'Some SC ' + str(len(self.composition)) + ' ' + str(self.composition[0])
