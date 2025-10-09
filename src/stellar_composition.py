#!/usr/bin/env python


class StellarComposition:

    def __init__(self, input_list_of_floats):
        self.composition = input_list_of_floats

    def __repr__(self):
        return f"Some SC {len(self.composition)} {self.composition[0]}"
