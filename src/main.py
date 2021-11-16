#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse as ap
import manager as mn
import pwd_utils as pu

def main():
    arguments = pu.parse_command_line_arguments()
    manager = mn.Manager(arguments)

    # Change arguments to the IDs of the system(s) to run:
    manager.run([360, 376]) # e.g. this is GD61 and WD0449

if __name__ == '__main__':
    main()
