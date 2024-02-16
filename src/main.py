#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse as ap
import manager as mn
import pwd_utils as pu

def main():
    arguments = pu.parse_command_line_arguments()
    manager = mn.Manager(arguments)
    manager.run()

if __name__ == '__main__':
    main()
