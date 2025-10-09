#!/usr/bin/env python

import manager as mn
import pwd_utils as pu


def main():
    pu.set_up_configuration()
    manager = mn.Manager()
    manager.run()


if __name__ == "__main__":
    main()
