""" Script for running vasp convergences """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import argparse
import casm.vaspwrapper

def main():
    """ Execute vasp.converge """
    parser = argparse.ArgumentParser(
        description='''Set up, run, and collect data on a VASP Convergence over
                       a specified property.''')

    parser.add_argument(
        'configdir',
        help='A path to a configuration in a CASM project')
    parser.add_argument(
        '-c', '--collect', action='store_true',
        help='Switch to collection mode instead of submission mode (default)')

    args = parser.parse_args()

    print("Begin vasp.converge")

    if args.collect:
        convergence = casm.vaspwrapper.Converge(args.configdir)

        convergence.collect()

        print("Finish vasp.converge\n\n")
    else:

        convergence = casm.vaspwrapper.Converge(args.configdir)

        convergence.submit()

        print("Finish vasp.converge\n\n")


if __name__ == '__main__':
    main()
