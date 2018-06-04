from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import casm.qewrapper
import sys

def main():
    print("Begin qe.relax")

    if len(sys.argv) != 2:
        print("Usage: qe.relax configdir")
        sys.exit()

    configdir = sys.argv[1]

    print("  Construct a Quantum Espresso relaxation")
    relaxation = casm.qewrapper.Relax(configdir)

    print("  Submitting...")
    relaxation.submit()

    print("Finish qe.relax")

if __name__ == "__main__":
    main()
