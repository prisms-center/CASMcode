from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import casm.questwrapper
import sys

def main():
    print("Begin quest.relax")

    if len(sys.argv) != 2:
        print("Usage: quest.relax configdir")
        sys.exit()

    configdir = sys.argv[1]

    relaxation = casm.questwrapper.Relax(configdir)

    relaxation.submit()

    print("Finish quest.relax\n\n")

if __name__ == "__main__":
    main()
