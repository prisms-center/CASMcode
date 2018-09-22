from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import sys

from casm.misc import compat, noindent
import casm.vaspwrapper

def main():
    print("Begin vasp.relax.report")

    if len(sys.argv) != 2:
        print("Usage: vasp.relax.report configdir")
        sys.exit()

    configdir = sys.argv[1]

    try:
      output = casm.vaspwrapper.Relax.properties(configdir)
    except:
      print(("Unable to report properties for directory %s. Please verify that it contains a completed VASP calculation."%configdir))
      raise

    compat.dump(json, output, 'properties.calc.json', 'w', cls=noindent.NoIndentEncoder, 
                indent=4, sort_keys=True)

    print("Finish vasp.relax.report\n\n")

if __name__ == "__main__":
    main()
