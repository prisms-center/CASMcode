#!/usr/bin/env python
from casm import noindent
import casm.vaspwrapper
import sys
import json

def main():
    print "Begin vasp.relax.report"

    if len(sys.argv) != 2:
        print "Usage: vasp.relax.report configdir"
        sys.exit()

    configdir = sys.argv[1]

    try:
      output = casm.vaspwrapper.Relax.properties(configdir)
    except:
      print("Unable to report properties for directory %s. Please verify that it contains a completed VASP calculation."%configdir)
      raise

    with open('properties.calc.json', 'w') as file:
      file.write(json.dumps(output, file, cls=noindent.NoIndentEncoder, indent=4, sort_keys=True))

    print "Finish vasp.relax.report\n\n"

if __name__ == "__main__":
    main()
