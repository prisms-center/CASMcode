from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import sys
import subprocess
import time

def main():
    print("Begin casm.view")

    if len(sys.argv) != 3:
        print("Usage: casm.view '...command to call visualization software...' /path/to/POSCAR")
        sys.exit()

    VESTAcommand = sys.argv[1]
    POSCARpath = sys.argv[2]

    print("Reading POSCAR:")
    print(open(POSCARpath).read())
    subprocess.Popen((VESTAcommand + " " + POSCARpath).split())

    time.sleep(1)

if __name__ == "__main__":
    main()
