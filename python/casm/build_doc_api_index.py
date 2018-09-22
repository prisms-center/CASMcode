#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import glob
import importlib
import os
import setuptools
import sys
import time

def write_index(f, package):
    if hasattr(package, '__all__'):
        f.write(package.__name__ + '\n')
        f.write('-'*len(package.__name__) + '\n\n')
        f.write('.. autosummary::\n')
        f.write('    :toctree:\n\n')
        for member in package.__all__:
            f.write('    ' + package.__name__ + '.' + member + '\n')
        f.write('\n')

def main():
    with open('doc/source/python/index.rst','w') as f:
        f.write('.. python/index.rst\n\n')

        f.write('CASM Python packages documentation\n')
        f.write('==================================\n\n')

        packages = sorted(setuptools.find_packages('casm'))
        for packagename in packages:
            print('casm' + '.' + packagename)
            #sys.stdout.flush()
            package = importlib.import_module('casm' + '.' + packagename)
            if hasattr(package, '__all__') and len(package.__all__):
                print(package.__all__)
                write_index(f, package)


if __name__ == "__main__":
    main()
