"""test_casm/test_vasp/misc.py"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import distutils.spawn
import os
from os.path import exists, dirname
import unittest
import shutil
import warnings

import test_casm

def cp_input(input, output):
    """Copy input directory tree to output directory tree
    
    Notes:
        - If output directory tree exists, rm existing first, then copy
        - Constructs any immediate directories necessary
    """
    if exists(output):
        shutil.rmtree(output)
    if not exists(dirname(output)):
        os.makedirs(dirname(output))
    shutil.copytree(input, output)

def casm_vasp_setup(self):
    """Implements common setup for casm.vasp tests
      - check for 'skip' or 'skip_MyTestCase' files
      - check for 'vasp' executable and set 'self.has_vasp'
      - check for 'CASM_VASP_POTCAR_DIR' and set 'self.has_potcars'
    
    Notes:
        Standalone implementation to allow easier use by subpackages
    """
    
    # First run common setup for 'casm'
    test_casm.casm_setup(self)
    
    # Now run common setup for 'casm.vasp'
    
    # Check for 'vasp' executable
    self.has_vasp = distutils.spawn.find_executable('vasp') is not None
    if not self.has_vasp:
        warnings.warn("\n'vasp' executable not detected: will test behaviour for system without VASP")

    # Check for 'CASM_VASP_POTCAR_DIR' environment variable
    self.has_potcars = 'CASM_VASP_POTCAR_DIR' in os.environ
    if not self.has_potcars:
        warnings.warn("\n'CASM_VASP_POTCAR_DIR' environment variable not found: will test behaviour that does not require POTCARs")


class CasmVaspTestCase(unittest.TestCase):
    """test_casm.test_vasp base unittest class
    
    Attributes:
        has_vasp (bool): True if 'vasp' executable exist. Some tests require
            that VASP is installed, but others testing input/output file io do not.
        has_potcars (bool): True if 'CASM_VASP_POTCAR_DIR' exists. Some tests require
            POTCARs, but others do not.
    """
    
    @classmethod
    def setUpClass(cls):
        """On inherited classes, run our `setUp` method"""
        if cls is not CasmVaspTestCase and cls.setUp is not CasmVaspTestCase.setUp:
            orig_setUp = cls.setUp
            def setUpOverride(self, *args, **kwargs):
                CasmVaspTestCase.setUp(self)
                return orig_setUp(self, *args, **kwargs)
            cls.setUp = setUpOverride

    def setUp(self):
        """Common Setup: 
          - check for 'skip' or 'skip_MyTestCase' files
          - check for 'vasp' executable and set 'self.has_vasp'
          - check for 'CASM_VASP_POTCAR_DIR' and set 'self.has_potcars'
        """
        casm_vasp_setup(self)
