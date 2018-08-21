"""test_casm/test_project/misc.py"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import unittest
import warnings

from distutils.spawn import find_executable
from os.path import join

import test_casm

def casm_project_setup(self):
    """Implements common setup for casm.project tests
      - check for 'skip' or 'skip_MyTestCase' files
      - check for 'CASM_TEST_PROJECTS_DIR' and set 'self.has_projects'

    Notes:
        Standalone implementation to allow easier use by subpackages
    """

    # First run common setup for 'casm'
    test_casm.casm_setup(self)

    # Check for 'casm' executable
    self.has_casm = find_executable('casm') is not None
    if not self.has_casm:
        warnings.warn("\n'casm' executable not detected: will test behaviour that does not require libcasm")

    # Check for 'CASM_TEST_PROJECTS_DIR' environment variable
    self.has_projects = 'CASM_TEST_PROJECTS_DIR' in os.environ and len(os.environ['CASM_TEST_PROJECTS_DIR']) > 0
    if not self.has_projects:
        warnings.warn("\n'CASM_TEST_PROJECTS_DIR' environment variable not found: will test behaviour that does not require a test project")
    else:
        self.test_projects_dir = os.environ['CASM_TEST_PROJECTS_DIR']
        self.ZrO_dir = join(self.test_projects_dir, '0.3.X', 'ZrO.0')

class CasmProjectTestCase(unittest.TestCase):
    """test_casm.test_project base unittest class

    Attributes:
        has_projects (bool): True if 'CASM_TEST_PROJECTS_DIR' exists and is non-zero length. Some
            tests require a test project, but others do not.
        test_projects_dir (str): Value of 'CASM_TEST_PROJECTS_DIR'; the location
            of test projects.
        ZrO_dir (str): Location of a ZrO test project that contains VASP calculations.

    """

    @classmethod
    def setUpClass(cls):
        """On inherited classes, run our `setUp` method"""
        if cls is not CasmProjectTestCase and cls.setUp is not CasmProjectTestCase.setUp:
            orig_setUp = cls.setUp
            def setUpOverride(self, *args, **kwargs):
                CasmProjectTestCase.setUp(self)
                return orig_setUp(self, *args, **kwargs)
            cls.setUp = setUpOverride

    def setUp(self):
        """Common Setup:
          - check for 'skip' or 'skip_MyTestCase' files
          - check for 'CASM_TEST_PROJECTS_DIR' and set 'self.has_projects'
        """
        casm_project_setup(self)
