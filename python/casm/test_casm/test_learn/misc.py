"""test_casm/test_learn/misc.py"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import unittest
import warnings

from distutils.spawn import find_executable
from os.path import join

import test_casm
from test_casm.test_project import casm_project_setup

def casm_learn_setup(self):
    """Implements common setup for casm.learn tests
      - check for 'skip' or 'skip_MyTestCase' files
      - check for 'CASM_TEST_PROJECTS_DIR' and set 'self.has_projects'

    Notes:
        Uses test_casm.test_project.casm_project_setup
    """

    # Use same setup as casm.project tests
    casm_project_setup(self)


class CasmLearnTestCase(unittest.TestCase):
    """test_casm.test_learn base unittest class

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
        if cls is not CasmLearnTestCase and cls.setUp is not CasmLearnTestCase.setUp:
            orig_setUp = cls.setUp
            def setUpOverride(self, *args, **kwargs):
                CasmLearnTestCase.setUp(self)
                return orig_setUp(self, *args, **kwargs)
            cls.setUp = setUpOverride

    def setUp(self):
        """Common Setup:
          - check for 'skip' or 'skip_MyTestCase' files
          - check for 'CASM_TEST_PROJECTS_DIR' and set 'self.has_projects'
        """
        casm_learn_setup(self)
