"""test_casm/misc.py"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import inspect
from os.path import abspath, basename, exists, dirname, join
import unittest

def check_skip_all(test):
    skipfile = join(test.classdir,'skip')
    if exists(skipfile):
        raise unittest.SkipTest("'" + basename(skipfile) + "' file found. Skipping " + test.__class__.__name__)

def check_skip_one(test):
    skipfile = join(test.classdir,'skip_' + test.__class__.__name__)
    if exists(skipfile):
        raise unittest.SkipTest("'" + basename(skipfile) + "' file found. Skipping " + test.__class__.__name__)

def casm_setup(self):
    """Implements common setup for casm tests
       - check for 'skip' or 'skip_MyTestCase' files
    
    Notes:
        Standalone implementation to allow easier use by subpackages
    """  
    self.classfile = abspath(inspect.getfile(self.__class__))
    self.classdir = dirname(self.classfile)
    
    check_skip_all(self)
    check_skip_one(self)

class CasmTestCase(unittest.TestCase):
    """CASM base unittest class
    
    Notes:
        From (9/15/2027): https://gist.github.com/twolfson/13f5f5784f67fd49b245
    """
    @classmethod
    def setUpClass(cls):
        """On inherited classes, run our `setUp` method"""
        # Inspired via http://stackoverflow.com/questions/1323455/python-unit-test-with-base-and-sub-class/17696807#17696807
        if cls is not CasmTestCase and cls.setUp is not CasmTestCase.setUp:
            orig_setUp = cls.setUp
            def setUpOverride(self, *args, **kwargs):
                CasmTestCase.setUp(self)
                return orig_setUp(self, *args, **kwargs)
            cls.setUp = setUpOverride

    def setUp(self):
        """Common Setup: 
            - check for 'skip' or 'skip_MyTestCase' files
        """
        casm_setup(self)
        
        
