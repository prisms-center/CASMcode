"""test_casm/test_vasp/test_vasp.py"""
import unittest
import os
from os.path import join, abspath, dirname
import inspect
import json

from casm import vasp
from casm.misc.contexts import working_dir, captured_output, print_stringIO
from misc import cp_input, before_all, has_vasp

class TestVasp(unittest.TestCase):
  
    def setUp(self):
        """
        Read test case data
        """
        before_all()
        self.currdir = dirname(abspath(inspect.getfile(inspect.currentframe()))) # script directory
        with open(join(self.currdir, 'test_cases.json'), 'r') as f:
            self.cases = json.load(f)["vasp"]
    
    def test_run(self):
        """
        Test vasp.run()
        """
        for i,case in enumerate(self.cases["run"]):
            input = join(self.currdir, 'input_data', case['input_data'])
            output = join(self.currdir, 'output_data','vasp','run','test_case_' + str(i))
            cp_input(input, output)
            with captured_output() as (sout, serr):
                if has_vasp:
                    vasp.run(output, **case["settings"])
                    self.assertTrue(True) # Todo: add tests
                else:
                    self.assertRaisesRegexp(OSError, "No such file or directory", 
                        vasp.run, output, **case["settings"])
            #print_stringIO(sout)
            #print_stringIO(serr)
            


