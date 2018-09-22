"""test_casm/test_scripts/test_casm_learn.py"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import unittest
import os
import json
import numpy as np
import shutil
import six
import sys
from os.path import join

from casm import project
from casm.misc import contexts
from casm.scripts import casm_learn

import test_casm
from test_casm.test_scripts import CasmScriptsTestCase

class TestCasmLearn(CasmScriptsTestCase):
    """Test casm-learn"""

    def setUp(self):
        pass

    def test_RFE(self):
        """Test casm-learn --exRFE"""
        if self.has_projects:

            # ZrO test project construction
            proj = project.Project(self.ZrO_dir, verbose=False)

            # files and directories
            fit_dir = join(proj.path, 'fit')
            fit_RFE = join(fit_dir, 'fit_RFE.json')
            default_bset = join(self.ZrO_dir, 'basis_sets', 'bset.default')
            default_bspecs = join(default_bset, 'bspecs.json')
            test_bset = join(self.ZrO_dir, 'basis_sets', 'bset.test')
            test_bspecs = join(test_bset, 'bspecs.json')
            test_clex = join(self.ZrO_dir, 'cluster_expansions', 'clex.formation_energy',
                            'calctype.default', 'ref.default', 'bset.test')
            test_eci = join(test_clex, 'eci.test', 'eci.json')

            # setup
            def _clean():
                if os.path.exists(fit_dir):
                    shutil.rmtree(fit_dir)
                if os.path.exists(test_bset):
                    shutil.rmtree(test_bset)
                if os.path.exists(test_clex):
                    shutil.rmtree(test_clex)

            _clean()
            os.mkdir(fit_dir)

            # test 'casm-learn --exREF'
            testargs = ["casm-learn", "--exRFE"]
            with contexts.captured_output(wd=fit_dir) as (sout, serr):
                with contexts.patch.object(sys, 'argv', testargs):
                    casm_learn.main()
            res = json.loads(sout.getvalue())
            self.assertEqual(res['feature_selection']['method'], 'RFE')
            self.assertEqual(res['estimator']['method'], 'LinearRegression')

            # set random_state for CV score
            res['problem_specs']['cv']['kwargs']['random_state'] = 0

            # have 'checkhull' only checking training configs
            res['checkhull'] = {
                'selection':'train',
                'write_results': True}

            # save 'fit_RFE.json'
            with open(fit_RFE, 'wb') as f:
                f.write(six.u(json.dumps(res, indent=2)).encode('utf-8'))

            # create 'test' bset and eci
            proj.capture('settings --new-bset test')
            proj.capture('settings --new-eci test')
            shutil.copyfile(default_bspecs, test_bspecs)
            stdout, stderr, returncode = proj.capture('bset -u')
            #print("OK 1")

            # create 'train' selection
            with contexts.captured_output(wd=fit_dir) as (sout, serr):
                proj.capture("select --set 'and(is_calculated,lt(comp(a),0.695))' -o train")
            #print("OK 2")

            # fit
            testargs = 'casm-learn -s fit_RFE.json --quiet'.split()
            #print testargs
            with contexts.captured_output(wd=fit_dir) as (sout, serr):
                with contexts.patch.object(sys, 'argv', testargs):
                    casm_learn.main()
            self.assertTrue(os.path.exists(join(fit_dir, 'fit_RFE_halloffame.pkl')))
            self.assertTrue(os.path.exists(join(fit_dir, 'fit_RFE_specs.pkl')))
            #print("OK 3")

            # view hall of fame
            testargs = 'casm-learn -s fit_RFE.json --hall --quiet'.split()
            #print testargs
            with contexts.captured_output(wd=fit_dir) as (sout, serr):
                with contexts.patch.object(sys, 'argv', testargs):
                    casm_learn.main()
            lines = sout.getvalue().splitlines()
            res = dict(zip(lines[0].split(), lines[2].split()))
            self.assertEqual(len(lines), 3)
            self.assertEqual(res['Selected'], '0111001101000000001101001001000100010010...')
            self.assertEqual(int(res['#Selected']), 25)
            self.assertAlmostEqual(float(res['CV']), 0.016540883) # depends on 'random_state'
            self.assertAlmostEqual(float(res['RMS']), 0.014697276)
            self.assertAlmostEqual(float(res['wRMS']), 0.014697276)
            self.assertEqual(res['Estimator'], 'LinearRegression')
            self.assertEqual(res['FeatureSelection'], 'RFE')
            #print("OK 4")

            # check hull
            testargs = 'casm-learn -s fit_RFE.json --checkhull --indiv 0 --quiet'.split()
            #print testargs
            with contexts.captured_output(wd=fit_dir) as (sout, serr):
                with contexts.patch.object(sys, 'argv', testargs):
                    casm_learn.main()
            self.assertTrue(os.path.exists(join(fit_dir, 'checkhull_fit_RFE_0_dft_gs')))
            self.assertTrue(os.path.exists(join(fit_dir, 'checkhull_fit_RFE_0_clex_gs')))
            #print("OK 5")

            # select ECI
            testargs = 'casm-learn -s fit_RFE.json --select 0'.split()
            #print testargs
            with contexts.captured_output(wd=fit_dir) as (sout, serr):
                with contexts.patch.object(sys, 'argv', testargs):
                    casm_learn.main()
            self.assertTrue(os.path.exists(test_eci))
            #print("OK 6")
