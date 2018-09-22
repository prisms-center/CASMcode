"""test_casm/test_project/test_selection.py"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import unittest
import os
from os.path import join
import json

import numpy as np
from casm import project

import test_casm
from test_casm.test_project import CasmProjectTestCase

class TestCasmSelection(CasmProjectTestCase):

    def setUp(self):
        pass

    def test_init(self):
        """Test project.Selection()"""
        if self.has_projects:
            # ZrO test project construction
            proj = project.Project(self.ZrO_dir, verbose=False)
            sel = project.Selection(proj)
            self.assertIs(sel.proj, proj)
            self.assertEqual(sel.path, "MASTER")
            self.assertEqual(sel.all, True)

    def test_data(self):
        """Test Selection.data"""
        if self.has_projects:
            # ZrO test project construction
            proj = project.Project(self.ZrO_dir, verbose=False)
            proj.capture("select --set-on")
            sel = project.Selection(proj)

            Nconfig = 336
            self.assertEqual(sel.all, True)
            self.assertEqual(sel.data.shape, (Nconfig,2))
            self.assertEqual(sel.data.columns[0], 'configname')
            self.assertEqual(sel.data.columns[1], 'selected')
            self.assertEqual(sel.data.dtypes[0], 'object')
            self.assertEqual(sel.data.dtypes[1], 'bool')
            self.assertEqual(list(sel.data['selected']), [True]*Nconfig)

            # selections are not automatically updated
            proj.capture("select --set-off")
            self.assertEqual(list(sel.data['selected']), [True]*Nconfig)

            # and query will not overwrite existing data
            sel.query(["selected"])
            self.assertEqual(list(sel.data['selected']), [True]*Nconfig)

            # unless you query w/force
            sel.query(["selected"], force=True)
            self.assertEqual(list(sel.data['selected']), [False]*Nconfig)

            # this will add a column
            sel.query(['comp'])
            self.assertEqual(sel.data.shape, (Nconfig,3))
            self.assertEqual(sel.data.columns[2], 'comp(a)')
            self.assertEqual(sel.data.dtypes[2], 'float64')

            # Check 'all' constructor arg
            proj.capture("select --set 'lt(scel_size,3)'")
            sel = project.Selection(proj, all=True)
            self.assertEqual(sel.data.shape, (Nconfig,2))
            sel = project.Selection(proj, all=False)
            self.assertEqual(sel.data.shape, (13,2))

    def test_save(self):
        """Test Selection.save, Selection.saveas"""
        if self.has_projects:
            # ZrO test project construction
            proj = project.Project(self.ZrO_dir, verbose=False)
            proj.capture("select --set-on")

            # we'll set 'selected' based on the scel_size
            Nconfig = 336
            sel = project.Selection(proj)
            self.assertEqual(sel.data.shape, (Nconfig,2))
            self.assertEqual(sel.data['selected'].sum(), Nconfig)
            sel.query(['scel_size'])
            sel.data.loc[:,'selected'] = sel.data['scel_size'] < 3
            self.assertEqual(sel.data['selected'].sum(), 13)

            # saveas new selection file
            test_select = join(proj.path, 'test_select')
            if os.path.isfile(test_select):
                os.remove(test_select)
            new_sel = sel.saveas(test_select)
            self.assertEqual(new_sel.data.shape, (336,3))
            self.assertEqual(new_sel.data['selected'].sum(), 13)

            # saveas an existing file should raise
            self.assertRaises(Exception, lambda x: sel.saveas(test_select))

            # then save the original selection back to the MASTER list
            sel.save()

            # then reload and check
            read_sel = project.Selection(proj, all=False)
            self.assertEqual(read_sel.data.shape, (13,2))

            # reload 'new_sel' and check
            read_new_sel = project.Selection(proj, test_select, all=False)
            self.assertEqual(read_new_sel.data.shape, (13,3))
            self.assertEqual(read_new_sel.data['selected'].sum(), 13)

            if os.path.isfile(test_select):
                os.remove(test_select)
