"""test_casm/test_project/test_project.py"""
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

class TestCasmProject(CasmProjectTestCase):
  
    def setUp(self):
        pass
    
    def test_init(self):
        """Test project.Project()"""
        if self.has_projects:
            # ZrO test project construction
            proj = project.Project(self.ZrO_dir, verbose=False)
            self.assertEqual(proj.path, self.ZrO_dir)
            self.assertEqual(proj.name, "ZrO")

    def test_prim(self):
        """Test Prim"""
        if self.has_projects:
            proj = project.Project(self.ZrO_dir, verbose=False)
            prim = proj.prim
            self.assertIs(prim.proj, proj)
            self.assertIs(proj.prim, prim)
            self.assertTrue(
                np.allclose(prim.lattice_matrix, 
                            np.array([[ 3.23398686, -1.61699343, -0.],
                                      [ 0.,          2.80071477,  0.        ],
                                      [ 0.,          0.,          5.16867834]])))
            self.assertAlmostEqual(prim.lattice_parameters['a'], 3.23398686, places=6)
            self.assertAlmostEqual(prim.lattice_parameters['b'], 3.23398686, places=6)
            self.assertAlmostEqual(prim.lattice_parameters['c'], 5.16867834, places=6)
            self.assertAlmostEqual(prim.lattice_parameters['alpha'], 90., places=6)
            self.assertAlmostEqual(prim.lattice_parameters['beta'], 90., places=6)
            self.assertAlmostEqual(prim.lattice_parameters['gamma'], 120., places=6)
            #print prim.basis
            self.assertEqual(prim.coordinate_mode, 'Fractional')
            self.assertEqual(prim.lattice_symmetry_s, 'D6h')
            self.assertEqual(prim.lattice_symmetry_hm, '6/mmm')
            self.assertEqual(prim.lattice_system, 'hexagonal')
            self.assertEqual(prim.crystal_symmetry_s, 'D6h')
            self.assertEqual(prim.crystal_symmetry_hm, '6/mmm')
            self.assertEqual(prim.crystal_system, 'hexagonal')
            self.assertEqual(prim.crystal_family, 'hexagonal')
            self.assertEqual(prim.space_group_number, '191:194')
            self.assertEqual(prim.components, ['Zr', 'Va', 'O']) 
            self.assertEqual(prim.elements, ['Zr', 'Va', 'O'])
            self.assertEqual(prim.n_independent_compositions, 1)
            self.assertEqual(prim.degrees_of_freedom, ['occupation'])
    
    def test_composition_axes(self):
        """Test CompositionAxes"""
        if self.has_projects:
            proj = project.Project(self.ZrO_dir, verbose=False)
            comp_axes = proj.composition_axes
            self.assertIs(proj.composition_axes, comp_axes)
            self.assertEqual(comp_axes.name, '0')
            self.assertEqual(comp_axes.components, ['Zr', 'Va', 'O'])
            self.assertEqual(comp_axes.n_independent_compositions, 1)
            self.assertEqual(comp_axes.mol_formula, 'Zr(2)Va(2-2a)O(2a)')
            self.assertEqual(comp_axes.param_formula, 'a(0.5-0.25Va+0.25O)')
            self.assertTrue(np.allclose(comp_axes.end_members['origin'],
                            np.array([ 2.,  2.,  0.])))
            self.assertTrue(np.allclose(comp_axes.end_members['a'],
                            np.array([ 2.,  0.,  2.])))
            
