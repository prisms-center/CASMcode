from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re
import six

# Note: This could use double-checking

# Schoenflies to Hermann-Mauguin
symmap = {
  'C1': '1',
  'Ci': '-1',
  'S2': '-1',
  'C2': '2',
  'Cs': 'm',
  'C1h': 'm',
  'C2h': '2/m',
  'D2': '222',
  'V': '222',
  'C2v': 'mm2',
  'D2h': 'mmm',
  'Vh': 'mmm',
  'C4': '4',
  'S4': '-4',
  'C4h': '4/m',
  'D4': '422',
  'C4v': '4mm',
  'D2d': '-42m',
  'Vd': '-42m',
  'D4h': '4/mmm',
  'C3': '3',
  'S6': '-3',
  'C3i': '-3',
  'D3': '32',
  'C3v': '3m',
  'D3d': '-3m',
  'C6': '6',
  'C3h': '-6',
  'C6h': '6/m',
  'D6': '622',
  'C6v': '6mm',
  'D3h': '-6m2',
  'D6h': '6/mmm',
  'T': '23',
  'Th': 'm-3',
  'O': '432',
  'Td': '-43m',
  'Oh': 'm-3m',
}

# lattice system
latsysmap = {
  'Ci': 'triclinic',
  'S2': 'triclinic',
  'C2h': 'monoclinic',
  'D2h': 'orthorhombic',
  'D4h': 'tetragonal',
  'D3d': 'rhombohedral',
  'D6h': 'hexagonal',
  'Oh': 'cubic',
}

# crystal system to Schoenflies
xtalsysmap = {
  'triclinic': ['C1', 'Ci', 'S2'],
  'monoclinic': ['C2', 'Cs', 'C1h', 'C2h'],
  'orthorhombic': ['D2', 'V', 'C2v', 'D2h', 'Vh'],
  'tetragonal': ['C4', 'S4', 'C4h', 'D4', 'C4v', 'D2d', 'Vd', 'D4h'],
  'trigonal': ['C3', 'S6', 'C3i', 'D3', 'C3v', 'D3d'],
  'hexagonal': ['C6', 'C3h', 'C6h', 'D6', 'C6v', 'D3h', 'D6h'],
  'cubic': ['T', 'Th', 'O', 'Td', 'Oh']
}

# crystal family to Schoenflies
xtalfamilymap = {
  'triclinic': ['C1', 'Ci', 'S2'],
  'monoclinic': ['C2', 'Cs', 'C1h', 'C2h'],
  'orthorhombic': ['D2', 'V', 'C2v', 'D2h', 'Vh'],
  'tetragonal': ['C4', 'S4', 'C4h', 'D4', 'C4v', 'D2d', 'Vd', 'D4h'],
  'hexagonal': ['C3', 'S6', 'C3i', 'D3', 'C3v', 'D3d', 'C6', 'C3h', 'C6h', 'D6', 'C6v', 'D3h', 'D6h'],
  'cubic': ['T', 'Th', 'O', 'Td', 'Oh']
}

# Schoenflies to International Tables space group number (range)
space_group_number_map = {
  'C1': '1',
  'Ci': '2',
  'S2': '2',
  'C2': '3:5',
  'Cs': '6:9',
  'C1h': '6:9',
  'C2h': '10:15',
  'D2': '16:24',
  'V': '16:24',
  'C2v': '25:46',
  'D2h': '47:74',
  'Vh': '47:74',
  'C4': '75:80',
  'S4': '81:82',
  'C4h': '83:88',
  'D4': '89:98',
  'C4v': '99:110',
  'D2d': '111:122',
  'Vd': '111:122',
  'D4h': '123:142',
  'C3': '143:146',
  'S6': '147:148',
  'C3i': '149:155',
  'D3': '149:155',
  'C3v': '156:161',
  'D3d': '162:167',
  'C6': '168:173',
  'C3h': '174',
  'C6h': '175:176',
  'D6': '177:182',
  'C6v': '183:186',
  'D3h': '187:190',
  'D6h': '191:194',
  'T': '195:199',
  'Th': '200:206',
  'O': '207:214',
  'Td': '215:220',
  'Oh': '221:230',
}

def lattice_symmetry(stdout):
    """
    Parse output from 'casm sym' for lattice point group, return Schoenflies form
    """
    for line in stdout.split('\n'):
        if re.search('Lattice point group is:', line):
            return line.split()[-1]
    return None

def hm_symmetry(schoenflies_symbol):
    """
    Determine Hermann-Mauguin symbol from Schoenflies symbol
    """
    return symmap[schoenflies_symbol]

def lattice_system(schoenflies_symbol):
    """
    Determine lattice system from lattice symmetry
    """
    return latsysmap[schoenflies_symbol]

def crystal_symmetry(stdout):
    """
    Parse output from 'casm sym' for crystal point group, return Schoenflies form
    """
    for line in stdout.split('\n'):
        if re.search('Crystal point group is:', line):
            return line.split()[-1]
    return None

def crystal_system(schoenflies_symbol):
    """ 
    Determine crystal system from crystal symmetry
    """
    for key, val in six.iteritems(xtalsysmap):
        if schoenflies_symbol in val:
            return key
    return None

def crystal_family(schoenflies_symbol):
    """ 
    Determine crystal family from crystal symmetry
    """
    for key, val in six.iteritems(xtalfamilymap):
        if schoenflies_symbol in val:
            return key
    return None

