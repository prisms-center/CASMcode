""" lcaoOUT class and associated functions, methods, and error class """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re
from .geom import Cell

class LcaoOUT(object):   #pylint: disable=too-few-public-methods
    """ Container object for reading and parsing lcao.out files """

    def __init__(self, filename="lcao.out"):
        self.filename = filename
        self.total_energy = None
        self.forces = []
        self.all_e_0 = []
        self.complete = False
        self.converged = True
        self.coord_mode = "cartesian"
        self.cell = None
        self.scale = 1.0

        self.read(filename)

    def read(self, filename="lcao.out"):
        """ Reads a lcao.out file and parses all the tags """
        ### Read some sort of file ###
        with open(filename) as lcao:
            while True:
                line = lcao.readline()
                if line == "":
                    break
                # if re.search(r"FINAL\s*RELAXED\s*ENERGY", line):
                #     self.total_energy = float(line.strip().split()[-1])
                if re.search(r"Binding", line):
                    self.total_energy = float(line.strip().split()[-2])
                elif re.search(r"total\s*force", line):
                    lcao.readline()
                    while True:
                        line = lcao.readline()
                        if line == "":
                            break
                        if re.search(r"STR\(Ry\)=", line):
                            break
                        elif re.search(r"f-defect", line):
                            break
                        self.forces.append([float(i) for i in line.strip().split()[1:]])
                elif re.search(r"\s*Binding\s*energy", line):
                    self.all_e_0.append(float(line.strip().split()[3]))
                elif re.search(r"\s*coordinate", line):
                    self.coord_mode = lcao.readline().strip().lower()
                elif re.search(r"^\s*primitive", line, re.IGNORECASE):
                    self.cell = Cell.seq(lcao)
                elif re.search(r"^\s*FINAL\s*RELAXED\s*ENERGY", line, re.IGNORECASE):
                    self.complete = True
                elif re.search(r"^\s*Normal\s*exit:\s*NO\s*CONVG", line, re.IGNORECASE):
                    self.converged = False
                elif re.search(r"^\s*scale\s*", line, re.IGNORECASE):
                    if line.split()[0] == "scale":
                        self.scale = float(lcao.readline().strip())
