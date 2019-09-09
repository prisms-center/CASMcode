from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, gzip,re
import numpy as np
from casm.quantumespresso.qeio import outfile, poscar

class QErunError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class QErun:
    """ An object containing values read from outfile

        The QErun class contains:
            self.total_energy: final total energy ('e_wo_entrp')
            self.forces: final forces on atoms (2d list of double)
            self.atom_type: list of atom types (list of str)
            self.atoms_per_type: list of number of atoms for each atom type (list of int)
            self.lattice:  final lattice (2d list of double)
            self.rec_lat: final reciprocal lattice (2d list of double)
            self.basis: final basis (2d list of double)
            self.coord_mode: coordinate mode for basis (always 'direct')
            self.is_complete: did the QE calculation run to completion? (bool)
            self.dos: Basic DOS output (3d list of double or bool)
            self.efermi: Fermi level for DOS use (double)
            self.dos_lm: l or lm-projected DOS output (4d list of double or bool)
            self.eigenvalues: eigenvalues and energies for doing band structure plots (list of 2d list of double)
            self.all_e_0: energy (e_0_energy) for each electronic step of each ionic step
            self.elec_conv: says whether the last scf steps converged within the maximum electron_maxstep
    """
    def __init__(self, outfilename, DOS=False, Band=False):
        """ Create a QErun object from an infile and outfile """

        #### Filename #########################
        self.outfilename = outfilename

        #### set default values ###############

        # final energy ('e_wo_entrp')
        self.total_energy = None

        # final forces on atoms (2d list)
        self.forces = []

        # list of atom types
        self.atom_type = []

        # list of number of atoms for each atom type
        self.atoms_per_type = []

        # final lattice (2d list)
        self.lattice = None

        # final rec lattice (2d list)
        self.rec_lat = None

        # final basis (2d list)
        self.basis = None

        # coordinate mode for basis
        self.coord_mode = ''

        # did the QE calculation run to completion?
        self.is_complete = False

        # Is DOS present?
        self.DOS = DOS                # Whether to even read DOS; grabbed from init args
        self.dos = None

        # What about band structure?
        self.Band = Band
        if self.Band:
            self.DOS = True        # Need to read DOS to do some band things


        # What is the fermi level?
        self.efermi = None

        # is l or lm-projected DOS present?
        self.dos_lm = None

        # Energy (e_0_energy) for each electronic step for each ionic step
        self.all_e_0 = []

        ##### read from infile and outfile ########
        self.iter_read()

    def is_complete(self):
        """ Return True if QE calculation ran to completion """
        return self.is_complete

    def iter_read(self):
        """ Create a QErun object from an  outfile """
        myoutfile = outfile.Outfile(self.outfilename)
        pos = poscar.Poscar(self.outfilename)
        self.is_complete = myoutfile.complete
        self.all_e_0 = myoutfile.E
        self.total_energy = self.all_e_0[-1]
        self.lattice = pos._lattice.tolist()
        self.rec_lat = pos._reciprocal_lattice
        self.coord_mode = pos.coord_mode
        self.atom_type = pos.type_atoms
        self.atoms_per_type = pos.num_atoms
        self.elec_conv = True 
        
        ###Read things that haven't been taken care of by other functions
        ### i.e. efermi, dos_lm, Band, DOS, dos, forces, nstep
        if os.path.isfile(self.outfilename):
            if self.outfilename.split(".")[-1].lower() == "gz":
                f = gzip.open(self.outfilename)
            else:
                f = open(self.outfilename)
        elif os.path.isfile(self.outfilename+".gz"):
            f = gzip.open(self.outfilename+".gz")
        else:
            raise QErunError("file not found: " + self.outfilename)

        line=f.readline()
        m = re.search("Final en.*",line)
        while not m and self.elec_conv:
               line = f.readline()
               if line=='':
                    self.elec_conv = False
               line = line.strip()
               m = re.search("Final en.*",line)
        f.seek(0)
        
        line=f.readline()
        m = re.search("highest occupied level.*|.*Fermi energy is.*",line)
        if not m:
            while not m:
               line = f.readline()
               if line=='':
                    raise QErunError("EOF reach without finding highest occupied level or Fermi energy")
               line = line.strip()
               m = re.search("highest occupied level.*|.*Fermi energy is.*",line)
        try:
            if line.split()[-1]=="ev":
                self.efermi = float(line.split()[-2])
            else:
                self.efermi = float(line.split()[-1])
        except ValueError:
            raise QErunError("could not convert efermi to float")

        line=f.readline()
        m = re.search("Forces acting on atoms.*",line)
        if not m:
            while not m:
               line = f.readline()
               if line=='':
                    raise QErunError("EOF reach without finding Forces acting on atoms")
               line = line.strip()
               m = re.search("Forces acting on atoms.*",line)
        forces=[]
        f.readline()
        for i in range(sum(self.atoms_per_type)):
            line=f.readline().split()
            forces+=[[float(line[-3]),float(line[-2]),float(line[-1])]]
        
        f.close()
        #Band & DOS stuff not implemented

        #parse basis from poscar object into 2d list
        names=[]
        for site in pos.basis:
            names+= [site.occ_alias]
        sortednames=sorted(names,key=lambda x: self.atom_type.index(x))
        indexmap=sorted(range(len(names)),key=lambda x: sortednames.index(names[x]))
        rawlist=[]
        for site in pos.basis:
            rawlist += [site.position]
        #make sure basis and forces correspond
        self.basis = map(lambda x:rawlist[x].tolist() ,indexmap)

        self.forces = map(lambda x:forces[x] ,indexmap)


