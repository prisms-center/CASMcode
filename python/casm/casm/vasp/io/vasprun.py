from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import xml.etree.cElementTree as etree
import os, gzip
import numpy as np

class VasprunError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class Vasprun:
    """ An object containing values read from vasprun.xml

        The Vasprun class contains:
            self.total_energy: final total energy ('e_wo_entrp')
            self.forces: final forces on atoms (2d list of double)
            self.atom_type: list of atom types (list of str)
            self.atoms_per_type: list of number of atoms for each atom type (list of int)
            self.lattice:  final lattice (2d list of double)
            self.rec_lat: final reciprocal lattice (2d list of double)
            self.basis: final basis (2d list of double)
            self.coord_mode: coordinate mode for basis (always 'direct')
            self.is_complete: did the VASP calculation run to completion? (bool)
            self.dos: Basic DOS output (3d list of double or bool)
            self.dos_efermi: Fermi level for DOS use (double)
            self.dos_lm: l or lm-projected DOS output (4d list of double or bool)
            self.eigenvalues: eigenvalues and energies for doing band structure plots (list of 2d list of double)
            self.all_e_0: energy (e_0_energy) for each electronic step of each ionic step
            self.nelm: NELM (max. number of electronic steps)
    """
    def __init__(self,filename, DOS=False, Band=False):
        """ Create a Vasprun object from a vasprun.xml file with name 'filename' """

        #### Filename #########################
        self.filename = filename

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
        self.coord_mode = 'direct'

        # did the VASP calculation run to completion?
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

        ##### read from vasprun.xml file ########
        self.iter_read()

    def is_complete(self):
        """ Return True if VASP calculation ran to completion """
        return self.is_complete

    def iter_read(self):
        """ Create a Vasprun object from a vasprun.xml file with name 'filename' """
        if os.path.isfile(self.filename):
            if self.filename.split(".")[-1].lower() == "gz":
                f = gzip.open(self.filename)
            else:
                f = open(self.filename)
        elif os.path.isfile(self.filename+".gz"):
            f = gzip.open(self.filename+".gz")
        else:
            raise VasprunError("file not found: " + self.filename)

        for event, elem in etree.iterparse(f):
                if elem.tag == 'parameters':
                    self.nelm = int(elem.find(".//i[@name='NELM']").text)

                if elem.tag == 'calculation':

                        #finding the energy of the relaxation
                        self.total_energy = None
                        energy = elem.findall('energy')
                        for child in energy[-1]:
                            if child.attrib['name']=='e_wo_entrp':
                                self.total_energy = float(child.text)
                                break
                            child.clear()

                        #finding the forces on the atoms
                        varray = elem.findall('varray')
                        forces = None
                        for v in varray:
                            if v.attrib['name'] == 'forces':
                                forces = v
                        self.forces = []
                        for child in forces:
                            self.forces.append([float(i) for i in child.text.strip().split()])
                            child.clear()
                        v.clear()

                        # find energy (e_0_energy) for each electronic step in this ionic step
                        self.all_e_0.append([])
                        for i in elem.findall("./scstep/energy/i[@name='e_0_energy']"):
                            self.all_e_0[-1].append(float(i.text))

                        if self.DOS:

                                # gather the DOS, if calculated and warranted
                                self.dos = None
                                self.dos_lm = None
                                self.efermi = None

                                is_dos = elem.find('dos')
                                if is_dos != None:
                                    for child in is_dos:
                                        if child.attrib['name'] == 'efermi':
                                            self.efermi = float(child.text)
                                            child.clear()
                                            break
                                        else:
                                            child.clear()

                                    total = is_dos.find('total').find('array')


                                    spins = len(total.find('set').findall('set'))
                                    points = len(total.find('set').findall('set')[0].findall('r'))
                                    fields = len(total.findall('field'))

                                    total_array = np.zeros([spins, points, fields])

                                    for s in range(spins):
                                        my_set = total.find('set').findall('set')[s].findall('r')
                                        for r in range(points):
                                             total_array[s,r] = map(float, my_set[r].text.split())

                                    self.dos = total_array


                                    partial = is_dos.find('partial')
                                    if partial != None:
                                        partial = partial.find('array')
                                    else:
                                        partial = None

                                    if partial:
                                        ions = len(partial.find('set').findall('set'))
                                        spins = len(partial.find('set').findall('set')[0].findall('set'))
                                        points = len(partial.find('set').findall('set')[0].findall('set')[0].findall('r'))
                                        fields = len(partial.findall('field'))

                                        partial_array = np.zeros([ions, spins, points, fields])

                                        for i in range(ions):
                                            i_set = partial.find('set').findall('set')[i].findall('set')
                                            for s in range(spins):
                                                my_set = i_set[s].findall('r')
                                                for r in range(points):
                                                     partial_array[i,s,r] = map(float, my_set[r].text.split())

                                        self.dos_lm = partial_array

                        if self.Band:
                                eig_parse_vals = elem.find('eigenvalues').find('array').find('set')
                                num_kpoints = len(eig_parse_vals.find('set').findall('set'))
                                num_bands = len(eig_parse_vals.find('set').find('set').findall('r'))
                                self.eigenvalues = [np.empty([num_kpoints,num_bands]),np.empty([num_kpoints,num_bands])]
                                for spin_ctr,spin in enumerate(eig_parse_vals):
                                    for kpoint_ctr,kpoint in enumerate(spin):
                                        for band_ctr,band in enumerate(kpoint):
                                            self.eigenvalues[spin_ctr][kpoint_ctr,band_ctr] = band.text.strip().split()[0]
                                eig_parse_vals.clear()

                elif elem.tag == 'kpoints':
                        if self.Band:
                                self.num_divisions = elem.find('generation').find('i')
                                if self.num_divisions != None:
                                        self.num_divisions = int(self.num_divisions.text.strip())
                                        kp_parse_vals = elem.find('varray') if elem.find('varray').attrib['name'] == 'kpointlist' else None
                                        if kp_parse_vals != None:
                                                num_paths = len(kp_parse_vals.findall('v'))/self.num_divisions
                                                kpoint_list = np.zeros([self.num_divisions*num_paths, 3])
                                                for i, row in enumerate(kp_parse_vals):
                                                        kpoint_list[i,:] = [float(val) for val in row.text.strip().split()]
                                                self.kpoint_divisions = np.empty([2*num_paths, 3])
                                                for i in range(num_paths):
                                                        self.kpoint_divisions[2*i] = kpoint_list[i*self.num_divisions,:]
                                                        self.kpoint_divisions[2*i+1] = kpoint_list[(i+1)*self.num_divisions-1, :]
                                                kp_parse_vals.clear()
                                else:
                                        self.kpoint_divisions = False

                elif elem.tag == 'atominfo':

                        #finding the atom_type and atoms_per_type
                        temp_node = elem.findall('array')
                        for child in temp_node:
                            if child.attrib['name'] == 'atomtypes':
                                child_node = child
                                break
                        t_atom_info = child_node.findall('set')[0].findall('rc')
                        self.atom_type = []
                        self.atoms_per_type = []
                        for child in t_atom_info:
                            self.atom_type.append(child[1].text)
                            self.atoms_per_type.append(int(child[0].text))


                elif elem.tag == 'structure':
                        self.is_complete = False
                        if 'name' in elem.attrib:
                                if elem.attrib['name'] == 'finalpos':
                                        #finding the final structure
                                        finalpos = elem
                                        self.is_complete = True

                                        # collect the final lattice
                                        self.lattice = []
                                        for i in range(3):
                                            self.lattice.append([float(x) for x in finalpos[0][0][i].text.strip().split()])

                                        # Collect the reciprocal lattice
                                        self.rec_lat = []
                                        for i in range(3):
                                            self.rec_lat.append([float(x) for x in finalpos[0][2][i].text.strip().split()])

                                        # collect the final basis
                                        self.basis = []
                                        for i,basis_atom in enumerate(finalpos[1]):
                                            self.basis.append([float(x) for x in basis_atom.text.strip().split()])
                                        finalpos.clear()

                else:
                        pass
        elem.clear()
        f.close()

