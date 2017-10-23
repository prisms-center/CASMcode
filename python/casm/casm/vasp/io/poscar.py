from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import numpy as np
import math

class PoscarError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Site:
    """ Site in a basis.

        Contains:
            self.cart = True if cartesian coordinate
            self.SD_FLAG = string for selective dynamics description, empty string by default (ex. "T F F").
            self.occupant = CASM specie name, empty string by default
            self.occ_alias = alias (POTCAR) name, empty string by default
            self.position = np.array coordinate
    """
    def __init__(self, cart, position, SD_FLAG = "", occupant = "", occ_alias = ""):
        """ Site constructor """
        self.cart = cart
        self.SD_FLAG = SD_FLAG
        self.occupant = occupant
        self.occ_alias = occ_alias
        if not isinstance(position, np.ndarray):
            raise PoscarError("Attempting to construct a Site and 'position' is not a numpy ndarray")
        self.position = position


    def write(self, file):
        """ Write this Site to a file """
        np.savetxt(file,self.position,fmt='%.8f',newline = '    ')
        file.write(self.SD_FLAG)
        file.write(self.occupant)
        file.write('\n')


class Poscar:
    """
    The Poscar class contains:
        self.header: Contains POSCAR header line with whitespace stripped from beginning and end
        self.scaling: lattice scale
        self._lattice: the lattice; lattice vectors are stored as rows in a numpy array
        self._reciprocal_lattice: the reciprocal lattice; lattice vectors are stored as rows in a numpy array
        self.SD_FLAG: True or False, Selective Dynamics Flag
        self.type_atoms: lists the specie names as in the POS (ex. [Mn3 Mn4])
        self.type_atoms_alias: lists the POTCAR names for the species
        self.num_atoms: lists the atoms as in the POS (ex. [1 1])
        self.coord_mode: Contains the coordinate mode text from POSCAR, with whitespace stripped from beginning and end
        self.basis: a list of Site objects

    """
    def __init__(self, filename, species=None, legacy_support=True):
        """ Construct a Poscar object from 'filename'

            Args:
                filename = POS/POSCAR file to read
                species = (default None) If given a Species dict, it is used to set self.type_atoms_alias for determining which POTCARs to use
        """
        self.read(filename, species, legacy_support)


    def read(self, filename, species=None, legacy_support=True):
        """ Reads a POS/POSCAR from 'filename'"""
        self.basis = []
        self._lattice = np.zeros((3,3))
        self._reciprocal_lattice = np.zeros((3,3))
        self.coord_mode = ''
        self.SD_FLAG = False
        self.num_atoms = []
        self.type_atoms = []
        self.type_atoms_alias = list(self.type_atoms)
        self.header = ""
        self.scaling = 0
        try:
            file = open(filename,'r')
        except IOError:
            raise PoscarError("Could not read file: " + filename)
        self.header = file.readline().strip()
        try:
            self._read_lattice(file)
            self._read_atominfo(file)
            self._read_flags(file)
            self._read_basis(file)
        except PoscarError as e:
            if legacy_support == False:
                raise e
            try:
                file.seek(0)
                file.readline()
                self._read_lattice(file)
                self._read_atominfo_legacy(file)
                self._read_flags(file)
                self._read_basis_legacy(file)
            except PoscarError as e:
                raise e
        file.close()

        # set type alias
        if species != None:
            self.update(species)


    def write(self, filename, sort=True):
        """ Write Poscar to 'filename'.

            If 'sort' = True, sort basis sites to minimize number of POTCARs
        """
        try:
            file = open(filename,'w')
        except IOError:
            raise PoscarError("Could not write: " + filename)
        file.write(str(self.header)+'\n')
        file.write(str(self.scaling)+'\n')
        np.savetxt(file,self._lattice, newline = '\n', fmt='%.8f' )

        if sort:
            type_line = ''
            num_line = ''
            d = self.basis_dict()
            for atom in sorted(d.keys()):
                type_line = type_line+'   '+atom
                num_line = num_line + '   ' +str(len(d[atom]))
            if(self.SD_FLAG):
                file.write('Selective Dynamics\n')
            file.write(type_line+'\n')
            file.write(num_line+'\n')
            file.write(self.coord_mode+'\n')
            for atom in sorted(d.keys()):
                for tb in d[atom]:
                    tb.write(file)
        else:
            file.write(' '.join(self.type_atoms_alias) + '\n')
            file.write(' '.join([str(x) for x in self.num_atoms]) + '\n')
            if(self.SD_FLAG):
                file.write('Selective Dynamics\n')
            file.write(self.coord_mode+'\n')
            for s in self.basis:
                s.write(file)

        file.close()


    def lattice(self, index = None):
        """ Returns the lattice, or lattice vector 'index', as numpy array"""
        if index != None:
            return self._lattice[index,:]
        return self._lattice


    def reciprocal_lattice(self, index = None):
        """ Returns the reciprocal lattice vector 'index', as numpy array"""
        if index != None:
            return self._reciprocal_lattice[index,:]
        return self._reciprocal_lattice


    def volume(self, lattice=None):
        """ Returns scalar triple product of lattice vectors """
        if lattice is None:
            lattice = self._lattice
        return np.inner(lattice[0,:], np.cross(lattice[1,:],lattice[2,:]))


    def reciprocal_volume(self, reciprocal_lattice=None):
        """ Returns scalar triple product of reciprocal lattice vector """
        if reciprocal_lattice is None:
            reciprocal_lattice = self._reciprocal_lattice
        return self.volume(reciprocal_lattice)


    def basis_dict(self):
        """ Return a dictionary where keys are unique specie 'alias' and values are lists of atoms."""
        basis_dict = dict()
        for i,atom in enumerate(self.type_atoms_alias):
            start = sum(self.num_atoms[0:(i)])
            end = start + self.num_atoms[i]
            if atom not in basis_dict.keys():
                basis_dict[atom] = self.basis[start:end]
            else:
                basis_dict[atom] += self.basis[start:end]
        return basis_dict


    def unsort_dict(self):
        """ Return a dict to unsort atom order.

            Returns 'unsort_dict', for which: unsorted_dict[orig_index] == sorted_index;

            unsorted_dict[sorted_index] == orig_index

            For example:
              'unsort_dict[0]' returns the index into the unsorted POSCAR of the first atom in the sorted POSCAR
        """

        # create basis_dict, but with initial position in POSCAR instead of coordinate
        basis_dict = dict()
        for i,atom in enumerate(self.type_atoms_alias):
            start = sum(self.num_atoms[0:(i)])
            end = start + self.num_atoms[i]
            if atom not in basis_dict.keys():
                basis_dict[atom] = range(start,end)
            else:
                basis_dict[atom] += range(start,end)

        print(basis_dict)

        orig_pos = []
        for atom in sorted(basis_dict.keys()):
            orig_pos += basis_dict[atom]

        print(orig_pos)

        new_pos = range(0,len(self.basis))

        return dict(zip(new_pos, orig_pos))


    def _read_lattice(self,file):
        """ Called by self.read() to read the lattice into self._lattice

            Args:
                file: an open POSCAR being read from

            self._lattice contains lattice vectors stored as rows in a numpy array (easy inversion)
        """
        try:
            line = file.readline().strip()
            self.scaling = float(line)
        except ValueError:
            raise PoscarError("Could not read lattice scaling: '" + line + "'")
        lat=[]
        for i in range(3):
            line = file.readline()
            try:
                lat.append([float(x) for x in line.split()])
            except ValueError:
                raise PoscarError("Could not read lattice vector: '" + line + "'")
        self._lattice = self.scaling*np.array(lat)
        if self._lattice.shape!=(3,3):
            raise PoscarError("Lattice shape error: " + np.array_str(self._lattice))
        self._reciprocal_lattice = 2.0*math.pi*np.linalg.inv(np.transpose(self._lattice))
        self.scaling=1.0
        return


    def _read_atominfo(self, file):
        """ Called by self.read() to read the 2 lines that correspond to the num_atoms and type_atoms
            in a VASP5 style POSCAR

            Args:
                file: an open POSCAR being read from

            self.type_atoms is a list of strings corresponding to atom names
            self.num_atoms is a list of int corresponding to the number of each atom type
            self.type_atoms_alias is set from self.type_atoms
        """
        typeline = file.readline().strip()
        numline = file.readline().strip()
        self.type_atoms = typeline.split()
        if len(self.type_atoms)==0:
            raise PoscarError("No atom names found: '" + typeline + "'")
        try:
            self.num_atoms = [int(n) for n in numline.split()]
        except ValueError:
            raise PoscarError("Could not read number of each atom type: '" + numline + "'")
        if len(self.num_atoms)!=len(self.type_atoms):
            raise PoscarError("Atom type and number lists are not the same length: \n'" + typeline + "'\n" + "'" + numline + "'")
        self.type_atoms_alias = list(self.type_atoms)


    def _read_atominfo_legacy(self, file):
        """ Called by self.read() to read the lines that correspond to the num_atoms and type_atoms
            in a VASP4 style POSCAR

            Args:
                file: an open POSCAR being read from

            self.type_atoms is a list of strings corresponding to atom names
            self.num_atoms is a list of int corresponding to the number of each atom type
            self.type_atoms_alias is set from self.type_atoms
        """
        numline = file.readline().strip()
        try:
            self.num_atoms = [int(n) for n in numline.split()]
        except ValueError:
            raise PoscarError("Could not read number of each atom type: '" + numline + "'")
        self.type_atoms = map(int, range(len(self.num_atoms)))
        self.type_atoms_alias = list(self.type_atoms)


    def _read_flags(self,file):
        """ Called by self.read() to read POSCAR flags: Selective Dynamics, and Coord Mode

            Args:
                file: an open POSCAR being read from

            self.SD_FLAG = True/False if selective dynamics on
            self.coord_mode = Contains the coordinate mode text, with whitespace stripped from beginning and end
        """
        line = file.readline().strip()
        if len(line)==0:
            raise PoscarError("Could not read Select Dynamics or Coord Mode line")
        if(line[0].lower()=='s'):
            self.SD_FLAG = True
            self.coord_mode = file.readline()
            if len(line)==0:
                raise PoscarError("Read Select Dynamics, but could not read Coord Mode line")
        else:
            self.SD_FLAG = False
            self.coord_mode = line

        if self.coord_mode[0].lower() not in ['c','d','k']:
            raise PoscarError("Read invalid coord mode: '" + self.coord_mode + "'")


    def _read_basis(self,file):
        """ Called by self.read() to read basis of POSCAR file into self.basis

            Args:
                file: an open POSCAR being read from

            self.basis is a list of Site objects
        """
        self.basis = []
        cart = not(self.coord_mode[0].lower()=='d')
        for i in range(len(self.num_atoms)):
            for j in range(self.num_atoms[i]):
                line = file.readline().strip()
                if len(line)==0:
                    raise PoscarError("Error reading basis: insufficient number of atoms")
                words = line.split()
                pos = []
                if (self.SD_FLAG):
                    if len(words) < 6:
                        raise PoscarError("Error reading basis: insufficient number of SD tags")
                    SD_FLAG = words[3]+' '+words[4]+' '+words[5]
                else:
                    SD_FLAG = ''
                try:
                    pos = [float(x) for x in words[0:3]]
                except ValueError:
                    raise PoscarError("Error reading basis coordinate: '" + line + "'")
                self.basis.append(Site(cart, np.array(pos), SD_FLAG, self.type_atoms[i], self.type_atoms[i]))


    def _read_basis_legacy(self,file):
        """ Called by self.read() to read basis of VASP4  POSCAR file into self.basis

            Args:
                file: an open POSCAR being read from

            self.basis is a list of Site objects
        """
        self.basis = []
        cart = not(self.coord_mode[0].lower()=='d')
        for i in range(len(self.num_atoms)):
            for j in range(self.num_atoms[i]):
                line = file.readline().strip()
                if len(line)==0:
                    raise PoscarError("Error reading basis: insufficient number of atoms")
                words = line.split()
                pos = []
                if (self.SD_FLAG):
                    if len(words) < 6:
                        raise PoscarError("Error reading basis: insufficient number of SD tags")
                    SD_FLAG = words[3]+' '+words[4]+' '+words[5]
                else:
                    SD_FLAG = ''
                try:
                    pos = [float(x) for x in words[0:3]]
                except ValueError:
                    raise PoscarError("Error reading basis coordinate: '" + line + "'")
                if len(words) == 7 or len(words) == 4:
                    atom_type = words[-1]
                    if j > 0:
                        if self.type_atoms[i] != atom_type:
                            raise PoscarError("Error reading basis: mismatch between type count in header and types in basis")
                    else:
                        self.type_atoms[i] = atom_type
                        self.type_atoms_alias[i] = atom_type
                else:
                    atom_type = self.type_atoms[i]
                self.basis.append(Site(cart, np.array(pos), SD_FLAG, atom_type, atom_type))



    def update(self, species):
        """ Set self.type_atoms_alias and self.basis[x].alias according to Species dict.

            Warning: Not sure we want to check equivalence regardless of atom name case
        """
        for i,atom in enumerate(self.type_atoms):
            for s in species:
                if s == atom:
                    self.type_atoms_alias[i] = species[s].alias
        index = 0
        for i in range(len(self.type_atoms_alias)):
            for j in range(self.num_atoms[i]):
                self.basis[index].occ_alias = self.type_atoms_alias[i]
                index += 1

        if 'MAGMOM' in species.tags:
            for base in self.basis:
                if base.occupant in species.keys():
                    base.mag = float(species[base.occupant].tags['MAGMOM'])

        return



