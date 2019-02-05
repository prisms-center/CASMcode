import os
import numpy as np

from casm.vasp.io.poscar import Poscar, Site
from casm.project import DirectoryStructure, ProjectSettings
import casm.aimswrapper


class GeometryError(Exception):
    def __init__(self,  msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Geometry:
    """
    The Geometry class contains:
        self._lattice: the lattice; lattice vectors are stored as rows in a numpy array
        self._reciprocal_lattice: the reciprocal lattice; lattice vectors are stored as rows in a numpy array
        self.sel_dyn: True or False, Selective Dynamics Flag
        self.type_atoms: lists the specie names as in the POS (ex. [Mn3 Mn4])
        self.type_atoms_alias: lists the POTCAR names for the species
        self.num_atoms: lists the atoms as in the POS (ex. [1 1])
        self.coord_mode: Contains the coordinate mode text from POSCAR, with whitespace stripped from beginning and end
        self.basis: a list of Site objects

    """
    def __init__(self, filename):
        """ Construct a Geometry object from 'filename'

            Args:
                filename = file to read
        """
        self.basis = []
        self._lattice = np.zeros((3, 3))
        self._reciprocal_lattice = np.zeros((3, 3))
        self.coord_mode = ''
        self.sel_dyn = False
        self.SD_FLAG = False
        self.num_atoms = []
        self.type_atoms = []
        self.type_atoms_alias = list(self.type_atoms)

        try:
            file = open(filename, 'rb')
        except IOError:
            raise GeometryError("Could not read file: " + filename)

        if filename.split('/')[-1] == 'POS':
            self.read_pos(filename)
        else:
            self._read_geo(file)

        file.close()

    @staticmethod
    def read_casm_settings():
        configdir = os.getcwd()
        _res = os.path.split(configdir)
        cfgname = os.path.split(_res[0])[1] + "/" + _res[1]
        casm_dirs = DirectoryStructure(configdir)
        casm_sets = ProjectSettings(configdir)
        clex = casm_sets.default_clex
        setfile = casm_dirs.settings_path_crawl("relax.json", cfgname, clex)
        settings = casm.aimswrapper.read_settings(setfile)
        return settings

    def write(self, filename, species_data):
        """ Write geometry.in to 'filename'.
        """
        try:
            file = open(filename, 'w')
        except IOError:
            raise GeometryError("Could not write: " + filename)

        settings = self.read_casm_settings()
    
        file.write('#auto-gemerated geometry.in by CASM\n')
        for i in range(3):
            file.write("lattice_vector %.8f %.8f %.8f\n" %
                       (self._lattice[i, 0], self._lattice[i, 1], self._lattice[i, 2]))

        if self.coord_mode[0].lower() == 'd':
            for s in self.basis:
                file.write('atom_frac %8.8f %8.8f %8.8f %s\n' % 
                           (s.position[0], s.position[1], s.position[2], s.occupant))
                file.write('  initial_moment %3.3f\n' % float(species_data[s.occupant].init_mom))
                if settings['is_slab']:
                    if s.position[2] < float(settings['fix_pos']) / self._lattice[2, 2]:
                        file.write('  constrain_relaxation .true.\n')
        else:
            for s in self.basis:
                file.write('atom_frac %8.8f %8.8f %8.8f %s\n' % 
                           (s.position[0], s.position[1], s.position[2], s.occupant))
                file.write('  initial_moment %3.3f\n' % float(species_data[s.occupant].init_mom))
                if settings['is_slab']:
                    if s.position[2] < float(settings['fix_pos']):
                        file.write('  constrain_relaxation .true.\n')

    def lattice(self, index=None):
        """ Returns the lattice, or lattice vector 'index', as numpy array"""
        if index is not None:
            return self._lattice[index, :]
        return self._lattice

    def reciprocal_lattice(self, index=None):
        """ Returns the reciprocal lattice vector 'index', as numpy array"""
        if index is not None:
            return self._reciprocal_lattice[index, :]
        return self._reciprocal_lattice

    def volume(self, lattice=None):
        """ Returns scalar triple product of lattice vectors """
        if lattice is None:
            lattice = self._lattice
        return np.inner(lattice[0, :], np.cross(lattice[1, :], lattice[2, :]))

    def reciprocal_volume(self, reciprocal_lattice=None):
        """ Returns scalar triple product of reciprocal lattice vector """
        if reciprocal_lattice is None:
            reciprocal_lattice = self._reciprocal_lattice
        return self.volume(reciprocal_lattice)

    def basis_dict(self):
        """ Return a dictionary where keys are unique specie 'alias' and values are lists of atoms."""
        basis_dict = dict()
        for i, atom in enumerate(self.type_atoms_alias):
            start = sum(self.num_atoms[0:i])
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
        for i, atom in enumerate(self.type_atoms_alias):
            start = sum(self.num_atoms[0:i])
            end = start + self.num_atoms[i]
            if atom not in basis_dict.keys():
                basis_dict[atom] = range(start, end)
            else:
                basis_dict[atom] += range(start, end)

        orig_pos = []
        for atom in sorted(basis_dict.keys()):
            orig_pos += basis_dict[atom]

        new_pos = range(0, len(self.basis))

        return dict(zip(new_pos, orig_pos))

    def _read_geo(self, file):
        """ Called by __init__ to read the lattice into self._lattice
            Args:
                file: an open geometry.skel being read from

            self._lattice contains lattice vectors stored as rows in a numpy array (easy inversion)
        """
        atom_read = False
        atom_name = None
        cart = None
        lat = []
        cont = file.readlines()
        for line in cont:
            pos = np.empty(3)
            if b'lattice_vector' in line:
                lat.append([float(x) for x in line.split()[1:4]])

            if b'atom ' in line:
                self.coord_mode = 'cartesian'
                cart = True
                atom_read = True
                word = line.split()
                try:
                    pos[0] = float(word[1])
                    pos[1] = float(word[2])
                    pos[2] = float(word[3])
                    atom_name = word[4].decode('UTF-8')
                except ValueError:
                    raise GeometryError("Error reading basis coordinate: '" + line + "'")

            if b'atom_frac ' in line:
                self.coord_mode = 'direct'
                cart = False
                atom_read = True
                word = line.split()
                try:
                    pos[0] = float(word[1])
                    pos[1] = float(word[2])
                    pos[2] = float(word[3])
                    atom_name = word[4].decode('UTF-8')
                except ValueError:
                    raise GeometryError("Error reading basis coordinate: '" + line + "'")

            if atom_read:
                sd_flags = 'T T T'
                self.basis.append(Site(cart=cart, position=np.array(pos),
                                       SD_FLAG=sd_flags, occupant=atom_name, occ_alias=atom_name))
                atom_read = False

#            ln += 1

        # done reading, analyze now
        all_names = []
        for a in self.basis:
            all_names.append(a.occupant)

        tmp = []
        symbols = [all_names[0]]
        for i in range(len(all_names)):
            tmp = []
            for j in range(len(symbols)):
                tmp.append(symbols[j])
            if tmp.count(all_names[i]) == 0:
                symbols.append(all_names[i])

        nums = {}
        for i in range(len(symbols)):
            k = 0
            for j in range(len(all_names)):
                if symbols[i] == all_names[j]:
                    k += 1
            nums[symbols[i]] = k

        self.type_atoms = symbols
        if len(self.type_atoms) == 0:
            raise GeometryError("No atom names found")
        try:
            self.num_atoms = [nums[sym] for sym in symbols]
        except ValueError:
            raise GeometryError("Could not read number of each atom type")
        if len(self.num_atoms) != len(self.type_atoms):
            raise GeometryError("Atom type and number lists are not the same length")
        self.type_atoms_alias = list(self.type_atoms)

        self._lattice = np.array(lat)
        if self._lattice.shape != (3, 3):
            raise GeometryError("Lattice shape error: " + np.array_str(self._lattice))
        self._reciprocal_lattice = 2 * np.pi * np.linalg.inv(np.transpose(self._lattice))
        self.scaling = 1.0

        del tmp, symbols, nums

    def check_constraints(self, cont, ln, total_lines):
        flags = ['T', 'T', 'T']

        if ln >= total_lines:
            return flags[0] + ' ' + flags[1] + ' ' + flags[2]

        next_line = cont[ln]
        k = 0

        if b'initial_moment' in next_line:
            raise ValueError('Initial Moments will be added automatically, please remove from geometry.skel.')

        if b'constrain_relaxation ' in next_line:
            limit = np.min([3, total_lines - ln])
            self.sel_dyn = True
            check_lines = cont[ln:ln+limit]
            check = check_lines[k].strip()
            while b'constrain_relaxation ' in check:
                try:
                    word = check.split()
                    if word[1] == b'.true.':
                        flags = ['F', 'F', 'F']
                    if word[1].lower() == b'x':
                        flags[0] = 'F'
                    if word[1].lower() == b'y':
                        flags[1] = 'F'
                    if word[1].lower() == b'z':
                        flags[2] = 'F'
                except ValueError:
                    raise GeometryError("Error reading constrints from line '" + check + "'")

                k += 1

                if k == len(check_lines):
                    break
                check = check_lines[k].strip()

        return flags[0] + ' ' + flags[1] + ' ' + flags[2]

    def read_pos(self, file):
        self.basis = Poscar(file).basis
        self._lattice = Poscar(file)._lattice
        self._reciprocal_lattice = Poscar(file).reciprocal_lattice()
        self.coord_mode = Poscar(file).coord_mode
        self.SD_FLAG = Poscar(file).SD_FLAG
        self.num_atoms = Poscar(file).num_atoms
        self.type_atoms = Poscar(file).type_atoms
        self.type_atoms_alias = Poscar(file).type_atoms_alias
