""" Object and functions for handling atom geometries """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

#pylint: disable=too-many-lines

import re
import os

from ..species import species_settings
from ..seq_exceptions import SiteError, CellError, GeomError, POSError

class Site(object): #pylint: disable=too-few-public-methods
    """ Site in a basis.

        Contains:
            self.occupant = CASM specie name, empty string by default
            self.occ_alias = alias (atom file) name, empty string by default
            self.position = vec of float
            self.charge = charge at this coordinate
    """
    def __init__(self, position, occupant="", occ_alias="", charge=None):
        """ Site constructor """
        self.occupant = occupant
        self.occ_alias = occ_alias
        self.charge = charge
        self.position = position
        try:
            _, _, _ = [float(x) for x in self.position]
        except:
            raise SiteError("Position %s could not be cast to a triplet of floats", self.position)

    def write(self, stream=None, index=None):
        """ Write this Site to a file or string """
        arg_string = ""
        if index is not None:
            arg_string += (self.occupant + "_" + str(index) + "    " + self.occ_alias + "    " +
                           ("  ").join([str(x) for x in self.position]))
            # stream.write(self.occupant + "_" + str(index) + "    " +
            #              self.occ_alias + "    " +
            #              ("  ").join([str(x) for x in self.position]))
        else:
            arg_string += (("  ").join([str(x) for x in self.position]) +
                           "    " +  self.occupant)
            # stream.write(("  ").join([str(x) for x in self.position]) +
            #              "    " +  self.occupant)
        arg_string += "\n"
        if stream is not None:
            stream.write(arg_string)
        return arg_string

    def write_charge(self, stream, total_charge=None):
        """ Write a charge at a site """
        if total_charge is None:
            if self.charge is not None:
                total_charge = float(self.charge)
            else:
                raise GeomError("No charge on site and no charge supplied, can't print charge!")
        else:
            total_charge = float(total_charge)
        stream.write("charge\n")
        stream.write("%g\n" % total_charge)
        stream.write("location of charge\n")
        stream.write(("  ").join([str(x) for x in self.position]) + "   " + self.occupant)

class Cell(object): #pylint: disable=too-many-instance-attributes
    """ Container  class for cell data

        Contains:
            self._lattice (matrix of float): primative lattice vectors
            self.coord_mode (str): LATTICE or CARTESIAN; value of "coordinate" keyword in lcao.in,
                default LATTICE
            self._scale (float): scale factor for cell, default 1.0
                default 1.0
    """

    @classmethod
    def POS(cls, POS):    #pylint: disable=invalid-name
        """ Shorthand for POS construction """
        return cls(POS=POS)

    @classmethod
    def seq(cls, seqfile):
        """ Shorthand for SeqQuest construction """
        return cls(seqfile=seqfile)

    @classmethod
    def mini(cls, mini_POS):    #pylint: disable=invalid-name
        """ Construction from partial POS """
        return cls(POS=mini_POS, mini=True)

    @property
    def lattice(self):
        """ Returns a copy of the lattice """
        return [x[:] for x in self._lattice[:]]

    @property
    def reciprocal_lattice(self):
        """ Returns a copy of the reciprocal lattice """
        return [x[:] for x in self._reciprocal_lattice[:]]

    def __init__(self, POS=None, seqfile=None, mini=False):
        """ Construct a Cell object from 'cell' file

            Args:
                POS (string): Path to a POS file
                seqfile (string or buffer): Path to a lcao.in, lcao.out, or lcao.geom file OR an
                    open file or buffer object. Info suppliments and overrides POS file """
        ### Initialize ###
        # Set to default values, if such things exist
        self._lattice = [[0 for _ in range(3)] for _ in range(3)]
        self._reciprocal_lattice = [[0 for _ in range(3)] for _ in range(3)]
        self._scale = 1.0

        if mini:
            if POS is None:
                raise POSError("Need a (mini)POS for mini mode")
            self.read_mini_POS(POS)
        else:
            if POS is not None:
                self.read_POS(POS)
        if seqfile is not None:
            self.read_seq(seqfile)

        if seqfile is None and POS is None:
            raise GeomError("Neither POS nor lcao.* provided, aborting!")

        ### Reciprocal lattice ###
        r_1 = _div(_cross(self._lattice[1], self._lattice[2]), _dot(self._lattice[0],
                                                                    _cross(self._lattice[1],
                                                                           self._lattice[2])))
        r_2 = _div(_cross(self._lattice[2], self._lattice[0]), _dot(self._lattice[1],
                                                                    _cross(self._lattice[2],
                                                                           self._lattice[0])))
        r_3 = _div(_cross(self._lattice[0], self._lattice[1]), _dot(self._lattice[2],
                                                                    _cross(self._lattice[0],
                                                                           self._lattice[1])))
        self._reciprocal_lattice = [r_1[:], r_2, r_3]

        # ### Ensuring 'scale' is always 1.0 ###
        # self.scale(scale=self._scale)
        # self._scale = 1.0

    def read_POS(self, POS):   #pylint: disable=invalid-name
        """ Read a POS and construct lattice coords in a Cell object """
        with open(POS) as stream:
            _ = [stream.readline() for _ in range(1)]
            self._read_lattice(stream)

    def read_mini_POS(self, mini_POS):   #pylint: disable=invalid-name
        """ Read a POS and construct lattice coords in a Cell object """
        with open(mini_POS) as stream:
            self._read_lattice(stream)

    def read_seq(self, seqfile):
        """ Read a seqquest file (or a buffer containing info from a file) and construct or amend
                a Cell object """
        # Check if we're a stream
        if hasattr(seqfile, 'read'):
            self._read_seq_stream(seqfile)
        else:
            self._read_seq_file(seqfile)

    def _read_seq_file(self, seqfile):    #pylint: disable=too-many-branches, too-many-statements, too-many-locals
        """ Read a seqquest file (or a buffer containing info from a file) and construct or amend
                a Cell object """
        ### Initialize ###
        lattice = [[None for _ in range(3)] for _ in range(3)]
        scale = None
        inner = 0
        with open(seqfile) as stream:
            while True:
                line = stream.readline()
                if line == "":
                    break
            # for line in stream:
                if re.search("cellfile", line):
                    inner = self._read_seq_file(stream.readline().strip())
                    break
                if re.search(r"scale\s+", line):
                    scale = float(stream.readline().strip())
                if re.search("primitive", line):
                    for i in range(3):
                        try:
                            line = stream.readline()
                            lattice[i][0], lattice[i][1], lattice[i][2] = [
                                float(x) for x in line.strip().split()]
                        except:
                            raise CellError("Error reading lattice vector %r" % line)

        ### Check validity ###
        if inner == 0:
            if self._lattice[0][0] is not None:
                # copy/deepcopy would preserve internal structure; we want an unshared true copy
                for i in range(3):
                    self._lattice[i] = lattice[i][:]
            else:
                raise CellError("Something went wrong, %s could not be read" % seqfile)

            # 'scale' keyword is optional
            if self._scale is not None:
                self._scale = scale

        return 1

    def _read_seq_stream(self, seqfile):    #pylint: disable=too-many-branches, too-many-statements
        """ Read a seqquest file (or a buffer containing info from a file) and construct or amend
                a Cell object """
        ### Initialize ###
        lattice = [[None for _ in range(3)] for _ in range(3)]
        scale = 1.0
        inner = 0
        # Check if we're a stream
        current_pos = seqfile.tell()
        while True:
            line = seqfile.readline()
            if line == "":
                break
        # for line in seqfile:
            if re.search("cellfile", line):
                inner = self._read_seq_file(seqfile.readline().strip())
                break
            elif re.search("grid dimensions", line):
                break
            elif re.search(r"\s*>", line):
                break
            elif re.search(r"scale\s+", line):
                scale = float(seqfile.readline().strip())
            elif re.search("primitive", line):
                pass
            # We were passed the stream with 'primitive' already consumed!
            else:
                i = 0
                lattice[i][0], lattice[i][1], lattice[i][2] = [
                    float(x) for x in line.strip().split()]
                for i in range(1, 3):
                    try:
                        line = seqfile.readline()
                        lattice[i][0], lattice[i][1], lattice[i][2] = [
                            float(x) for x in line.strip().split()]
                    except:
                        raise CellError("Error reading lattice vector %r" % line)
            current_pos = seqfile.tell()

        # Rewind the stream before passing control back
        seqfile.seek(current_pos)

        ### Check validity ###
        if inner == 0:
            if self._lattice[0][0] is not None:
                # copy/deepcopy would preserve internal structure; we want an unshared true copy
                for i in range(3):
                    self._lattice[i] = lattice[i][:]
            else:
                raise CellError("Something went wrong, %s could not be read" % seqfile)

            # 'scale' keyword is optional
            if self._scale is not None:
                self._scale = scale

        return 1

    def _read_lattice(self, stream):
        """ Called by self.read_POS() to read a lattice into self._lattice

            Args:
                stream (buffer): an open POS being read from

            self._lattice contains lattice vectors stored as rows in a numpy array (easy inversion)
        """
        try:
            line = stream.readline().strip()
            self._scale = float(line)
        except ValueError:
            raise POSError("Could not read lattice scaling: '" + line + "'")
        lat = []
        for _ in range(3):
            line = stream.readline()
            try:
                lat.append([float(x) for x in line.split()])
            except ValueError:
                raise POSError("Could not read lattice vector: '" + line + "'")
        self._lattice = lat
        if len(self._lattice) != 3:
            raise POSError("Lattice shape error: " + self._lattice)
        for line in self._lattice:
            if len(line) != 3:
                raise POSError("Lattice shape error: " + self._lattice)

    def write_POS(self, filename):  #pylint: disable=invalid-name
        """ Write POS to filename """
        # If a "Geom" has already written here:
        if os.path.isfile(filename):
            mini_geom = Geom.mini(filename)
            write_POS(self, mini_geom, filename)
        else:
            with open(filename, "w") as stream:
                stream.write("  %.8f" % (self._scale) + "\n")
                for i in range(3):
                    stream.write("     " + ("     ").join(["%.8f" % x
                                                           for x in self._lattice[i]]) + "\n")

    def get_lattice(self, index=None):
        """ Returns the lattice, or lattice vector 'index', as listed list"""
        if index != None:
            return self._lattice[index][:]
        else:
            return [x[:] for x in self.lattice[:]]

    def get_reciprocal_lattice(self, index=None):
        """ Returns the reciprocal lattice vector 'index', as numpy array"""
        if index != None:
            return self._reciprocal_lattice[index][:]
        return self._reciprocal_lattice

    def volume(self):
        """ Returns scalar triple product of lattice vectors """
        return volume(self._lattice)
        # return abs(_dot(self._lattice[0], _cross(self._lattice[1], self._lattice[2])))

    def reciprocal_volume(self):
        """ Returns scalar triple product of reciprocal lattice vector """
        return volume(self._reciprocal_lattice)
        # return abs(_dot(self._reciprocal_lattice[0], _cross(self._reciprocal_lattice[1],
        #                                                     self._reciprocal_lattice[2])))

    def scale(self, scale=1.0, scalex=1.0, scaley=1.0, scalez=1.0, ndim=3):    #pylint: disable=too-many-arguments
        """ Applies scale settings """
        # Point structures don't get any cell dimensions scaled!
        if ndim == 0:
            pass
        elif ndim == 1 or ndim == 2:
            raise CellError("Handling of ndim != 0 or 3 is not yet implemented!")
        elif ndim == 3:
            for i in range(3):
                for j in range(3):
                    self._lattice[i][j] *= scale * [scalex, scaley, scalez][i]
        else:
            raise CellError("Invalid ndim value %r encountered!" % ndim)
        # Now our scale is 1.0
        self._scale = 1.0

    def scalex(self, scalex=1.0):
        """ Alias for global scale """
        self.scale(scalex=scalex)

    def scaley(self, scaley=1.0):
        """ Alias for global scale """
        self.scale(scaley=scaley)

    def scalez(self, scalez=1.0):
        """ Alias for global scale """
        self.scale(scalez=scalez)

class Geom(object): #pylint: disable=too-many-instance-attributes
    """ Container class for geometry data

        Contains:
            self.natoms (int): no of atoms in the unit cell
            self.basis (list of Sites): atoms in the crystal/molecule/etc
            self.type_atoms (list of str): lists the specie names as in the POS (ex. [Mn3 Mn4])
            self.type_atoms_alias (list of str): lists the *atm file names for the species
            self.num_atoms (list of int): lists the atoms as in the POS (ex. [1 1])
            *self.charge (float): net charge of the geom, if any
            *self.charge_loc (vec of float): location of the point charge, if any
            self.coord_mode (string): coord mode of POS/geom
    """

    @classmethod
    def POS(cls, POS, species=None):    #pylint: disable=invalid-name
        """ Shorthand for POS construction """
        return cls(POS=POS, species=species)

    @classmethod
    def seq(cls, seqfile, species=None):
        """ Shorthand for SeqQuest construction """
        return cls(seqfile=seqfile, species=species)

    @classmethod
    def geom(cls, seqfile, species=None):
        """ Shorthand for geom construction """
        return cls(seqfile=seqfile, species=species)

    @classmethod
    def mini(cls, mini_POS):    #pylint: disable=invalid-name
        """ Construction from partial POS """
        return cls(POS=mini_POS, mini=True)

    def __init__(self, POS=None, seqfile=None, species=None, mini=False):
        """ Construct a Geom object from 'geom' file

            Args:
                POS (string): Path to a POS file
                seqfile (string or buffer): Path to a lcao.in, lcao.out, or lcao.geom file OR an
                    open file or buffer object. Info suppliments and overrides POS file
                species (dict of IndividualSpecies): Species info for atoms in the
                    crystal/molecule/etc
        """
        ### Initialize ###
        # Set to default values, if such things exist
        self.natoms = 0
        self.basis = []
        # self.lattice = [[0 for _ in range(3)] for _ in range(3)]
        self.type_atoms = []
        self.type_atoms_alias = []
        self.num_atoms = []
        self.charge_loc = None
        self.charge = None
        self.title = ""
        self.coord_mode = None

        if mini:
            if POS is None:
                raise POSError("Need a (mini)POS for mini mode")
            self.read_mini_POS(POS)
        else:
            if POS is not None:
                self.read_POS(POS)
        if seqfile is not None:
            self.read_seq(seqfile)

        if seqfile is None and POS is None:
            raise GeomError("Neither POS nor lcao.* provided, aborting!")

        if species is not None:
            if isinstance(species, dict):
                self.update(species)
            else:
                species = species_settings(species)
                self.update(species)

        self.natoms = len(self.basis)

    def read_POS(self, POS):   #pylint: disable=invalid-name
        """ Read a POS and construct atom coords in a Geom object """
        with open(POS) as stream:
            self.title = stream.readline().strip()
            _ = [stream.readline() for _ in range(4)]
            # self._read_lattice(stream)
            self._read_atominfo(stream)
            self._read_basis(stream)

    def read_mini_POS(self, mini_POS):   #pylint: disable=invalid-name
        """ Read a POS and construct atom coords in a Geom object """
        with open(mini_POS) as stream:
            self.title = stream.readline().strip()
            # _ = [stream.readline() for _ in range(4)]
            # self._read_lattice(stream)
            self._read_atominfo(stream)
            self._read_basis(stream)

    def _read_atominfo(self, stream):
        """ Called by self.read_POS() to read the lines that correspond to the num_atoms and
            type_atoms in a CASM style POS

            Args:
                stream (buffer): an open POS being read from

            self.type_atoms is a list of strings corresponding to atom names
            self.num_atoms is a list of int corresponding to the number of each atom type
            self.type_atoms_alias is set from self.type_atoms
        """
        self.type_atoms = stream.readline().strip().split()
        numline = stream.readline().strip()
        try:
            self.num_atoms = [int(n) for n in numline.split()]
        except ValueError:
            raise POSError("Could not read number of each atom type: '" + numline + "'")
        # self.type_atoms = [int(x) for x in range(len(self.num_atoms))]
        # self.type_atoms = map(int, range(len(self.num_atoms)))
        self.type_atoms_alias = list(self.type_atoms)

    def _read_basis(self, stream):
        """ Called by self.read_POS() to read basis of POSCAR file into self.basis

            Args:
                stream (buffer): an open POSCAR being read from

            self.basis is a list of Site objects
        """
        self.basis = []
        coord_mode = stream.readline().strip()
        if coord_mode.lower() == "cartesian":
            self.coord_mode = "CARTESIAN"
        elif coord_mode.lower() == "direct":
            self.coord_mode = "LATTICE"
        else:
            raise POSError("Invalid coordinate mode %s found in POS!" % coord_mode)
        for i in range(len(self.num_atoms)):
            for _ in range(self.num_atoms[i]):
                line = stream.readline().strip()
                if len(line) == 0:
                    raise POSError("Error reading basis: insufficient number of atoms")
                words = line.split()
                pos = []
                try:
                    pos = [float(x) for x in words[0:3]]
                except ValueError:
                    raise POSError("Error reading basis coordinate: '" + line + "'")
                self.basis.append(Site(pos[:], self.type_atoms[i], self.type_atoms[i]))

    def read_seq(self, seqfile, title=None):
        """ Read a seqquest file (or a buffer containing info from a file) and construct or amend
                a Geom object
        """

        if title is not None:
            self.title = title
        # Check if we're a stream
        if hasattr(seqfile, 'read'):
            self._read_seq_stream(seqfile)
        else:
            self._read_seq_file(seqfile)

    def _read_seq_file(self, seqfile):    #pylint: disable=too-many-branches, too-many-statements
        """ Read a seqquest file (or a buffer containing info from a file) and construct or amend
                a Geom object """
        ### Initialize ###
        num_atoms = {}
        basis = []
        type_atoms = []
        type_atoms_alias = []
        pos = [0, 0, 0]
        inner = 0
        with open(seqfile) as stream:
            while True:
                line = stream.readline()
                if line == "":
                    break
            # for line in stream:
                if re.search("geomfile", line):
                    # inner = self._read_seq_file(stream.readline().strip())
                    inner = self._read_seq_file("lcao.geom_in")
                    break
                if re.search("atom, type, position", line):
                    self.title = line.strip().split(";")[-1]
                    while True:
                        line = stream.readline()
                        if line == "":
                            break
                    # for line in stream:
                        if re.search("kgrid", line):
                            break
                        try:
                            atom_id, atom_alias, pos[0], pos[1], pos[2] = line.split()
                        except:
                            raise GeomError("Error reading basis site: %r" % line)
                        atom_type = "_".join(atom_id.split("_")[:-1])
                        pos = [float(x) for x in pos]

                        basis.append(Site(pos[:], atom_type, atom_alias))
                        if not atom_type in type_atoms:
                            type_atoms += [atom_type]
                        if not atom_alias in type_atoms_alias:
                            type_atoms_alias += [atom_alias]
                        if atom_type in num_atoms.keys():
                            num_atoms[atom_type] += 1
                        else:
                            num_atoms[atom_type] = 1
                if re.search("kgrid", line):
                    break

        ### Check validity ###
        if inner == 0:
            if len(basis) == 0:
                raise GeomError("Something went wrong, %s could not be read" % seqfile)

        # Want a sorted num_atoms list
        self.num_atoms = [num_atoms[x] for x in type_atoms]
        self.basis = basis
        self.type_atoms = type_atoms
        self.type_atoms_alias = type_atoms_alias

        return 1

    def _read_seq_stream(self, seqfile):    #pylint: disable=too-many-branches, too-many-statements
        """ Read a seqquest file (or a buffer containing info from a file) and construct or amend
                a Geom object """
        ### Initialize ###
        num_atoms = {}
        basis = []
        type_atoms = []
        type_atoms_alias = []
        pos = [0, 0, 0]
        inner = 0

        current_pos = seqfile.tell()
        while True:
            line = seqfile.readline()
            if line == "":
                break
        # for line in seqfile:
            if re.search("geomfile", line):
                # inner = self._read_seq_file(seqfile.readline().strip())
                inner = self._read_seq_file("lcao.geom_in")
                break
            if re.search("kgrid", line):
                seqfile.seek(current_pos)
                break
            else:
                if re.search("atom, type, position", line):
                    self.title = line.strip().split(";")[-1]
                else:
                    try:
                        atom_id, atom_alias, pos[0], pos[1], pos[2] = line.split()
                    except:
                        raise GeomError("Error reading basis site: %r" % line)
                    atom_type = atom_id.split("_")[0]
                    pos = [float(x) for x in pos]

                    basis.append(Site(pos[:], atom_type, atom_alias))
                    if not atom_type in type_atoms:
                        type_atoms += [atom_type]
                    if not atom_alias in type_atoms_alias:
                        type_atoms_alias += [atom_alias]
                    if atom_type in num_atoms.keys():
                        num_atoms[atom_type] += 1
                    else:
                        num_atoms[atom_type] = 1
                # Rewind the stream before passing control back
                current_pos = seqfile.tell()

        ### Check validity ###
        if inner == 0:
            if len(basis) == 0:
                raise GeomError("Something went wrong, %s could not be read" % seqfile)

        # Want a sorted num_atoms list
        self.num_atoms = [num_atoms[x] for x in type_atoms]
        self.basis = basis
        self.type_atoms = type_atoms
        self.type_atoms_alias = type_atoms_alias
        self.natoms = len(self.basis)

        return 1

    def update(self, species):
        """ Set self.type_atoms_alias and self.basis[x].alias according to Species dict.

            Warning: Not sure we want to check equivalence regardless of atom name case
        """
        for i, atom in enumerate(self.type_atoms):
            for s in species:   #pylint: disable=invalid-name
                if s == atom:
                    self.type_atoms_alias[i] = species[s].alias
        index = 0
        charge_loc = False
        self.charge = None  # We need to RESET charge before accumulating in it
        for i in range(len(self.type_atoms_alias)):
            for _ in range(self.num_atoms[i]):
                self.basis[index].occ_alias = self.type_atoms_alias[i]
                if species[self.basis[index].occupant].charge is not None:
                    if self.charge is None:
                        self.charge = 0.0
                    self.charge += float(
                        species[self.basis[index].occupant].charge)
                    if species[self.basis[index].occupant].charge_loc:
                        if not charge_loc:
                            self.charge_loc = self.basis[index].position
                            charge_loc = True
                            old_charge_loc = index
                        else:
                            raise POSError("Multiple locations for the charge center supplied:\n"
                                           "     site %i: %s at %s in the POS\n"
                                           "     site %i: %s at %s in the POS"
                                           % (old_charge_loc, self.basis[old_charge_loc].occupant,
                                              self.basis[old_charge_loc].position, index,
                                              self.basis[index].occupant,
                                              self.basis[index].position))
                index += 1

    def write_geom(self, filename=None):
        """ Write geom to filename or string"""
        # with open(filename, "w") as stream:
        #     stream.write("atom, type, position; %s\n" % self.title)
        #     for i, site in enumerate(self.basis):
        #         site.write(stream, i)
        arg_string = "atom, type, position; %s\n" % self.title
        for i, site in enumerate(self.basis):
            arg_string += "  " + site.write(index=i)
        if filename is not None:
            with open(filename, "w") as stream:
                stream.write(arg_string)
        return arg_string

    def write_POS(self, filename):  #pylint: disable=invalid-name
        """ Write POS to filename """
        # If a "Cell" has already written here:
        if os.path.isfile(filename):
            mini_cell = Cell.mini(filename)
            write_POS(mini_cell, self, filename)
        else:
            with open(filename, "w") as stream:
                stream.write(self.title + "\n")
                stream.write((" ").join(self.type_atoms) + "\n")
                stream.write((" ").join([str(x) for x in self.num_atoms]) + "\n")
                if self.coord_mode is "lattice":
                    stream.write("Direct" + "\n")
                else:
                    stream.write(self.coord_mode + "\n")
                for site in self.basis:
                    site.write(stream)

    def basis_dict(self):
        """ Return a dictionary where keys are unique specie 'alias' and values are lists of
            atoms.
        """
        basis_dict = dict()
        for i, atom in enumerate(self.type_atoms_alias):
            start = sum(self.num_atoms[0:(i)])
            end = start + self.num_atoms[i]
            if atom not in basis_dict.keys():
                basis_dict[atom] = self.basis[start:end]
            else:
                basis_dict[atom] += self.basis[start:end]
        return basis_dict

    def unsort_dict(self):
        """
             Returns 'unsort_dict', for which: unsorted_dict[orig_index] == sorted_index;

             unsorted_dict[sorted_index] == orig_index

             For example, 'unsort_dict[0]' returns the index into the unsorted POSCAR of the first
                atom in the sorted POSCAR
         """

        # create basis_dict, but with initial position in POSCAR instead of coordinate
        basis_dict = dict()
        for i, atom in enumerate(self.type_atoms_alias):
            start = sum(self.num_atoms[0:(i)])
            end = start + self.num_atoms[i]
            if atom not in basis_dict.keys():
                basis_dict[atom] = range(start, end)
            else:
                basis_dict[atom] += range(start, end)

        print(basis_dict)

        orig_pos = []
        for atom in sorted(basis_dict.keys()):
            orig_pos += basis_dict[atom]

        print(orig_pos)

        new_pos = range(0, len(self.basis))

        return dict(zip(new_pos, orig_pos))

    def scale(self, scale=1.0, scalex=1.0, scaley=1.0, scalez=1.0, ndim=3):    #pylint: disable=too-many-arguments
        """ Applies scale settings """
        # Lattice coord-mode doesn't need scaling!
        if self.coord_mode == "lattice":
            return
        # Point structures don't get any cell dimensions scaled!
        if ndim == 0:
            pass
        elif ndim == 1 or ndim == 2:
            raise CellError("Handling of ndim != 0 or 3 is not yet implemented!")
        elif ndim == 3:
            for site in self.basis:
                for i in range(3):
                    site.position[i] *= scale * [scalex, scaley, scalez][i]
        else:
            raise CellError("Invalid ndim value %r encountered!" % ndim)

    def scalex(self, scalex=1.0):
        """ Alias for global scale """
        self.scale(scalex=scalex)

    def scaley(self, scaley=1.0):
        """ Alias for global scale """
        self.scale(scaley=scaley)

    def scalez(self, scalez=1.0):
        """ Alias for global scale """
        self.scale(scalez=scalez)

    def to_cart(self, cell):
        """ Switch the mode to 'cartesian', taking an input cell for scale """
        if self.coord_mode == "CARTESIAN":
            return
        for site in self.basis:
            site.position = [sum([site.position[j] * cell.lattice[j][i] for j in range(3)])
                             for i in range(3)][:]
        self.coord_mode = "CARTESIAN"

def write_POS(cell, geom, filename):    #pylint: disable=invalid-name
    """ Writes a POS file using a Cell and Geom object together """
    with open(filename, "w") as stream:
        stream.write(geom.title + "\n")
        stream.write("  %.8f" % (cell._scale) + "\n")
        for i in range(3):
            stream.write("     " + ("     ").join(["%.8f" % x for x in cell.lattice[i]]) + "\n")
        stream.write((" ").join(geom.type_atoms) + "\n")
        stream.write((" ").join([str(x) for x in geom.num_atoms]) + "\n")
        if geom.coord_mode == "LATTICE":
            stream.write("Direct" + "\n")
        else:
            stream.write(geom.coord_mode + "\n")
        for site in geom.basis:
            site.write(stream)

def volume(lattice):
    """ Volume of a generic list-of-vecs """
    return abs(_dot(lattice[0], _cross(lattice[1], lattice[2])))

def _cross(a, b):   #pylint: disable=invalid-name
    """ Calculates the cross product of 2 vector-lists to avoid importing numpy """
    x = a[1]*b[2] - a[2]*b[1]   #pylint: disable=invalid-name
    y = a[2]*b[0] - a[0]*b[2]   #pylint: disable=invalid-name
    z = a[0]*b[1] - a[1]*b[0]   #pylint: disable=invalid-name
    return [x, y, z]

def _dot(a, b):   #pylint: disable=invalid-name
    """ Calculates the dot product of 2 vector-lists to avoid importing numpy """
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def _div(a, b): #pylint: disable=invalid-name
    """ Calculates the division of a vector-list by a scalar to avoid importing numpy """
    return [float(x)/b for x in a]
