"""Setup and helpers for Setup block from lcao.in"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re
import os
from ..geom import Cell, Geom
from ..seq_exceptions import LcaoError

class Setup(dict):  #pylint: disable=too-many-public-methods
    """ Special structure for setup options commands in lcao.in """
    # __*_keys are tuples because ORDER IS IMPORTANT and these are
    # iterated over later
    __start_str = "setup data"
    __end_str = "end setup phase data"

    __functional_keys = (
        "functional",
        )

    __functional_values = (
        "CAPZ", "CAPZSP",
        "LDA", "LDA-SP",
        "PBE", "PBE-SP",
        "GGA", "GGA-SP",
        "PW91", "PW91SP",
        "BLYP", "BLYPSP",
        "AM05", "AM05SP",
        None,
        )

    __vdw_keys = (
        "vdw_potential",
        )

    __vdw_values = (
        "OFF",
        "DFT-D2",
        "D2",
        "ULG12",
        "ULG",
        None,
        )

    __ndim_keys = (
        "dimension",
        )

    __ndim_values = (
        0,
        1,
        2,
        3,
        )

    __coord_keys = (
        "coordinate",
        )

    __coord_values = (
        "LATTICE",
        "CARTESIAN",
        None,
        )

    __ionopt_keys = (
        "ionopt",
        )

    __ionopt_values = (
        None,
        -2,
        0,
        1,
        2,
        3,
        )

    __geom_keys = (
        "atom, type",
        )

    __cell_keys = (
        "primitive",
        )

    __notes_keys = (
        "notes",
        )

    __atom_file_keys = (
        "atom file",
        )

    __vdw_data_keys = (
        "vdw_data",
        )

    __int_num_keys = (
        "atom types",
        "number",
        )

    __float_num_keys = (
        "spin",
        "dielectric",
        "scale", "scaleu",
        "scalex", "scaley", "scalez",
        "strfac",
        "charge",
        )

    __int_vec_keys = (
        "grid",
        "kgrid",
        )

    __float_list_keys = (
        "masses",
        "energies",
        )

    __float_vec_keys = (
        "efield",
        "location",
        )

    __float_mat_keys = (
        "strain",
        )

    __all_keys = set(__functional_keys +
                     __vdw_keys +
                     __ndim_keys +
                     __coord_keys +
                     __ionopt_keys +
                     __geom_keys +
                     __cell_keys +
                     __int_num_keys +
                     __float_num_keys +
                     __int_vec_keys +
                     __float_list_keys +
                     __float_vec_keys +
                     __float_mat_keys +
                     __notes_keys +
                     __atom_file_keys +
                     __vdw_data_keys)

    @property
    def start_str(self):
        """ view into protected member """
        return self.__start_str

    @property
    def end_str(self):
        """ view into protected member """
        return self.__end_str

    @property
    def functional_keys(self):
        """ view into protected member """
        return self.__functional_keys

    @property
    def functional_values(self):
        """ view into protected member """
        return self.__functional_values

    @property
    def vdw_keys(self):
        """ view into protected member """
        return self.__vdw_keys

    @property
    def vdw_values(self):
        """ view into protected member """
        return self.__vdw_values

    @property
    def ndim_keys(self):
        """ view into protected member """
        return self.__ndim_keys

    @property
    def ndim_values(self):
        """ view into protected member """
        return self.__ndim_values

    @property
    def coord_keys(self):
        """ view into protected member """
        return self.__coord_keys

    @property
    def coord_values(self):
        """ view into protected member """
        return self.__coord_values

    @property
    def ionopt_keys(self):
        """ view into protected member """
        return self.__ionopt_keys

    @property
    def ionopt_values(self):
        """ view into protected member """
        return self.__ionopt_values

    @property
    def geom_keys(self):
        """ view into protected member """
        return self.__geom_keys

    @property
    def notes_keys(self):
        """ view into protected member """
        return self.__notes_keys

    @property
    def atom_file_keys(self):
        """ view into protected member """
        return self.__atom_file_keys

    @property
    def vdw_data_keys(self):
        """ view into protected member """
        return self.__vdw_data_keys

    @property
    def cell_keys(self):
        """ view into protected member """
        return self.__cell_keys

    @property
    def int_num_keys(self):
        """ view into protected member """
        return self.__int_num_keys

    @property
    def float_num_keys(self):
        """ view into protected member """
        return self.__float_num_keys

    @property
    def int_vec_keys(self):
        """ view into protected member """
        return self.__int_vec_keys

    @property
    def float_list_keys(self):
        """ view into protected member """
        return self.__float_list_keys

    @property
    def float_vec_keys(self):
        """ view into protected member """
        return self.__float_vec_keys

    @property
    def float_mat_keys(self):
        """ view into protected member """
        return self.__float_mat_keys

    def __setitem__(self, key, value):  #pylint: disable=too-many-branches, too-many-statements
        """ Overriding default 'setitem' to prevent method writing after init """
        # check if we're 1. frozen, and 2. if the attr already exists
        if self.__is_frozen and not self.__contains__(key):
            raise TypeError(
                "%r could not be added: %r is frozen to adding new keys" % (key, self))
        else:
            if key in self.__functional_keys:
                if not value in self.__functional_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__functional_values))
            elif key in self.__vdw_keys:
                if not value in self.__vdw_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__vdw_values))
            elif key in self.__ndim_keys:
                if not value in self.__ndim_values and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__ndim_values))
            elif key in self.__coord_keys:
                if not value in self.__coord_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__coord_values))
            elif key in self.__ionopt_keys:
                if not value in self.__ionopt_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__ionopt_values))
            elif key in self.__geom_keys:
                if not isinstance(value, Geom) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "Geom objects"))
            elif key in self.__cell_keys:
                if not isinstance(value, Cell) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "cell objects"))
            elif key in self.__int_num_keys:
                if not isinstance(value, int) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "int objects"))
            elif key in self.__float_num_keys:
                if not isinstance(value, float) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float objects"))
            elif key in self.__int_vec_keys:
                if not isinstance(value, list) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "int 3-vecs"))
                if value is None:
                    pass
                elif len(value) != 3:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "int 3-vecs"))
                elif (not isinstance(value[0], int) or
                      not isinstance(value[1], int) or
                      not isinstance(value[2], int)):
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "int 3-vecs"))
            elif key in self.__float_vec_keys:
                if not isinstance(value, list) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float 3-vecs"))
                elif value is None:
                    pass
                elif len(value) != 3:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float 3-vecs"))
                elif (not isinstance(value[0], float) or
                      not isinstance(value[1], float) or
                      not isinstance(value[2], float)):
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float 3-vecs"))
            elif key in self.__float_list_keys:
                if not isinstance(value, list) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float lists"))
                elif value is None:
                    pass
                else:
                    for i_value in value:
                        if not isinstance(i_value, float):
                            raise TypeError(
                                "Key %r in %r object could not be set to %r: valid values are %r"
                                % (key, self.__class__, value, "float lists"))
            elif key in self.__float_mat_keys:
                if not isinstance(value, list) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float 3x3-matrices"))
                if value is None:
                    pass
                elif len(value) != 3:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float 3x3-matrices"))
                elif (not isinstance(value[0], list) or
                      not isinstance(value[1], list) or
                      not isinstance(value[2], list)):
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float 3x3-matrices"))
                elif (not len(value[0]) != 3 or
                      not len(value[1]) != 3 or
                      not len(value[2]) != 3):
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float 3x3-matrices"))
                elif (not isinstance(value[0][0], float) or not isinstance(value[0][1], float) or   #pylint: disable=too-many-boolean-expressions
                      not isinstance(value[0][2], float) or not isinstance(value[1][0], float) or
                      not isinstance(value[1][1], float) or not isinstance(value[1][2], float) or
                      not isinstance(value[2][0], float) or not isinstance(value[2][1], float) or
                      not isinstance(value[2][2], float)):
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float 3x3-matrices"))
            # Everything else is "bonus strings" and we are agnostic to type here
            else:
                pass
            super(Setup, self).__setitem__(key, value)

    def __init__(self):
        super(Setup, self).__init__()
        self.__is_frozen = False
        for key in self.__all_keys:
            self.__setitem__(key, None)
        self.__setitem__("unparsed", [])
        self.geom = None
        self.cell = None
        self.__is_frozen = False

    def read_stream(self, stream):  #pylint: disable=too-many-branches, too-many-statements
        """ Parse StringIO from lcao.in to _Setup object """
        while True: #pylint: disable=too-many-nested-blocks
            unparsed = True
            line = stream.readline()
            if line == "":
                return
        # for line in stream:
            # End of setup block
            elif re.search(r"end\s*setup", line, re.IGNORECASE):
                break
            else:
                for key in self.functional_keys + self.vdw_keys + self.coord_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = stream.readline().strip()
                        unparsed = False
                for key in self.ndim_keys + self.int_num_keys + self.ionopt_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = int(stream.readline().strip())
                        unparsed = False
                for key in self.float_num_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = float(stream.readline().strip())
                        unparsed = False
                for key in self.int_vec_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = [int(x) for x in stream.readline().strip().split()[:3]]
                        unparsed = False
                for key in self.float_vec_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = [float(x) for x in stream.readline().strip().split()[:3]]
                        unparsed = False
                for key in self.float_list_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        if self["atom types"] is None:
                            raise LcaoError("'atom types' must be specified before %s!" % key)
                        self[key] = []
                        for _ in range(self["atom types"]):
                            self[key].append(float(stream.readline().strip().split()[0]))
                        unparsed = False
                for key in self.float_mat_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        for _ in range(3):
                            self[key].append(float(stream.readline().strip().split()[:3]))
                        unparsed = False
                for key in self.geom_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self.geom = Geom.seq(stream)
                        self[key] = self.geom
                        unparsed = False
                    if re.search(r"^\s*geomfile", line, re.IGNORECASE):
                        # We need to figure out where lcao.in came from
                        if hasattr(stream, 'name'):
                            geom_loc = os.path.join(os.path.dirname(stream.name), "lcao.geom_in")
                        # Otherwise, guess
                        else:
                            geom_loc = "lcao.geom_in"
                        self.geom = Geom.geom(geom_loc)
                        self[key] = self.geom
                        unparsed = False
                for key in self.cell_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self.cell = Cell.seq(stream)
                        self[key] = self.cell
                        unparsed = False
                for key in self.notes_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = []
                        while True:
                            line = stream.readline().strip()
                            if re.search(r"end_notes", line):
                                break
                            else:
                                self[key].append(line)
                        unparsed = False
                for key in self.atom_file_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        if self["atom types"] is None:
                            raise LcaoError("'atom types' must be specified before %s!" % key)
                        self[key] = [stream.readline().strip()]
                        for _ in range(self["atom types"] - 1):
                            stream.readline()
                            self[key].append(stream.readline().strip())
                        unparsed = False
                for key in self.vdw_data_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = []
                        while True:
                            line = stream.readline().strip()
                            if re.search(r"end_vdw", line):
                                break
                            else:
                                self[key].append(line)
                        unparsed = False

                    # There may be other commands we don't know what to do with
                if unparsed:
                    self["unparsed"] += [line]

        # Update default values, scales, etc
        self._update_defaults()

    def _update_defaults(self): #pylint: disable=too-many-branches, too-many-statements
        """ Applies lcao.in defaults to settings, and applies scale parameters """
        # Functional defaults to PBE
        if self["functional"] is None:
            self["functional"] = 'PBE'
        # Dimension is REQUIRED
        if self["dimension"] is None:
            raise LcaoError("'dimension' must be specified in lcao.in file!")
        # Coordinate defaults to cartesian
        if self["coordinate"] is None:
            if self.geom is not None:
                if self.geom.coord_mode is not None:
                    self["coordinate"] = self.geom.coord_mode.upper()
        if self["coordinate"] is None:
            self["coordinate"] = 'cartesian'.upper()
        if self.geom is not None:
            self.geom.coord_mode = self["coordinate"]
        # Fix 'scale' things
        if self["scale"] == 1.0:
            self["scale"] = None
        if self["scale"] is not None:
            self.geom.scale(self["scale"])
            self.cell.scale(self["scale"])
            self["scale"] = None
        if self["scalex"] == 1.0:
            self["scalex"] = None
        if self["scalex"] is not None:
            self.geom.scale(self["scalex"])
            self.cell.scale(self["scalex"])
            self["scalex"] = None
        if self["scaley"] is not None:
            self.geom.scale(self["scaley"])
            self.cell.scale(self["scaley"])
            self["scaley"] = None
        if self["scalez"] is not None:
            self.geom.scale(self["scalez"])
            self.cell.scale(self["scalez"])
            self["scalez"] = None
        # Fix coord mode to always be CARTESIAN
        if self["coordinate"] == "LATTICE":
            self.geom.to_cart(self.cell)
        # Grid dimensions
        if self["grid"] is None:
            raise LcaoError("'grid' must be specified in lcao.in file! \
                (automatic grid scaling not yet supported...)")
        # atom types
        if self["atom types"] is None:
            if self["atom file"] is not None:
                self["atom types"] = len(self["atom file"])
            elif self.geom is not None:
                self["atom types"] = len(self.geom.type_atoms_alias)
        # charge
        if self["charge"] is not None:
            self.geom.charge = self["charge"]
        if self["location"] is not None:
            self.geom.charge_loc = self["location"]

        # natoms
        if self["number"] is None:
            if self.geom is not None:
                self["number"] = self.geom.natoms

        # kgrid
        if self["dimension"] >= 1:
            if self["kgrid"] is None:
                raise LcaoError("'kgrid' must be specified in lcao.in file!\
                    (automatic grid scaling not yet supported...)")

    def construct_args(self, geom_in_file=False):   #pylint: disable=too-many-branches, too-many-statements
        """ Constructs and returns the setup block in a lcao.in file """
        arg_string = "setup data:\n"
        if self['notes'] is not None:
            arg_string += "notes\n"
            for line in self['notes']:
                arg_string += "  " + line + "\n"
            arg_string += "end_notes\n"
        if self['functional'] is not None:
            arg_string += "functional\n  "
            arg_string += self['functional'] + "\n"
        if self['spin'] is not None:
            arg_string += "spin polarization\n  "
            arg_string += "%.8f\n" % self['spin']
        if self['vdw_potential'] is not None:
            arg_string += "vdw_potential\n  "
            arg_string += self['vdw_potential'] + "\n"
        if self['efield'] is not None:
            arg_string += "efield\n  "
            arg_string += "%.8f %.8f %.8f \n" % tuple(self['efield'])
        if self['dielectric'] is not None:
            arg_string += "dielectric constant\n  "
            arg_string += "%.8f\n" % self['dielectric']
        if self['dimension'] is not None:
            arg_string += "dimension of system\n  "
            arg_string += "%i\n" % self['dimension']
        if self['coordinate'] is not None:
            arg_string += "coordinate units\n  "
            arg_string += self['coordinate'].upper() + "\n"
        if self['scale'] is not None:
            arg_string += "scale\n  "
            arg_string += "%.8f\n" % self['scale']
        if self['scaleu'] is not None:
            arg_string += "scaleu\n  "
            arg_string += "%.8f\n" % self['scaleu']
        if self['scalex'] is not None:
            arg_string += "scalex\n  "
            arg_string += "%.8f\n" % self['scalex']
        if self['scaley'] is not None:
            arg_string += "scaley\n  "
            arg_string += "%.8f\n" % self['scaley']
        if self['scalez'] is not None:
            arg_string += "scalez\n  "
            arg_string += "%.8f\n" % self['scalez']
        if self['strain'] is not None:
            arg_string += "strain\n"
            for i in range(3):
                arg_string += "  %.8f %.8f %.8f\n" % tuple(self['strain'][i])
            if self['strfac'] is not None:
                arg_string += "strfac\n  "
                arg_string += "%.8f\n" % self['strfac']
        if self.cell is not None:
            arg_string += "primitive lattice vectors\n"
            for i in range(3):
                arg_string += "  %.8f %.8f %.8f\n" % tuple(self.cell.lattice[i])
        if self['grid'] is not None:
            arg_string += "grid dimensions\n  "
            arg_string += "%i %i %i\n" % tuple(self['grid'])
        if self['atom types'] is not None:
            arg_string += "atom types\n  "
            arg_string += "%i\n" % self['atom types']
        if self['atom file'] is not None:
            for line in self['atom file']:
                arg_string += "atom file\n"
                arg_string += "  " + line + "\n"
        if self['vdw_data'] is not None:
            arg_string += "vdw_data\n"
            for line in self['vdw_data']:
                arg_string += line + "\n"
            arg_string += "end_vdw"
        if self['masses'] is not None:
            arg_string += "masses\n"
            for line in self['masses']:
                arg_string += "  %.8f\n" % line
        if self['energies'] is not None:
            arg_string += "energies\n"
            for line in self['energies']:
                arg_string += "  %.8f\n" % line
        if self['ionopt'] is not None:
            arg_string += "ionopt\n  "
            arg_string += "%i\n" % self['ionopt']
        if self['charge'] is not None:
            arg_string += "charge\n  "
            arg_string += "%.8f\n" % self['charge']
            if self['location'] is not None:
                arg_string += "location of LMCC charge\n  "
                arg_string += "%.8f %.8f %.8f\n" % tuple(self['location'])
        if self['number'] is not None:
            arg_string += "number of atoms\n  "
            arg_string += "%i\n" % self['number']
        if self.geom is not None:
            if geom_in_file:
                arg_string += self.geom.write_geom()
            else:
                arg_string += "geomfile\n"
        if self['kgrid'] is not None:
            arg_string += "kgrid\n  "
            arg_string += "%i %i %i\n" % tuple(self['kgrid'])
        arg_string += "end setup phase data\n"
        return arg_string
