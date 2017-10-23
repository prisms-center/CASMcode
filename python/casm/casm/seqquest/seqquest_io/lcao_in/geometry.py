"""Setup and helpers for Geometry block from lcao.in"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re

class Geometry(dict): #pylint: disable=too-few-public-methods
    """ Container for Geometry part of Lcao.in """

    __gmethod_keys = (
        "gmethod",
        )

    __gmethod_values = (
        "BROYDEN",
        "ASD",
        "DMDAT",
        "DAMPED",
        "STEEPEST",
        None)

    __atom_list_list_keys = (
        "gfixed",
        "grelax",
        )

    __atom_three_list_keys = (
        "frame",
        )

    __vgfix_keys = (
        "vgfixa",
        "vgfixd",
        "vgfixp",
        )

    __int_num_keys = (
        "gsteps",
        "ghistory",
        )

    __float_num_keys = (
        "gconv",
        "gblend",
        "timestep",
        )

    __all_keys = set(__gmethod_keys +
                     __atom_list_list_keys +
                     __atom_three_list_keys +
                     __vgfix_keys +
                     __int_num_keys +
                     __float_num_keys)

    @property
    def gmethod_keys(self):
        """ view into protected member """
        return self.__gmethod_keys

    @property
    def gmethod_values(self):
        """ view into protected member """
        return self.__gmethod_values


    @property
    def atom_list_list_keys(self):
        """ view into protected member """
        return self.__atom_list_list_keys

    @property
    def atom_three_list_keys(self):
        """ view into protected member """
        return self.__atom_three_list_keys

    @property
    def vgfix_keys(self):
        """ view into protected member """
        return self.__vgfix_keys

    @property
    def float_num_keys(self):
        """ view into protected member """
        return self.__float_num_keys

    @property
    def int_num_keys(self):
        """ view into protected member """
        return self.__int_num_keys

    def __setitem__(self, key, value):  #pylint: disable=too-many-branches, too-many-statements, too-many-nested-blocks
        """ Overriding default 'setitem' to prevent method writing after init """
        # check if we're 1. frozen, and 2. if the attr already exists
        if self.__is_frozen and not self.__contains__(key): #pylint: disable=too-many-nested-blocks
            raise TypeError(
                "%r could not be added: %r is frozen to adding new keys" % (key, self))
        else:
            if key in self.__gmethod_keys:
                if not value in self.__gmethod_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__gmethod_values))
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
            elif key in self.__atom_list_list_keys:
                if not isinstance(value, list) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float objects"))
                if value is None:
                    pass
                else:
                    for v in value: #pylint: disable=invalid-name
                        if len(v) != 2:
                            raise TypeError(
                                "Key %r in %r object could not be set to %r: valid values are %r"
                                % (key, self.__class__, value, "pairs of atom numbers/IDs"))
                        for vv in v: #pylint: disable=invalid-name
                            if not isinstance(vv, int):
                                raise TypeError(
                                    "Key %r in %r object could not be set to %r: valid values are %r" #pylint: disable=line-too-long
                                    % (key, self.__class__, value, "pairs of atom numbers/IDs"))
            elif key in self.__atom_three_list_keys:
                if not isinstance(value, list) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float objects"))
                if value is None:
                    pass
                else:
                    for pair in value:
                        if len(pair.split()) != 3:
                            raise TypeError(
                                "Key %r in %r object could not be set to %r: valid values are %r"
                                % (key, self.__class__, value, "triplets of atom numbers/IDs"))
            elif key in self.__vgfix_keys:
                if not isinstance(value, list) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float objects"))
                if value is None:
                    pass
                else:
                    if len(value) != 2:
                        raise TypeError(
                            "Key %r in %r object could not be set to %r: valid values are %r"
                            % (key, self.__class__, value,
                               "n_fixed x_axis y_axis z_axis\n list of n_fixed atoms"))
                    if len(value[0].split()) != 4:
                        raise TypeError(
                            "Key %r in %r object could not be set to %r: valid values are %r"
                            % (key, self.__class__, value,
                               "n_fixed x_axis y_axis z_axis\n list of n_fixed atoms"))
                    if len(value[1].split()) != int(value[0].split()[0]):
                        raise TypeError(
                            "Key %r in %r object could not be set to %r: valid values are %r"
                            % (key, self.__class__, value,
                               "n_fixed x_axis y_axis z_axis\n list of n_fixed atoms"))

            # Everything else is "bonus strings" and we are agnostic to type here
            else:
                pass
            super(Geometry, self).__setitem__(key, value)


    def __init__(self):
        super(Geometry, self).__init__()
        self.__is_frozen = False
        for key in self.__all_keys:
            self.__setitem__(key, None)
        self.__setitem__("unparsed", [])
        self.__is_frozen = False

    def read_stream(self, stream):  #pylint: disable=too-many-branches, too-many-statements
        """ Parse StringIO from lcao.in to geometry object """
        while True: #pylint: disable=too-many-nested-blocks
            unparsed = True
            line = stream.readline()
            if line == "":
                return
        # for line in stream:
            # End of setup block
            elif re.search(r"end\s*geometry", line, re.IGNORECASE):
                break
            else:
                for key in self.gmethod_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = stream.readline().strip()
                        unparsed = False
                for key in self.int_num_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = int(stream.readline().strip())
                        unparsed = False
                for key in self.float_num_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = float(stream.readline().strip())
                        unparsed = False
                for key in self.atom_list_list_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        if self[key] is None:
                            self[key] = [[int(x) for x in stream.readline().strip().split()]]
                        else:
                            self[key].append([int(x) for x in stream.readline().strip().split()])
                        unparsed = False
                for key in self.__atom_three_list_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = stream.readline().strip()
                        unparsed = False
                for key in self.vgfix_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = [stream.readline().strip(), stream.readline().strip()]
                        unparsed = False
                # There may be other commands we don't know what to do with
                if unparsed:
                    self["unparsed"] += [line]

        # Update default values, scales, etc
        self._update_defaults()

    def _update_defaults(self): #pylint: disable=no-self-use
        """ Applies lcao.in defaults to settings """
        # all keywords are optional
        return

    def construct_args(self):   #pylint: disable=too-many-branches, too-many-statements
        """ Constructs and returns the setup block in a lcao.in file """
        arg_string = "geometry optimization\n"
        if self['gmethod'] is not None:
            arg_string += "gmethod\n  "
            arg_string += self['gmethod'] + "\n"
        for key in self.int_num_keys:
            if self[key] is not None:
                arg_string += key + "\n  "
                arg_string += "%i\n" % self[key]
        for key in self.float_num_keys:
            if self[key] is not None:
                arg_string += key + "\n  "
                arg_string += "%.8f\n" % self[key]
        for key in self.atom_list_list_keys + self.atom_three_list_keys:
            if self[key] is not None:
                for line in self[key]:
                    arg_string += key + "\n  "
                    arg_string += "  ".join([str(x) for x in line]) + "\n"
        for key in self.vgfix_keys:
            if self[key] is not None:
                arg_string += key + "\n  "
                arg_string += self[key][0] + "\n  "
                arg_string += self[key][1] + "\n"
        arg_string += "end geometry\n"
        return arg_string
