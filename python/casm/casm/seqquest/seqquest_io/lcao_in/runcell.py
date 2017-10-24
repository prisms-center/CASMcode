"""Setup and helpers for Runcell block from lcao.in"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re

class Runcell(dict): #pylint: disable=too-few-public-methods
    """ Container for Runcell part of Lcao.in """
    # __*_keys are tuples because ORDER IS IMPORTANT and these are
    # iterated over later
    __start_str = "cell optimization data"
    __end_str = "end cell optimization"

    __ucmethod_keys = (
        "ucmethod",
        )

    __ucmethod_values = (
        "BROYDEN",
        "STEEPEST",
        None,
        )

    __constraint_keys = (
        "constraint",
        )

    __constraint_values = (
        "isotropic",
        "vfixed",
        "xrelax",
        "xfixed",
        "yrelax",
        "yfixed",
        "zrelax",
        "zfixed",
        None,
        )

    __float_num_keys = (
        "pressure",
        "str_broyden",
        "bulk_modulus",
        "uc_convergence",
        )

    __int_num_keys = (
        "ucsteps",
        "uchistory",
        )

    __float_vec_keys = (
        "uniaxial_pressure",
        )

    __float_mat_keys = (
        "stress",
        )

    __all_keys = set(__ucmethod_keys +
                     __constraint_keys +
                     __float_num_keys +
                     __int_num_keys +
                     __float_vec_keys +
                     __float_mat_keys)


    @property
    def ucmethod_keys(self):
        """ view into protected member """
        return self.__ucmethod_keys


    @property
    def ucmethod_values(self):
        """ view into protected member """
        return self.__ucmethod_values

    @property
    def constraint_keys(self):
        """ view into protected member """
        return self.__constraint_keys

    @property
    def constraint_values(self):
        """ view into protected member """
        return self.__constraint_values

    @property
    def float_num_keys(self):
        """ view into protected member """
        return self.__float_num_keys

    @property
    def int_num_keys(self):
        """ view into protected member """
        return self.__int_num_keys

    @property
    def float_vec_keys(self):
        """ view into protected member """
        return self.__float_vec_keys

    @property
    def float_mat_keys(self):
        """ view into protected member """
        return self.__int_num_keys

    def __setitem__(self, key, value):  #pylint: disable=too-many-branches, too-many-statements
        """ Overriding default 'setitem' to prevent method writing after init """
        # check if we're 1. frozen, and 2. if the attr already exists
        if self.__is_frozen and not self.__contains__(key):
            raise TypeError(
                "%r could not be added: %r is frozen to adding new keys" % (key, self))
        else:
            if key in self.__ucmethod_keys:
                if not value in self.__ucmethod_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__ucmethod_values))
            elif key in self.__constraint_keys:
                if not value in self.__constraint_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__constraint_values))   #pylint: disable=line-too-long
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
            super(Runcell, self).__setitem__(key, value)

    def __init__(self):
        super(Runcell, self).__init__()
        self.__is_frozen = False
        for key in self.__all_keys:
            self.__setitem__(key, None)
        self.__setitem__("unparsed", [])
        self.__is_frozen = False

    def read_stream(self, stream):  #pylint: disable=too-many-branches, too-many-statements
        """ Parse StringIO from lcao.in to Runcell object """
        while True: #pylint: disable=too-many-nested-blocks
            unparsed = True
            line = stream.readline()
            if line == "":
                return
        # for line in stream:
            # End of setup block
            elif re.search(r"end\s*cell", line, re.IGNORECASE):
                break
            else:
                for key in self.ucmethod_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = stream.readline().strip()
                        unparsed = False
                for key in self.constraint_keys:
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
                for key in self.float_vec_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = [float(x) for x in stream.readline().strip().split()[:3]]
                        unparsed = False
                for key in self.float_mat_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        for _ in range(3):
                            self[key].append(float(stream.readline().strip().split()[:3]))
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

    def construct_args(self):
        """ Constructs and returns the cell block in a lcao.in file """
        arg_string = "cell optimization data\n"
        if self['ucmethod'] is not None:
            arg_string += "ucmethod\n  "
            arg_string += self['ucmethod'] + "\n"
        # for key in self.constraint_keys:
        #     if self[key] is True:
        #         arg_string += key + "\n  "
        if self['constraint'] is not None:
            arg_string += "constraint\n  "
            arg_string += self['constraint'] + "\n"
        for key in self.int_num_keys:
            if self[key] is not None:
                arg_string += key + "\n  "
                arg_string += "%i\n" % self[key]
        for key in self.float_num_keys:
            if self[key] is not None:
                arg_string += key + "\n  "
                arg_string += "%.8f\n" % self[key]
        if self['uniaxial_pressure'] is not None:
            arg_string += "uniaxial_pressure\n  "
            arg_string += "%.8f %.8f %.8f \n" % tuple(self['uniaxial_pressure'])
        if self['stress'] is not None:
            arg_string += "stress\n"
            for i in range(3):
                arg_string += "  %.8f %.8f %.8f\n" % tuple(self['stress'][i])
        arg_string += "end cell optimization\n"
        return arg_string
