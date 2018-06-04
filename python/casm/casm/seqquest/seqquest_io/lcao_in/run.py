"""Setup and helpers for Run block from lcao.in"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re
from .bandstructure import Bandstructure
from .dynamics import Dynamics
from .geometry import Geometry
from .runcell import Runcell

class Run(dict):  #pylint: disable=too-many-public-methods
    """ Special structure for setup options commands in lcao.in """
    # __*_keys are tuples because ORDER IS IMPORTANT and these are
    # iterated over later
    __start_str = "run phase data"
    __end_str = "end run phase data"

    __spmeth_keys = (
        "spmeth",
        )

    __spmeth_values = (
        "SIMPLE",
        "LINEAR",
        None
        )

    __float_num_keys = (
        "temperature",
        "blend",
        "scfbl2",
        "convergence",
        "cutfrc",
        "cutfac",
        "spconv",
        "spblend",
        )

    __sci_num_keys = (
        "cutgrd",
        )

    __int_num_keys = (
        "states",
        "iterations",
        "history",
        "spsteps",
        )

    __bool_keys = (
        "no ges",
        "closed",
        "kclosed",
        )

    __dynamics_keys = (
        "dynamics",
        )

    __geometry_keys = (
        "geometry",
        )

    __runcell_keys = (
        "cell",
        )

    __bandstructure_keys = (
        "bandstructure",
        )

    __all_keys = set(__spmeth_keys +
                     __float_num_keys +
                     __sci_num_keys +
                     __int_num_keys +
                     __bool_keys +
                     __dynamics_keys +
                     __geometry_keys +
                     __runcell_keys +
                     __bandstructure_keys)

    @property
    def spmeth_keys(self):
        """ view into protected member """
        return self.__spmeth_keys

    @property
    def spmeth_values(self):
        """ view into protected member """
        return self.__spmeth_values

    @property
    def float_num_keys(self):
        """ view into protected member """
        return self.__float_num_keys

    @property
    def sci_num_keys(self):
        """ view into protected member """
        return self.__sci_num_keys

    @property
    def int_num_keys(self):
        """ view into protected member """
        return self.__int_num_keys

    @property
    def bool_keys(self):
        """ view into protected member """
        return self.__bool_keys

    @property
    def dynamics_keys(self):
        """ view into protected member """
        return self.__dynamics_keys

    @property
    def geometry_keys(self):
        """ view into protected member """
        return self.__geometry_keys

    @property
    def runcell_keys(self):
        """ view into protected member """
        return self.__runcell_keys

    @property
    def bandstructure_keys(self):
        """ view into protected member """
        return self.__bandstructure_keys

    @property
    def start_str(self):
        """ view into protected member """
        return self.__start_str

    @property
    def end_str(self):
        """ view into protected member """
        return self.__end_str


    def __setitem__(self, key, value):  #pylint: disable=too-many-branches, too-many-statements
        """ Overriding default 'setitem' to prevent method writing after init """
        # check if we're 1. frozen, and 2. if the attr already exists
        if self.__is_frozen and not self.__contains__(key):
            raise TypeError(
                "%r could not be added: %r is frozen to adding new keys" % (key, self))
        else:
            if key in self.__spmeth_keys:
                if not value in self.__spmeth_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__spmeth_values))
                     # __bool_keys +
            elif key in self.__dynamics_keys:
                if not isinstance(value, Dynamics) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "dynamics objects"))
            elif key in self.__geometry_keys:
                if not isinstance(value, Geometry) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "geometry objects"))
            elif key in self.__runcell_keys:
                if not isinstance(value, Runcell) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "runcell objects"))
            elif key in self.__bandstructure_keys:
                if not isinstance(value, Bandstructure) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "bandstructure objects"))
            elif key in self.__int_num_keys:
                if not isinstance(value, int) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "int objects"))
            elif key in self.__float_num_keys + self.__sci_num_keys:
                if not isinstance(value, float) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "float objects"))
            elif key in self.__bool_keys:
                if value is not True and value is not False and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, "True, False, None"))
            # Everything else is "bonus strings" and we are agnostic to type here
            else:
                pass
            super(Run, self).__setitem__(key, value)

    def __init__(self):
        super(Run, self).__init__()
        self.__is_frozen = False
        for key in self.__all_keys:
            self.__setitem__(key, None)
        self.__setitem__("unparsed", [])
        self.bandstructure = None
        self.dynamics = None
        self.geometry = None
        self.runcell = None
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
            elif re.search(r"end\s*run", line, re.IGNORECASE):
                break
            else:
                for key in self.spmeth_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = stream.readline().strip()
                        unparsed = False
                for key in self.int_num_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = int(stream.readline().strip())
                        unparsed = False
                for key in self.float_num_keys + self.sci_num_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = float(stream.readline().strip().replace("d", "e"))
                        unparsed = False
                for key in self.bool_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = True
                        unparsed = False
                    elif re.search(r"^\s*~"+key, line, re.IGNORECASE):
                        self[key] = False
                        unparsed = False
                for key in self.geometry_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self.geometry = Geometry()
                        self.geometry.read_stream(stream)
                        self[key] = self.geometry
                        unparsed = False
                for key in self.dynamics_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self.dynamics = Dynamics()
                        self.dynamics.read_stream(stream)
                        self[key] = self.dynamics
                        unparsed = False
                for key in self.runcell_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self.runcell = Runcell()
                        self.runcell.read_stream(stream)
                        self[key] = self.runcell
                        unparsed = False
                for key in self.bandstructure_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self.bandstructure = Bandstructure()
                        self.bandstructure.read_stream(stream)
                        self[key] = self.bandstructure
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
        arg_string = "run phase input data\n"
        if self['spmeth'] is not None:
            arg_string += "spmeth\n  "
            arg_string += self['spmeth'] + "\n"
        for key in self.int_num_keys:
            if self[key] is not None:
                arg_string += key + "\n  "
                arg_string += "%i\n" % self[key]
        for key in self.float_num_keys:
            if self[key] is not None:
                arg_string += key + "\n  "
                arg_string += ("%.8f\n" % self[key])
        for key in self.sci_num_keys:
            if self[key] is not None:
                arg_string += key + "\n  "
                arg_string += ("%.3g\n" % self[key]).replace("e", ".d")
        for key in self.bool_keys:
            if self[key] is True:
                arg_string += key + "\n"
        if self.dynamics is not None:
            arg_string += self.dynamics.construct_args()
        if self.geometry is not None:
            arg_string += self.geometry.construct_args()
        if self.runcell is not None:
            arg_string += self.runcell.construct_args()
        if self.bandstructure is not None:
            arg_string += self.bandstructure.construct_args()
        arg_string += "end run phase data\n"
        return arg_string
