"""Setup and helpers for Commands block from lcao.in"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re

class Commands(dict):
    """ Special structure for command options commands in lcao.in """
    # __*_keys are tuples because ORDER IS IMPORTANT and these are
    # iterated over later
    __output_keys = (
        "output",
        )
    __output_values = (
        0,
        1,
        2,
        3,
        4,
        5,
        None,
        )

    __record_keys = (
        "set",
        )

    __do_keys = (
        "setup",
        "iters",
        "force",
        "tests",
        "dynamics",
        "relax",
        "cell",
        "neb",
        "spinopt",
        "bands",
        )

    __do_values = (
        "do",
        "no",
        None,
        )

    __bool_keys = (
        "keepsym",
        "redusym",
        )

    __all_keys = set(__output_keys + __record_keys + __do_keys + __bool_keys)

    @property
    def output_keys(self):
        """ View into protected member """
        return self.__output_keys

    @property
    def output_values(self):
        """ View into protected member """
        return self.__output_values

    @property
    def record_keys(self):
        """ View into protected member """
        return self.__record_keys

    @property
    def do_keys(self):
        """ View into protected member """
        return self.__do_keys

    @property
    def do_values(self):
        """ View into protected member """
        return self.__do_values

    @property
    def bool_keys(self):
        """ View into protected member """
        return self.__bool_keys

    def __setitem__(self, key, value):
        """ Overriding default 'setitem' to prevent method writing after init """
        # check if we're 1. frozen, and 2. if the attr already exists
        if self.__is_frozen and not self.__contains__(key):
            raise TypeError(
                "%r could not be added: %r is frozen to adding new keys" % (key, self))
        else:
            if key in self.__do_keys:
                if not value in self.__do_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__do_values))
            elif key in self.__bool_keys:
                if not value in [True, False, None]:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, [True, False, None]))
            elif key in self.__output_keys:
                if not value in self.__output_values:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, self.__do_values))
            elif key in self.__record_keys:
                if not isinstance(value, int) and value is not None:
                    raise TypeError(
                        "Key %r in %r object could not be set to %r: valid values are %r"
                        % (key, self.__class__, value, ("int or None")))
            # Everything else is "bonus strings" and we are agnostic to type here
            else:
                pass
            super(Commands, self).__setitem__(key, value)

    def __init__(self):
        super(Commands, self).__init__()
        self.__is_frozen = False
        for key in self.__all_keys:
            self.__setitem__(key, None)
        self.__setitem__("unparsed", [])
        self.__is_frozen = False

    def read_stream(self, stream):
        """ Parse StringIO from lcao.in to Commands object """
        while True:
            line = stream.readline()
            if line == "":
                return
        # for line in stream:
            # End of commands block
            if re.search(r"setup\s*data", line, re.IGNORECASE):
                return
            # Parse the line
            else:
                unparsed = True
                s_args = line.split()
                if s_args[0] in self.do_values:
                    self[s_args[1]] = s_args[0]
                    unparsed = False
                for key in self.bool_keys:
                    if re.search(r"^\s*"+key, line, re.IGNORECASE):
                        self[key] = True
                        unparsed = False
                    elif re.search(r"^\s*~"+key, line, re.IGNORECASE):
                        self[key] = False
                        unparsed = False
                if s_args[0] in self.output_keys:
                    self[s_args[0]] = int(s_args[-1])
                    unparsed = False
                if s_args[0] in self.record_keys:
                    self[s_args[0]] = int(s_args[-1])
                    unparsed = False
                if unparsed:
                    # There may be other commands we don't know what to do with
                    self["unparsed"] += [line]

    def construct_args(self):
        """ Constructs and returns the commands block in a lcao.in file """
        arg_string = ""
        for key in self.__output_keys:
            if self[key] is not None:
                arg_string += (key + " level " + str(self[key]) + "\n")
        for key in self.__do_keys:
            if self[key] is not None:
                arg_string += (self[key] + " " + key + "\n")
        for key in self.__bool_keys:
            if self[key] is not None:
                arg_string += (self[key] + "\n")
        for key in self.__record_keys:
            if self[key] is not None:
                arg_string += (key + " record_length " + str(self[key]) + "\n")
        for line in self["unparsed"]:
            arg_string += line

        return arg_string
