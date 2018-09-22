""" lcaoMOD class and associated functions, methods, and error class """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re

from StringIO import StringIO

from .lcao_in.commands import Commands

class LcaoMOD(object):   #pylint: disable=too-few-public-methods
    """ Container object for reading and parsing lcao.mod files

        lcao.mod files provide a means of modifying a lcao.in file
            for use with the 'initial' and 'final' tags in relax.json.
            The schema is identical to that of lcao.in, except that empty
            sections are not specified (e.g., don't put 'run phase data
            \n end run phase data' if there are no tags to be modified),
            and care needs to be taken to not set contradictory tags (e.g.,
            if you set 'kclosed', you must unset 'closed'!)."""

    def __init__(self, filename="lcao.mod"):
        self.filename = filename
        self.tags = dict()
        self.commands = None
        self.setup = None
        self.run = None

        self.read(filename)

    def read(self, filename):
        """ Reads a lcao.mod file and parses all the tags """
        with open(filename) as stream:
            commands_in = StringIO()
            run_in = StringIO()
            neb_in = StringIO()
            # Break up sections into buffers
            while True:
                line = stream.readline()
                if line == "":
                    return
                if re.search(r"setup\s*data", line, re.IGNORECASE):
                    pass
                if re.search(r"setup\s*data", line, re.IGNORECASE):
                    pass
                if re.search(r"setup\s*data", line, re.IGNORECASE):
                    pass
                if re.search(r"setup\s*data", line, re.IGNORECASE):
                    pass



    def _read_command(self, stream):
        """ todo """
        pass

    def _read_setup(self, stream):
        """ todo """
        pass

    def _read_run(self, stream):
        """ todo """
        if False:     # Check for bandstructure block
            self._read_bandstructure(stream)
        elif False:       # Check for dynamics block
            self._read_dynamics(stream)
        elif False:     # Check for geometry block
            self._read_geometry(stream)
        elif False:     # Check for cell block
            self._read_cell(stream)
        pass

    def _read_bandstructure(self, stream):
        """ todo """
        pass

    def _read_dynamics(self, stream):
        """ todo """
        pass

    def _read_geometry(self, stream):
        """ todo """
        pass

    def _read_cell(self, stream):
        """ todo """
        pass

    def _read_neb(self, stream):
        """ todo """
        pass
