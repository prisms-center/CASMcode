"""Setup and helpers for Bandstructure block from lcao.in"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

class Bandstructure(dict): #pylint: disable=too-few-public-methods
    """ Container for Bandstructure part of Lcao.in """

    def read_stream(self, stream):
        """ Parse StringIO from lcao.in to Bandstructure object """
        pass

    def construct_args(self):
        """ Constructs and returns the bandstructure block in a lcao.in file """
        pass
