"""Setup and helpers for Dynamics block from lcao.in"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

class Dynamics(dict): #pylint: disable=too-few-public-methods
    """ Container for Dynamics part of Lcao.in """

    def read_stream(self, stream):
        """ Parse StringIO from lcao.in to Dynamics object """
        pass

    def construct_args(self):
        """ Constructs and returns the dynamics block in a lcao.in file """
        pass
