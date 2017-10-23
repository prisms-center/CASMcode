""" Class for creating seqquest i/o files """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import shutil
from .lcao_in import LcaoIN
from .lcao_out import LcaoOUT
from .species import species_settings

QUEST_INPUT_FILE_LIST = ["lcao.in", "lcao.in_geom"]
DEFAULT_QUEST_MOVE_LIST = []
DEFAULT_QUEST_COPY_LIST = ["lcao.in"]
DEFAULT_QUEST_REMOVE_LIST = []

def job_complete(jobdir=None):
    """return True if the seqquest job at path 'jobdir' is complete"""
    if jobdir is None:
        jobdir = os.getcwd()
    outfile = os.path.join(jobdir, "lcao.out")
    if (not os.path.isfile(outfile)) and (not os.path.isfile(outfile+".gz")):
        return False
    if LcaoOUT(outfile).complete:
        return True
    return False

def get_lcao_tag(keys, jobdir=None):
    """ Opens lcao.in in 'jobdir' and returns 'key' value """
    if jobdir is None:
        jobdir = os.getcwd()
    tlcaoin = LcaoIN(os.path.join(jobdir, "lcao.in"))
    return tlcaoin.get(keys)

class SeqquestIO(object):
    """ Generate a set of SeqQuest input files from settings files

        Contains:
            self.geom: Geometry object
            self.lcao_in: LcaoIN object
            self.species: Species dict
    """

    def __init__(self, lcao_in_file, super_poscarfile, speciesfile, extra_input_files=None):
        """ Construct a Seqquest IO object

            Args:
                lcao_in_file:  path to lcao.in template file
                super_poscarfile: path to POS for this configuration
                speciesfile: path to SPECIES file
        """
        self.species = species_settings(speciesfile)
        self.lcao_in = LcaoIN(lcao_in_file, speciesfile=speciesfile, POS=super_poscarfile)
        if extra_input_files is not None:
            self.extra_input_files = extra_input_files

    def write_atomfile(self, dirpath):
        """ Write symlinks to the atom files """
        links = []
        # Get all the unique *.atm files
        for indiv_spec in self.species.values():
            if not indiv_spec.atm_location in links:
                links.append(indiv_spec.atm_location)
        # Create symlinks
        for link in links:
            my_link = os.path.join(self.species.values()[0].atm_dir, link)
            os.symlink(my_link, os.path.join(dirpath, link))

    def write(self, dirpath):
        """ Write SeqQuest input files in directory 'dirpath' """
        self.write_atomfile(dirpath)
        self.lcao_in.write(filename=os.path.join(dirpath, "lcao.in"),
                           geom_filename=os.path.join(dirpath, "lcao.geom_in"))


        # copy extra input files
        for e_file in self.extra_input_files:
            shutil.copy(e_file, dirpath)
