from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import numpy as np
import os, re

"""
Sample Species File:
POTCAR_DIR_PATH = /absolute/path/to/vasp_potentials
SPECIES    ALIAS    POTCAR  POTCAR_location    MAGMOM    NBANDS
Mn3        Mn       0       -                  3         5
Mn4        Mn       1       PAW_PBE/Mn         4         6

"""


class SpeciesError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class SpeciesDict(dict):
    """
        SpeciesDict subclasses dict so it can hold additional metadata without disrupting routines that
            rely on species_setting containing *only* key/value pairs corresponding to IndividualSpecies
    """
    def set_available_tags(self, tags):
        """ Set a metadata member that lists the available tags listed in a SPECIES file """
        self.tags = tags

class IndividualSpecies:
    """
        The IndividualSpecies class contains:
            self.name: the name as listed in the POS file
            self.alias: the species file lists the name, POSCAR printing uses this string
            self.tags: (dict) the tags that need to be modified in the INCAR for this specie (ex. MAGMOM = -1)
                        All values are stored as strings.
            self.write_potcar: decides if a POTCAR file needs to be written for this specie
                     (for example, you have Mn3 and Mn4 but want to write only one POTCAR, set
                      one of them to true and the other to false)
            self.potcardir_base: common directory for all POTCARs
            self.potcar_location: location of POTCAR directory relative to self.potcardir_base
            self.potcardir: directory containing particular POTCAR (self.potcardir_base joined with self.potcar_location)
    """

    def __init__(self, values, tags, potcardir_base):
        """ Construct an IndividualSpecies.

            Args:
                values: (str list) entries in SPECIES file row
                tags: (str list) column names 4+ in SPECIES file, INCAR tags that need to be modified
                potcardir_base: (str) common directory for all POTCARs
        """
        if len(values) != (len(tags)+4):
            raise SpeciesError("Length of values != length of tags + 4.\nvalues = " + str(values) + "\ntags = " + str(tags))
        self.name = values[0]
        self.alias = values[1]
        try:
            self.write_potcar = not (int(values[2]) == 0)
        except ValueError:
            raise SpeciesError("Could not read POTCAR: " + str(values))
        self.potcardir_base = potcardir_base
        self.potcar_location = values[3]
        self.potcardir = os.path.join(self.potcardir_base, self.potcar_location.lstrip("/"))
        self.tags = dict()
        for i,key in enumerate(tags):
            self.tags[key] = values[i+4]


    def write_header(self,file):
        """ Write header to a file """
        file.write("POTCAR_DIR_PATH = " + self.potcardir_base)
        headers = "{0:<12} {1:<12} {2:<12} {3:<36} ".format("SPECIES","ALIAS","POTCAR","POTCAR_location")
        for key in sorted(self.tags.keys()):
            headers += "{0:<12}".format(key)
        file.write(headers)


    def write(self,file):
        """ Write IndividualSpecies line"""
        values = "{0:<12} {1:<12} {2:<12} {3:<36} ".format(self.name, self.alias, self.write_potcar, self.potcar_location)
        for key in sorted(self.tags.keys()):
            values += "{0:<12}".format(self.tags[key])
        file.write(values)


    def print_all(self):
        print(self.name)
        print(self.alias)
        print(self.tags)
        print(self.write_potcar)
        print(self.potcardir_base)
        print(self.potcar_location)
        print(self.potcardir)

# replacing Species class with 'species_settings' function that returns a dict of IndividualSpecies

#class Species:
#    """
#    The Species class:
#        -specieList: a list of IndividualSpecies stored as key value pairs
#                     (specieList['Mn3'] would be the Mn3 specie)
#    """
#
#    def __init__(self,filename):
#        self.POTCAR_DIR_PATH = None
#        self.specieList = []
#        self.read(filename)
#
#    def read(self,filename):
#        try:
#            file = open(filename)
#        except IOError:
#            raise SpecieError('IOError',filename)
#
#        # Read POTCAR_DIR_PATH from first line
#        firstline = file.readline().strip().split()
#        if len(firstline) != 3 or firstline[0] != "POTCAR_DIR_PATH" or firstline[1] != "=" or not os.path.isdir(firstline[2].strip("\'\"")):
#            raise SpecieError(firstline,"expected: POTCAR_DIR_PATH = /path/to/POTCAR_DIR")
#        self.POTCAR_DIR_PATH = firstline[2].strip("\'\"")
#
#        # Parsing the header
#        header = file.readline().strip()
#        columnNames = header.split()
#        if len(columnNames)<4:
#            raise SpecieError(header,'insufficient number of columns')
#        tagList = columnNames[4:]
#        self.specieList = []
#        for line in file:
#            self.specieList.append(IndividualSpecie(line.strip(),tagList))
#
#
#    def print_all(self):
#        for specie in self.specieList:
#            specie.print_all()


def species_settings(filename):
    """ Returns a SpeciesDict of IndividualSpecies objects, with keys equal to their names. """
    try:
        file = open(filename)
    except IOError:
        raise SpeciesError("Could not open: '" + filename + "'")

    # Read POTCAR_DIR_PATH from first line
    line = file.readline()
    m = re.match("POTCAR_DIR_PATH\s*=\s*(.*)",line)
    if not m:
        raise SpeciesError("Could not read POTCAR_DIR_PATH.\nExpected: POTCAR_DIR_PATH = /path/to/POTCAR_DIR\nFound: '" + line + "'")
    POTCAR_DIR_PATH = m.group(1)

    # Parsing the header
    header = file.readline().strip()
    column_names = header.split()
    if len(column_names) < 4:
        raise SpeciesError("Insufficient number of columns in SPECIES file")
    tags = column_names[4:]
    species_settings = SpeciesDict()
    species_settings.set_available_tags(tags)
    for line in file:
        if line.strip():
            values = line.strip().split()
            species_settings[values[0]] = IndividualSpecies(values, tags, POTCAR_DIR_PATH)

    file.close()

    return species_settings


def write_species_settings(species, filename):
    """ Write a SPECIES file from a species dict """
    try:
        file = open(filename,'w')
    except:
        raise SpeciesError("Could not open file for writing: '" + filename + "'")
    species.keys()[0].write_header(file)
    for s in sorted(species.keys()):
        species[s].write(file)
    file.close()


