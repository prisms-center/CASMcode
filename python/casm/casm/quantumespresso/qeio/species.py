from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import numpy as np
import os, re

"""
Sample Species File:
PSEUDO_DIR_PATH = /absolute/path/to/quantum_espresso_potentials
SPECIES    ALIAS    UPF     UPF_location
Mn3        Mn       0       -         
Mn4        Mn       1       PAW_PBE/Mn

"""


class SpeciesError(Exception):
    def __init__(self,msg):
        self.msg = msg
    
    def __str__(self):
        return self.msg


class IndividualSpecies:
    """
        The IndividualSpecies class contains:
            self.name: the name as listed in the POS file
            self.alias: the species file lists the name, POSCAR printing uses this string
            self.tags: (dict) the tags that need to be modified in the Infile for this specie (ex. MAGMOM = -1)
                        All values are stored as strings. 
            self.use_upf: decides if a UPF file needs to be used for this specie
                     (for example, you have Mn3 and Mn4 but want to use only one UPF, set
                      one of them to true and the other to false)
            self.pseudo_base: common directory for all pseudopotentials
            self.pseudo_location: location of UPF directory relative to self.pseudo_base
            self.pseudodir: directory containing particular UPF (self.pseudo_base joined with self.pseudo_location)
    """

    def __init__(self, values, tags, pseudo_base):
        """ Construct an IndividualSpecies.
            
            Args:
                values: (str list) entries in SPECIES file row
                tags: (str list) column names 4+ in SPECIES file, INCAR tags that need to be modified
                pseudo_base: (str) common directory for all UPFs
        """
        if len(values) != (len(tags)+4):
            raise SpeciesError("Length of values != length of tags + 4.\nvalues = " + str(values) + "\ntags = " + str(tags))
        self.name = values[0]
        self.alias = values[1]
        try:
            self.use_upf = not (int(values[2]) == 0)
        except ValueError:
            raise SpeciesError("Could not read UPF: " + str(values))
        self.pseudo_base = pseudo_base
        self.pseudo_location = values[3]
        self.pseudodir = os.path.join(self.pseudo_base, self.pseudo_location.lstrip("/"))
        self.tags = dict()
        for i,key in enumerate(tags):
            self.tags[key] = values[i+4]
    
    
    def write_header(self,file):
        """ Write header to a file """
        file.write("PSEUDO_DIR_PATH = " + self.pseudo_base + '\n')
        headers = "{0:<12} {1:<12} {2:<12} {3:<36} ".format("SPECIES","ALIAS","UPF","UPF_location")
        for key in sorted(self.tags.keys()):
            headers += "{0:<12}".format(key)
        headers= headers + '\n'
        file.write(headers)
    
    
    def write(self,file):
        """ Write IndividualSpecies line"""
        values = "{0:<12} {1:<12} {2:<12} {3:<36} ".format(self.name, self.alias, self.use_upf, self.pseudo_location)
        for key in sorted(self.tags.keys()):
            values += "{0:<12}".format(self.tags[key])
        values = values + '\n'
        file.write(values)
        

    def print_all(self):
        print(self.name)
        print(self.alias)
        print(self.tags)
        print(self.use_upf)
        print(self.pseudo_base)
        print(self.pseudo_location)
        print(self.pseudodir)

def species_settings(filename):
    """ Returns a dict of IndividualSpecies objects, with keys equal to their names. """
    try:
        file = open(filename)
    except IOError:
        raise SpeciesError("Could not open: '" + filename + "'")
    
    # Read PSEUDO_DIR_PATH from first line
    line = file.readline()
    m = re.match("PSEUDO_DIR_PATH\s*=\s*(.*)",line)
    if not m:
        raise SpeciesError("Could not read PSEUDO_DIR_PATH.\nExpected: PSEUDO_DIR_PATH = /path/to/PSEUDO_DIR\nFound: '" + line + "'")
    PSEUDO_DIR_PATH = m.group(1)
    
    # Parsing the header
    header = file.readline().strip()
    column_names = header.split()
    if len(column_names) < 4:
        raise SpeciesError("Insufficient number of columns in SPECIES file")
    tags = column_names[4:]
    species_settings = dict()
    for line in file:
        if line.strip():
            values = line.strip().split()
            species_settings[values[0]] = IndividualSpecies(values, tags, PSEUDO_DIR_PATH)
    
    file.close()
    
    return species_settings


def write_species_settings(species, filename):
    """ Write a SPECIES file from a species dict """
    try:
        file = open(filename,'w')
    except:
        raise SpeciesError("Could not open file for writing: '" + filename + "'")
    species[species.keys()[0]].write_header(file)
    for s in sorted(species.keys()):
        species[s].write(file)
    file.close()
    
        
