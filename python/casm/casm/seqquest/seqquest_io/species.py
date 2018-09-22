"""
Sample Species File:
atm_dir_path = /absolute/path/to/atm_files
SPECIES   ALIAS    atm_location       CHARGE    CHARGE_LOC ENERGY_REF  gfixed
In        In       In_d0.atm          0           False     -15.875    False
Ga        Ga       Ga_d0.atm          0           False     -9.433     False
As        As       As.atm             0           False     -16.20     False
Va        Va       vac.atm            0           False     0          True
Va_p      Va       vac.atm            1           True      0          True

Other valid tags:
    MASS, SPIN_POL

"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import re
from .seq_exceptions import SpeciesError

class IndividualSpecies(object):    #pylint: disable=too-many-instance-attributes
    """
        The IndividualSpecies class contains:
            self.name(str): the name as listed in the POS file
            self.alias (str): the species file lists the name, geom printing uses this string
            self.atm_dir (str): common directory for all *.atm files
            self.atm_location (str): location of specific *.atm file relative to self.atm_dir
            self.atm_file (str): abspath of particular *.atm (atm_dir/atm_file)
            self.charge (float): charge of this site, if any
            self.charge_loc (bool): whether or not this site can be a charge center
            self.mass (float): atomic masses
            self.spin (float): excess electrons in majority spin
            self.gfixed (bool): whether species should be gfixed in space
            self.tags (dict): the misc tags that need to be modified in the lcao.in for this specie
                                (ex. MAGMOM = 2) All values are stored as strings.
    """
    def __init__(self, values, tags, tags_idx, atm_dir, #pylint: disable=too-many-arguments
                 charge_idx=None, charge_loc_idx=None,
                 nrg_idx=None, mass_idx=None, spin_idx=None,
                 gfixed_idx=None):
        """ Construct an IndividualSpecies.

            Args:
                values (str list): entries in SPECIES file row
                tags (str list): column names not specifically handled in SPECIES file, lcao.in
                    tags that need to be modified
                tags_idx (int list): indices into columsn given by tags
                atm_dir (str):common directory for all *.atm files
                charge_idx (int): index of column with 'CHARGE'
                charge_loc_idx (int): index of column with 'CHAGE_LOC'
                nrg_idx (int): index of column with 'ENERGY_REF'
                mass_idx (int): index of column with 'MASS'
                spin_idx (int): index of column with 'SPIN_POL'
                gfixed_idx (int): index of column with 'gfixed'
        """
        tag_len = (len(tags) + 3
                   + bool(charge_idx) + bool(charge_loc_idx)
                   + bool(nrg_idx) + bool(mass_idx) + bool(spin_idx) + bool(gfixed_idx))
        if len(values) != tag_len:
            raise SpeciesError("Length of values != length of tags + %i.\nvalues = " % tag_len
                               + str(values) + "\ntags = " + str(tags))
        self.name = values[0]
        self.alias = values[1]
        self.atm_dir = atm_dir
        self.atm_location = values[2]
        self.atm_file = os.path.join(self.atm_dir, self.atm_location.lstrip("/"))
        self.tags = dict()
        if charge_idx is not None:
            self.charge = float(values[charge_idx])
            self.charge_loc = (values[charge_loc_idx].lower() == "true")
        else:
            self.charge = None
            self.charge_loc = None
        if nrg_idx is not None:
            self.nrg = float(values[nrg_idx])
        else:
            self.nrg = None
        if mass_idx is not None:
            self.mass = float(values[mass_idx])
        else:
            self.mass = None
        if spin_idx is not None:
            self.spin = float(values[spin_idx])
        else:
            self.spin = None
        if gfixed_idx is not None:
            self.gfixed = (values[gfixed_idx].lower() == 'true')
        else:
            self.gfixed = False

        for tag, idx in zip(tags, tags_idx):
            self.tags[tag] = values[idx]

    def write_header(self, stream):
        """ Write header to a file """
        stream.write("atm_dir_path = " + self.atm_dir + "\n")
        headers = "{0:<12} {1:<12} {2:<36} {3:<12} {4:<12} ".format(
            "SPECIES", "ALIAS", "atm_location", "CHARGE", "CHARGE_LOC")
        for key in sorted(self.tags.keys()):
            headers += "{0:<12}".format(key)
        stream.write(headers + "\n")


    def write(self, stream):
        """ Write IndividualSpecies line"""
        values = "{0:<12} {1:<12} {2:<36} {3:<12} {4:<12} ".format(
            self.name, self.alias, self.atm_location, self.charge, self.charge_loc)
        for key in sorted(self.tags.keys()):
            values += "{0:<12}".format(self.tags[key])
        stream.write(values + "\n")


    def print_all(self):
        """ Print contents of SPECIES """
        print(self.name)
        print(self.alias)
        print(self.tags)
        print(self.atm_dir)
        print(self.atm_location)
        print(self.atm_file)

def species_settings(filename): #pylint: disable=too-many-locals, too-many-branches
    """ Returns a dict of IndividualSpecies objects, with keys equal to their names. """
    with open(filename) as stream:
        # Read atm_dir_path from first line
        line = stream.readline()
        match = re.match(r"atm_dir_path\s*=\s*(.*)", line)
        if not match:
            raise SpeciesError("""Could not read atm_dir_path.
Expected: atm_dir_path = /path/to/atm_dir
Found: '" + line + "'""")
        atm_dir_path = match.group(1)

        # Allocating
        tags = []
        tags_idx = []
        charge_idx = None
        charge_loc_idx = None
        nrg_idx = None
        mass_idx = None
        spin_idx = None
        gfixed_idx = None

        ### Parsing the header ###
        header = stream.readline().strip()
        column_names = header.split()
        # Check if the header is even valid
        if len(column_names) < 3:
            raise SpeciesError("Insufficient number of columns in SPECIES file")
        for i, name in enumerate(column_names):
            if name in ['SPECIES', 'ALIAS', 'atm_location']:
                pass
            # Look for 'CHARGE'
            elif name == 'CHARGE':
                charge_idx = i
                try:
                    charge_loc_idx = column_names.index("CHARGE_LOC")
                except ValueError:
                    raise SpeciesError("Keywords CHARGE and CHARGE_LOC must both be specified!")
            elif name == 'CHARGE_LOC':
                pass

            # Look for 'ENERGY_REF'
            elif name == 'ENERGY_REF':
                nrg_idx = i

            # Look for 'MASS'
            elif name == 'MASS':
                mass_idx = i
            # Look for 'SPIN_POL'
            elif name == 'SPIN_POL':
                spin_idx = i
            # Look for 'gfixed'
            elif name == 'gfixed':
                gfixed_idx = i

            else:
                tags += [name]
                tags_idx += [i]
        my_species_settings = dict()
        for line in stream:
            if line.strip():
                values = line.strip().split()
                my_species_settings[values[0]] = IndividualSpecies(values, tags, tags_idx,
                                                                   atm_dir_path, charge_idx,
                                                                   charge_loc_idx, nrg_idx,
                                                                   mass_idx, spin_idx, gfixed_idx)

    return my_species_settings


def write_species_settings(species, filename):
    """ Write a SPECIES file from a species dict """
    with open(filename, 'w') as stream:
        species.values()[0].write_header(stream)
        for spec in sorted(species.keys()):
            species[spec].write(stream)
