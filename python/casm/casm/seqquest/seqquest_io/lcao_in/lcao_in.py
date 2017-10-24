""" lcaoIN class and associated functions, methods, and error class """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re

from .setup import Setup
from .run import Run
from .commands import Commands
from ..species import species_settings
from ..geom import Geom, Cell

class Neb(object):  #pylint: disable=too-few-public-methods
    """ Container for Neb part of Lcao.in """
    pass

class LcaoIN(object):   #pylint: disable=too-many-instance-attributes
    """ Container object for reading, parsing, and writing lcao.in files """

    def __getitem__(self, key):
        """ Giving dict-like access to props """
        if key not in ["commands", "setup", "run", "neb"]:
            raise KeyError(key)
        else:
            return self.__dict__[key]

    def __init__(self, filename="lcao.in", speciesfile=None, POS=None):
        self.commands = Commands()
        self.setup = Setup()
        self.run = None
        self.neb = None
        self.read(filename)
        if POS is not None and self.setup.geom is None:
            self.setup.geom = Geom.POS(POS, speciesfile)
            self.setup.cell = Cell.POS(POS)
            self.setup.geom.to_cart(self.setup.cell)
            self.setup["coordinate"] = "CARTESIAN"
            if self.setup["atom types"] is None:
                self.setup["atom types"] = len(self.setup.geom.num_atoms)
            if self.setup["number"] is None:
                self.setup["number"] = self.setup.geom.natoms
        if speciesfile is not None:
            self.update(speciesfile)

    def read(self, filename="lcao.in"):
        """ Reads a lcao.in file and parses all the tags """
        ### Read some sort of file ###
        with open(filename) as lcao:
            # Parse command block, always comes first
            self.commands.read_stream(lcao)
            # Enter setup block, always comes second
            self.setup.read_stream(lcao)
            # Now we have optional sections
            while True:
                line = lcao.readline()
                if line == "":
                    break
            # for line in lcao:
                # Enter run block
                if re.search(r"run\s*phase", line, re.IGNORECASE):
                    self.run = Run()
                    self.run.read_stream(lcao)

    def write(self, filename="lcao.in", geom_in_file=False, geom_filename="lcao.geom_in"):
        """ Writes the lcaoIN object as filename """
        with open(filename, "w") as stream:
            stream.write(self.commands.construct_args())
            stream.write(self.setup.construct_args(geom_in_file))
            if self.run is not None:
                stream.write(self.run.construct_args())
            if self.neb is not None:
                stream.write(self.neb.construct_args())
        if geom_in_file is False and self.setup.geom is not None:
            self.setup.geom.write_geom(geom_filename)

    def update(self, speciesfile):
        """ Updates lcao.in and geom with info from a speciesfile """
        species = species_settings(speciesfile)
        # Update atom files, if not present
        if self.setup['atom file'] is None:
            atom_files = []
            energies = []
            for indiv_spec in species.values():
                if not indiv_spec.alias + " = " + indiv_spec.atm_location in atom_files:
                    atom_files.append(indiv_spec.alias + " = " + indiv_spec.atm_location)
                    if indiv_spec.nrg is not None:
                        energies.append(indiv_spec.nrg)
            self.setup['atom file'] = atom_files[:]
            if len(energies) == len(atom_files):
                self.setup['energies'] = energies[:]
            self.setup["atom types"] = len(self.setup["atom file"])
        # Update the geom, if available
        if self.setup.geom is not None:
            self.setup.geom.update(species)
            # Update gfixed, if now available
            site_idx = 1
            gfixed_sites = []
            for n, occ in zip(self.setup.geom.num_atoms, self.setup.geom.type_atoms):   #pylint: disable=invalid-name
                if species[occ].gfixed:
                    gfixed_sites.append([site_idx, site_idx + n - 1])
                site_idx += n
            if gfixed_sites != [] and self.run.geometry["gfixed"] is None:
                self.run.geometry["gfixed"] = [x[:] for x in gfixed_sites] # Ensure copy via slices


        # Update charge, if now available
        if self.setup['charge'] is None:
            if self.setup.geom.charge != 0.0 and self.setup.geom.charge is not None:
                self.setup['charge'] = self.setup.geom.charge
                self.setup['location'] = self.setup.geom.charge_loc
                if self.setup['ionopt'] is None:
                    self.setup['ionopt'] = -2

    def get(self, keys):
        """ Crawl down object tree to get self[key[0]][key[1]][key[2]]... etc from lcao.in """
        obj = self
        for ikey in keys:
            obj = obj[ikey]
        return obj

    def set(self, keys, value):
        """ Crawl down object tree to set self[key[0]][key[1]][key[2]]... etc from lcao.in """
        obj = self
        for ikey in keys[:-1]:
            obj = obj[ikey]
        obj[keys[-1]] = value
