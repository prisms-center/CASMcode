from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import numpy as np
from casm.vasp.io import poscar, kpoints, species, incar

class VaspIO:
    """ Generate a set of VASP input files from settings files
       
        Contains:
            self.incar: Incar object
            self.poscar: Poscar object
            self.kpoints: Kpoints object
            self.species: Species dict
    """
    def __init__(self, incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, sort=False):
        """ Construct a VaspIO object
           
            Args:
                incarfile:  path to INCAR file
                prim_kpointsfile: path to primitive KPOINTS file
                prim_poscarfile: path to primitive POSCAR file
                super_poscarfile: path to POSCAR file for this configuration
                speciesfile: path to SPECIES file
             
            This functions reads the input files and generates self.kpoints appropriate for self.poscar 
            given that 'prim_kpointsfile' is for 'prim_poscarfile'.
        """
        prim = poscar.Poscar(prim_poscarfile)
        prim_kpoints = kpoints.Kpoints(prim_kpointsfile)
        
        
        self.species = species.species_settings(speciesfile)
        self.poscar = poscar.Poscar(super_poscarfile, self.species)
        self.incar = incar.Incar(incarfile, self.species, self.poscar, sort)
        self.kpoints = prim_kpoints.super_kpoints(prim, self.poscar)
    

    def write_potcar(self, filename, sort=False):
        """ Write an appropriate POTCAR """
        if sort == False:
            try:
                file = open(filename,'w')
            except:
                raise VaspIOError("Could not open file for writing: '" + filename + "'")
            for name in self.poscar.type_atoms:
                potcar = open( os.path.join(self.species[name].potcardir,'POTCAR'))
                file.write( potcar.read())
                potcar.close()
            file.close()
        else:
            
            if True:
                # new code
                
                # dict: key = alias, value = list of Sites
                pos = self.poscar.basis_dict()
                
                with open(filename,'w') as file:
                    # for each alias
                    for alias in sorted(pos.keys()):
                        # find matching IndividualSpecies with write_potcar == True
                        for name in self.species:
                            if self.species[name].alias == alias and self.species[name].write_potcar:
                                # add to POTCAR file
                                with open( os.path.join(self.species[name].potcardir,'POTCAR')) as potcar:
                                    file.write( potcar.read())
                                break
                
            else:
                # original code:
            
                # dict: key = alias, value = list of Sites
                pos = self.poscar.basis_dict()
                
                # list of unique basis occupants
                basis = []
                
                # dict: key = alias, value = list of Site occupants (? why not just use pos ?)
                alias = dict()
                for k in sorted(pos.keys()):
                    alias[k] = []
                    for atom in pos[k]:
                        if atom.occupant not in basis:
                            basis.append(atom.occupant)
                        if atom.occupant not in alias[k]:
                            alias[k].append(atom.occupant)
                
                # list of IndividualSpecies names
                specie = [species[x].name for x in self.species]
                
                # list of potcardir
                potcar_list = []
                
                # loop over unique basis occupants
                for s in basis:
                    try:
                        if (self.species[s].write_potcar):
                            potcar_list.append(self.species[s].potcardir)
                        else:
                            # ? don't understand this ?
                            if len(alias[self.species[s].alias])<=1:
                                potcar_list.append(self.species[s].potcardir)
                    except ValueError:
                        continue
                
                # write POTCAR
                with open(filename,'w') as outFile:
                    for f in potcar_list:
                        with open(os.path.join(f, 'POTCAR')) as inFile:
                            outFile.write(inFile.read())
    
            
    def write(self,dirpath):
        """ Write VASP input files in directory 'dirpath' """
        self.poscar.write(os.path.join(dirpath,'POSCAR'), self.sort)
        self.incar.write(os.path.join(dirpath,'INCAR'))
        self.kpoints.write(os.path.join(dirpath,'KPOINTS'))
        self.write_potcar(os.path.join(dirpath,'POTCAR'), self.sort)
        

