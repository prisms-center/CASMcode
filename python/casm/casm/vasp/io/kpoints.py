from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import numpy as np
import math
import copy

from casm.wrapper.misc import remove_chars

class KpointsError(Exception):
    def __init__(self,msg):
        self.msg = msg
    
    def __str__(self):
        return self.msg


class Kpoints:
    """
    The Kpoints class contains:
        self.header: (str) the first line from the KPOINTS file being read
        self.num_points: (int) contains the value in the second line (0=>automatic)
        self.subdivisions: (list of int) the number of kpoints along each of the vectors in reciprocal space
                               or the kpoint total "length" if the mode is Automatic
        self.automode: (str) Gamma/Monkhorst-Pack/Automatic
        self.shift: (list of float) the shifts that are added to the automatically generated points
    """
    def __init__(self,filename):
        """ Constructs a Kpoints object from a KPOINTS file """
        self.read(filename)
    
    
    def read(self,filename):
        """ Reads a KPOINTS file """
        try:
            file = open(filename)
        except IOError:
            raise KpointsError('IOError',filename)
            
        # read header
        self.header = file.readline().strip()
        
        # read second line
        line = file.readline().strip()
        if len(line)<=0:
            raise KpointsError("Could not read number of points: '" + line + "'")
        words = line.split()
        if words[0]!='0':
            raise KpointsError("Non-automatic kpoint generation not implemented")
        try:
            self.num_points = [int(words[0])]
        except ValueError:
            raise KpointsError("Illegal line for the number of points")
        
        # read third line
        self.automode = file.readline().strip()
        if self.automode[0].lower() not in ['m','g', 'a']:
            raise KpointsError("Illegal mode: '" + self.automode + "'")
        
        # read subdivisions line
        line = file.readline()
        if self.automode[0].lower() == 'a':
            nEntries = 1
        else:
            nEntries = 3
        try:
            self.subdivisions = [int(word) for word in line.split()[0:nEntries]]
        except ValueError:
            raise KpointsError("The subdivisions line could not be read: '" + line + "'")
        
        # read shift line
        self.shift = [0.0, 0.0, 0.0]
        line = file.readline()
        if not line == '':
            if self.automode[0].lower() == 'a':
                raise KpointsError("Fully automatic k-point mesh generation doesn't support shifts! \
                                    \n Please remove shift line: '" + line + "'")
            if len(line.split()) != 0:
                if len(line.split()) < 3:
                    raise KpointsError("The shift line could not be understood: '" + line + "'")
                try:
                    self.shift = [float(word) for word in line.split()[0:3]]
                except ValueError:
                    raise KpointsError("The shift line could not be read: '" + line + "'")
        
        file.close()
    
    
    def super_kpoints(self, prim, super):
        """ Assuming 'self' is the kpoints associated with a PRIM, it uses a scaling method to calculate
                  the kpoint-mesh for a supercell, such that it has a equal or greater kpoint
                  density than the prim. If the kpoints associated with a PRIM are mode 'a' then this
                  process is bypassed: VASP will correctly scale the kpoints by the supercell reciprocal
                  lattice at runtime. 
            
            Returns:
                super_kpoints: a Kpoints object for the supercell
            
            Args:
                prim: Poscar object for the prim OR None (if self.automode = 'a')
                super: a Poscar object for the supercell (not used if self.automode = 'a')
        """
        
        super_kpoints = copy.deepcopy(self)
        
        if self.automode[0].lower() == 'a':
           # Do nothing if we're using a fully-automatic k-point mesh - VASP will deal with this for us
           pass

        elif True:

            if prim is None:
                raise KpointsError("No POSCAR was provided for the PRIM, so the PRIM KPOINTS could not be scaled!")
            
            super_kpoints.subdivisions = [1, 1, 1]
            
            # calculate prim volumetric kpoint densities
            prim_density = self.density(prim)
            
            # calculate recip lattice vector lengths
            super_recip_vec_lengths = [np.linalg.norm(super.reciprocal_lattice(x)) for x in range(3)]
            
            # while supercell kpoint density is less than prim kpoint density
            while super_kpoints.density(super) < prim_density:
                
                # increase the number of subdivisions along the least dense super recip vector
                linear_density = [super_kpoints.subdivisions[x]/super_recip_vec_lengths[x] for x in range(3)]
                min_index = linear_density.index(min(linear_density))
                super_kpoints.subdivisions[min_index] += 1
                
                # set all subdivisions to be at similar linear density 
                scale = super_kpoints.subdivisions[min_index] / super_recip_vec_lengths[min_index]
                for i in range(3):
                    super_kpoints.subdivisions[i] = int(math.ceil(scale * super_recip_vec_lengths[i]-0.1))
            
            # end while
                
            
        else:
            # calculate recip lattice vector lengths
            super_recip_vec_lengths = [np.linalg.norm(super.reciprocal_lattice(x)) for x in range(3)]  
            prim_recip_vec_lengths = [np.linalg.norm(prim.reciprocal_lattice(x)) for x in range(3)]
            
            # set estimated super_kpoints subdivisions, using prim subdivisions/recip length along shortest recip vector
            short_ind = np.argmin(np.array(prim_recip_vec_lengths))
            shortest = prim_recip_vec_lengths[short_ind]
            effective_subdivisions = self.subdivisions[short_index]
            scale = effective_subdivisions / shortest
            for i in range(3):
                super_kpoints.subdivisions[i] = int(math.ceil(scale * super_recip_vec_lengths[i]))
            
            # calculate kpoint densities
            prim_density = self.density(prim)
            super_density = super_kpoints.density(super)
            
            # increase effective prim subdivisions until super_density >= prim_density
            while(super_density < prim_density):
                effective_subdivisions += 1
                scale = effective_subdivisions / shortest
                for i in range(3):
                    super_kpoints.subdivisions[i] = int(math.ceil(scale * super_recip_vec_lengths[i]))
                super_density = super_kpoints.density(super)
        
        return super_kpoints
    
    
    def density(self, poscar):
        """ Return the kpoint density with respect to a Poscar.
            
            Args:
                poscar: a Poscar object
        """
        return (self.subdivisions[0] * self.subdivisions[1] * self.subdivisions[2]) / poscar.reciprocal_volume()
    
    
    def write(self,filename):
        """ Write a KPOINTS file """
        try:
            file = open(filename,'w')
        except IOError:
            raise KpointsError("Write failed")
        file.write(self.header+'\n')
        file.write(remove_chars(self.num_points, '[\[\],]')+'\n')
        file.write(self.automode+'\n')
        file.write(remove_chars(self.subdivisions, '[\[\],]')+'\n')
        if self.automode[0].lower() != 'a':
            file.write(remove_chars(self.shift, '[\[\],]')+'\n')
        file.close()
        return

