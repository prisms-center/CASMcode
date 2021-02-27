import os
import copy

import numpy as np

from casm.aims.io.geometry import Geometry
import casm.aimswrapper
from casm.project import DirectoryStructure, ProjectSettings


class KpointsError(Exception):
    def __init__(self, msg):
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
    def __init__(self, filename):
        """ Constructs a Kpoints object"""
        try:
            file = open(filename)
        except IOError:
            raise KpointsError('IOError' + filename)

        for line in file:
            try:
                if "k_grid" in line:
                    self.subdivisions = [float(x) for x in line.split()[1:4]]
                else:
                    pass
            except IOError:
                raise KpointsError("Error reading k_grid from line: '" + line + "'\nIn file: '" + filename + "'")

            try:
                if "k_offset" in line:
                    self.shift = [float(x) for x in line.split()[1:4]]
                else:
                    pass
            except IOError:
                raise KpointsError("Error reading k_offset from line: '" + line + "'\nIn file: '" + filename + "'")
        file.close()

    def super_kpoints(self, prim, super_cell):
        """ Assuming 'self' is the kpoints associated with a PRIM, it uses a scaling method to calculate
                  the kpoint-mesh for a supercell, such that it has a equal or greater kpoint
                  density than the prim. 
            
            Returns:
                super_kpoints: a Kpoints object for the supercell
            
            Args:
                prim: Poscar object for the prim
                super_cell: a Kpoints object for the supercell
        """
        configdir = os.getcwd()
        _res = os.path.split(configdir)
        cfgname = os.path.split(_res[0])[1] + "/" + _res[1]
        casm_dirs = DirectoryStructure(configdir)
        casm_sets = ProjectSettings(configdir)
        clex = casm_sets.default_clex
        setfile = casm_dirs.settings_path_crawl("relax.json", cfgname, clex)
        settings = casm.aimswrapper.read_settings(setfile)

        super_kpoints = copy.deepcopy(self)

        if prim is None:
            raise KpointsError("This error means there is no geometry.skel...")

        super_kpoints.subdivisions = [1, 1, 1]

        # calculate prim volumetric kpoint densities
        prim_density = self.density(prim)

        # calculate recip lattice vector lengths
        super_recip_vec_lengths = [np.linalg.norm(super_cell.reciprocal_lattice(x)) for x in range(3)]

        # while supercell kpoint density is less than prim kpoint density
        while super_kpoints.density(super_cell) < prim_density:
            # increase the number of subdivisions along the least dense super recip vector
            linear_density = [super_kpoints.subdivisions[x] / super_recip_vec_lengths[x] for x in range(3)]
            min_index = linear_density.index(min(linear_density))
            super_kpoints.subdivisions[min_index] += 1
            
            # set all subdivisions to be at similar linear density 
            scale = super_kpoints.subdivisions[min_index] / super_recip_vec_lengths[min_index]
            for i in range(2):
                super_kpoints.subdivisions[i] = int(np.ceil(scale * super_recip_vec_lengths[i] - 0.1))

        if settings["is_slab"]:  # slab override for FHI-aims ONLY
            ox = prim.lattice(0)
            oy = prim.lattice(1)

            sx = super_cell.lattice(0)
            sy = super_cell.lattice(1)

            norm_prim_x = np.linalg.norm(ox)
            norm_prim_y = np.linalg.norm(oy)

            norm_supr_x = np.linalg.norm(sx)
            norm_supr_y = np.linalg.norm(sy)

            super_kpoints.subdivisions[0] = int(np.ceil(self.subdivisions[0] * (norm_prim_x / norm_supr_x)))
            super_kpoints.subdivisions[1] = int(np.ceil(self.subdivisions[1] * (norm_prim_y / norm_supr_y)))
            super_kpoints.subdivisions[2] = 1

        return super_kpoints

    def density(self, kpts_obj):
        """ Return the kpoint density with respect to a Kpoints object. """
        return ((self.subdivisions[0] * self.subdivisions[1] * self.subdivisions[2]) /
                Geometry.reciprocal_volume(kpts_obj))

    def parse(self, filename):
        """ Parse k_grid and k_offset  """
        try:
            file = open(filename, 'r')
        except IOError:
            raise KpointsError("Write failed")

        tmp = file.readlines()
        file.close()

        ln = 0
        for line in tmp:
            if "k_grid" in line:
                tmp[ln] = "k_grid " + str(self.subdivisions[0]) + " " + str(self.subdivisions[1]) + \
                          " " + str(self.subdivisions[2]) + "\n"
            if "k_offset" in line:
                tmp[ln] = "k_offset " + str(self.shift[0]) + " " + str(self.shift[1]) + \
                          " " + str(self.shift[2]) + "\n"
            ln += 1

        file = open(filename, 'w')
        for line in tmp:
            file.write(line)
        file.close()

        del tmp

        return
