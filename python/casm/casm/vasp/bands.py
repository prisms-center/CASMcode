import os
import sys
import shutil as sh

from pymatgen.io.vasp import Kpoints, Procar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath


class BandsError:
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Bands(object):
    """
    Computes band structure with pymatgen standard settings
     - k-point density is 1000 / atom
     - high symmetry path is determined from cell geometry
    """

    def __init__(self, rundir=None):
        print("Constructing a VASP BandRun object")
        sys.stdout.flush()

        # store path to .../relaxdir, and create if not existing
        if rundir is None:
            raise BandsError('Can not create from nothing directory')

        self.band_dir = os.path.abspath(os.path.join(rundir, 'band_structure'))
        self.contcar_dir = os.path.abspath(os.path.join(rundir, 'run.final'))
        self.new_incar = []

        print("  Run directory: %s" % self.band_dir)
        sys.stdout.flush()

    def setup(self):
        """ Create VASP input files and take CONTCAR from last converged run """

        print('Getting VASP input files for band structure computation in directory: ')
        print('  %s' % self.band_dir)

        if not os.path.isdir(os.path.join(self.band_dir)):
            os.mkdir(os.path.join(self.band_dir))

        s = Structure.from_file(os.path.join(self.contcar_dir, 'CONTCAR'))
        irr_bri_zone = HighSymmKpath(s)

        s.to(os.path.join(self.band_dir, 'POSCAR'), fmt='POSCAR')
        Kpoints.automatic_linemode(100, irr_bri_zone).write_file(os.path.join(self.band_dir, 'KPOINTS'))
        sh.copyfile(os.path.join(self.contcar_dir, 'POTCAR'), os.path.join(self.band_dir, 'POTCAR'))
        self.manage_tags(os.path.join(self.contcar_dir, 'INCAR'))
        with open(os.path.join(self.band_dir, 'INCAR'), 'w') as f:
            for line in self.new_incar:
                f.write(line)

    def manage_tags(self, incar_file):
        remove_tags = ['NSW', 'EDIFFG', 'IBRION', 'ISIF']
        with open(os.path.abspath(incar_file)) as f:
            for line in f:
                if any(remove_tags) not in line:
                    self.new_incar.append(line)
        nbands = int(Procar(os.path.join(self.contcar_dir, 'PROCAR')).nbands)
        self.new_incar.append('NBANDS = %i' % (nbands * 1.5))  # use 50% more bands just to make sure
        self.new_incar.append('NEDOS = 15001')
        self.new_incar.append('EMIN = -15')
        self.new_incar.append('EMAX =  10')
        self.new_incar.append('ICORELEVEL = 1')
