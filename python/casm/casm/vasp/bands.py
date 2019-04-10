import re
import os
import sys
import time
import subprocess

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
            raise BandsError('Can not create from nothing-directory')

        self.band_dir = os.path.abspath(os.path.join(rundir, 'calctype.default', 'band_structure'))
        self.contcar_dir = os.path.abspath(os.path.join(rundir, 'calctype.default', 'run.final'))
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

        s.to(filename=os.path.join(self.band_dir, 'POSCAR'), fmt='POSCAR')
        Kpoints.automatic_linemode(100, irr_bri_zone).write_file(os.path.join(self.band_dir, 'KPOINTS'))
        sh.copyfile(os.path.join(self.contcar_dir, 'POTCAR'), os.path.join(self.band_dir, 'POTCAR'))
        self.manage_tags(os.path.join(self.contcar_dir, 'INCAR'))
        with open(os.path.join(self.band_dir, 'INCAR'), 'w') as f:
            for line in self.new_incar:
                f.write(line)

    @staticmethod
    def exec_dft(jobdir=None, stdout="std.out", stderr="std.err", command=None, ncpus=None,
                 poll_check_time=5.0, err_check_time=60.0, settings=None):
        """ Run selected DFT software using subprocess.

            The 'command' is executed in the directory 'jobdir'.

            Args:
                jobdir:     directory to run.  If jobdir is None, the current directory is used.
                stdout:     filename to write to.  If stdout is None, "std.out" is used.
                stderr:     filename to write to.  If stderr is None, "std.err" is used.
                command:    (str or None) FHI-aims execution command
                            If command != None: then 'command' is run in a subprocess
                            Else, if ncpus == 1, then command = "aims"
                            Else, command = "mpirun -np {NCPUS} aims"
                ncpus:      (int) if '{NCPUS}' is in 'command' string, then 'ncpus' is substituted in the command.
                            if ncpus==None, $PBS_NP is used if it exists, else 1
                settings:   The global CASM settings parsed from relax.json
                poll_check_time: how frequently to check if the vasp job is completed
                err_check_time: how frequently to parse vasp output to check for errors
        """
        print("Begin band strcture run:")
        sys.stdout.flush()

        if jobdir is None:
            jobdir = os.getcwd()

        currdir = os.getcwd()
        os.chdir(jobdir)

        if ncpus is None:
            if "PBS_NP" in os.environ:
                ncpus = os.environ["PBS_NP"]
            else:
                ncpus = 1

        if command is None:
            if ncpus == 1:
                command = settings['software']
            else:
                command = "mpirun -np {NCPUS} " + settings['software']

        if re.search("NCPUS", command):
            command = command.format(NCPUS=str(ncpus))

        print("  jobdir:", jobdir)
        print("  exec:", command)
        sys.stdout.flush()

        err = None
        sout = open(os.path.join(jobdir, stdout), 'w')
        serr = open(os.path.join(jobdir, stderr), 'w')

        p = subprocess.Popen(command.split(), stdout=sout, stderr=serr)

        # wait for process to end, and periodically check for errors
        last_check = time.time()
        while p.poll() is None:
            time.sleep(poll_check_time)
            if time.time() - last_check > err_check_time:
                last_check = time.time()
                err = error_check(jobdir, os.path.join(jobdir, stdout))
                if err is not None:
                    # FreezeErrors are fatal and usually not helped with abort_scf
                    if "FreezeError" in err.keys():
                        print("  DFT run seems frozen, killing job")
                        sys.stdout.flush()
                        p.kill()

        # close output files
        sout.close()
        serr.close()

        os.chdir(currdir)

        print("Run ended")
        sys.stdout.flush()

        # check finished job for errors
        if err is None:
            err = error_check(jobdir, os.path.join(jobdir, stdout))
            if err is not None:
                print("  Found errors:", end='')
                for e in err:
                    print(e, end='')
            print("\n")

        return err

    def submit(self, settings):
        pass

    def manage_tags(self, incar_file):
        remove_tags = ['NSW', 'EDIFFG', 'IBRION', 'ISIF', 'ISMEAR']
        with open(os.path.abspath(incar_file)) as f:
            for line in f:
                if not any(tag in line for tag in remove_tags):
                    self.new_incar.append(line)
        nbands = int(Procar(os.path.join(self.contcar_dir, 'PROCAR')).nbands)
        self.new_incar.append('ISMEAR = 2\n')
        self.new_incar.append('NBANDS = %i\n' % int(nbands * 1.5))  # use 50% more bands just to make sure
        self.new_incar.append('NEDOS = 5001\n')
        self.new_incar.append('EMIN = -15\n')
        self.new_incar.append('EMAX =  15\n')
        self.new_incar.append('ICORELEVEL = 1')


class FreezeError(object):
    """DFT code appears frozen"""
    def __init__(self):
        self.pattern = None

    def __str__(self):
        return "DFT code appears to be frozen"

    @staticmethod
    def error(jobdir=None):
        """ Check if code is frozen

        Args:
            jobdir: job directory
        Returns:
            True if:
                1) no file has been modified for 5 minutes
        """

        # Check if any files modified in last 300s
        most_recent = None
        most_recent_file = None
        for f in os.listdir(jobdir):
            t = time.time() - os.path.getmtime(os.path.join(jobdir, f))
            if most_recent is None:
                most_recent = t
                most_recent_file = f
            elif t < most_recent:
                most_recent = t
                most_recent_file = f

        print("Most recent file output (" + most_recent_file + "):", most_recent, " seconds ago.")
        sys.stdout.flush()
        if most_recent < 300:
            return False

    @staticmethod
    def fix():
        """ Fix by killing the job and resubmitting."""
        print("  Kill job and continue...")
        raise BandsError('DFT code was frozen [no outpot 5 Minutes], killed. FIXME if needed...')


def error_check(jobdir, stdoutfile):
    """ Check vasp stdout for errors """
    err = dict()

    # Error to check line by line, only look for first of each type
    sout = open(stdoutfile, 'r')

    # Error to check for once
    possible = [FreezeError()]
    for p in possible:
        if p.error(jobdir=jobdir):
            err[p.__class__.__name__] = p

    sout.close()
    if len(err) == 0:
        return None
    else:
        return err
