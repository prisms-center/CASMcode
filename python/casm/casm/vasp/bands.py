import re
import os
import sys
import time
import json
import matplotlib
import subprocess

import shutil as sh

from casm.project import DirectoryStructure, ProjectSettings
from casm.project.io import read_project_settings, read_band_settings

from pymatgen.io.vasp import Kpoints, Procar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.plotter import BSDOSPlotter

from prisms_jobs import Job, JobDB

matplotlib.use('Agg')


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
        if rundir is None:
            raise BandsError('Can not create from nothing-directory')

        print("Constructing a VASP BandRun object")
        sys.stdout.flush()

        print("  Reading CASM settings")
        self.casm_directories = DirectoryStructure(rundir)
        self.casm_settings = ProjectSettings(rundir)
        if self.casm_settings is None:
            raise BandsError("Not in a CASM project. The '.casm' directory was not found.")

        # fixed to default_clex for now
        self.clex = self.casm_settings.default_clex

        _res = os.path.split(rundir)
        self.configname = os.path.split(_res[0])[1] + "/" + _res[1]

        print("  Reading DFT and plot settings for configuration: " + self.configname)
        sys.stdout.flush()

        setfile = self.casm_directories.settings_path_crawl("relax.json", self.configname, self.clex)
        if setfile is None:
            raise BandsError("Could not find relax.json in settings directory")
        else:
            print("  Read DFT settings from:" + setfile)
        self.settings = read_project_settings(setfile)

        bandfile = self.casm_directories.settings_path_crawl("bands.json", self.configname, self.clex)
        if bandfile is None:
            raise BandsError("Could not find bands.json in settings directory")
        else:
            print("  Read band settings from:" + bandfile)
        self.band_settings = read_band_settings(bandfile)

        self.submit_dir = rundir
        self.band_dir = os.path.abspath(os.path.join(rundir, 'calctype.default', 'band_structure'))
        self.contcar_dir = os.path.abspath(os.path.join(rundir, 'calctype.default', 'run.final'))
        self.new_incar = []
        self.status_file = os.path.abspath(os.path.join(rundir, 'calctype.default', 'status.json'))

        print("  Run directory: %s" % self.band_dir)
        sys.stdout.flush()

    def setup_chg(self):
        """ Create VASP input files and take CONTCAR from last converged run """

        print('Getting VASP input files for charge computation in directory: ')
        print('  %s' % self.band_dir)

        if not os.path.isdir(os.path.join(self.band_dir)):
            os.mkdir(os.path.join(self.band_dir))

        s = Structure.from_file(os.path.join(self.contcar_dir, 'CONTCAR'))

        s.to(filename=os.path.join(self.band_dir, 'POSCAR'), fmt='POSCAR')
        Kpoints.automatic_density(s, 1000, force_gamma=True).write_file(os.path.join(self.band_dir, 'KPOINTS'))
        sh.copyfile(os.path.join(self.contcar_dir, 'POTCAR'), os.path.join(self.band_dir, 'POTCAR'))
        self.manage_tags_chg(os.path.join(self.contcar_dir, 'INCAR'))
        with open(os.path.join(self.band_dir, 'INCAR'), 'w') as f:
            for line in self.new_incar:
                f.write(line)

    def setup_band(self):
        """ Create VASP input files and take CONTCAR from last converged run """

        print('Getting VASP input files for band structure computation in directory: ')
        print('  %s' % self.band_dir)

        if not os.path.isdir(os.path.join(self.band_dir)):
            os.mkdir(os.path.join(self.band_dir))

        s = Structure.from_file(os.path.join(self.contcar_dir, 'CONTCAR'))
        irr_bri_zone = HighSymmKpath(s)

        s.to(filename=os.path.join(self.band_dir, 'POSCAR'), fmt='POSCAR')
        Kpoints.automatic_linemode(self.band_settings['band_subdiv'],
                                   irr_bri_zone).write_file(os.path.join(self.band_dir, 'KPOINTS'))
        sh.copyfile(os.path.join(self.contcar_dir, 'POTCAR'), os.path.join(self.band_dir, 'POTCAR'))
        self.manage_tags_bands(os.path.join(self.contcar_dir, 'INCAR'))
        with open(os.path.join(self.band_dir, 'INCAR'), 'w') as f:
            for line in self.new_incar:
                f.write(line)

    def exec_chg(self, jobdir=None, stdout="chg.out", stderr="chg.err", command=None, ncpus=None,
                 poll_check_time=5.0, err_check_time=60.0):
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
                poll_check_time: how frequently to check if the vasp job is completed
                err_check_time: how frequently to parse vasp output to check for errors
        """
        print("Begin charge density run:")
        sys.stdout.flush()

        if jobdir is None:
            jobdir = os.getcwd()

        currdir = os.getcwd()
        os.chdir(jobdir)

        if ncpus is None:
            if "PBS_NP" in os.environ:
                ncpus = os.environ["PBS_NP"]
            elif "SLURM_NPROCS" in os.environ:
                ncpus = os.environ["SLURM_NPROCS"]
            else:
                ncpus = 1

        if command is None:
            if ncpus == 1:
                command = self.settings['run_cmd']
            else:
                command = "mpirun -np {NCPUS} " + self.settings['run_cmd']

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
                        print("  DFT for CHG run seems frozen, killing job")
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

    def exec_band(self, jobdir=None, stdout="band.out", stderr="band.err", command=None, ncpus=None,
                  poll_check_time=5.0, err_check_time=60.0):
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
                poll_check_time: how frequently to check if the vasp job is completed
                err_check_time: how frequently to parse vasp output to check for errors
        """
        print("Begin band structure run:")
        sys.stdout.flush()

        if jobdir is None:
            jobdir = os.getcwd()

        currdir = os.getcwd()
        os.chdir(jobdir)

        if ncpus is None:
            if "PBS_NP" in os.environ:
                ncpus = os.environ["PBS_NP"]
            elif "SLURM_NPROCS" in os.environ:
                ncpus = os.environ["SLURM_NPROCS"]
            else:
                ncpus = 1

        if command is None:
            if ncpus == 1:
                command = self.settings['run_cmd']
            else:
                command = "mpirun -np {NCPUS} " + self.settings['run_cmd']

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

    def run(self):
        self.setup_chg()
        result = self.exec_chg(jobdir=self.band_dir)
        if result is not None:
            raise BandsError('Self consistent charge computation did not complete, check what happended.')
        self.setup_band()
        result = self.exec_band(jobdir=self.band_dir)
        if result is None:
            self.plot_bandos(plot_dir=self.band_dir)

    def submit(self):
        """Submit a job for this band structure computation"""
        print("Submitting configuration: " + self.configname)
        # first, check if the job has already been submitted and is not completed
        db = JobDB()
        print("Calculation directory: ", self.band_dir)
        sub_id = db.select_regex_id("rundir", self.band_dir)

        if sub_id is not []:
            for j in sub_id:
                job = db.select_job(j)
                if job["jobstatus"] != "?":
                    print("JobID: " + job["jobid"],
                          "  Jobstatus:" + job["jobstatus"] + "  Not submitting.")
                    sys.stdout.flush()
                    return

        with open(self.status_file, 'r') as f:
            state = json.load(f)

        if 'bands' in state:
            if state['bands'] == 'plotted':
                print('Band Structure in ' + self.configname + ' is already plotted, skipping submit.')
                return
        status = state['status']

        # check the current status
        if status == "complete":
            print("Preparing to submit the band structure job")
            sys.stdout.flush()

            # cd to configdir, submit jobs from configdir, then cd back to currdir
            currdir = os.getcwd()
            os.chdir(self.submit_dir)

            # construct command to be run
            cmd = ""
            if self.settings["prerun"] is not None:
                cmd += self.settings["prerun"] + "\n"
            cmd += "python -c \"from casm.vasp.bands import Bands; Bands('" + self.submit_dir + "').run()\"\n"
            if self.settings["postrun"] is not None:
                cmd += self.settings["postrun"] + "\n"

            print("Constructing the job")
            sys.stdout.flush()

            # construct a pbs.Job
            job = Job(name=self.configname + '_BANDS',
                      account=self.settings["account"],
                      nodes=int(self.settings["nodes"]),
                      ppn=int(self.settings['ppn']),
                      walltime=self.settings["walltime"],
                      pmem=self.settings["pmem"],
                      qos=self.settings["qos"],
                      queue=self.settings["queue"],
                      message=self.settings["message"],
                      email=self.settings["email"],
                      priority=self.settings["priority"],
                      command=cmd)

            print("Submitting")
            sys.stdout.flush()
            # submit the job
            job.submit()

            # return to current directory
            os.chdir(currdir)

            print("CASM band structure job submission complete\n")
            sys.stdout.flush()

            return

        elif status == "not_converging":
            print("Status:" + status + "  Not submitting.")
            sys.stdout.flush()
            return

        elif status != "incomplete":
            raise BandsError("unexpected relaxation status: '" + status)

    def manage_tags_chg(self, incar_file):
        remove_tags = ['NSW', 'EDIFFG', 'IBRION', 'ISIF', 'ISMEAR', 'LCHARG']
        with open(os.path.abspath(incar_file)) as f:
            for line in f:
                if not any(tag in line for tag in remove_tags):
                    self.new_incar.append(line)
        nbands = int(Procar(os.path.join(self.contcar_dir, 'PROCAR')).nbands)
        self.new_incar.append('ISMEAR = 2\n')
        self.new_incar.append('NBANDS = %i\n' % int(nbands * 1.5))  # use 50% more bands just to make sure
        self.new_incar.append('LCHARG = .TRUE.\n')

    def manage_tags_bands(self, incar_file):
        remove_tags = ['NSW', 'EDIFFG', 'IBRION', 'ISIF', 'ISMEAR', 'LCHARG']
        with open(os.path.abspath(incar_file)) as f:
            for line in f:
                if not any(tag in line for tag in remove_tags):
                    self.new_incar.append(line)
        nbands = int(Procar(os.path.join(self.contcar_dir, 'PROCAR')).nbands)
        self.new_incar.append('ISMEAR = 2\n')
        self.new_incar.append('NBANDS = %i\n' % int(nbands * 1.5))  # use 80% more bands just to make sure
        self.new_incar.append('ICHARG = 11\n')
        self.new_incar.append('NEDOS = 5001\n')
        self.new_incar.append('EMIN = -15\n')
        self.new_incar.append('EMAX =  15\n')
        self.new_incar.append('ICORELEVEL = 1')

    def plot_bandos(self, plot_dir=None):
        if plot_dir is None:
            raise BandsError('Need to supply a diretory to plot in')

        currdir = os.getcwd()
        os.chdir(plot_dir)

        run = Vasprun('vasprun.xml')

        if not run.converged:
            raise BandsError('The band structure computation of VASP did not converge for ' + self.configname + '.\n'
                             'Not plotting bands / DoS from it. Check what happend and rerun...')

        dos = run.complete_dos
        bst = run.get_band_structure(kpoints_filename='KPOINTS', efermi=run.efermi, line_mode=True)

        bsplot = BSDOSPlotter(bs_projection=self.band_settings['band_bs_projection'],
                              dos_projection=self.band_settings['band_dos_projection'],
                              vb_energy_range=self.band_settings['band_vb_energy_range'],
                              cb_energy_range=self.band_settings['band_cb_energy_range'],
                              fixed_cb_energy=self.band_settings['band_fixed_cb_energy'],
                              egrid_interval=self.band_settings['band_egrid_interval'],
                              font=self.band_settings['band_font'],
                              axis_fontsize=self.band_settings['band_axis_fontsize'],
                              tick_fontsize=self.band_settings['band_tick_fontsize'],
                              legend_fontsize=self.band_settings['band_legend_fontsize'],
                              bs_legend=self.band_settings['band_bs_legend'],
                              dos_legend=self.band_settings['band_dos_legend'],
                              rgb_legend=self.band_settings['band_rgb_legend'],
                              fig_size=self.band_settings['band_fig_size']).get_plot(bst, dos)

        bsplot.savefig(self.band_settings['band_plot_name'], dpi=self.band_settings['band_plot_dpi'])

        os.chdir(currdir)

        state = {'status': 'complete', 'bands': 'plotted'}

        with open(self.status_file, 'w') as f:
            json.dump(state, f)


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
