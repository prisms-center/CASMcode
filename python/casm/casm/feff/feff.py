import matplotlib
matplotlib.use('Agg')

import re
import os
import sys
import time
import json
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as sci_int

from casm.project import DirectoryStructure, ProjectSettings
from casm.project.io import read_project_settings, read_feff_settings

from pymatgen.io.feff.sets import MPXANESSet
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from prisms_jobs import Job, JobDB


class FeffError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Feff(object):
    """
    Computes band structure with pymatgen standard settings
     - k-point density is 1000 / atom
     - high symmetry path is determined from cell geometry
    """

    def __init__(self, rundir=None):
        if rundir is None:
            raise FeffError('Can not create from nothing-directory')

        print("Constructing a FEFF XAS object")
        sys.stdout.flush()

        print("  Reading CASM settings")
        self.casm_directories = DirectoryStructure(rundir)
        self.casm_settings = ProjectSettings(rundir)
        if self.casm_settings is None:
            raise FeffError("Not in a CASM project. The '.casm' directory was not found.")

        # fixed to default_clex for now
        self.clex = self.casm_settings.default_clex

        _res = os.path.split(rundir)
        self.configname = os.path.split(_res[0])[1] + "/" + _res[1]

        print("  Reading DFT and plot settings for configuration: " + self.configname)
        sys.stdout.flush()

        setfile = self.casm_directories.settings_path_crawl("relax.json", self.configname, self.clex)
        if setfile is None:
            raise FeffError("Could not find relax.json in settings directory")
        else:
            print("  Read DFT settings from:" + setfile)
        self.settings = read_project_settings(setfile)

        fefffile = self.casm_directories.settings_path_crawl("feff.json", self.configname, self.clex)
        if fefffile is None:
            raise FeffError("Could not find feff.json in settings directory")
        else:
            print("  Read FEFF settings from:" + fefffile)
        self.feff_settings = read_feff_settings(fefffile)

        self.submit_dir = rundir
        self.feff_dir = os.path.abspath(os.path.join(rundir, 'calctype.default', 'xanes'))
        self.contcar_dir = os.path.abspath(os.path.join(rundir, 'calctype.default', 'run.final'))
        self.new_incar = []
        self.status_file = os.path.abspath(os.path.join(rundir, 'calctype.default', 'status.json'))

        self.feff_program_sequence = ['rdinp', 'atomic', 'dmdw', 'pot', 'ldos', 'screen', 'opconsat',
                                      'xsph', 'fms', 'mkgtr', 'path', 'genfmt', 'ff2x', 'sfconv',
                                      'compton', 'eels']

        print("  Run directory: %s" % self.feff_dir)
        sys.stdout.flush()

    def setup(self):
        print('Getting structure for FEFF computation in directory: ')
        print('  %s' % self.feff_dir)

        st = Structure.from_file(os.path.join(self.contcar_dir, 'CONTCAR'))

        to_compute = []
        sp_compute = []

        for s in SpacegroupAnalyzer(st).get_symmetry_dataset()['equivalent_atoms']:
            if to_compute.count(str(s)) == 0:
                to_compute.append(str(s))
                sp_compute.append(str(st.species[int(s)].symbol))

        for ab in range(len(to_compute)):
            edge = ['K']

            for e in edge:
                self.feff_settings['feff_user_tags']['EDGE'] = e
                mp_xanes = MPXANESSet(absorbing_atom=int(to_compute[ab]), structure=st,
                                      nkpts=self.feff_settings['feff_nkpts'],
                                      radius=self.feff_settings['feff_radius'],
                                      user_tag_settings=self.feff_settings['feff_user_tags'])

                mp_xanes.write_input(os.path.join(self.feff_dir, str(sp_compute[ab]) + '_' +
                                                  str(e) + '_' + str(to_compute[ab])),
                                                  make_dir_if_not_present=True)

    def exec_feff(self, jobdir=None, command=None, ncpus=None, poll_check_time=5.0, err_check_time=60.0):
        """ Run FEFF program sequence

            The 'command' is executed in the directory 'jobdir'.

            Args:
                jobdir:     directory to run.  If jobdir is None, the current directory is used.
                command:    (str or None) FHI-aims execution command
                            If command != None: then 'command' is run in a subprocess
                            Else, if ncpus == 1, then command = "aims"
                            Else, command = "mpirun -np {NCPUS} aims"
                ncpus:      (int) if '{NCPUS}' is in 'command' string, then 'ncpus' is substituted in the command.
                            if ncpus==None, $PBS_NP is used if it exists, else 1
                poll_check_time: how frequently to check if the vasp job is completed
                err_check_time: how frequently to parse vasp output to check for errors
        """
        print("Begin FEFF run:")
        sys.stdout.flush()

        if jobdir is None:
            raise FeffError('Can NOT run FEFF in higher directories')

        if command is not None:
            raise FeffError('You can not specify a single program to run from FEFF.')

        currdir = os.getcwd()
        os.chdir(jobdir)

        if ncpus is None:
            if "PBS_NP" in os.environ:
                ncpus = os.environ["PBS_NP"]
            elif "SLURM_NPROCS" in os.environ:
                ncpus = os.environ["SLURM_NPROCS"]
            else:
                ncpus = 1

        err = None

        for p in self.feff_program_sequence:

            if ncpus == 1:
                command = os.path.join(self.feff_settings['feff_base_dir'], p)
            else:
                command = "mpirun -np {NCPUS} " + os.path.join(self.feff_settings['feff_base_mpi'], p)

            if re.search("NCPUS", command):
                command = command.format(NCPUS=str(ncpus))

            print("  jobdir:", jobdir)
            print("  exec:", command)
            sys.stdout.flush()

            stdout = p + '.log'
            stderr = p + '.err'

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
                            print("  FEFF run seems frozen, killing job")
                            sys.stdout.flush()
                            p.kill()

            # close output files
            sout.close()
            serr.close()

        os.chdir(currdir)

        print("Run ended")
        sys.stdout.flush()

        return err

    def run(self):
        self.setup()
        result = self.exec_feff(jobdir=self.feff_dir)
        if result is None:
            self.plot_feff(plot_dir=self.feff_dir)

    def submit(self):
        """Submit a job for this band structure computation"""
        print("Submitting configuration: " + self.configname)
        # first, check if the job has already been submitted and is not completed
        db = JobDB()
        print("Calculation directory: ", self.feff_dir)
        sub_id = db.select_regex_id("rundir", self.feff_dir)

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

        if 'feff' in state:
            if state['feff'] == 'plotted':
                print('XAS in ' + self.configname + ' is already plotted, skipping submit.')
                return
        status = state['status']

        # check the current status
        if status == "complete":
            print("Preparing to submit the FEFF jobs")
            sys.stdout.flush()

            # cd to configdir, submit jobs from configdir, then cd back to currdir
            currdir = os.getcwd()
            os.chdir(self.submit_dir)

            st = Structure.from_file(os.path.join(self.contcar_dir, 'CONTCAR'))

            to_compute = []
            sp_compute = []

            for s in SpacegroupAnalyzer(st).get_symmetry_dataset()['equivalent_atoms']:
                if to_compute.count(str(s)) == 0:
                    to_compute.append(str(s))
                    sp_compute.append(str(st.species[int(s)].symbol))

            for ab in range(len(to_compute)):
                edge = ['K']

                for e in edge:
                    rdir = os.path.join(self.feff_dir, str(sp_compute[ab]) + '_' +
                                        str(e) + '_' + str(to_compute[ab]))

                    # construct command to be run
                    cmd = ""
                    if self.settings["prerun"] is not None:
                        cmd += self.settings["prerun"] + "\n"
                    cmd += "python -c \"from casm.feff.feff import Feff; Feff('" + rdir + "').run()\"\n"
                    if self.settings["postrun"] is not None:
                        cmd += self.settings["postrun"] + "\n"

                    print("Constructing the job")
                    sys.stdout.flush()

                    # construct the Job
                    job = Job(name=self.configname + '_FEFF',
                              account=self.settings["account"],
                              nodes=int(self.settings["nodes"]),
                              ppn=int(self.settings['ppn']),
                              walltime=self.settings["walltime"],
                              pmem=self.settings["pmem"],
                              queue=self.settings["queue"],
                              message=self.settings["message"],
                              email=self.settings["email"],
                              command=cmd)

                    print("Submitting: " + rdir)
                    sys.stdout.flush()
                    # submit the job
                    job.submit()

            # return to current directory
            os.chdir(currdir)

            print("CASM FEFF job submission complete\n")
            sys.stdout.flush()

            return

        elif status == "not_converging":
            print("Status:" + status + "  Not submitting.")
            sys.stdout.flush()
            return

        elif status != "incomplete":
            raise FeffError("unexpected relaxation status: '" + status)

    def plot_feff(self, plot_dir=None):
        if plot_dir is None:
            raise FeffError('Need to supply a diretory to plot in')

        currdir = os.getcwd()
        os.chdir(plot_dir)

        st = Structure.from_file(os.path.join(self.contcar_dir, 'CONTCAR'))

        to_compute = []
        sp_compute = []
        for s in SpacegroupAnalyzer(st).get_symmetry_dataset()['equivalent_atoms']:
            if to_compute.count(str(s)) == 0:
                to_compute.append(str(s))
                sp_compute.append(st.species[int(s)].symbol)

        wt, ct = np.unique(SpacegroupAnalyzer(st).get_symmetry_dataset()['equivalent_atoms'], return_counts=True)
        weights = dict(zip(wt, ct))

        print('Plotting FEFF data in ' + plot_dir)
        print('  Atom indices in Structure to compute: ', to_compute)
        print('    Species are: ', sp_compute)
        print('    Weights are: ', weights)

        num_grid = 2000
        data = dict()

        for ab in range(len(to_compute)):
            edge = ['K']

            for e in edge:
                emin = 1E31
                emax = -1E9

                mufile = os.path.join(plot_dir,
                                      str(sp_compute[ab]) + '_' + str(e) + '_' + str(to_compute[ab]), 'xmu.dat')

                if not os.path.isfile(mufile):
                    print('  ***  WARNING: File not found:   %s' % mufile)
                    print(' You might want to rerun this for complete averages...')
                    continue

                w, e, k, mu0, chi, p = np.genfromtxt(mufile, unpack=True)

                for i in range(len(mu0)):
                    mu0 *= weights[int(to_compute[ab])]

                if self.feff_settings['feff_use_omega']:
                    if np.amin(w) < emin:
                        emin = np.amin(w)
                    if np.amax(w) > emax:
                        emax = np.amax(w)
                    data[str(sp_compute[ab]) + str(e) + str(to_compute[ab])] = {'data': [w, mu0], 'mins': [emin, emax]}
                else:
                    if np.amin(e) < emin:
                        emin = np.amin(e)
                    if np.amax(e) > emax:
                        emax = np.amax(e)
                    data[str(sp_compute[ab]) + str(e) + str(to_compute[ab])] = {'data': [e, mu0], 'mins': [emin, emax]}

                fname = os.path.join(plot_dir, str(sp_compute[ab]) + '_'
                                     + str(e) + '_' + str(to_compute[ab]), 'xanes_single.png')
                self.plot_single(data, str(sp_compute[ab]) + str(e) + str(to_compute[ab]), fname)
                print(fname)

        # find the min / max x values to interpolate
        absvals = dict()
        for ab in range(len(to_compute)):
            edge = ['K']
            for e in edge:
                absmin = 1E31
                absmax = 0
                if not str(sp_compute[ab]) + str(e) + str(to_compute[ab]) in data:
                    print('Check: ', str(sp_compute[ab]) + str(e) + str(to_compute[ab]))
                    continue
                if data[str(sp_compute[ab]) + str(e) + str(to_compute[ab])]['mins'][0] < absmin:
                    absmin = data[str(sp_compute[ab]) + str(e) + str(to_compute[ab])]['mins'][0]
                if data[str(sp_compute[ab]) + str(e) + str(to_compute[ab])]['mins'][1] > absmax:
                    absmax = data[str(sp_compute[ab]) + str(e) + str(to_compute[ab])]['mins'][1]
                absvals[str(sp_compute[ab]) + str(e)] = [absmin, absmax]

        # interpolate and average now
        avg = np.empty((num_grid, 2))
        avg.fill(0)

        for ab in range(len(to_compute)):
            edge = ['K']
            for e in edge:
                if not str(sp_compute[ab]) + str(e) + str(to_compute[ab]) in data:
                    print('Check: ', str(sp_compute[ab]) + str(e) + str(to_compute[ab]))
                    continue
                xi = np.linspace(absvals[str(sp_compute[ab]) + str(e)][0],
                                 absvals[str(sp_compute[ab]) + str(e)][1], num_grid)
                yi = sci_int.griddata(data[str(sp_compute[ab]) + str(e) + str(to_compute[ab])]['data'][0],
                                      data[str(sp_compute[ab]) + str(e) + str(to_compute[ab])]['data'][1],
                                      xi[:], method='linear', fill_value=0)
                for i in range(num_grid):
                    avg[i, 0] = xi[i]
                    avg[i, 1] += yi[i]

                avg[:, 1] /= float(len(to_compute[ab]))

        counts = []
        for s in sp_compute:
            counts.append(str(s))

        fname = 'average_xanes_O_K.png'  # + str(sp_compute[ab]) + '_' + str(e) + '.png'
        print(fname)
        lbl = r'Average O K-Edge (' + str(counts.count('O')) + ' total, $\sigma=$ 1 eV)'
        # lbl = r'Average ' + str(e) + '-Edge ' + str(sp_compute[ab]) + ' (' + str(len(to_compute[ab])) + ' total)'
        self.plot_avgs(avg, lbl, fname)

        os.chdir(currdir)

        state = {'status': 'complete', 'feff': 'plotted'}

        with open(self.status_file, 'w') as f:
            json.dump(state, f)

    def plot_single(self, in_data, loc_index, oname):
        plt.rcParams['xtick.major.size'] = 8
        plt.rcParams['xtick.major.width'] = 3
        plt.rcParams['xtick.minor.size'] = 4
        plt.rcParams['xtick.minor.width'] = 3
        plt.rcParams['xtick.labelsize'] = 16
        plt.rcParams['ytick.major.size'] = 8
        plt.rcParams['ytick.major.width'] = 3
        plt.rcParams['ytick.minor.size'] = 4
        plt.rcParams['ytick.minor.width'] = 3
        plt.rcParams['ytick.labelsize'] = 16
        plt.rcParams['axes.linewidth'] = 3

        g_y = self.gauss_broad(in_data[loc_index]['data'][0], in_data[loc_index]['data'][1],
                               self.fwhm2sigma(float(self.feff_settings['feff_plot_sigma'])))
        scale = 1 / np.nanmax(g_y)
        plt.plot(in_data[loc_index]['data'][0], g_y * scale, '-',
                 color='navy', label=r'O ($\sigma$ = 1 eV)', linewidth=1.0)

        plt.title(r'XANES', fontsize=16)
        plt.xlabel(r'Energy [eV]', fontsize=16)
        plt.ylabel(r'Absorption $\mu_{0}$', fontsize=16)
        plt.legend(loc='upper right', ncol=1, fontsize=8)
        plt.tight_layout()
        plt.savefig(oname, dpi=170)
        plt.clf()

    def plot_avgs(self, in_data, lbls, oname):
        plt.rcParams['xtick.major.size'] = 8
        plt.rcParams['xtick.major.width'] = 3
        plt.rcParams['xtick.minor.size'] = 4
        plt.rcParams['xtick.minor.width'] = 3
        plt.rcParams['xtick.labelsize'] = 16
        plt.rcParams['ytick.major.size'] = 8
        plt.rcParams['ytick.major.width'] = 3
        plt.rcParams['ytick.minor.size'] = 4
        plt.rcParams['ytick.minor.width'] = 3
        plt.rcParams['ytick.labelsize'] = 16
        plt.rcParams['axes.linewidth'] = 3

        plt.title(r'Average XANES', fontsize=16)
        plt.xlabel(r'Energy [eV]', fontsize=16)
        plt.ylabel(r'Absorption $\mu_{0}$', fontsize=16)

        g_y = self.gauss_broad(in_data[:, 0], in_data[:, 1],
                               float(self.fwhm2sigma(self.feff_settings['feff_plot_sigma'])))
        plt.plot(in_data[:, 0], g_y / np.nanmax(g_y), '-', color='navy', label=lbls, linewidth=1.0)

        plt.tight_layout()
        plt.savefig(oname, dpi=170)
        plt.close()

    @staticmethod
    def fwhm2sigma(fwhm):
        return fwhm / np.sqrt(8 * np.log(2))

    @staticmethod
    def gauss_broad(x_data, y_data, sigma):
        yvals = np.zeros(x_data.shape)
        for ii in range(len(x_data)):
            kernel = np.exp(-(x_data - x_data[ii]) ** 2 / (2 * sigma ** 2))
            kernel = kernel / sum(kernel)
            yvals[ii] = sum(y_data * kernel)
        return yvals


class FreezeError(object):
    def __init__(self):
        self.pattern = None

    def __str__(self):
        return "FEFF appears to be frozen"

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
        if most_recent < 1200:
            return False

    @staticmethod
    def fix():
        """ Fix by killing the job and resubmitting."""
        print("  Kill job and continue...")
        raise FeffError('FEFF code was frozen [no output 20 Minutes], killed. FIXME if needed...')


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
