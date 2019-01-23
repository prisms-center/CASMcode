import os
import sys

import casm
import casm.project
import casm.aimswrapper

from casm.aims.aims import continue_job, run, complete_job
from casm.aims.io.basis import basis_settings
from casm.aims.io.geometry import Geometry
from casm.aims.io.io import AIMS_INPUT_FILE_LIST, job_complete
from casm.aims.io.parser import Parser


class AimsRelaxError:
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class AimsRelax(object):
    """The Relax class contains functions for setting up, executing, and parsing a VASP relaxation.

        The relaxation is initialized in a directory containing VASP input files, called 'relaxdir'.
        It then creates the following directory structure:
        .../relaxdir/
            run.0/
            run.1/
            ...
            run.final/

        'run.i' directories are only created when ready.
        'run.final' is a final constant volume run {"ISIF":2, "ISMEAR":-5, "NSW":0, "IBRION":-1}.

        Contains:
            self.relaxdir  (.../relax)
            self.rundir    (list of .../relax/run.i)
            self.finaldir (.../relax/run.final)
    """
    def __init__(self, relaxdir=None, settings=None):
        """
        Construct a VASP relaxation job object.

        Args:
            relaxdir:  path to vasp relaxation directory
            settings:   dictionary-like object containing settings, or if None, it reads
                        the json file: .../relaxdir/relax.json

                possible settings keys are:
                    used by vasp.run() function:
                        "ncpus": number of ncpus to run mpi on
                        "aims_cmd": (default "aims") shell command to execute FHI-aims
                                    or None to use default mpirun
                        "strict_kpoint": force strict copying of kpoints in control.skel file,
                                         otherwise kpoints are scaled based on supercell size
                    used by not_converging():
                        "run_limit": (default 10) maximum number of runs to allow
                                     before setting status to "not_converging"
        """
        print("Constructing a FHI-aims Relax object")
        sys.stdout.flush()

        # store path to .../relaxdir, and create if not existing
        if relaxdir is None:
            relaxdir = os.getcwd()
        self.relaxdir = os.path.abspath(relaxdir)

        print("  Relax directory:", self.relaxdir)
        sys.stdout.flush()

        # find existing .../relaxdir/run.run_index directories, store paths in self.rundir list
        self.rundir = []
        self.errdir = []
        self.update_rundir()
        self.update_errdir()

        if settings is None:
            self.settings = dict()
        else:
            self.settings = settings

        # set default settings:
        if "aims_cmd" not in self.settings:
            self.settings["aims_cmd"] = None
        if "compress" not in self.settings:
            self.settings["compress"] = []

        storedir = 'run.' + str(settings["basis"])
        self.finaldir = os.path.join(self.relaxdir, storedir)

        print("FHI-aims Relax object constructed\n")
        sys.stdout.flush()

    def add_rundir(self):
        """Make a new run.i directory"""
        os.mkdir(os.path.join(self.relaxdir, "run." + str(len(self.rundir))))
        self.update_rundir()
        self.update_errdir()

    def update_rundir(self):
        """Find all .../config/vasp/relax/run.i directories, store paths in self.rundir list"""
        self.rundir = []
        run_index = len(self.rundir)
        while os.path.isdir(os.path.join(self.relaxdir, "run." + str(run_index))):
                self.rundir.append(os.path.join(self.relaxdir, "run." + str(run_index)))
                run_index += 1

    def add_errdir(self):
        """Move run.i to run.i_err.j directory"""
        os.rename(self.rundir[-1], self.rundir[-1] + "_err." + str(len(self.errdir)))
        self.update_errdir()

    def update_errdir(self):
        """Find all .../config/vasp/relax/run.i_err.j directories, store paths in self.errdir list"""
        self.errdir = []
        if len(self.rundir) == 0:
            pass
        else:
            err_index = len(self.errdir)
            while os.path.isdir(self.rundir[-1] + "_err." + str(err_index)):
                    self.errdir.append(self.rundir[-1] + "_err." + str(err_index))
                    err_index += 1

    def setup(self, initdir, settings):
        """ mv all files and directories (besides initdir) into initdir """
        print("Moving files into initial run directory:", initdir)

        if settings["basis"] == 'light':
            initdir = os.path.abspath(initdir)
            for p in os.listdir(self.relaxdir):
                if p in AIMS_INPUT_FILE_LIST and os.path.join(self.relaxdir, p) != initdir:
                    os.rename(os.path.join(self.relaxdir, p), os.path.join(initdir, p))
            print("\n")

        if settings["basis"] == 'tight':
            configdir = os.getcwd()
            _res = os.path.split(configdir)
            cfgname = os.path.split(_res[0])[1] + "/" + _res[1]
            casm_dirs = casm.project.DirectoryStructure(configdir)
            casm_sets = casm.project.ProjectSettings(configdir)
            clex = casm_sets.default_clex
            aimsfiles = casm.aimswrapper.aims_input_file_names(casm_dirs, cfgname, clex)
            bin_0, bin_1, bin_2, basisfile = aimsfiles

            print('Taking geometry.in.next_step from previous light run...')

            old_light_dir = os.path.join(self.relaxdir, 'run.light')
            old_geo_file = os.path.join(old_light_dir, 'geometry.in.next_step')

            if not os.path.isdir(old_light_dir):
                if not os.path.isfile(old_geo_file):
                    err_str = 'The directory ', old_light_dir, ' exists which points to a previously converged run.'
                    err_str += 'However, I cannot find ', old_geo_file, ' which usually means that the' \
                                                                        ' starting structure'
                    err_str += 'was already converged. This is very unusual and you need to check what went on.'
                    err_str += 'When you are sure that the structure converged properly, ' \
                               'copy the correct geometry file to ', old_geo_file, ' and restart...'
                    raise AimsRelaxError(err_str)
                err_str = 'Previous converged light run does not exist, not running tight on unrelaxed structures...'
                err_str += 'Stopping, please run light first'
                raise AimsRelaxError(err_str)

            my_basis = basis_settings(basisfile)
            newgeo = Geometry(os.path.join(old_light_dir, "geometry.in.next_step"))
            newgeo.write(os.path.join(self.relaxdir, "geometry.in"), my_basis)
            initdir = os.path.abspath(initdir)
            for p in os.listdir(self.relaxdir):
                if p in AIMS_INPUT_FILE_LIST and os.path.join(self.relaxdir, p) != initdir:
                    os.rename(os.path.join(self.relaxdir, p), os.path.join(initdir, p))

        print('')
        sys.stdout.flush()

    def complete(self):
        """Check if the relaxation is complete.

           Completion criteria: "Have a nice day." in std.out
        """
        outfile = os.path.join(self.finaldir, "std.out")
        if not os.path.isfile(outfile):
            return False
        if not Parser(outfile).complete():
            return False
        return True

    def converged(self):
        """Check if configuration is relaxed.

           This is called when self.rundir[-1] is complete.

           Convergence criteria: FHI-aims returns "Have a nice day."
        """
        outcarfile = os.path.join(self.rundir[-1], "std.out")
        if not os.path.isfile(outcarfile):
            return False
        if not Parser(outcarfile).complete():
            return False
        return True

    def not_converging(self):
        """Check if configuration is not converging.

           This is called when self.rundir[-1] is complete and not a constant volume job and self.converged() == False.

           Not converging criteria: >= 10 runs without completion
        """
        if len(self.rundir) >= int(self.settings["run_limit"]):
            return True
        return False

    def run(self):
        """ Perform a series of FHI-aims jobs to relax a structure. """

        print("Begin FHI-aims relaxation run")
        sys.stdout.flush()

        # get current status of the relaxation:
        status, task = self.status()
        print("\n++  status:" + status + "  next task:", task)
        sys.stdout.flush()

        while status == "incomplete":
            if task == "setup":
                self.add_rundir()
                self.setup(self.rundir[-1], self.settings)
            elif task == "relax":
                self.add_rundir()
                continue_job(self.rundir[-2], self.rundir[-1], self.settings)
            elif task == "constant":
                self.add_rundir()
                continue_job(self.rundir[-2], self.rundir[-1], self.settings)
            else:
                # probably hit walltime
                print('Assuming walltime has been hit last run, continuing there...')
                self.add_rundir()
                continue_job(self.rundir[-2], self.rundir[-1], self.settings)

            while True:
                # run FHI-aims
                result = run(self.rundir[-1], ncpus=self.settings["ncpus"])
                # if no errors, continue
                if result is None or result == self.not_converging:
                    break
                # else, attempt to fix first error
                self.add_errdir()
                os.mkdir(self.rundir[-1])
                # self.add_rundir()
#                err = result.itervalues().next()

#                print "\n++  status:", "error", "  next task:", "fix_error"
#                sys.stdout.flush()

#                print "Attempting to fix error:", str(err)
#                err.fix(self.errdir[-1],self.rundir[-1], self.settings)
#                print ""
#                sys.stdout.flush()

#                if (self.settings["backup"] != None) and len(self.rundir) > 1:
#                    print "Restoring from backups:"
#                    for f in self.settings["backup"]:
#                        if os.path.isfile(os.path.join(self.rundir[-2], f + "_BACKUP.gz")):
#                            f_in = gzip.open(os.path.join(self.rundir[-2], f + "_BACKUP.gz", 'rb'))
#                            f_out = open(os.path.join(self.rundir[-1], f, 'wb'))
#                            f_out.write(f_in.read())
#                            f_in.close()
#                            f_out.close()
#                            print f, " restored!"
#                    sys.stdout.flush()

            status, task = self.status()
            print("\n++  status:" + status + "  next task:", task)
            sys.stdout.flush()

        if status == "complete":
            if not os.path.isdir(self.finaldir):
                # mv final results to relax.final
                print("mv" + os.path.basename(self.rundir[-1]) + os.path.basename(self.finaldir))
                sys.stdout.flush()
                os.rename(self.rundir[-1], self.finaldir)
                self.rundir.pop()
                complete_job(self.finaldir, self.settings)

        return status, task

    def status(self):
        """ Determine the status of a vasp relaxation series of runs. Individual runs in the series
            are in directories labeled "run.0", "run.1", etc.

            Returns a tuple: (status = "incomplete" or "complete" or "not_converging",
                                task = continuedir or "relax" or "constant" or None)

            The first value is the status of the entire relaxation.

            The second value is the current task, where 'continuedir' is the path to a
            vasp job directory that is not yet completed, "relax" indicates another
            volume relaxation job is required, and "constant" that a constant volume run is required.
        """
        # if not yet started
        if len(self.rundir) == 0:
            return "incomplete", "setup"

        # check status of relaxation runs
        self.update_rundir()

        # if the latest run is complete:
        if job_complete(self.rundir[-1]):
            return "complete", None

        # elif not converging, return 'not_converging' error
        elif self.not_converging():
            return "not_converging", None
        else:
            return "incomplete", "relax"
