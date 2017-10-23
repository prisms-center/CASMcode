from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, math, sys, shutil, gzip
import casm.vasp
from casm.vasp import io

class Relax(object):
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
                    used by casm.vasp.run() function:
                        "ncpus": number of ncpus to run mpi on
			"npar" or "ncore": number of ways to parallelize
                        "kpar": number of ways to parallelize k-points
                        "vasp_cmd": (default, see casm.vasp.run) shell command to execute vasp, or None to use default mpirun
                        "strict_kpoint": force strict copying of KPOINTS file, otherwise kpoints are scaled based on supercell size
                    used by not_converging():
                        "run_limit": (default 10) maximum number of runs to allow before setting status to "not_converging"

        """

        print("Constructing a VASP Relax object")
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

        self.finaldir = os.path.join(self.relaxdir, "run.final")

        if settings is None:
            self.settings = dict()
        else:
            self.settings = settings

        # set default settings:
        if not "npar" in self.settings:
            self.settings["npar"] = None
        if not "kpar" in self.settings:
            self.settings["kpar"] = None
        if not "ncore" in self.settings:
            self.settings["ncore"] = None
        if not "vasp_cmd" in self.settings:
            self.settings["vasp_cmd"] = None
        if not "ncpus" in self.settings:
            self.settings["ncpus"] = None
        if not "run_limit" in self.settings:
            self.settings["run_limit"] = 10
        if not "nrg_convergence" in self.settings:
            self.settings["nrg_convergence"] = None
        if not "compress" in self.settings:
            self.settings["compress"] = []
        if not "move" in self.settings:
            self.settings["move"] = []
        if not "copy" in self.settings:
            self.settings["copy"] = []
        if not "backup" in self.settings:
            self.settings["backup"] = []
        if not "remove" in self.settings:
            self.settings["remove"] = []
        if not "extra_input_files" in self.settings:
            self.settings["extra_input_files"] = []
        if not "initial" in self.settings:
            self.settings["initial"] = None
        if not "final" in self.settings:
            self.settings["final"] = None
        if not "err_types" in self.settings:
            self.settings["err_types"] = ['SubSpaceMatrixError']

        print("VASP Relax object constructed\n")
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
        while os.path.isdir( os.path.join(self.relaxdir, "run." + str(run_index))):
                self.rundir.append( os.path.join(self.relaxdir, "run." + str(run_index)) )
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
        initdir = os.path.abspath(initdir)
        for p in os.listdir(self.relaxdir):
            if (p in (io.VASP_INPUT_FILE_LIST + self.settings["extra_input_files"])) and (os.path.join(self.relaxdir, p) != initdir):
                os.rename(os.path.join(self.relaxdir,p), os.path.join(initdir,p))
        print("")
        sys.stdout.flush()

        # Keep a backup copy of the base INCAR
        shutil.copyfile(os.path.join(initdir,"INCAR"),os.path.join(self.relaxdir,"INCAR.base"))

        # If an initial incar is called for, copy it in and set the appropriate flag
        if (self.settings["initial"] != None) and (os.path.isfile(os.path.join(self.relaxdir,self.settings["initial"]))):
            new_values = io.Incar(os.path.join(self.relaxdir,self.settings["initial"])).tags
            io.set_incar_tag(new_values, initdir)
            print("  Set INCAR tags:", new_values, "\n")
            sys.stdout.flush()

    def complete(self):
        """Check if the VASP relaxation is complete.

           Completion criteria: .../config/vasp/relax/run.final/OUTCAR exists and is complete
        """
        outcarfile = os.path.join(self.finaldir,"OUTCAR")
        if not os.path.isfile(outcarfile):
            return False
        if not io.Outcar(outcarfile).complete():
            return False
        return True


    def converged(self):
        """Check if configuration is relaxed.

           This is called when self.rundir[-1] is complete and not a constant volume job.

           Convergence criteria is: at least 2 relaxation jobs are complete, and:
                                    1) the last job completed with <= 3 ionic steps
                                    or 2) the last two jobs had final E0 differ by less than
                                          self.settings["nrg_convergence"]
        """
        if len(self.rundir) >= 2:
            if io.ionic_steps(self.rundir[-1]) <= 3:
                return True
            if self.settings["nrg_convergence"] != None:
                if io.job_complete(self.rundir[-1]) and io.job_complete(self.rundir[-2]):
                    o1 = io.Oszicar(os.path.join(self.rundir[-1],"OSZICAR"))
                    o2 = io.Oszicar(os.path.join(self.rundir[-2],"OSZICAR"))
                    if abs( o1.E[-1] - o2.E[-1]) < self.settings["nrg_convergence"]:
                        return True

        return False


    def not_converging(self):
        """Check if configuration is not converging.

           This is called when self.rundir[-1] is complete and not a constant volume job and self.converged() == False.

           Not converging criteria: >= 10 runs without completion
        """
        if len(self.rundir) >= int(self.settings["run_limit"]):
            return True
        return False


    def run(self):
        """ Perform a series of vasp jobs to relax a structure. Performs a series of vasp calculations until
            convergence is reached according to the criteria in 'status()'. Then performs a final constant volume run
            {"ISIF":2, "ISMEAR":-5, "NSW":0, "IBRION":-1}.
        """

        print("Begin VASP relaxation run")
        sys.stdout.flush()

        # get current status of the relaxation:
        (status, task) = self.status()
        print("\n++  status:", status, "  next task:", task)
        sys.stdout.flush()

        while status == "incomplete":
            if task == "setup":
                self.add_rundir()
                self.setup(self.rundir[-1], self.settings)

            elif task == "relax":
                self.add_rundir()
                casm.vasp.continue_job(self.rundir[-2], self.rundir[-1], self.settings)
                shutil.copyfile(os.path.join(self.relaxdir,"INCAR.base"),os.path.join(self.rundir[-1],"INCAR"))

            elif task == "constant":
                self.add_rundir()
                casm.vasp.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

                # set INCAR to ISIF = 2, ISMEAR = -5, NSW = 0, IBRION = -1
                if (self.settings["final"] != None) and (os.path.isfile(os.path.join(self.relaxdir,self.settings["final"]))):
                    new_values = io.Incar(os.path.join(self.relaxdir, self.settings["final"])).tags
                else:
                    new_values = {"ISIF":2, "ISMEAR":-5, "NSW":0, "IBRION":-1}

                # set INCAR system tag to denote 'final'
                if io.get_incar_tag("SYSTEM", self.rundir[-1]) is None:
                    new_values["SYSTEM"] = "final"
                else:
                    new_values["SYSTEM"] = io.get_incar_tag("SYSTEM", self.rundir[-1]) + " final"

                io.set_incar_tag( new_values, self.rundir[-1])
                print("  Set INCAR tags:", new_values, "\n")
                sys.stdout.flush()

            else:
                # probably hit walltime
                self.add_rundir()
                casm.vasp.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

            while True:
                # run vasp
                result = casm.vasp.run(self.rundir[-1], npar=self.settings["npar"],ncore=self.settings["ncore"],command=self.settings["vasp_cmd"],ncpus=self.settings["ncpus"],kpar=self.settings["kpar"],err_types=self.settings["err_types"])

                # if no errors, continue
                if result is None or self.not_converging():
                    # Check for actions that should be taken after the initial run
                    if len(self.rundir) == 1:
                        if self.settings["fine_ngx"]:
                            outcarfile = os.path.join(self.rundir[-1], "OUTCAR")
                            if not os.path.isfile(outcarfile):
                                # This is an error but I'm not sure what to do about it
                                pass
                            else:
                                init_outcar = io.Outcar(outcarfile)
                                if not init_outcar.complete:
                                    # This is an error but I'm not sure what to do about it
                                    pass
                                elif (init_outcar.ngx is None or
                                      init_outcar.ngy is None or
                                      init_outcar.ngz is None):
                                    # This is an error but I'm not sure what to do about it
                                    pass
                                else:
                                    ng_tags = {
                                        "ngx" : init_outcar.ngx*2,
                                        "ngy" : init_outcar.ngy*2,
                                        "ngz" : init_outcar.ngz*2}
                                    print(ng_tags)
                                    io.set_incar_tag(ng_tags, self.relaxdir, "INCAR.base")
                    break

                # else, attempt to fix first error
                self.add_errdir()
                os.mkdir(self.rundir[-1])
                # self.add_rundir()
                err = result.itervalues().next()

                print("\n++  status:", "error", "  next task:", "fix_error")
                sys.stdout.flush()

                print("Attempting to fix error:", str(err))
                err.fix(self.errdir[-1],self.rundir[-1], self.settings)
                print("")
                sys.stdout.flush()

                if (self.settings["backup"] != None) and len(self.rundir) > 1:
                    print("Restoring from backups:")
                    for f in self.settings["backup"]:
                        if os.path.isfile(os.path.join(self.rundir[-2], f + "_BACKUP.gz")):
                            f_in = gzip.open(os.path.join(self.rundir[-2], f + "_BACKUP.gz", 'rb'))
                            f_out = open(os.path.join(self.rundir[-1], f, 'wb'))
                            f_out.write(f_in.read())
                            f_in.close()
                            f_out.close()
                            print(f, " restored!")
                    sys.stdout.flush()

            (status, task) = self.status()
            print("\n++  status:", status, "  next task:", task)
            sys.stdout.flush()

        if status == "complete":
            if not os.path.isdir(self.finaldir):
                # mv final results to relax.final
                print("mv", os.path.basename(self.rundir[-1]), os.path.basename(self.finaldir))
                sys.stdout.flush()
                os.rename(self.rundir[-1], self.finaldir)
                self.rundir.pop()
                casm.vasp.complete_job(self.finaldir, self.settings)

        return (status, task)


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

        # check if all complete
        if io.job_complete(self.finaldir):
            return ("complete",None)

        # check status of relaxation runs
        self.update_rundir()

        # if not yet started
        if len(self.rundir) == 0:
            return ("incomplete", "setup")

        # if the latest run is complete:
        if io.job_complete(self.rundir[-1]):

            # if it is a final constant volume run
            if io.get_incar_tag("SYSTEM", self.rundir[-1]) != None:
                if io.get_incar_tag("SYSTEM", self.rundir[-1]).split()[-1].strip().lower() == "final":
                # if io.get_incar_tag("ISIF", self.rundir[-1]) == 2 and \
                #    io.get_incar_tag("NSW", self.rundir[-1]) == 0 and \
                #    io.get_incar_tag("ISMEAR", self.rundir[-1]) == -5:
                    return ("complete", None)

            # elif constant volume run (but not the final one)
            if io.get_incar_tag("ISIF", self.rundir[-1]) in [0,1,2]:
                if io.get_incar_tag("NSW", self.rundir[-1]) == len(io.Oszicar(os.path.join(self.rundir[-1],"OSZICAR")).E):
                    return ("incomplete", "relax")      # static run hit NSW limit and so isn't "done"
                else:
                    return ("incomplete", "constant")

            # elif convergence criteria met
            if self.converged():
                return ("incomplete", "constant")

            # elif not converging, return 'not_converging' error
            if self.not_converging():
                return ("not_converging", None)

            # else continue relaxing
            else:
                return ("incomplete", "relax")

        # elif not converging, return 'not_converging' error
        elif self.not_converging():
            return ("not_converging", None)

        # else if the latest run is not complete, continue running it
        return ("incomplete", self.rundir[-1])
