from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, math, sys, shutil, gzip
import casm.quantumespresso
from casm.quantumespresso import qeio

class Relax(object):
    """The Relax class contains functions for setting up, executing, and parsing a Quantum Espresso relaxation.

        The relaxation is initialized in a directory containing Quantum Espresso input files, called 'relaxdir'.
        It then creates the following directory structure:
        .../relaxdir/
            run.0/
            run.1/
            ...
            run.final/



        'run.i' directories are only created when ready.
        'run.final' is a final constant volume run 

        Contains:
            self.relaxdir  (.../relax)
            self.rundir    (list of .../relax/run.i)
            self.finaldir (.../relax/run.final)
    """
    def __init__(self, relaxdir=None, settings=None):
        """
        Construct a Quantum Espresso relaxation job object.

        Args:
            relaxdir:  path to quantum espresso relaxation directory
            settings:   dictionary-like object containing settings, or if None, it reads
                        the json file: .../relaxdir/relax.json

                possible settings keys are:
                    used by casm.quantumespresso.run() function:
                        "ncpus": number of ncpus to run mpi on
			"npar" or "ncore": number of ways to parallelize
                        "kpar": number of ways to parallelize k-points
                        "qe_cmd": (default, see quantumespresso.run) shell command to execute quantumespresso, or None to use default mpirun
                        "strict_kpoint": force strict copying of KPOINTS file, otherwise kpoints are scaled based on supercell size
                    used by not_converging():
                        "run_limit": (default 10) maximum number of runs to allow before setting status to "not_converging"

        """

        print("Constructing a Quantum Espresso Relax object")
        sys.stdout.flush()

        # store path to .../relaxdir, and create if not existing
        if relaxdir == None:
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

        if settings == None:
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
        if not "qe_cmd" in self.settings:
            self.settings["qe_cmd"] = None
        if not "ncpus" in self.settings:
            self.settings["ncpus"] = None
        if not "run_limit" in self.settings:
            self.settings["run_limit"] = 10
        if not "nrg_convergence" in self.settings:
            self.settings["nrg_convergence"] = None
        if not "compress" in self.settings:
            self.settings["compress"] = []
        if not "err_types" in self.settings:
            self.settings["err_types"] = ['SubSpaceMatrixError']
        ## added these because of key errors
        if not "extra_input_files" in self.settings:
            self.settings["extra_input_files"] = []
        if not "move" in self.settings:
            self.settings["move"] = []
        if not "copy" in self.settings:
            self.settings["copy"] = []
        if not "remove" in self.settings:
            self.settings["remove"] = []
        if not "compress" in self.settings:
            self.settings["compress"] = []
        if not "backup" in self.settings:
            self.settings["backup"] = []
        if not "initial" in self.settings:
            self.settings["initial"] = None
        if not "final" in self.settings:
            self.settings["final"] = None

        print("Quantum Espresso Relax object constructed\n")
        sys.stdout.flush()


    def add_rundir(self):
        """Make a new run.i directory"""
        os.mkdir(os.path.join(self.relaxdir, "run." + str(len(self.rundir))))
        self.update_rundir()
        self.update_errdir()


    def update_rundir(self):
        """Find all .../config/qe/relax/run.i directories, store paths in self.rundir list"""
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
        """Find all .../config/qe/relax/run.i_err.j directories, store paths in self.errdir list"""
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

        infilename=settings["infilename"]

        print("Moving files into initial run directory:", initdir)
        initdir = os.path.abspath(initdir)
        for p in os.listdir(self.relaxdir):
            if (p in ([infilename] + self.settings["extra_input_files"])) and (os.path.join(self.relaxdir, p) != initdir):
                os.rename(os.path.join(self.relaxdir,p), os.path.join(initdir,p))
        print("")
        sys.stdout.flush()

        # Keep a backup copy of the base Infile
        shutil.copyfile(os.path.join(initdir,infilename),os.path.join(self.relaxdir,infilename + ".base"))

        # If an initial infile is called for, copy it in and set the appropriate flag
        if (self.settings["initial"] != None) and (os.path.isfile(os.path.join(self.relaxdir,self.settings["initial"]))):
            new_values = qeio.Infile(os.path.join(self.relaxdir,self.settings["initial"])).tags
            qeio.set_infile_tag(new_values,infilename,initdir)
            print("  Set Infile tags:", new_values, "\n")
            sys.stdout.flush()

    def complete(self):
        """Check if the Quantum Espresso relaxation is complete.

           Completion criteria: .../config/qe/relax/run.final/outfilename exists and is complete
        """
        outfilename= self.settings["outfilename"]
        myoutfile = os.path.join(self.finaldir,outfilename)
        if not os.path.isfile(myoutfile):
            return False
        if not qeio.Outfile(myoutfile).complete():
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
        outfilename=self.settings["outfilename"]
        if len(self.rundir) >= 2:
            if qeio.ionic_steps(outfilename,self.rundir[-1]) <= 3:
                return True
            if self.settings["nrg_convergence"] != None:
                if qeio.job_complete(outfilename,self.rundir[-1]) and qeio.job_complete(outfilename,self.rundir[-2]):
                    o1 = qeio.Outfile(os.path.join(self.rundir[-1],outfilename))
                    o2 = qeio.Outfile(os.path.join(self.rundir[-2],outfilename))
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
        """ Perform a series of quantum espresso jobs to relax a structure. Performs a series of quantum espresso calculations until
            convergence is reached according to the criteria in 'status()'. Then performs a final constant volume run
        """

        print("Begin Quantum Espresso relaxation run")
        sys.stdout.flush()

        # get current status of the relaxation:
        (status, task) = self.status()
        print("\n++  status:", status, "  next task:", task)
        sys.stdout.flush()

        infilename=self.settings["infilename"]
        outfilename=self.settings["outfilename"]

        while status == "incomplete":
            if task == "setup":
                self.add_rundir()
                self.setup(self.rundir[-1], self.settings)

            elif task == "relax":
                self.add_rundir()
                casm.quantumespresso.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

            elif task == "constant":
                self.add_rundir()
                casm.quantumespresso.continue_job(self.rundir[-2], self.rundir[-1],self.settings)

                # set Infile to calculation = relax
                if (self.settings["final"] != None) and (os.path.isfile(os.path.join(self.relaxdir,self.settings["final"]))):
                    new_values = qeio.Infile(os.path.join(self.relaxdir, self.settings["final"])).tags
                else:
                    new_values = {"calculation":"'relax'"}

                # set Infile title tag to denote 'final'
                if qeio.get_infile_tag("title",infilename, self.rundir[-1]) is None:
                    new_values["title"] = "'final'"
                else:
                    new_values["title"] = "'{}'".format(qeio.get_infile_tag("title", self.rundir[-1]) + " final")

                qeio.set_infile_tag(new_values,infilename,self.rundir[-1])
                print("  Set Infile tags:", new_values, "\n")
                sys.stdout.flush()

            else:
                # probably hit walltime
                self.add_rundir()
                casm.quantumespresso.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

            while True:
                # run quantum espresso
                result = casm.quantumespresso.run(infilename,outfilename,self.rundir[-1],command=self.settings["qe_cmd"],ncpus=self.settings["ncpus"],err_types=self.settings["err_types"])

                # if no errors, continue
                if result == None or self.not_converging():
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
                casm.quantumespresso.complete_job(self.finaldir, self.settings)

        return (status, task)


    def status(self):
        """ Determine the status of a quantum espresso relaxation series of runs. Individual runs in the series
            are in directories labeled "run.0", "run.1", etc.

            Returns a tuple: (status = "incomplete" or "complete" or "not_converging",
                                task = continuedir or "relax" or "constant" or None)

            The first value is the status of the entire relaxation.

            The second value is the current task, where 'continuedir' is the path to a
            quantum espresso job directory that is not yet completed, "relax" indicates another
            volume relaxation job is required, and "constant" that a constant volume run is required.
        """
        infilename=self.settings["infilename"]
        outfilename=self.settings["outfilename"]
        # check if all complete
        if qeio.job_complete(outfilename,self.finaldir):
            return ("complete",None)

        # check status of relaxation runs
        self.update_rundir()
        # if not yet started
        if len(self.rundir) == 0:
            return ("incomplete", "setup")

        # if the latest run is complete:
        if qeio.job_complete(outfilename,self.rundir[-1]):
            # if it is a final constant volume run
            if qeio.get_infile_tag("title", infilename, self.rundir[-1]) != None:
                print(qeio.get_infile_tag("title", infilename, self.rundir[-1])[1:-1].split()[-1].strip().lower())
                if qeio.get_infile_tag("title", infilename, self.rundir[-1])[1:-1].split()[-1].strip().lower() == "final":
                    return ("complete", None)

            # elif constant volume run (but not the final one)
            if qeio.get_infile_tag("calculation", infilename, self.rundir[-1]) in ['relax','scf','nscf']:
                if qeio.get_infile_tag("nsteps", infilename, self.rundir[-1]) == len(qeio.Outfile(os.path.join(self.rundir[-1],outfilename)).E):
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
