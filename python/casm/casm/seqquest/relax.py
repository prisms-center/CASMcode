""" Relax class """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import re
# import math
import sys
import shutil
from . import seqquest
from . import seqquest_io

class RelaxError(Exception):
    """ Relax-specific errors """
    pass

class Relax(object):
    """The Relax class contains functions for setting up, executing, and parsing a seqquest calc.

        The relaxation is initialized in a directory containing seqquest files, called 'relaxdir'.
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
        Construct a seqquest relaxation job object.

        Args:
            configdir:  path to seqquest relaxation directory
            settings:   dictionary-like object containing settings, or if None, it reads
                        the json file: .../relaxdir/relax.json

                possible settings keys are:
                    used by seqquest.run() function:
                        "ncpus": number of ncpus to run mpi on
                        "run_cmd": (default, see seqquest.run) shell command to execute seqquest,
                            or None to use default mpirun
                    used by not_converging():
                        "run_limit": (default 10) maximum number of runs to allow before setting
                            status to "not_converging"

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
        if "run_cmd" not in self.settings:
            self.settings["vasp_cmd"] = None
        if "ncpus" not in self.settings:
            self.settings["ncpus"] = None
        if "run_limit" not in self.settings:
            self.settings["run_limit"] = 10
        if "compress" not in self.settings:
            self.settings["compress"] = []
        if "err_types" not in self.settings:
            self.settings["err_types"] = ['SubSpaceMatrixError']

        print("SeqQuest Relax object constructed\n")
        sys.stdout.flush()


    def add_rundir(self):
        """Make a new run.i directory"""
        os.mkdir(os.path.join(self.relaxdir, "run." + str(len(self.rundir))))
        self.update_rundir()
        self.update_errdir()


    def update_rundir(self):
        """Find all .../config/calctype.CALC/relax/run.i directories, store paths in
            self.rundir list"""
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
        """Find all .../config/calctype.CALC/relax/run.i_err.j directories, store paths in
            self.errdir list"""
        self.errdir = []
        if len(self.rundir) == 0:
            pass
        else:
            err_index = len(self.errdir)
            while os.path.isdir(self.rundir[-1] + "_err." + str(err_index)):
                self.errdir.append(self.rundir[-1] + "_err." + str(err_index))
                err_index += 1


    def setup(self, initdir, _):
        """ mv all files and directories (besides initdir) into initdir """

        print("Moving files into initial run directory:", initdir)
        initdir = os.path.abspath(initdir)
        for ifile in os.listdir(self.relaxdir):
            if ((ifile in seqquest_io.QUEST_INPUT_FILE_LIST + self.settings["extra_input_files"])
                    and (os.path.join(self.relaxdir, ifile) != initdir)):
                os.rename(os.path.join(self.relaxdir, ifile), os.path.join(initdir, ifile))
            elif re.search(r".*\.atm", ifile):
                if os.path.islink(os.path.join(self.relaxdir, ifile)):
                    os.symlink(os.readlink(os.path.join(self.relaxdir, ifile)), os.path.join(initdir, ifile))
                else:
                    shutil.copy(os.path.join(self.relaxdir, ifile), os.path.join(initdir, ifile))
            elif ifile == "lcao.geom_in":
                os.rename(os.path.join(self.relaxdir, ifile), os.path.join(initdir, ifile))

        print("")
        sys.stdout.flush()

        # Keep a backup copy of the base lcao.in
        shutil.copyfile(os.path.join(initdir, "lcao.in"),
                        os.path.join(self.relaxdir, "lcao.in.base"))

        # Initial lcao.in files are NOT yet supported!
        if self.settings["initial"] != None:
            raise RelaxError("the 'initial' keyword in relax/converge.json is not yet supported for\
                             SeqQuest!")

        # # If an initial incar is called for, copy it in and set the appropriate flag
        # if ((self.settings["initial"] != None) and
        #       (os.path.isfile(os.path.join(self.relaxdir,self.settings["initial"])))):
        #     new_values = io.Incar(os.path.join(self.relaxdir,self.settings["initial"])).tags
        #     io.set_incar_tag(new_values, initdir)
        #     print "  Set INCAR tags:", new_values, "\n"
        #     sys.stdout.flush()

    def complete(self):
        """Check if the seqquest relaxation is complete.

           Completion criteria: .../config/calctype.CALC/relax/run.final/lcao.out exists
            and is complete
        """
        outfile = os.path.join(self.finaldir, "lcao.out")
        if not os.path.isfile(outfile):
            return False
        if not seqquest_io.LcaoOUT(outfile).complete:
            return False
        return True

    def converged(self):
    # def converged(self):
        """Check if configuration is relaxed.

           For now, this always returns True, as SeqQuest doesn't have the anisotropic basis
           set problem that VASP has.
        """
        outfile = os.path.join(self.rundir[-1], "lcao.out")
        return seqquest_io.LcaoOUT(outfile).complete

           # This is called when self.rundir[-1] is complete and not a constant volume job.

           # Convergence criteria is: at least 2 relaxation jobs are complete, and:
           #                          1) the last job completed with <= 3 ionic steps
           #                          or 2) the last two jobs had final E0 differ by less than
           #                                self.settings["nrg_convergence"]
        #"""
        # if len(self.rundir) >= 2:
        #     if io.ionic_steps(self.rundir[-1]) <= 3:
        #         return True
        #     if self.settings["nrg_convergence"] != None:
        #         if io.job_complete(self.rundir[-1]) and io.job_complete(self.rundir[-2]):
        #             o1 = io.Oszicar(os.path.join(self.rundir[-1],"OSZICAR"))
        #             o2 = io.Oszicar(os.path.join(self.rundir[-2],"OSZICAR"))
        #             if abs( o1.E[-1] - o2.E[-1]) < self.settings["nrg_convergence"]:
        #                 return True

        # return False


    def not_converging(self):
        """Check if configuration is not converging.

           This is called when self.rundir[-1] is complete and not a constant volume job and
            self.converged() == False.

           Not converging criteria: >= 10 runs without completion
        """
        if len(self.rundir) >= int(self.settings["run_limit"]):
            return True
        outfile = os.path.join(self.rundir[-1], "lcao.out")
        if not seqquest_io.LcaoOUT(outfile).converged:
            # This indicates that we have an electronic step problem, which, for now, we do NOT know how to solve
            return True
        return False


    def run(self):  #pylint: disable=too-many-branches, too-many-statements
        """ Perform a seqquest jobs to relax a structure. Then performs a final run, using settings
            in lcao.in file pointed to by "final" if present.
        """

        print("Begin SeqQuest relaxation run")
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
                seqquest.continue_job(self.rundir[-2], self.rundir[-1], self.settings)
                shutil.copyfile(os.path.join(self.relaxdir, "lcao.in.base"),
                                os.path.join(self.rundir[-1], "lcao.in"))

            elif task == "constant":
                self.add_rundir()
                seqquest.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

                # set INCAR to ISIF = 2, ISMEAR = -5, NSW = 0, IBRION = -1
                # Implement final settings, if present
                # Final settings not yet implemented for seqquest!!!
                if self.settings["final"] != None:
                    raise RelaxError("the 'final' keyword in relax/converge.json is not yet\
                                     supported for SeqQuest!")
                else:
                    new_values = dict()
                # if (self.settings["final"] != None) and (os.path.isfile(os.path.join(self.relaxdir
                #                                                           ,self.settings["final"]
                #                                                           ))):
                #     new_values = io.Incar(os.path.join(self.relaxdir,
                #                                        self.settings["final"])).tags
                # else:
                #     new_values = {"ISIF":2, "ISMEAR":-5, "NSW":0, "IBRION":-1}

                # set lcao.in notes system tag to denote 'final'
                final_lcao = seqquest_io.LcaoIN(os.path.join(self.rundir[-1], "lcao.in"))
                # notes_value = seqquest_io.get_lcao_tag(['setup', 'notes'], self.rundir[-1])
                if final_lcao.get(['setup', 'notes']) is None:
                    final_lcao.set(['setup', 'notes'], ["final"])
                else:
                    final_lcao.set(['setup', 'notes'], final_lcao.get(['setup', 'notes'])
                                   + ["final"])
                new_values[('setup', 'notes')] = ["final"]

                final_lcao.write(os.path.join(self.rundir[-1], "lcao.in"))

                print("  Set lcao.in tags:", new_values, "\n")
                sys.stdout.flush()

            else:
                # probably hit walltime
                self.add_rundir()
                seqquest.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

            while True:
                # run vasp
                # result = seqquest.run(self.rundir[-1], command=self.settings["run_cmd"],
                #           ncpus=self.settings["ncpus"],err_types=self.settings["err_types"])
                result = seqquest.run(self.rundir[-1], command=self.settings["run_cmd"],    #pylint: disable=assignment-from-none
                                      ncpus=self.settings["ncpus"],
                                      err_types=self.settings["err_types"])

                # if no errors, continue
                if result is None or self.not_converging():
                    break

                # else, attempt to fix first error
                self.add_errdir()
                os.mkdir(self.rundir[-1])
                # self.add_rundir()
                err = result.itervalues().next()

                print("\n++  status:", "error", "  next task:", "fix_error")
                sys.stdout.flush()

                print("Attempting to fix error:", str(err))
                err.fix(self.errdir[-1], self.rundir[-1], self.settings)
                print("")
                sys.stdout.flush()

                # backup files are NOT yet supported!
                if self.settings["backup"] != None:
                    raise RelaxError("the 'backup' keyword in relax/converge.json is not yet\
                                      supported for SeqQuest!")
                # if (self.settings["backup"] != None) and len(self.rundir) > 1:
                #     print "Restoring from backups:"
                #     for f in self.settings["backup"]:
                #         if os.path.isfile(os.path.join(self.rundir[-2], f + "_BACKUP.gz")):
                #            f_in = gzip.open(os.path.join(self.rundir[-2], f + "_BACKUP.gz", 'rb'))
                #             f_out = open(os.path.join(self.rundir[-1], f, 'wb'))
                #             f_out.write(f_in.read())
                #             f_in.close()
                #             f_out.close()
                #             print f, " restored!"
                #     sys.stdout.flush()

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
                seqquest.complete_job(self.finaldir, self.settings)

        return (status, task)


    def status(self):   #pylint: disable=too-many-return-statements
        """ Determine the status of a quest relaxation series of runs. Individual runs in the series
            are in directories labeled "run.0", "run.1", etc.

            Returns a tuple: (status = "incomplete" or "complete" or "not_converging",
                                task = continuedir or "relax" or "constant" or None)

            The first value is the status of the entire relaxation.

            The second value is the current task, where 'continuedir' is the path to a
            quest job directory that is not yet completed, "relax" indicates another
            volume relaxation job is required, and "constant" that a constant volume run is required
        """

        # check if all complete
        if seqquest_io.job_complete(self.finaldir):
            return ("complete", None)

        # check status of relaxation runs
        self.update_rundir()

        # if not yet started
        if len(self.rundir) == 0:
            return ("incomplete", "setup")

        # if the latest run is complete:
        if seqquest_io.job_complete(self.rundir[-1]):

            # if it is a final constant volume run
            final_lcao = seqquest_io.LcaoIN(os.path.join(self.rundir[-1], "lcao.in"))
            # notes_value = seqquest_io.get_lcao_tag(['setup', 'notes'], self.rundir[-1])
            if final_lcao.get(['setup', 'notes']) is None:
                pass
            else:
                if "final" in final_lcao.get(['setup', 'notes']):
                    return ("complete", None)

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

