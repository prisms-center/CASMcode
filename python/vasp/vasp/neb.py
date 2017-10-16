import os, math, sys, shutil
import vasp
import io

class Neb(object):
    """The class contains all the funtions for performing Nudged Elastic Band(NEB) Calcutations with VASP
        
       Also supports VTST tools to perform Climbing Image NEB."""

    ##write more intro

    def __init__(self, calcdir=None, settings=None):
        """
           Construct a VASP NEB calculation object
        
           Reads the settings from a neb.json """
        ## write more intro
        
        print "Constructing a VASP NEB  calculation object"
        sys.stdout.flush()

        if calcdir is None:
            calcdir = os.getcwd()
        self.calcdir = os.path.abspath(calcdir)
        ## calcdir example : $ROOT/training_data/$(calc_subdir)/$(config_info)/calctype.default/N_images_X

        print "  NEB directory:", self.calcdir
        sys.stdout.flush()

        # find existing .../calcdir/run.run_index directories, store paths in self.rundir list
        self.rundir = []
        self.errdir = [] ## keeps track of the error directories of present run only
        self.update_rundir()
        self.update_errdir() ## redundent but keep it for the flow

        self.finaldir = os.path.join(self.calcdir, "run.final")

        ## setting up the settings options
        if settings is None:
            self.settings = dict()
        else:
            self.settings = settings

        ##set defult values ## To be done ## check them
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
        if not "err_types" in self.settings:
            self.settings["err_types"] = ['SubSpaceMatrixError']

        ## temporary for testing
        if settings is None:
            if not "extra_input_files" in self.settings:
                self.settings["extra_input_files"] = []
            if not "initial" in self.settings:
                self.settings["initial"] = None
            if not "final" in self.settings:
                self.settings["final"] = None
            if not "n_images" in self.settings:
                self.settings["n_images"] = 3

        ##make sure all the defaults are set in self.settings

        print "VASP NEB object constructed\n"
        sys.stdout.flush()


    def add_rundir(self):
        """Make a new run.i directory"""
        os.mkdir(os.path.join(self.calcdir, "run." + str(len(self.rundir))))
        self.update_rundir()
        self.update_errdir()


    def update_rundir(self):
        """Find all .../config/calctype.default/run.i directories, store paths in self.rundir list"""
        self.rundir = []
        run_index = len(self.rundir)
        while os.path.isdir(os.path.join(self.calcdir, "run." + str(run_index))):
            self.rundir.append(os.path.join(self.calcdir, "run." + str(run_index)))
            run_index += 1


    def add_errdir(self):
        """Move run.i to run.i_err.j directory"""
        os.rename(self.rundir[-1], self.rundir[-1] + "_err." + str(len(self.errdir)))
        self.update_errdir()


    def update_errdir(self):
        """Find all .../config/calctype.default/run.i_err.j directories, store paths in self.errdir list"""
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
        ### still has to handle the images

        print "Moving files and directories into initial run directory:", initdir
        initdir = os.path.abspath(initdir)
        for p in os.listdir(self.calcdir):
            if (p in (io.VASP_INPUT_FILE_LIST + self.settings["extra_input_files"])) and (os.path.join(self.calcdir, p) != initdir):
                os.rename(os.path.join(self.calcdir, p), os.path.join(initdir, p))
        ## copying the folders with image poscars into initdir
        for i in range(settings["n_images"]+2):
            folder_name = str(i).zfill(2) #max(2, len(str(settings["n_images"]))+1 )) ##too fancy!!!
            shutil.move(os.path.join(self.calcdir,folder_name),initdir)

        print ""
        sys.stdout.flush()

        # Keep a backup copy of the base INCAR
        shutil.copyfile(os.path.join(initdir, "INCAR"), os.path.join(self.calcdir, "INCAR.base"))
        os.rename(os.path.join(initdir, "POSCAR"), os.path.join(self.calcdir, "POSCAR.start"))
        # If an initial incar is called for, copy it in and set the appropriate flag
        if (self.settings["initial"] != None) and (os.path.isfile(os.path.join(self.calcdir, self.settings["initial"]))):
            new_values = io.Incar(os.path.join(self.calcdir, self.settings["initial"])).tags
            io.set_incar_tag(new_values, initdir)
            print "  Set INCAR tags:", new_values, "\n"
            sys.stdout.flush()


    def complete(self): 
        """Check if the VASP calculation is complete.

           Completion criteria: self.calcdir/self.rundir[-1]/OUTCAR exists and is complete
        """
        ## test if this gives write output when job is still running
        for img in range(1, self.settings["n_images"]+1):
            outcarfile = os.path.join(self.rundir[-1], str(img).zfill(2), "OUTCAR")
            if not os.path.isfile(outcarfile):
                return False
            if not io.Outcar(outcarfile).complete(): ## check this
                return False
        return True


    def converged(self):
        ### might be redundent. presently return complete()
        return self.complete()


    def not_converging(self):
        """Check if configuration is not converging.

           This is called when self.rundir[-1] is not complete

           Not converging criteria: #runs >= self.settings["run_limit"] runs without completion
        """
        if len(self.rundir) >= int(self.settings["run_limit"]):
            return True
        return False


    def run(self):
        """ Perform a series of vasp jobs to run calculations. Performs a series of vasp calculations until
            convergence is reached according to the criteria in 'status()'.
        """

        print "Begin VASP relaxation run"
        sys.stdout.flush()

        # get current status of the relaxation:
        (status, task) = self.status()
        print "\n++  status:", status, "  next task:", task
        sys.stdout.flush()

        while status == "incomplete":
            if task == "setup":
                self.add_rundir()
                self.setup(self.rundir[-1], self.settings)

            elif task == "new_run":
                self.add_rundir()
                # if "n_images" in settings then image CONTCARs will be copied
                vasp.continue_job(self.rundir[-2], self.rundir[-1], self.settings)
                shutil.copyfile(os.path.join(self.calcdir, "INCAR.base"),
                                os.path.join(self.rundir[-1], "INCAR")) ## should it be enforced??
            elif task == "constant":
                self.add_rundir()
                # if "n_images" in settings then image CONTCARs will be copied
                vasp.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

                # set INCAR to ISIF = 2, ISMEAR = -5, NSW = 0, IBRION = -1
                if (self.settings["final"] != None) and (os.path.isfile(os.path.join(self.calcdir, self.settings["final"]))):
                    new_values = io.Incar(os.path.join(self.calcdir, self.settings["final"])).tags
                else:
                    new_values = {"ISIF":2, "ISMEAR":-5, "NSW":0, "IBRION":-1}

                # set INCAR system tag to denote 'final'
                if io.get_incar_tag("SYSTEM", self.rundir[-1]) is None:
                    new_values["SYSTEM"] = "final"
                else:
                    new_values["SYSTEM"] = io.get_incar_tag("SYSTEM", self.rundir[-1]) + " final"

                io.set_incar_tag(new_values, self.rundir[-1])
                print "  Set INCAR tags:", new_values, "\n"
                sys.stdout.flush()

            else: ## redundent
                self.add_rundir()
                vasp.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

            while True:
                # run vasp
                result = vasp.run(self.rundir[-1], stdout="stdout",
                                  npar=self.settings["npar"], ncore=self.settings["ncore"],
                                  command=self.settings["vasp_cmd"], ncpus=self.settings["ncpus"],
                                  kpar=self.settings["kpar"], err_types=self.settings["err_types"],
                                  is_neb=True)

                # if no errors, continue
                if result is None or self.not_converging():
                    break

                # else, attempt to fix first error
                self.add_errdir()
                os.mkdir(self.rundir[-1])
                # self.add_rundir()
                err = result.itervalues().next()

                print "\n++  status:", "error", "  next task:", "fix_error"
                sys.stdout.flush()

                print "Attempting to fix error:", str(err)
                err.fix(self.errdir[-1], self.rundir[-1], self.settings)
                print ""
                sys.stdout.flush()

            (status, task) = self.status()
            print "\n++  status:", status, "  next task:", task
            sys.stdout.flush()

        if status == "complete":
            if not os.path.isdir(self.finaldir):
                # mv final results to relax.final
                print "mv", os.path.basename(self.rundir[-1]), os.path.basename(self.finaldir)
                sys.stdout.flush()
                os.rename(self.rundir[-1], self.finaldir)
                self.rundir.pop()
                vasp.complete_job(self.finaldir, self.settings)

        return (status, task)


    def status(self):
        """ Determine the status of a vasp relaxation series of runs. Individual runs in the series
            are in directories labeled "run.0", "run.1", etc.

            Returns a tuple: (status = "incomplete" or "complete" or "not_converging",
                                task = "setup" or "new_run" or None)

            The first value is the status of the entire relaxation.

            The second value is the current task, where 
            "setup" means to start running calculations in run.0 in the current folder
            "new_run" indicates another calculation job is required as present run is incomplete and need to start a next run.
        """

        # check if all complete
        if io.job_complete(os.path.join(self.finaldir, "01")):
            return ("complete", None)

        # check status of relaxation runs
        self.update_rundir()

        # if not yet started
        if len(self.rundir) == 0:
            return ("incomplete", "setup")

        # check if all complete
        if io.job_complete(os.path.join(self.rundir[-1], "01")):
            # if it is a final constant volume run
            if io.get_incar_tag("SYSTEM", self.rundir[-1]) != None:
                if io.get_incar_tag("SYSTEM", self.rundir[-1]).split()[-1].strip().lower() == "final":
                    return ("complete", None)

            # elif constant volume run (but not the final one)
            if io.get_incar_tag("ISIF", self.rundir[-1]) in [0, 1, 2]:
                if io.get_incar_tag("NSW", self.rundir[-1]) == len(io.Oszicar(os.path.join(self.rundir[-1], "01",  "OSZICAR")).E):
                    return ("incomplete", "new_run")    # static run hit NSW limit and so isn't "done"
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
                return ("incomplete", "new_run")

        # elif not converging, return 'not_converging' error
        elif self.not_converging():
            return ("not_converging", None)

        # else if the latest run is not complete, continue with a new run
        return ("incomplete", 'new_run')
