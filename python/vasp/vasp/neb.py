import os, math, sys, shutil, gzip
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) ##sets unbuffered output
import vasp, io

class Neb(object):
    """The class contains all the funtions for performing Nudged Elastic Band(NEB) Calcutations with VASP
        
       Also supports VTST tools to perform Climbing Image NEB."""

    def __init(self, calcdir = None, settings = None):
        """
           Construct a VASP NEB calculation object
        
           Reads the settings from a neb.json """
           
        print "Constructing a VASP NEB  calculation object"
        sys.stdout.flush()

        if calcdir is None:
            calcdir = os.getcwd()
        self.calcdir = os.path.abspath(calcdir)
        ## calcdir example : $ROOT/training_data/$(config_info)/calctype.default/

        print "  NEB directory:", self.calcdir
        sys.stdout.flush()

        # find existing .../calcdir/run.run_index directories, store paths in self.rundir list
        self.rundir = []
        self.errdir = [] ## keeps track of the error directories of present run
        self.update_rundir() 
        self.update_errdir() ## redundent but keep it for the flow
        
        ## setting up the settings options
        if settings is None:
            self.settings = dict()
        else:
            self.settings = settings

        ##set defult values ## To be done     
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
        while os.path.isdir( os.path.join(self.calcdir, "run." + str(run_index))):
                self.rundir.append( os.path.join(self.calcdir, "run." + str(run_index)) )
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

        print "Moving files into initial run directory:", initdir
        initdir = os.path.abspath(initdir)
        for p in os.listdir(self.calcdir):
            if (p in (io.VASP_INPUT_FILE_LIST + self.settings["extra_input_files"])) and (os.path.join(self.calcdir, p) != initdir):
                os.rename(os.path.join(self.calcdir,p), os.path.join(initdir,p))
        print ""
        sys.stdout.flush()

        # Keep a backup copy of the base INCAR
        shutil.copyfile(os.path.join(initdir,"INCAR"),os.path.join(self.calcdir,"INCAR.base"))

        # If an initial incar is called for, copy it in and set the appropriate flag
        if (self.settings["initial"] != None) and (os.path.isfile(os.path.join(self.calcdir,self.settings["initial"]))):
            new_values = io.Incar(os.path.join(self.calcdir,self.settings["initial"])).tags
            io.set_incar_tag(new_values, initdir)
            print "  Set INCAR tags:", new_values, "\n"
            sys.stdout.flush()


    def complete(self): 
        """Check if the VASP calculation is complete.

           Completion criteria: self.calcdir/self.rundir[-1]/OUTCAR exists and is complete
        """
        ## test if this gives write output when job is still running 
        outcarfile = os.path.join(self.rundir[-1],"OUTCAR")
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
                vasp.continue_job(self.rundir[-2], self.rundir[-1], self.settings)
                shutil.copyfile(os.path.join(self.relaxdir,"INCAR.base"),os.path.join(self.rundir[-1],"INCAR")) ## should it be enforced??

            else: ## redundent
                # probably hit walltime
                self.add_rundir()
                vasp.continue_job(self.rundir[-2], self.rundir[-1], self.settings)

            while True:
                # run vasp
                result = vasp.run(self.rundir[-1], npar=self.settings["npar"],ncore=self.settings["ncore"],command=self.settings["vasp_cmd"],
                                  ncpus=self.settings["ncpus"],kpar=self.settings["kpar"],err_types=self.settings["err_types"])

                # if no errors, continue
                if result is None or self.not_converging():
                    # # Check for actions that should be taken after the initial run ## is it required
                    # if len(self.rundir) == 1:
                    #     if self.settings["fine_ngx"]:
                    #         outcarfile = os.path.join(self.rundir[-1], "OUTCAR")
                    #         if not os.path.isfile(outcarfile):
                    #             # This is an error but I'm not sure what to do about it
                    #             pass
                    #         else:
                    #             init_outcar = io.Outcar(outcarfile)
                    #             if not init_outcar.complete:
                    #                 # This is an error but I'm not sure what to do about it
                    #                 pass
                    #             elif (init_outcar.ngx is None or
                    #                   init_outcar.ngy is None or
                    #                   init_outcar.ngz is None):
                    #                 # This is an error but I'm not sure what to do about it
                    #                 pass
                    #             else:
                    #                 ng_tags = {
                    #                     "ngx" : init_outcar.ngx*2,
                    #                     "ngy" : init_outcar.ngy*2,
                    #                     "ngz" : init_outcar.ngz*2}
                    #                 print ng_tags
                    #                 io.set_incar_tag(ng_tags, self.relaxdir, "INCAR.base")
                    break

                # else, attempt to fix first error
                self.add_errdir()
                os.mkdir(self.rundir[-1]) ##why??
                # self.add_rundir()
                err = result.itervalues().next()

                print "\n++  status:", "error", "  next task:", "fix_error"
                sys.stdout.flush()

                print "Attempting to fix error:", str(err)
                err.fix(self.errdir[-1],self.rundir[-1], self.settings)
                print ""
                sys.stdout.flush()

                # if (self.settings["backup"] != None) and len(self.rundir) > 1: ##is it required?
                #     print "Restoring from backups:"
                #     for f in self.settings["backup"]:
                #         if os.path.isfile(os.path.join(self.rundir[-2], f + "_BACKUP.gz")):
                #             f_in = gzip.open(os.path.join(self.rundir[-2], f + "_BACKUP.gz", 'rb'))
                #             f_out = open(os.path.join(self.rundir[-1], f, 'wb'))
                #             f_out.write(f_in.read())
                #             f_in.close()
                #             f_out.close()
                #             print f, " restored!"
                #     sys.stdout.flush()

            (status, task) = self.status()
            print "\n++  status:", status, "  next task:", task
            sys.stdout.flush()

        # if status == "complete":
        #     if not os.path.isdir(self.finaldir):
        #         # mv final results to relax.final
        #         print "mv", os.path.basename(self.rundir[-1]), os.path.basename(self.finaldir)
        #         sys.stdout.flush()
        #         os.rename(self.rundir[-1], self.finaldir)
        #         self.rundir.pop()
        #         vasp.complete_job(self.finaldir, self.settings)

        return (status, task)


    def status(self):
        """ Determine the status of a vasp relaxation series of runs. Individual runs in the series
            are in directories labeled "run.0", "run.1", etc.

            Returns a tuple: (status = "incomplete" or "complete" or "not_converging",
                                task = "setup" or "new_run" or None)

            The first value is the status of the entire relaxation.

            The second value is the current task, where 
            "setup" means to start running calculations in run.0 in the current folder
            "new_run" indicates another calculation job is required as present run is incomplete and start a next run.
        """

        # check status of relaxation runs
        self.update_rundir()

        # check if all complete
        if io.job_complete(self.rundir[-1]):
            return ("complete",None)

        # if not yet started
        if len(self.rundir) == 0:
            return ("incomplete", "setup")

        # elif not converging, return 'not_converging' error
        if self.not_converging():
            return ("not_converging", None)

        # else if the latest run is not complete, continue with a new run
        return ("incomplete", 'new_run')



    """ Required funtions in the class
    
    --> **setup(self,calcdir,settings) #still has to set up the folders for the images 
    --> **run(self): ## clean up some feature. fine_nxg and backup
    ** - Important to do

    other functions
    --> status(self): ## checks the status of the job ## may be more status options
        
    To Do:
    1) Edit io/kpoints.py file for new INCAR options
    2) Edit io/vaspio.py and make a write_neb_vasp_input(INPUT)
    3) dig up all new errors possible for NEB calcs
    4) start editing neb.py at casm/vaspwrapper/ 
    5) make vasp.neb, vasp.neb_setup executables along with casm-calc support to run the calculation

    Things to discuss
      1) n_images vs intermediate_pos
      2) position of atoms in pos or information on which atoms moved.
      3) make it executable independently ina any config folder. so where should the pos be? and should report already present non neb calc??
