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
        self.errdir = []
        self.update_rundir() ##define these functions
        self.update_errdir()
        
        ## setting up the settings options
        if settings is None:
            self.settings = dict()
        else:
            self.settings = settings

        ##set defult values ## To be done     
        ##make sure all the defau;ts are set in self.settings

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

    """ Required funtions in the class
    
    --> **setup(self,calcdir,settings,diffusion_properties??) #sets up the folders for the images ## Also handles the errordirs if implemented
    --> **run(self): ## setup the job and runs the VASP NEB calculation
    ** - Important to do

    other functions
    --> status(self): ## checks the status of the job
    --> funtions to keep track of status 1) complete() 2)converged()??
        
    do we need error directories ## do it when working on errors
    --> def add_errdir(self):
    --> def update_errdir(self):

    To Do:
    1) Edit io/kpoints.py file for new INCAR options
    2) Edit io/vaspio.py and make a write_neb_vasp_input(INPUT)
    3) dig up all new errors possible for NEB calcs
    4) start editing neb.py at casm/vaspwrapper/ 
    5) make vasp.neb, vasp.neb_setup executables along with casm-calc support to run the calculation

    Things to discuss
      1) possible input formats for the difussion hop settings
      2) possible input option for various VASP options. Should we link to calctype and use a "settings/calctype.neb/" to read all the files. """
