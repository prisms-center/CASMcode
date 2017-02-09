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

        print "  NEB directory:", self.calcdir
        sys.stdout.flush()

        ## setting up the settings options
        if settings is None:
            self.settings = dict()
        else:
            self.settings = settings

        ##set defult values ## To be done
        
        print "VASP NEB object constructed\n"
        sys.stdout.flush()


    """ Required funtions in the class
    
    --> setup(self,calcdir,settings,diffusion_properties??) #sets up the folders for the images ## Also handles the errordirs if implemented
    --> run(self): ## setup the job and runs the VASP NEB calculation

    other functions
    --> status(self): ## checks the status of the job
    --> funtions to keep track of status 1) complete() 2)converged()??
        
    do we need error directories ## do it when working on errors
    --> def add_errdir(self):
    --> def update_errdir(self):

    To Do:
    1) Edit io/* files for new INCAR options
    2) dig up all new error s possible for NEB calcs
    3) make and vasp.neb, vasp.neb_setup executables along with casm-calc support to run the calculation

    Things to discuss
      1) possible input formats for the difussion hop settings
      2) possible input option for various VASP options. Should we link to calctype and use a "settings/calctype.neb/" to read all the files. """
