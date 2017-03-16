import os, math, sys, json, re, warnings
import pbs
import vasp
import casm
import casm.project
import vaspwrapper

class Relax(object):
    """The Relax class contains functions for setting up, executing, and parsing a VASP relaxation.

        The relaxation creates the following directory structure:
        config/
          calctype.name/
              run.0/
              ....

        'run.i' directories are only created when ready.

        This automatically looks for VASP settings files using:
          casm.project.DirectoryStructure.settings_path_crawl

    Attributes
    ----------

      casm_settings: casm.project.ProjectSettings instance
        CASM project settings

      casm_directories: casm.project.DirectoryStructure instance
        CASM project directory hierarchy

      settings: dict
        Settings for pbs and the relaxation, see vaspwrapper.read_settings

      configdir: str
        Directory where configuration results are stored. The result of:
          casm.project.DirectoryStructure.configuration_dir(self.configname)

      configname: str
        The name of the configuration to be calculated

      auto: boolean
        True if using pbs module's JobDB to manage pbs jobs

      sort: boolean
        True if sorting atoms in POSCAR by type

      clex: casm.project.ClexDescription instance
        The cluster expansion being worked on. Used for the 'calctype' settings.
        Currently, fixed to self.casm_settings.default_clex.

    """
    def __init__(self, configdir=None, auto=True, sort=True):
        """
        Construct a VASP relaxation job object.

        Arguments
        ----------

            configdir: str, optional, default=None
              Path to configuration directory. If None, uses the current working
              directory

            auto: boolean, optional, default=True,
              Use True to use the pbs module's JobDB to manage pbs jobs

            sort: boolean, optional, default=True,
              Use True to sort atoms in POSCAR by type

        """
        print "Construct a casm.vaspwrapper.Relax instance:"

        if configdir is None:
            configdir = os.getcwd()
        print "  Input directory:", configdir

        # get the configname from the configdir path
        #_res = os.path.split(configdir)
        _res = configdir.strip().split('/') 
        self.configname = _res[-2] + "/" + _res[-1]
        tmp_index = _res.index("training_data")
        self.run_subdir = ""
        ## anything between "training_data" and "configname" is dumped into run_subdir
        if temp_index < len(_res)-3:
            for i in range(tmp_index+1,len(_res)-2)]:
                self.run_subdir += _res[i] + "/"
            self.run_subdir = run_subdir[:-1] ##remove the trailing "/"
                
        print "  Configuration:", self.configname

        print "  Reading CASM settings"
        self.casm_directories=casm.project.DirectoryStructure(configdir)
        self.casm_settings = casm.project.ProjectSettings(configdir)
        if self.casm_settings is None:
            raise vaspwrapper.VaspWrapperError("Not in a CASM project. The file '.casm' directory was not found.")

        if os.path.abspath(configdir) != self.configdir:
            print ""
            print "input configdir:", configdir
            print "determined configname:", self.configname
            print "expected configdir given configname:", self.configdir
            raise vaspwrapper.VaspWrapperError("Mismatch between configname and configdir")

        # fixed to default_clex for now
        self.clex = self.casm_settings.default_clex

        # store path to .../config/calctype.name, and create if not existing
        self.calcdir = self.casm_directories.calctype_dir(self.configname, self.clex)
        try:
            os.mkdir(self.calcdir)
        except:
            pass
        print "  Calculations directory:", self.calcdir

        # read the settings json file
        print "  Reading relax.json settings file"
        sys.stdout.flush()
        setfile = self.casm_directories.settings_path_crawl("relax.json", self.configname, self.clex)

        if setfile is None:
            raise vaspwrapper.VaspWrapperError("Could not find \"relax.json\" in an appropriate \"settings\" directory")
            sys.stdout.flush()

        else:
            print "  Read settings from:", setfile
        self.settings = vaspwrapper.read_settings(setfile)

        # set default settings if not present
        if not "ncore" in self.settings:
            self.settings["ncore"] = None
        if not "npar" in self.settings:
            self.settings["npar"] = None
        if not "kpar" in self.settings:
            self.settings["kpar"] = None
        if not "vasp_cmd" in self.settings:
            self.settings["vasp_cmd"] = None
        if not "ncpus" in self.settings:
            self.settings["ncpus"] = None
        if not "run_limit" in self.settings:
            self.settings["run_limit"] = None
        if not "prerun" in self.settings:
            self.settings["prerun"] = None
        if not "postrun" in self.settings:
            self.settings["postrun"] = None

        self.auto = auto
        self.sort = sort
        print "  DONE\n"
        sys.stdout.flush()


        """
        To Do:
        1) Find the best way to edit DirectoryStructure to perform calculations in KMC directory struc


        """
