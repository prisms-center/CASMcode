import os, math, sys, json, re, warnings, shutil
import pbs
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper


class ConfigPropertiesBase(object):
    """ The Object holds all the properties that relate to a configuration

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

    def __init__(self, configname, calctype=None):
        self.configname = configname
        self.casm_directories = casm.project.DirectoryStructure(self.configdir)
        self.casm_settings = casm.project.ProjectSettings(self.configdir)
        if self.casm_settings is None:
            raise vaspwrapper.VaspWrapperError("Not in a CASM project. The file '.casm' directory was not found.")

        # fixed to default_clex for now
        self.clex = self.casm_settings.default_clex
        if calctype:
            self.clex.calctype = calctype

        # store path to .../config/calctype.name, and create if not existing
        # will be appended by n_images at the end after reading the settings file
        self.calcdir = self.casm_directories.calctype_dir(self.configname, self.clex, self.calc_subdir)
        self.results_subdir = '' #everything between $(calcdir)/run.*/ and OSZICAR and OUTCAR files
        setfile = self.casm_directories.settings_path_crawl("calc.json", self.configname, self.clex, self.calc_subdir)
        if setfile is None:
            raise vaspwrapper.VaspWrapperError("Could not find \"calc.json\" in an appropriate \"settings\" directory")
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

    @property
    def configdir(self): ## check the path
        return self.casm_directories.configuration_dir(self.configname, self.calc_subdir)

    @property
    def calc_subdir(self):
        # everything between training_data and configname is dumpded in calc_subdir #redundent
        return ""
