"""Defines the neb module methods"""

import os
import sys
import json
import shutil
import vasp
import casm
import pandas
from casm.project import Project, Selection
import vaspwrapper
from casm.vaspwrapper import VaspCalculatorBase
import vasp.Relax

class Relax(VaspCalculatorBase):
    """
    The Neb class contains functions for setting up, executing, and parsing a VASP neb calculation.

    Attributes
    ----------
    selection : casm.project.Selection
        selection of configuration
    calctype : string
        calctype to setup and run the neb calculations
    auto : bool
    sort : bool

    Methods
    -------
    from_configuration_dir(configuration_dir='string', calctype='string', bool, bool)
        returns a instance of the Neb class instantited with a single configuration
    config_properties(config_data=dict/Panda.DataFrame)
        return a dict of the properties required to setup a configuration
    pre_setup
        creates folder and makes POS files for each image
    setup
        sets up the input vasp files for the selection
    config_setup
        sets up the input vasp files for a single configuration
    get_vasp_input_files(config_data=dict/Pandas.DataFrame, settings=dict)
        returns filenames of a vasp neb calculation
    submit
        submit a job for each configuration
    run
        runs the neb calcutation on the selection
    report
        reports results for the selection
    run_cmd(configdir='string', calctype='string')
        return a string of command to run a single configuration
    finalize(config_data=dict/pandas_data, super_poscarfile='string')
        checks convergnce and write a properties file for the selection
    properties(calcdir='string', super_poscarfile='string', speciesfile='string')
        return a dict containing all the relaxed properties for a configuration

    Notes
    -----
    The Relax class contains functions for setting up, executing, and parsing a VASP relaxation.

    The relaxation creates the following directory structure:
    config/
        calctype.name/
            run.0/
            ....

    'run.i' directories are only created when ready.

    This automatically looks for VASP settings files using:
        casm.project.DirectoryStructure.settings_path_crawl

    """
    def __init__(self, selection, calctype=None, auto=True, sort=True):
        """Construct a VASP neb job object"""
        print "Construct a casm.vaspwrapper.Relax instance:"
        VaspCalculatorBase.__init__(self, selection, calctype, auto, sort)
        self.calculator = vasp.relax

    @classmethod
    def from_configuration_dir(cls, configuration_dir, calctype, auto=True, sort=True):
        """returns a instance of the Neb class instantited with a single configuration"""
        # change config_dir to configuration_dir all over
        proj = Project(configuration_dir)
        sel = Selection(proj, "EMPTY", "config", False)
        split_path = configuration_dir.split(os.path.sep)
        index = split_path.index("training_data")
        configname = '/'.join(split_path[index+1:])
        sel.data = pandas.DataFrame({"configname":configname, "selected":1}, index=range(2))
        sel_config = sel.saveas(os.path.join(proj.path, ".casm/tmp", configname.replace('/', '.')), True)
        obj = cls(sel_config, calctype, auto, sort)
        return obj

    def pre_setup(self):
        """Setus up folders and writes POS files"""
        self.selection.write_pos()

    def setup(self):
        """Setup initial relaxation run for the selection"""
        super(Relax, self).setup()

    def submit(self):
        """Submit a job for each configuration"""
        super(Relax, self).submit()

    def run(self):
        """Runs the neb calcutation on the selection"""
        super(Relax, self).run()

    def report(self):
        """Reports results for the selection"""
        super(Relax, self).report()

    @staticmethod
    def run_cmd(configdir, calctype):
        """Return a string of command to run a single configuration"""
        return "python -c \"import casm.vaspwrapper; obj = casm.vaspwrapper.Relax.from_configuration_dir('{0}', '{1}'); obj.run()\"\n".format(configdir,
                                                                                                                                              calctype)

    def finalize(self, config_data, super_poscarfile=None):
        """Checks convergnce and write a properties file for the selection"""
        if super_poscarfile is None:
            super_poscarfile = os.path.join(config_data["configdir"], "POS")
        super(Relax, self).finalize(config_data, super_poscarfile)
        sys.stdout.flush()

    def properties(self, calcdir, super_poscarfile=None, speciesfile=None):
        """Make properties output as a list of dict of each image properties"""
        output = super(Relax, self).properties(calcdir, super_poscarfile, speciesfile)
        return output
