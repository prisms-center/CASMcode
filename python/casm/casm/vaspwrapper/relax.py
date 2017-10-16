import os, math, sys, json, re, warnings, shutil
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
from casm.vaspwrapper import VaspCalculatorBase
import vasp.Relax

class Relax(VaspCalculatorBase):
    """The Relax class contains functions for setting up, executing, and parsing a VASP relaxation.

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
        """
        Construct a VASP neb job object.
        """
        print "Construct a casm.vaspwrapper.Relax instance:"
        VaspCalculatorBase.__init__(selection, calctype, auto, sort)
        self.calculator = vasp.relax

    @classmethod
    def relax(cls, configdir, calctype, auto=True, sort=True):
        sel = Selection.selection_from_configs([configdir])
        obj = cls(sel, calctype, auto, sort)
        return obj

    def pre_setup(self):
        self.selection.write_pos()

    def setup(self):
        """ Setup initial relaxation run for the selection
        """
        super(Relax, self).setup()

    def submit(self):
        super(Relax, self).submit()

    def run(self):
        super(Relax, self).run()

    def report(self):
        super(Relax, self).report()

    @staticmethod
    def run_cmd(configdir, calctype):
        return "python -c \"import casm.vaspwrapper; obj = casm.vaspwrapper.Relax.relax('{0}', '{1}'); obj.run()\"\n".format(configdir, calctype)

    def finalize(self, config_data, super_poscarfile=None):
        if super_poscarfile is None:
            super_poscarfile = os.path.join(config_data["configdir"], "POS")
        super(Relax, self).finalize(config_data, super_poscarfile)
        sys.stdout.flush()

    def properties(self, calcdir, super_poscarfile=None, speciesfile=None):
        """Make properties output as a list of dict of each image properties"""
        output = super(Relax, self).properties(calcdir, super_poscarfile, speciesfile)
        return output
