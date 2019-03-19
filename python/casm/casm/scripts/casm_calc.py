from __future__ import absolute_import, division, print_function, unicode_literals

import argparse
import json
from os import getcwd
from os.path import join, abspath
import sys
import six

from casm.misc import compat, noindent
from casm.project import Project, Selection
import casm.project.io

from casm.qewrapper import Relax as QERelax
from casm.vaspwrapper import Relax as VaspRelax
from casm.aimswrapper.relax import Relax as AimsRelax


class CasmCalcError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


"""
casm-calc --configs selection
          --software "quantumespresso" "vasp" "aims"
          --scheduler "pbs"
          --run / --submit / --setup / --report
"""

configs_help = """
CASM selection file or one of 'CALCULATED', 'ALL', or 'MASTER' (Default)
"""

path_help = """
Path to CASM project. Default=current working directory.
"""

run_help = """
Run calculation for all selected configurations.
"""
method_help = """
Choose what method to use to calculate training data options are vasp or quantumespresso (default="").
Overrides the calculator tag in relax.json. If calculator tag in relax.json is empty then VASP will be used.
"""

submit_help = """
Submit calculation for all selected configurations.
"""

setup_help = """
Setup calculation for all selected configurations.
"""

report_help = """
Report calculation results (print calc.properties.json file) for all selected configurations.
"""


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser(description='Submit calculations for CASM')
    parser.add_argument('-c', '--configs', help=configs_help, type=str, default="MASTER")
    parser.add_argument('--path', help=path_help, type=str, default=None)
    parser.add_argument('-m', '--method', help=method_help, type=str, default=None)
    parser.add_argument('--run', help=run_help, action="store_true", default=False)
    parser.add_argument('--submit', help=submit_help, action="store_true", default=False)
    parser.add_argument('--setup', help=setup_help, action="store_true", default=False)
    parser.add_argument('--report', help=report_help, action="store_true", default=False)
    args = parser.parse_args(argv)
  
    if args.path is None:
        args.path = getcwd()
  
    try:
        proj = Project(abspath(args.path))
        sel = Selection(proj, args.configs, all=False)
        if sel.data["configname"] is not None:
            configname = sel.data["configname"][0]
            casm_settings = proj.settings
            if casm_settings is None:
                raise CasmCalcError("Not in a CASM project. The file '.casm' directory was not found.")
        else:
            raise CasmCalcError("Not in CASM project.")

        casm_directories = proj.dir
        print("  Reading relax.json settings file")
        sys.stdout.flush()
        setfile = casm_directories.settings_path_crawl("relax.json", configname, casm_settings.default_clex)

        if setfile is None:
            raise casm.qewrapper.QEWrapperError("Could not find \"relax.json\" in \"settings\" directory")
        else:
            print("Using " + str(setfile) + " as settings...")

        settings = casm.project.io.read_project_settings(setfile)
        print("DFT software is:", settings['software'])

        relaxation = None

        if args.setup:
            sel.write_pos()
            for configname in sel.data["configname"]:
                if settings['software'] == "quantumespresso":
                    relaxation = QERelax(proj.dir.configuration_dir(configname))
                elif settings['software'] == "aims":
                    relaxation = AimsRelax(proj.dir.configuration_dir(configname))
                elif settings['software'] == 'vasp':
                    relaxation = VaspRelax(proj.dir.configuration_dir(configname))
                relaxation.setup()
    
        elif args.submit:
            sel.write_pos()
            for configname in sel.data["configname"]:
                if settings['software'] == "quantumespresso":
                    relaxation = QERelax(proj.dir.configuration_dir(configname))
                elif settings['software'] == "aims":
                    relaxation = AimsRelax(proj.dir.configuration_dir(configname))
                elif settings['software'] == 'vasp':
                    relaxation = VaspRelax(proj.dir.configuration_dir(configname))
                relaxation.submit()

        elif args.run:
            sel.write_pos()
            for configname in sel.data["configname"]:
                if settings['software'] == "quantumespresso":
                    relaxation = QERelax(proj.dir.configuration_dir(configname))
                elif settings['software'] == "aims":
                    relaxation = AimsRelax(proj.dir.configuration_dir(configname))
                elif settings['software'] == 'vasp':
                    relaxation = VaspRelax(proj.dir.configuration_dir(configname))
                relaxation.run()
    
        elif args.report:
            for configname in sel.data["configname"]:
                configdir = proj.dir.configuration_dir(configname)
                clex = proj.settings.default_clex
                calcdir = proj.dir.calctype_dir(configname, clex)
                finaldir = join(calcdir, "run.final")
                try:
                    if settings['software'] == "quantumespresso":
                        if settings["outfilename"] is None:
                            print("WARNING: No output file specified in relax.json using std.out")
                            settings["outfilename"] = "std.out"
                            outfilename = settings["outfilename"]
                            output = casm.qewrapper.Relax.properties(finaldir, outfilename)
                        else:
                            output = VaspRelax.properties(finaldir)
                except IOError:
                    print(("Unable to report properties for directory {}.\n"
                           "Please verify that it contains a completed calculation.".format(configdir)))

                calc_props = proj.dir.calculated_properties(configname, clex)
                print("writing:", calc_props)
                # TODO: Fix this line to work properly
                compat.dump(json, output, calc_props, 'w', cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)

    except ValueError as e:
        raise CasmCalcError(e)

if __name__ == "__main__":
    main()
