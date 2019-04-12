from __future__ import absolute_import, division, print_function, unicode_literals

import os
import sys
import argparse

from casm.project.io import read_project_settings
from casm.project import Project, Selection
from casm.feff import Feff


class CasmFeffError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


"""
casm-feff --run / --submit / --setup / --plot
"""

configs_help = """
CASM selection file or one of 'CALCULATED', 'ALL', or 'MASTER' (Default)
"""

run_help = """
Run calculation for all selected configurations.
"""

submit_help = """
Submit calculation for all selected configurations.
"""

setup_help = """
Setup calculation for all selected configurations.
"""

plot_help = """
Plots XANES for selected configurations.
"""


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(description='Submit calculations for CASM')
    parser.add_argument('-c', '--configs', help=configs_help, type=str, default="MASTER")
    parser.add_argument('--run', help=run_help, action="store_true", default=False)
    parser.add_argument('--submit', help=submit_help, action="store_true", default=False)
    parser.add_argument('--setup', help=setup_help, action="store_true", default=False)
    parser.add_argument('--plot', help=plot_help, action="store_true", default=False)
    args = parser.parse_args(argv)

    args.path = os.getcwd()

    proj = Project(os.path.abspath(args.path))
    sel = Selection(proj, args.configs, all=False)
    if sel.data["configname"] is not None:
        configname = sel.data["configname"][0]
        casm_settings = proj.settings
        if casm_settings is None:
            raise CasmFeffError("Not in a CASM project. The file '.casm' directory was not found.")
    else:
        raise CasmFeffError("Not in CASM project.")

    casm_directories = proj.dir
    print("  Reading relax.json settings file")
    sys.stdout.flush()

    setfile = casm_directories.settings_path_crawl("relax.json", configname, casm_settings.default_clex)
    if setfile is None:
        raise CasmFeffError("Could not find \"relax.json\" in \"settings\" directory")
    else:
        print("Using " + str(setfile) + " as computation settings...")

    fefffile = casm_directories.settings_path_crawl("bands.json", configname, casm_settings.default_clex)
    if fefffile is None:
        raise CasmFeffError("Could not find \"feff.json\" in \"settings\" directory")
    else:
        print("Using " + str(fefffile) + " for FEFF settings...")

    settings = read_project_settings(setfile)
    print("DFT software is:", settings['software'])

    if settings['software'] != 'vasp':
        raise CasmFeffError('This is currently ONLY VASP capable.')

    feff_calculator = None

    if args.setup:
        sel.write_pos()
        for configname in sel.data["configname"]:
            if settings['software'] == "quantumespresso":
                raise CasmFeffError('QE not implemented, use VASP')
            elif settings['software'] == "aims":
                raise CasmFeffError('FHI-aims not implemented, use VASP')
            elif settings['software'] == 'vasp':
                feff_calculator = Feff(proj.dir.configuration_dir(configname))
            feff_calculator.setup()

    elif args.submit:
        sel.write_pos()
        for configname in sel.data["configname"]:
            if settings['software'] == "quantumespresso":
                raise CasmFeffError('QE not implemented, use VASP')
            elif settings['software'] == "aims":
                raise CasmFeffError('FHI-aims not implemented, use VASP')
            elif settings['software'] == 'vasp':
                feff_calculator = Feff(proj.dir.configuration_dir(configname))
            feff_calculator.submit()

    elif args.run:
        sel.write_pos()
        for configname in sel.data["configname"]:
            if settings['software'] == "quantumespresso":
                raise CasmFeffError('QE not implemented, use VASP')
            elif settings['software'] == "aims":
                raise CasmFeffError('FHI-aims not implemented, use VASP')
            elif settings['software'] == 'vasp':
                feff_calculator = Feff(proj.dir.configuration_dir(configname))
            feff_calculator.run()

    elif args.plot:
        sel.write_pos()
        for configname in sel.data["configname"]:
            if settings['software'] == "quantumespresso":
                raise CasmFeffError('QE not implemented, use VASP')
            elif settings['software'] == "aims":
                raise CasmFeffError('FHI-aims not implemented, use VASP')
            elif settings['software'] == 'vasp':
                feff_calculator = Feff(proj.dir.configuration_dir(configname))
            feff_calculator.plot_bfeff(plot_dir=os.path.abspath(os.path.join(proj.dir.configuration_dir(configname),
                                                                              'calctype.default', 'band_structure')))


if __name__ == "__main__":
    main()
