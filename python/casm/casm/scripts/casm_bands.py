from __future__ import absolute_import, division, print_function, unicode_literals

import os
import sys
import argparse

from casm.project.io import read_project_settings
from casm.project import Project, Selection
from casm.vasp.bands import Bands as VaspBand


class CasmBandsError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


"""
casm-bands --run / --submit / --setup / --report
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


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(description='Submit calculations for CASM')
    parser.add_argument('-c', '--configs', help=configs_help, type=str, default="MASTER")
    parser.add_argument('--run', help=run_help, action="store_true", default=False)
    parser.add_argument('--submit', help=submit_help, action="store_true", default=False)
    parser.add_argument('--setup', help=setup_help, action="store_true", default=False)
    args = parser.parse_args(argv)

    args.path = os.getcwd()

    proj = Project(os.path.abspath(args.path))
    sel = Selection(proj, args.configs, all=False)
    if sel.data["configname"] is not None:
        configname = sel.data["configname"][0]
        casm_settings = proj.settings
        if casm_settings is None:
            raise CasmBandsError("Not in a CASM project. The file '.casm' directory was not found.")
    else:
        raise CasmBandsError("Not in CASM project.")

    casm_directories = proj.dir
    print("  Reading relax.json settings file")
    sys.stdout.flush()
    setfile = casm_directories.settings_path_crawl("relax.json", configname, casm_settings.default_clex)

    if setfile is None:
        raise CasmBandsError("Could not find \"relax.json\" in \"settings\" directory")
    else:
        print("Using " + str(setfile) + " as settings...")

    settings = read_project_settings(setfile)
    print("DFT software is:", settings['software'])

    if settings['software'] != 'vasp':
        raise CasmBandsError('This is currently ONLY VASP capable.')

    band_calculator = None

    if args.setup:
        sel.write_pos()
        for configname in sel.data["configname"]:
            if settings['software'] == "quantumespresso":
                raise CasmBandsError('QE not implemented, use VASP')
            elif settings['software'] == "aims":
                raise CasmBandsError('QE not implemented, use VASP')
            elif settings['software'] == 'vasp':
                band_calculator = VaspBand(proj.dir.configuration_dir(configname))
            band_calculator.setup()

    elif args.submit:
        sel.write_pos()
        for configname in sel.data["configname"]:
            if settings['software'] == "quantumespresso":
                raise CasmBandsError('QE not implemented, use VASP')
            elif settings['software'] == "aims":
                raise CasmBandsError('QE not implemented, use VASP')
            elif settings['software'] == 'vasp':
                band_calculator = VaspBand(proj.dir.configuration_dir(configname))
            band_calculator.submit()

    elif args.run:
        sel.write_pos()
        for configname in sel.data["configname"]:
            if settings['software'] == "quantumespresso":
                raise CasmBandsError('QE not implemented, use VASP')
            elif settings['software'] == "aims":
                raise CasmBandsError('QE not implemented, use VASP')
            elif settings['software'] == 'vasp':
                band_calculator = VaspBand(proj.dir.configuration_dir(configname))
            band_calculator.run()


if __name__ == "__main__":
    main()
