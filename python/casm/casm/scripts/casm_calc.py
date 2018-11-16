from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import argparse
import json
from os import getcwd
from os.path import join, abspath
import sys

from casm import qewrapper, questwrapper, vaspwrapper
from casm.misc import compat, noindent
from casm.project import Project, Selection

# casm-calc --configs selection
#           --type "config", "diff_trans", etc.
#           --calctype "default"
#           --setup /  --run / --submit / --report

configs_help = """
CASM selection file or one of 'CALCULATED', 'ALL', or 'MASTER' (Default)
"""

configtype_help = """
Type of configuartions 'config' (Default), 'diff_trans', 'diff_trans_config' or scel
"""

calctype_help = """
calctype on which to run the calculations
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

available_calculators = {
  "vasp":{
    "relax": vaspwrapper.Relax,
    "neb": vaspwrapper.Neb
  },
  "quantumexpresso":{
    "relax": qewrapper.Relax
  },
  "seqquest":{
    "relax": questwrapper.Relax
  }
}

def main(argv = None):
    if argv is None:
        argv = sys.argv[1:]
      
    parser = argparse.ArgumentParser(description = 'Submit calculations for CASM')
    parser.add_argument('-c', '--configs', help=configs_help, type=str, default="MASTER")
    parser.add_argument('-t', '--type', help=configtype_help, type=str, default="config")
    parser.add_argument('--calctype', help=calctype_help, type=str, default="")
    parser.add_argument('--path', help=path_help, type=str, default=None)
    parser.add_argument('--run', help=run_help, action="store_true", default=False)
    parser.add_argument('--submit', help=submit_help, action="store_true", default=False)
    parser.add_argument('--setup', help=setup_help, action="store_true", default=False)
    parser.add_argument('--report', help=report_help, action="store_true", default=False)
    args = parser.parse_args(argv)
    
    if args.path is None:
        args.path = getcwd()
    
    try:
        proj = Project(abspath(args.path))
        sel = Selection(proj, args.configs, args.type, all=False)
        if args.calctype == "":
            #get default calctype
            args.calctype = proj.settings.default_clex.calctype
    
        global_settings = json.load(open(join(proj.dir.calctype_settings_dir(args.calctype), "calc.json")))
        software = global_settings["software"]
        method = global_settings["method"]
        # Construct with Selection:
        # - This provides access to the Project, via sel.proj
        # - From the project you can make calls to run interpolation and query lattice relaxations
        calculator = available_calculators[software][method](sel, args.calctype)
    
        if args.setup:
            calculator.setup()
    
        elif args.submit:
            calculator.submit()
    
        elif args.run:
            calculator.run()
    
    except Exception as e:
      print(e)
      sys.exit(1)

if __name__ == "__main__":
  main()
