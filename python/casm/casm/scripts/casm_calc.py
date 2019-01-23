from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import argparse
import json
from os import getcwd
from os.path import join, abspath
import sys
import six

from casm.misc import compat, noindent
from casm.project import Project, Selection
import casm.qewrapper
from casm.vaspwrapper import Relax

# casm-calc --configs selection
#           --software "quantumespresso" "vasp"
#           --scheduler "pbs"
#           --run / --submit / --setup / --report

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

def main(argv = None):
  if argv is None:
    argv = sys.argv[1:]
    
  parser = argparse.ArgumentParser(description = 'Submit calculations for CASM')
  parser.add_argument('-c', '--configs', help=configs_help, type=str, default="MASTER")
  parser.add_argument('--path', help=path_help, type=str, default=None)
  parser.add_argument('-m','--method', help=method_help, type=str, default=None)
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
      configname=sel.data["configname"][0]
      casm_settings=proj.settings
      if casm_settings == None:
        raise casm.qewrapper.QEWrapperError("Not in a CASM project. The file '.casm' directory was not found.")
      casm_directories=proj.dir
      print("  Reading relax.json settings file")
      sys.stdout.flush()
      setfile = casm_directories.settings_path_crawl("relax.json",configname,casm_settings.default_clex)
      if setfile == None:
          raise casm.qewrapper.QEWrapperError("Could not find \"relax.json\" in an appropriate \"settings\" directory")
          sys.stdout.flush()
      else:
          print("Using "+str(setfile)+" as settings...")
      settings = casm.qewrapper.read_settings(setfile)
      if settings["software"] is None:
        settings["software"]="vasp"
      software=settings["software"]
      print("Relevant software is:", software)
    if args.setup:
      sel.write_pos()
      for configname in sel.data["configname"]:
        if software == "quantumespresso":
          relaxation = casm.qewrapper.Relax(proj.dir.configuration_dir(configname))
        else:
          relaxation = Relax(proj.dir.configuration_dir(configname))
        relaxation.setup()
    
    elif args.submit:
      sel.write_pos()
      for configname in sel.data["configname"]:
        if software == "quantumespresso":
          relaxation = casm.qewrapper.Relax(proj.dir.configuration_dir(configname))
        else:
          relaxation = Relax(proj.dir.configuration_dir(configname))
        relaxation.submit()
    
    elif args.run:
      sel.write_pos()
      for configname in sel.data["configname"]:
        if software == "quantumespresso":
          relaxation = casm.qewrapper.Relax(proj.dir.configuration_dir(configname))
        else:
          relaxation = Relax(proj.dir.configuration_dir(configname))
        relaxation.run()
    
    elif args.report:
      for configname in sel.data["configname"]:
        configdir = proj.dir.configuration_dir(configname)
        clex = proj.settings.default_clex
        calcdir = proj.dir.calctype_dir(configname, clex)
        finaldir = join(calcdir, "run.final")
        try:
          if software == "quantumespresso":
            if settings["outfilename"] is None:
                print("WARNING: No output file specified in relax.json using default outfilename of std.out")
                settings["outfilename"]="std.out"
            outfilename = settings["outfilename"]
            output = casm.qewrapper.Relax.properties(finaldir,outfilename)
          else:
            output = Relax.properties(finaldir)
          calc_props = proj.dir.calculated_properties(configname, clex)
          print("writing:", calc_props)
          with open(calc_props, 'wb') as f:
              f.write(six.u(json.dumps(output, cls=noindent.NoIndentEncoder, indent=2)).encode('utf-8'))
        #compat.dump(json, output, calc_props, 'w', cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)
        except:
          print(("Unable to report properties for directory {}.\n" 
                "Please verify that it contains a completed calculation.".format(configdir)))
  except Exception as e:
    print(e)
    sys.exit(1)

if __name__ == "__main__":
  main()
