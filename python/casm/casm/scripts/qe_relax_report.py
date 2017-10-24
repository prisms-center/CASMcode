from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import os
import sys

from casm.misc import compat, noindent
import casm.project
import casm.qewrapper

def main():
    print("Begin qe.relax.report")

    if len(sys.argv) != 2:
        print("Usage: qe.relax.report configdir")
        sys.exit()

    configdir = sys.argv[1]

    if configdir == None:
        configdir = os.getcwd()

    configdir=os.path.abspath(configdir)
    _res = os.path.split(configdir)
    configname = os.path.split(_res[0])[1] + "/" + _res[1]

    casm_settings=casm.project.ProjectSettings()
    if casm_settings == None:
        raise casm.qewrapper.QEWrapperError("Not in a CASM project. The file '.casm' directory was not found.")

    casm_directories=casm.project.DirectoryStructure()

    calcdir=casm_directories.calctype_dir(configname,casm_settings.default_clex)

    print("  Reading relax.json settings file")
    sys.stdout.flush()
    setfile = casm_directories.settings_path_crawl("relax.json",configname,casm_settings.default_clex)

    if setfile == None:
        raise casm.qewrapper.QEWrapperError("Could not find \"relax.json\" in an appropriate \"settings\" directory")
        sys.stdout.flush()

    else:
        print("Using "+str(setfile)+" as settings...")

    settings = casm.qewrapper.read_settings(setfile)

    if not "outfilename" in settings:
        print("WARNING: No output file specified in relax.json using default outfilename of std.out")
        settings["outfilename"]="std.out"

    outfilename = settings["outfilename"]

    qedir=os.path.join(configdir,calcdir,"run.final")

    try:
        output = casm.qewrapper.Relax.properties(qedir,outfilename)
    except:
        print(("Unable to report properties for directory %s. Please verify that it contains a completed Quantum Espresso calculation with name %s" %(qedir, outfilename)))
        raise

    compat.dump(json, output, 'properties.calc.json', 'w', cls=noindent.NoIndentEncoder, 
                indent=4, sort_keys=True)

    print("Finish qe.relax.report")

if __name__ == "__main__":
    main()
