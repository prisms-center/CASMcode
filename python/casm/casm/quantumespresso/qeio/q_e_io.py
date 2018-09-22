from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import re
import shutil
import six
import sys
from casm.quantumespresso.qeio import infile, outfile, species, poscar

#VASP_INPUT_FILE_LIST = ["INCAR", "STOPCAR", "POTCAR", "KPOINTS", "POSCAR",\
#                        "EXHCAR", "CHGCAR", "WAVECAR", "TMPCAR"]

DEFAULT_QE_MOVE_LIST = []

DEFAULT_QE_COPY_LIST = []

DEFAULT_QE_REMOVE_LIST = []


class QuantumEspressoIOError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


def job_complete(outfilename,jobdir=None):
    """Return True if Quantum Espresso job at path 'jobdir' is complete"""
    if jobdir == None:
        jobdir = os.getcwd()

    myoutfile = os.path.join(jobdir, outfilename)
    if (not os.path.isfile(myoutfile)) and (not os.path.isfile(myoutfile+".gz")):
        return False
    if outfile.Outfile(myoutfile).complete:
        return True
    return False


def get_infile_tag(key, infilename,jobdir=None):
    """Opens Infile in 'jobdir' and returns 'key' value."""
    if jobdir == None:
        jobdir = os.getcwd()
    tinfile = infile.Infile(os.path.join(jobdir,infilename))
    if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['CONTROL']:
        if ("CONTROL" in tinfile.namelists.keys()) and (key in tinfile.namelists["CONTROL"].tags.keys()):
            return tinfile.namelists["CONTROL"].tags[key]
    if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['SYSTEM']:
        if ("SYSTEM" in tinfile.namelists.keys()) and (key in tinfile.namelists["SYSTEM"].tags.keys()):
            return tinfile.namelists["SYSTEM"].tags[key]
    if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['ELECTRONS']:
        if ("ELECTRONS" in tinfile.namelists.keys()) and (key in tinfile.namelists["ELECTRONS"].tags.keys()):
            return tinfile.namelists["ELECTRONS"].tags[key]
    if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['IONS']:
        if ("IONS" in tinfile.namelists.keys()) and (key in tinfile.namelists["IONS"].tags.keys()):
            return tinfile.namelists["IONS"].tags[key]
    if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['CELL']:
        if ("CELL" in tinfile.namelists.keys()) and (key in tinfile.namelists["CELL"].tags.keys()):
            return tinfile.namelists["CELL"].tags[key]    
    return None


def set_infile_tag(tag_dict,infilename,jobdir=None):
    """Opens Infile in 'jobdir', sets 'key' value, and writes Infile
        If 'val' == None, the tag is removed from the Infile.
        DOES NOT TYPE CHECK THE TAGS IN TAG_DICT
    """
    if jobdir == None:
        jobdir = os.getcwd()
    myinfile = os.path.join(jobdir,infilename)
    tinfile = infile.Infile(myinfile)

    for key, val in six.iteritems(tag_dict):
        if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['CONTROL']:
            if ("CONTROL" in tinfile.namelists.keys()):
                if (val == None) or (str(val).strip() == ""):
                    if (key in tinfile.namelists["CONTROL"].tags.keys()):
                        del tinfile.namelists["CONTROL"].tags[key]
                else:
                    tinfile.namelists["CONTROL"].tags[key] = val
        if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['SYSTEM']:
            if ("SYSTEM" in tinfile.namelists.keys()):
                if (val == None) or (str(val).strip() == ""):
                    if (key in tinfile.namelists["SYSTEM"].tags.keys()):
                        del tinfile.namelists["SYSTEM"].tags[key]
                else:
                    tinfile.namelists["SYSTEM"].tags[key] = val
        if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['ELECTRONS']:
            if ("ELECTRONS" in tinfile.namelists.keys()):
                if (val == None) or (str(val).strip() == ""):
                    if (key in tinfile.namelists["ELECTRONS"].tags.keys()):
                        del tinfile.namelists["ELECTRONS"].tags[key]
                else:
                    tinfile.namelists["ELECTRONS"].tags[key] = val
        if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['IONS']:
            if ("IONS" in tinfile.namelists.keys()):
                if (val == None) or (str(val).strip() == ""):
                    if (key in tinfile.namelists["IONS"].tags.keys()):
                        del tinfile.namelists["IONS"].tags[key]
                else:
                    tinfile.namelists["IONS"].tags[key] = val
        if key in infile.QUANTUM_ESPRESSO_TOTAL_LIST['CELL']:
            if ("CELL" in tinfile.namelists.keys()):
                if (val == None) or (str(val).strip() == ""):
                    if (key in tinfile.namelists["CELL"].tags.keys()):
                        del tinfile.namelists["CELL"].tags[key]
                else:
                    tinfile.namelists["CELL"].tags[key] = val
    tinfile.write(myinfile)


def ionic_steps(outfilename,jobdir=None):
    """Find the number of ionic steps completed in 'jobdir'"""
    if jobdir == None:
        jobdir = os.getcwd()
    try:
        toutfile = outfile.Outfile(os.path.join(jobdir,outfilename))
        return len(toutfile.E)
    except:
        raise QuantumEspressoIOError("Could not read number of ionic steps from " + os.path.join(jobdir,outfilename))


def write_quantum_espresso_input(dirpath, infilename, super_poscarfile, speciesfile, sort=True, extra_input_files=[], strict_kpoints=False):
    """ Write Quantum Espresso input files in directory 'dirpath' """
    print("Setting up Quantum Espresso input files:", dirpath)


    # read Infile for K_POINTS, CELL_PARAMETERS, and generate super kpoints
    print("  Reading Infile:", infilename)
    myinfile = infile.Infile(infilename, None, None, sort)
    print("  Reading SPECIES:", speciesfile)
    species_settings = species.species_settings(speciesfile)
    print("  Reading K_POINTS from:", infilename)
    prim_kpoints = myinfile.cards["K_POINTS"]
    print("  Reading K_POINTS reference positions from CELL_PARAMETERS:", infilename)
    prim = poscar.Poscar(infilename,species_settings)
    print("  Reading supercell POS:", super_poscarfile)
    super = poscar.Poscar(super_poscarfile, species_settings)
    print("  Generating supercell KPOINTS")
    if strict_kpoints:
        super_kpoints = prim_kpoints
    else:
        super_kpoints = prim_kpoints.super_kpoints(prim, super)

    #get raw infile name (without directories)
    directories=re.split("/",infilename)
    if len(directories) > 1:
        infilename=directories[-1]

    # write main input file
    print("  Writing supercell positions to Infile")
    myinfile.rewrite_poscar_info(super,species_settings)
    print("  Writing supercell K_POINTS to Infile")
    myinfile.cards["K_POINTS"]=super_kpoints
    print("  Writing Infile:", os.path.join(dirpath,infilename))
    myinfile.write(os.path.join(dirpath,infilename))

    # copy extra input files
    print("  Copying extra input files", end=' ')
    for s in extra_input_files:
        print(s, end=' ')
        shutil.copy(s,dirpath)
    print("")

    print("Quantum Espresso input files complete\n")
    sys.stdout.flush()


