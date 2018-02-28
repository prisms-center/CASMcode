from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import re
import shutil
import six
import sys
from casm.vasp.io import incar, kpoints, oszicar, outcar, species, poscar

VASP_INPUT_FILE_LIST = ["INCAR", "STOPCAR", "POTCAR", "KPOINTS", "POSCAR",\
                        "EXHCAR", "CHGCAR", "WAVECAR", "TMPCAR"]

DEFAULT_VASP_MOVE_LIST = ["POTCAR"]

DEFAULT_VASP_COPY_LIST = ["INCAR", "KPOINTS"]

DEFAULT_VASP_REMOVE_LIST = ["IBZKPT", "CHG", "CHGCAR", "WAVECAR", "TMPCAR", "EIGENVAL",\
                            "DOSCAR", "PROCAR", "PCDAT", "XDATCAR", "LOCPOT", "ELFCAR", "PROOUT"]


class VaspIOError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def job_complete(jobdir=None):
    """Return True if vasp job at path 'jobdir' is complete"""
    if jobdir is None:
        jobdir = os.getcwd()
    outcarfile = os.path.join(jobdir, "OUTCAR")
    if (not os.path.isfile(outcarfile)) and (not os.path.isfile(outcarfile+".gz")):
        return False
    if outcar.Outcar(outcarfile).complete:
        return True
    return False


def get_incar_tag(key, jobdir=None):
    """Opens INCAR in 'jobdir' and returns 'key' value."""
    if jobdir is None:
        jobdir = os.getcwd()
    tincar = incar.Incar(os.path.join(jobdir,"INCAR"))
    for k in tincar.tags:
        if key.lower() == k.lower():
            return tincar.tags[k]
    return None


def set_incar_tag(tag_dict,jobdir=None, name=None):
    """Opens INCAR in 'jobdir', sets 'key' value, and writes INCAR
        If 'val' is None, the tag is removed from the INCAR.
    """
    if name is None:
        name = "INCAR"
    if jobdir is None:
        jobdir = os.getcwd()
    incarfile = os.path.join(jobdir,name)
    tincar = incar.Incar(incarfile)

    for key, val in six.iteritems(tag_dict):
        for k in tincar.tags:
            if key.lower() == k.lower():
                if (val is None) or (str(val).strip() == ""):
                    del tincar.tags[k]
                else:
                    tincar.tags[k] = val
                break

        if val != None and str(val).strip() != "":
            tincar.tags[key] = val

    tincar.write(incarfile)


def ionic_steps(jobdir=None):
    """Find the number of ionic steps completed in 'jobdir'"""
    try:
        toszicar = oszicar.Oszicar(os.path.join(jobdir,"OSZICAR"))
        return len(toszicar.E)
    except:
        raise VaspIOError("Could not read number of ionic steps from " + os.path.join(jobdir,"OSZICAR"))


def write_potcar(filename, poscar, species, sort=True):
    """ Write an appropriate POTCAR """
    if sort == False:
        with open(filename,'w') as file:
            for name in poscar.type_atoms:
                with open( os.path.join(species[name].potcardir,'POTCAR')) as potcar:
                    file.write( potcar.read())
    else:
        # dict: key = alias, value = list of Sites
        pos = poscar.basis_dict()

        with open(filename,'w') as file:
            # for each alias
            for alias in sorted(pos.keys()):
                # find matching IndividualSpecies with write_potcar == True
                for name in species:
                    if species[name].alias == alias and species[name].write_potcar:
                        # add to POTCAR file
                        with open( os.path.join(species[name].potcardir,'POTCAR')) as potcar:
                            file.write( potcar.read())
                        break

def write_stopcar(mode='e', jobdir=None):
    """ Write STOPCAR file with two modes:
        mode = 'e' for 'VASP stops at the next electronic step'
        mode = 'i' for 'VASP stops at the next ionic step' """
    if jobdir is None:
        jobdir = os.getcwd()
    if mode.lower()[0] == 'e':
        stop_string = "LABORT = .TRUE."
    elif mode.lower()[0] == 'i':
        stop_string = "LSTOP = .TRUE."
    else:
        raise VaspIOError("Invalid STOPCAR mode specified: " + str(mode))

    filename = os.path.join(jobdir,'STOPCAR')

    try:
        stopcar_write = open(filename,'w')
    except IOError as e:
        raise e

    stopcar_write.write(stop_string)
    stopcar_write.close()



def write_vasp_input(dirpath, incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, sort=True, extra_input_files=[], strict_kpoints=False):
    """ Write VASP input files in directory 'dirpath' """
    print("Setting up VASP input files:", dirpath)

    # read prim and prim kpoints
    print("  Reading KPOINTS:", prim_kpointsfile)
    prim_kpoints = kpoints.Kpoints(prim_kpointsfile)
    if prim_poscarfile != None:
        print("  Reading KPOINTS reference POSCAR:", prim_poscarfile)
        prim = poscar.Poscar(prim_poscarfile)
    else:
        prim = None

    # read species, super poscar, incar, and generate super kpoints
    print("  Reading SPECIES:", speciesfile)
    species_settings = species.species_settings(speciesfile)
    print("  Reading supercell POS:", super_poscarfile)
    super = poscar.Poscar(super_poscarfile, species_settings)
    print("  Reading INCAR:", incarfile)
    super_incar = incar.Incar(incarfile, species_settings, super, sort)
    print("  Generating supercell KPOINTS")
    if strict_kpoints:
        super_kpoints = prim_kpoints
    else:
        super_kpoints = prim_kpoints.super_kpoints(prim, super)


    # write main input files
    print("  Writing supercell POSCAR:", os.path.join(dirpath,'POSCAR'))
    super.write(os.path.join(dirpath,'POSCAR'), sort)
    print("  Writing INCAR:", os.path.join(dirpath,'INCAR'))
    super_incar.write(os.path.join(dirpath,'INCAR'))
    print("  Writing supercell KPOINTS:", os.path.join(dirpath,'KPOINTS'))
    super_kpoints.write(os.path.join(dirpath,'KPOINTS'))
    print("  Writing POTCAR:", os.path.join(dirpath,'POTCAR'))
    write_potcar(os.path.join(dirpath,'POTCAR'), super, species_settings, sort)

    # copy extra input files
    if len(extra_input_files):
        print("  Copying extra input files", end=' ')
    for s in extra_input_files:
        print("    ", s)
        shutil.copy(s,dirpath)

    print("  DONE\n")
    sys.stdout.flush()


