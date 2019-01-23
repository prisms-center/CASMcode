import os
import sys

from casm.aims.io.steps import Steps
from casm.aims.io.basis import basis_settings, write_basis
from casm.aims.io.parser import Parser
from casm.aims.io.control import Control
from casm.aims.io.kpoints import Kpoints
from casm.aims.io.geometry import Geometry

AIMS_INPUT_FILE_LIST = ["control.in", "geometry.in"]

DEFAULT_AIMS_COPY_LIST = []
DEFAULT_AIMS_MOVE_LIST = []
DEFAULT_AIMS_GZIP_LIST = ["std.out"]
DEFAULT_AIMS_REMOVE_LIST = []


class AimsIOError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


def job_complete(jobdir=None):
    """Return True if aims job at path 'jobdir' is complete"""
    if jobdir is None:
        jobdir = os.getcwd()
    runfile = os.path.join(jobdir, "std.out")
    if not os.path.isfile(runfile) and not os.path.isfile(runfile + ".gz"):
        return False
    if Parser(runfile).complete:
        return True
    return False


def ionic_steps(jobdir=None):
    """Find the number of ionic steps completed in 'jobdir'"""
    try:
        steps = Steps(os.path.join(jobdir, "std.out"))
        return len(steps.E)
    except IOError:
        raise AimsIOError("Could not read number of ionic steps from " + os.path.join(jobdir, "std.out"))


def write_aims_input(dirpath, controlfile, prim_posfile, super_posfile, basisfile, strict_kpoints):
    """ Update FHI-aims input files in directory 'dirpath' """
    print("Setting up FHI-aims input files:", dirpath)

    print("  Reading basis:", basisfile)
    basis_set = basis_settings(basisfile)

    print("  Reading supercell POS:", super_posfile)
    super_cell = Geometry(super_posfile)
    prim = Geometry(prim_posfile)

    print("  Reading control.in:", controlfile)
    super_ctrl = Control(controlfile)

    print("  Parsing k_grid:", controlfile)
    prim_kpoints = Kpoints(controlfile)

    print("  Generating supercell KPOINTS")
    if strict_kpoints:
        super_kpts = prim_kpoints
    else:
        super_kpts = prim_kpoints.super_kpoints(prim, super_cell)

    # write main input files
    print("  Writing supercell geometry.in:", os.path.join(dirpath, 'geometry.in'))
    super_cell.write(os.path.join(dirpath, 'geometry.in'), basis_set)

    print("  Writing control.in:", os.path.join(dirpath, 'control.in'))
    super_ctrl.write(os.path.join(dirpath, 'control.in'))

    print("  Parsing supercell k_grid and k_offset into control.in:", os.path.join(dirpath, 'control.in'))
    super_kpts.parse(os.path.join(dirpath, 'control.in'))

    print("  Adding basis set data to control.in:", os.path.join(dirpath, 'control.in'))
    write_basis(os.path.join(dirpath, 'control.in'), basis_set)

    print("  DONE\n")
    sys.stdout.flush()


def write_abort_scf(mode='e', jobdir=None):
    """ Write abort_scf file with two modes:
        mode = 'e' for stop at the next electronic step
        mode = 'i' for stop at the next ionic step
    """
    if jobdir is None:
        jobdir = os.getcwd()
    if mode.lower()[0] == 'e':
        filename = os.path.join(jobdir, 'abort_scf')
        stop_string = " "
    elif mode.lower()[0] == 'i':
        filename = os.path.join(jobdir, 'abort_opt')
        stop_string = " "
    else:
        raise AimsIOError("Invalid abort mode specified: " + str(mode))

    try:
        stop_write = open(filename, 'w')
    except IOError as e:
        raise e

    stop_write.write(stop_string)
    stop_write.close()
