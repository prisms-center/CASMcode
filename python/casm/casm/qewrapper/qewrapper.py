from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, shutil, six, re, subprocess, json
import warnings
import casm.quantumespresso.qeio as qeio
from casm import quantumespresso

class QEWrapperError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


def read_settings(filename):
    """Returns a JSON object reading JSON files containing settings for Quantum Espresso PBS jobs.

   Returns:
        settings = a JSON object containing the settings file contents
                     This can be accessed like a dict: settings["account"], etc.
                     ** All values are expected to be 'str' type. **

    The required keys are:
        "queue": queue to submit job in
        "ppn": processors (cores) per node to request
        "atom_per_proc": max number of atoms per processor (core)   <--- is this a vasp parsing field or pbs?
        "walltime": walltime to request (ex. "48:00:00")

    The optional keys are:
        "account": account to submit job under (default None)
        "pmem": string for requested memory (default None)
        "priority": requested job priority (default "0")
        "message": when to send messages about jobs (ex. "abe", default "a")
        "email": where to send messages (ex. "me@fake.com", default None)
        "qos": quality of service, 'qos' option (ex. "fluxoe")
        "quantumespresso_cmd": quantum espresso execution command (default is "pw.x < infilename > outfilename" (ncpus=1) or "mpirun -np {NCPUS} pw.x < infilename > outfilename" (ncpus!=1))
        "ncpus": number of cpus (cores) to run on (default $PBS_NP)
        "run_limit": number of quantum espresso runs until "not_converging" (default 10)
        "nrg_convergence": converged if last two runs complete and differ in energy by less than this amount (default None)
        "move": files to move at the end of a run (ex. [ ".wfc"], default [])
        "copy": files to copy from run to run ( default [infilename])
        "remove": files to remove at the end of a run ( default [".wfc",".igk",".save"]
        "compress": files to compress at the end of a run (ex. [outfilename], default [])
        "backup": files to compress to backups at the end of a run, used in conjunction with move (ex. [".wfc"])

        "encut": [START, STOP, STEP] values for converging ecutwfc to within nrg_convergence (ex. ["450", "Auto", "10"],
                 default ["Auto", "Auto", "10"] where "Auto" is either the largest ENMAX in all UPFS called in SPECIES for START,
                 or 2.0 * largest ENMAX for STOP)
        "kpoints": [start, stop, step] values for converging KPOINTS to within nrg_convergence (ex. ["5", "50", "1"],
                 default ["5", "Auto", "1"] <---- Needs to be adjusted for grid convergence
        "extra_input_files": extra input files to be copied from the settings directory, e.g., OCCUPATIONS file.
        "initial" : location of infile with tags for the initial run, if desired 
        "final" : location of infile with tags for the final run, if desired 
        "err_types" : list of errors to check for. Allowed entries are "IbzkptError" and "SubSpaceMatrixError". Default: ["SubSpaceMatrixError"] <---- STILL NEED TO IMPLEMENT
    """
    try:
        with open(filename, 'rb') as file:
            settings = json.loads(file.read().decode('utf-8'))
    except (IOError, ValueError) as e:
        print("Error reading settings file:", filename)
        raise e

    required = ["queue", "ppn", "atom_per_proc", "walltime"]


    optional = ["account","pmem","priority","message","email","qos","npar","ncore", "kpar", "ncpus","vasp_cmd","qe_cmd","run_limit","nrg_convergence", \
                "encut", "kpoints","extra_input_files", "move", "copy", "remove", "compress", "backup", "initial", "final", "strict_kpoints", "err_types", \
                "infilename","outfilename","software"]

    for key in required:
        if not key in settings:
            raise QEWrapperError( key + "' missing from: '" + filename + "'")

    for key in optional:
        if not key in settings:
            if key.lower() in ["extra_input_files", "remove", "compress", "backup","move","copy"]:
                settings[key] = []
            else:
                settings[key] = None

    if type(settings["remove"]) == list:
        if 'default' in settings["remove"]:
            settings["remove"] += qeio.DEFAULT_QE_REMOVE_LIST
    elif type(settings["remove"]) == str:
        if settings["remove"].lower() == 'default':
            settings["remove"] = qeio.DEFAULT_QE_REMOVE_LIST
        else:
            settings["remove"] = [settings["remove"]]
    if settings["priority"] == None:
        settings["priority"] = 0
    if settings["extra_input_files"] == None:
        settings["extra_input_files"] = []
    if settings["strict_kpoints"] == None:
        settings["strict_kpoints"] = False
    for k in settings.keys():
        if k not in required:
            if k not in optional:
                raise QEWrapperError("unknown key '" + k + "' found in: '" + filename + "'")

    return settings


def write_settings(settings, filename):
    """ Write 'settings' as json file, 'filename' """
    with open(filename, 'wb') as file:
        file.write(six.u(json.dumps( settings, file, indent=4)).encode('utf-8'))

      
def qe_input_file_names(dirstruc, configname, clex, infilename):
    # Find required input files in CASM project directory tree

    myinfile = dirstruc.settings_path_crawl(infilename,configname, clex)
    super_poscarfile = dirstruc.POS(configname)
    speciesfile = dirstruc.settings_path_crawl("SPECIES", configname, clex)

    # Verify that required input files exist
    if myinfile is None:
        raise quantumespresso.QuantumEspressoError("qe_input_file_names failed. No file found in CASM project matching: " + infilename )
    if super_poscarfile is None:
        raise quantumespresso.QuantumEspressoError("qe_input_file_names failed. No POS file found for this configuration.")
    if speciesfile is None:
        raise quantumespresso.QuantumEspressoError("qe_input_file_names failed. No SPECIES file found in CASM project.")

    return (myinfile, super_poscarfile, speciesfile)

