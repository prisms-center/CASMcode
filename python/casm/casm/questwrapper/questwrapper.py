""" Wrapper for handling seqquest/casm integration """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import six

from casm.seqquest import seqquest_io

class QuestWrapperError(Exception):
    """ Errors related to QuestWrapper """
    pass


def read_settings(filename):
    """Returns a JSON object reading JSON files containing settings for seqquest cluster jobs.

   Returns:
        settings = a JSON object containing the settings file contents
                     This can be accessed like a dict: settings["account"], etc.
                     ** All values are expected to be 'str' type. **

    The required keys are:
        "queue": queue to submit job in
        "ppn": processors (cores) per node to request
        "atom_per_proc": max number of atoms per processor (core)
        "walltime": walltime to request (ex. "48:00:00")

    The optional keys are:
        "account": account to submit job under (default None)
        "pmem": string for requested memory (default None)
        "priority": requested job priority (default "0")
        "message": when to send messages about jobs (ex. "abe", default "a")
        "email": where to send messages (ex. "me@fake.com", default None)
        "qos": quality of service, 'qos' option (ex. "fluxoe")
        "run_cmd": quest execution command (default is "vasp" (ncpus=1) or "mpirun -np {NCPUS} vasp" (ncpus!=1))
        "ncpus": number of cpus (cores) to run on (default $PBS_NP)
        "run_limit": number of vasp runs until "not_converging" (default 10)
        "nrg_convergence": converged if last two runs complete and differ in energy by less than this amount (default None)
        "move": files to move at the end of a run (ex. ["POTCAR", "WAVECAR"], default ["POTCAR"])
        "copy": files to copy from run to run (ex. ["INCAR", "KPOINTS"], default ["INCAR, KPOINTS"])
        "remove": files to remove at the end of a run (ex. ["IBZKPT", "CHGCAR"], default ["IBKZPT", "CHG", "CHGCAR", "WAVECAR", "TMPCAR", "EIGENVAL", "DOSCAR", "PROCAR", "PCDAT", "XDATCAR", "LOCPOT", "ELFCAR", "PROOUT"]
        "compress": files to compress at the end of a run (ex. ["OUTCAR", "vasprun.xml"], default [])
        "backup": files to compress to backups at the end of a run, used in conjunction with move (ex. ["WAVECAR"])
        "cont_relax": name of a calctype to initialize a calculation from, using the lcao.geom from ../calctype.cont_relax/run.final (ex. "default")
        "extra_input_files": extra input files to be copied from the settings directory, e.g., a vdW kernel file.
        "initial" : location of INCAR with tags for the initial run, if desired (e.g. to generate a PBE WAVECAR for use with M06-L)
        "final" : location of INCAR with tags for the final run, if desired (e.g. "ISMEAR = -5", etc). Otherwise, the settings enforced are ("ISMEAR = -5", "NSW = 0", "IBRION = -1", "ISIF = 2")
        "err_types" : list of errors to check for. Allowed entries are "IbzkptError" and "SubSpaceMatrixError". Default: ["SubSpaceMatrixError"]
        "preamble" : a text file containing anything that MUST be run before python is invoked (e.g. module.txt which contains "module load python", or "source foo")
        "prop": USED IN vasp.converge ONLY. Property to converge with respect to (current options are "KPOINT" and "ENCUT")
        "prerun" : bash commands to run before vasp.Relax.run (default None)
        "postrun" : bash commands to run after vasp.Relax.run completes (default None)
    """
    try:
        with open(filename, 'rb') as stream:
            settings = json.loads(stream.read().decode('utf-8'))
    except (IOError, ValueError) as e:
        print("Error reading settings file:", filename)
        raise e

    required = ["queue", "ppn", "atom_per_proc", "walltime"]

    optional = ["account","pmem","priority","message","email","qos",
                "ncpus","run_cmd","run_limit","nrg_convergence", "encut", "kpoints",
                "extra_input_files", "move", "copy", "remove", "compress", "backup", "initial",
                "final", "strict_kpoints", "err_types", "preamble", "prop", "prop_start",
                "prop_stop", "prop_step", "tol", "tol_amount", "name", "fine_ngx",
                "prerun", "postrun", "cont_relax"]
    for key in required:
        if not key in settings:
            raise QuestWrapperError( key + "' missing from: '" + filename + "'")

    for key in optional:
        if not key in settings:
            if key.lower() in ["extra_input_files", "remove", "compress", "backup"]:
                settings[key] = []
            elif key.lower() in ["move"]:
                settings[key] = seqquest_io.DEFAULT_QUEST_MOVE_LIST
            elif key.lower() in ["copy"]:
                settings[key] = seqquest_io.DEFAULT_QUEST_COPY_LIST
            # elif key.lower() in ["remove"]:
            #     settings[key] = vasp.io.DEFAULT_VASP_REMOVE_LIST
            else:
                settings[key] = None

    if isinstance(settings["remove"], list):
        if 'default' in settings["remove"]:
            settings["remove"] += seqquest_io.DEFAULT_QUEST_REMOVE_LIST
    elif isinstance(settings["remove"], str):
        if settings["remove"].lower() == 'default':
            settings["remove"] = seqquest_io.DEFAULT_QUEST_REMOVE_LIST
        else:
            settings["remove"] = [settings["remove"]]
    if settings["priority"] is None:
        settings["priority"] = 0
    if settings["extra_input_files"] is None:
        settings["extra_input_files"] = []
    if settings["strict_kpoints"] is None:
        settings["strict_kpoints"] = False
    if settings["fine_ngx"] is None:
        settings["fine_ngx"] = False
    for k in settings.keys():
        if k not in required:
            if k not in optional:
                raise QuestWrapperError("unknown key '" + k + "' found in: '" + filename + "'")

    return settings


def write_settings(settings, filename):
    """ Write 'settings' as json file, 'filename' """
    with open(filename, 'wb') as stream:
        stream.write(six.u(json.dumps(settings, stream, indent=4)).encode('utf-8'))

def read_properties(filename):
    """ Read a properties.calc.json"""
    required = ["atom_type", "atoms_per_type", "coord_mode", "relaxed_basis", "relaxed_energy", "relaxed_forces", "relaxed_lattice"]
    optional = ["relaxed_magmom", "relaxed_mag_basis"]

    with open(filename, 'rb') as myfile:
        properties = json.load(myfile.read().decode('utf-8'))

    for key in required:
        if not key in properties:
            raise QuestWrapperError(key + "' missing from: '" + filename + "'")

    for key in optional:
        if not key in properties:
            properties[key] = None

    return properties

def quest_input_file_names(dir, configname, clex):
    """
    Collect casm.questwrapper input files from the CASM project hierarchy

    Looks for:

      lcao.in:
        The base lcao.in file used for calculations. Found via:
          DirectoryStructure.settings_path_crawl

      POS:
        The CASM-generated POS file giving the initial structure to be calculated.

      SPECIES:
        The SPECIES file specifying Vasp settings for each species in the structure.


    Arguments
    ---------

      dir: casm.project.DirectoryStructure instance
        CASM project directory hierarchy

      configname: str
        The name of the configuration to be calculated

      clex: casm.project.ClexDescription instance
        The cluster expansion being worked on. Used for the 'calctype' settings.


    Returns
    -------

      filepaths: tuple(lcao.in, POS, SPECIES)
        A tuple containing the paths to the questwrapper input files


    Raises
    ------
      If any required file is not found.

    """
    # Find required input files in CASM project directory tree
    lcao_in = dir.settings_path_crawl("lcao.in", configname, clex)
    super_poscarfile = dir.POS(configname)
    speciesfile = dir.settings_path_crawl("SPECIES", configname, clex)

    # Verify that required input files exist
    if lcao_in is None:
        raise seqquest.SeqQuestError("Relax.setup failed. No lcao.in file found in CASM\
                                    project.")
    if super_poscarfile is None:
        raise seqquest.SeqQuestError("Relax.setup failed. No POS file found for this\
                                    configuration.")
    if speciesfile is None:
        raise seqquest.SeqQuestError("Relax.setup failed. No SPECIES file found in CASM\
                                    project.")

    return (lcao_in, super_poscarfile, speciesfile)

