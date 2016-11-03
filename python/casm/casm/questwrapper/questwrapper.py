""" FIXME """

import json
import seqquest.seqquest_io

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
        "run_cmd": vasp execution command (default is "vasp" (ncpus=1) or "mpirun -np {NCPUS} vasp" (ncpus!=1))
        "ncpus": number of cpus (cores) to run on (default $PBS_NP)
        "run_limit": number of vasp runs until "not_converging" (default 10)
        "nrg_convergence": converged if last two runs complete and differ in energy by less than this amount (default None)
        "move": files to move at the end of a run (ex. ["POTCAR", "WAVECAR"], default ["POTCAR"])
        "copy": files to copy from run to run (ex. ["INCAR", "KPOINTS"], default ["INCAR, KPOINTS"])
        "remove": files to remove at the end of a run (ex. ["IBZKPT", "CHGCAR"], default ["IBKZPT", "CHG", "CHGCAR", "WAVECAR", "TMPCAR", "EIGENVAL", "DOSCAR", "PROCAR", "PCDAT", "XDATCAR", "LOCPOT", "ELFCAR", "PROOUT"]
        "compress": files to compress at the end of a run (ex. ["OUTCAR", "vasprun.xml"], default [])
        "backup": files to compress to backups at the end of a run, used in conjunction with move (ex. ["WAVECAR"])
        "extra_input_files": extra input files to be copied from the settings directory, e.g., a vdW kernel file.
        "initial" : location of INCAR with tags for the initial run, if desired (e.g. to generate a PBE WAVECAR for use with M06-L)
        "final" : location of INCAR with tags for the final run, if desired (e.g. "ISMEAR = -5", etc). Otherwise, the settings enforced are ("ISMEAR = -5", "NSW = 0", "IBRION = -1", "ISIF = 2")
        "err_types" : list of errors to check for. Allowed entries are "IbzkptError" and "SubSpaceMatrixError". Default: ["SubSpaceMatrixError"]
        "preamble" : a text file containing anything that MUST be run before python is invoked (e.g. module.txt which contains "module load python", or "source foo")
        "prop": USED IN vasp.converge ONLY. Property to converge with respect to (current options are "KPOINT" and "ENCUT")
        "prop_start": USED IN vasp.converge ONLY. Starting value of "prop", e.g. 450 (for ENCUT) or 5 (for KPOINTS) or [4 4 4] (for KPOINTS)
        "prop_stop": USED IN vasp.converge ONLY. Ending value of "prop", e.g. 550 (for ENCUT) or 20 (for KPOINTS).
        "prop_step": USED IN vasp.converge ONLY. Delta value of "prop", e.g. 10 (for ENCUT) or 2 (for KPOINTS) or [1 1 2] (for KPOINTS)
        "tol" : USED IN vasp.converge ONLY. Tolerance type for convergence, e.g. relaxed_energy. Optional
        "tol_amount" : USED IN vasp.converge ONLY. Tolerance criteria convergence, e.g. 0.001. If the abs difference between two runs in their "tol" is smaller than "tol_amount", the "prop" is considered converged.
        "name" : USED IN vasp.converge ONLY. Name used in the .../config/calctype.calc/NAME/property_i directory scheme, where, if not specified, "prop"_converge is used as NAME
    """
    try:
        stream = open(filename)
        settings = json.load(stream)
        stream.close()
    except (IOError, ValueError) as e:
        print "Error reading settings file:", filename
        raise e

    required = ["queue", "ppn", "atom_per_proc", "walltime"]

    optional = ["account","pmem","priority","message","email","qos",
                "ncpus","run_cmd","run_limit","nrg_convergence", "encut", "kpoints",
                "extra_input_files", "move", "copy", "remove", "compress", "backup", "initial",
                "final", "strict_kpoints", "err_types", "preamble", "prop", "prop_start",
                "prop_stop", "prop_step", "tol", "tol_amount", "name", "fine_ngx",
                "cont_relax"]
    for key in required:
        if not key in settings:
            raise QuestWrapperError( key + "' missing from: '" + filename + "'")

    for key in optional:
        if not key in settings:
            if key.lower() in ["extra_input_files", "remove", "compress", "backup"]:
                settings[key] = []
            elif key.lower() in ["move"]:
                settings[key] = seqquest.seqquest_io.DEFAULT_QUEST_MOVE_LIST
            elif key.lower() in ["copy"]:
                settings[key] = seqquest.seqquest_io.DEFAULT_QUEST_COPY_LIST
            # elif key.lower() in ["remove"]:
            #     settings[key] = vasp.io.DEFAULT_VASP_REMOVE_LIST
            else:
                settings[key] = None

    if isinstance(settings["remove"], list):
        if 'default' in settings["remove"]:
            settings["remove"] += seqquest.seqquest_io.DEFAULT_QUEST_REMOVE_LIST
    elif isinstance(settings["remove"], str):
        if settings["remove"].lower() == 'default':
            settings["remove"] = seqquest.seqquest_io.DEFAULT_QUEST_REMOVE_LIST
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
    stream = open(filename, 'w')
    json.dump(settings, stream, indent=4)
    stream.close()

def read_properties(filename):
    """ Read a properties.calc.json"""
    required = ["atom_type", "atoms_per_type", "coord_mode", "relaxed_basis", "relaxed_energy", "relaxed_forces", "relaxed_lattice"]
    optional = ["relaxed_magmom", "relaxed_mag_basis"]

    with open(filename, 'r') as myfile:
        properties = json.load(myfile)

    for key in required:
        if not key in properties:
            raise QuestWrapperError(key + "' missing from: '" + filename + "'")

    for key in optional:
        if not key in properties:
            properties[key] = None

    return properties
