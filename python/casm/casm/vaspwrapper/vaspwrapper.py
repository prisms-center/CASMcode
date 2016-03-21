import os, shutil, re, subprocess, json
import vasp.io

class VaspWrapperError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


def read_settings(filename):
    """Returns a JSON object reading JSON files containing settings for VASP PBS jobs.

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
        "npar": vasp incar setting (default None)
        "ncore": vasp incar setting (default None)
        "kpar": vasp incar setting (default None)
        "vasp_cmd": vasp execution command (default is "vasp" (ncpus=1) or "mpirun -np {NCPUS} vasp" (ncpus!=1))
        "ncpus": number of cpus (cores) to run on (default $PBS_NP)
        "run_limit": number of vasp runs until "not_converging" (default 10)
        "nrg_convergence": converged if last two runs complete and differ in energy by less than this amount (default None)
        "move": files to move at the end of a run (ex. ["POTCAR", "WAVECAR"], default ["POTCAR"])
        "copy": files to copy from run to run (ex. ["INCAR", "KPOINTS"], default ["INCAR, KPOINTS"])
        "remove": files to remove at the end of a run (ex. ["IBZKPT", "CHGCAR"], default ["IBKZPT", "CHG", "CHGCAR", "WAVECAR", "TMPCAR", "EIGENVAL", "DOSCAR", "PROCAR", "PCDAT", "XDATCAR", "LOCPOT", "ELFCAR", "PROOUT"]
        "compress": files to compress at the end of a run (ex. ["OUTCAR", "vasprun.xml"], default [])
        "backup": files to compress to backups at the end of a run, used in conjunction with move (ex. ["WAVECAR"])

        "encut": [START, STOP, STEP] values for converging ENCUT to within nrg_convergence (ex. ["450", "Auto", "10"],
                 default ["Auto", "Auto", "10"] where "Auto" is either the largest ENMAX in all POTCARS called in SPECIES for START,
                 or 2.0 * largest ENMAX for STOP)
        "kpoints": [start, stop, step] values for converging KPOINTS to within nrg_convergence (ex. ["5", "50", "1"],
                 default ["5", "Auto", "1"] where "Auto" can only be used for STOP and means to keep increasing the KPOINTS length
                 by STEP until either nrg_convergence or walltime is reached). For meaning of the KPOINTS length parameter see
                 the VASP documentation at http://cms.mpi.univie.ac.at/vasp/vasp/Automatic_k_mesh_generation.html
        "extra_input_files": extra input files to be copied from the settings directory, e.g., a vdW kernel file.
        "initial" : location of INCAR with tags for the initial run, if desired (e.g. to generate a PBE WAVECAR for use with M06-L)
        "final" : location of INCAR with tags for the final run, if desired (e.g. "ISMEAR = -5", etc). Otherwise, the settings enforced are ("ISMEAR = -5", "NSW = 0", "IBRION = -1", "ISIF = 2")
        "err_types" : list of errors to check for. Allowed entries are "IbzkptError" and "SubSpaceMatrixError". Default: ["SubSpaceMatrixError"]
    """
    try:
        file = open(filename)
        settings = json.load(file)
        file.close()
    except (IOError, ValueError) as e:
        print "Error reading settings file:", filename
        raise e

    required = ["queue", "ppn", "atom_per_proc", "walltime"]

    optional = ["account","pmem","priority","message","email","qos","npar","ncore", "kpar", "ncpus","vasp_cmd","run_limit","nrg_convergence", \
                "encut", "kpoints","extra_input_files", "move", "copy", "remove", "compress", "backup", "initial", "final","err_types"]
    for key in required:
        if not key in settings:
            raise VaspWrapperError( key + "' missing from: '" + filename + "'")

    for key in optional:
        if not key in settings:
            if key.lower() in ["extra_input_files", "remove", "compress", "backup"]:
                settings[key] = []
            elif key.lower() in ["move"]:
                settings[key] = vasp.io.DEFAULT_VASP_MOVE_LIST
            elif key.lower() in ["copy"]:
                settings[key] = vasp.io.DEFAULT_VASP_COPY_LIST
            # elif key.lower() in ["remove"]:
            #     settings[key] = vasp.io.DEFAULT_VASP_REMOVE_LIST
            else:
                settings[key] = None

    if type(settings["remove"]) == list:
        if 'default' in settings["remove"]:
            settings["remove"] += vasp.io.DEFAULT_VASP_REMOVE_LIST
    elif type(settings["remove"]) == str:
        if settings["remove"].lower() == 'default':
            settings["remove"] = vasp.io.DEFAULT_VASP_REMOVE_LIST
        else:
            settings["remove"] = [settings["remove"]]
    if settings["priority"] == None:
        settings["priority"] = 0
    if settings["extra_input_files"] == None:
        settings["extra_input_files"] = []

    for k in settings.keys():
        if k not in required:
            if k not in optional:
                raise VaspWrapperError("unknown key '" + k + "' found in: '" + filename + "'")

    return settings


def write_settings(settings, filename):
    """ Write 'settings' as json file, 'filename' """
    file = open(filename,'w')
    json.dump( settings, file, indent=4)
    file.close()


