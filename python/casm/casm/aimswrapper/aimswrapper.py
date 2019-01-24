import json

from casm.project.io import read_project_settings

from casm.aims.aims import AimsError
from casm.aims.io.io import DEFAULT_AIMS_COPY_LIST, DEFAULT_AIMS_MOVE_LIST


class AimsWrapperError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


def read_settings(filename):
    """Returns a JSON object reading JSON files containing settings for FHI-aims PBS jobs.

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
        "aims_cmd": FHI-aims execution command (default is "aims" (ncpus=1) or "mpirun -np {NCPUS} aims" (ncpus>1))
        "ncpus": number of cpus (cores) to run on (default $PBS_NP)
        "run_limit": number of vasp runs until "not_converging" (default 10)
        "err_types" : list of errors to check for, currently ony FrozenError
        "prerun" : bash commands to run before aims.Relax.run (default None)
        "postrun" : bash commands to run after aims.Relax.run completes (default None)
    """
    try:
        file = open(filename)
        settings = read_project_settings(filename)
        file.close()
    except IOError as e:
        print("Error reading settings file:", filename)
        raise e

    required = ["is_slab", "basis"]
    optional = ["fix_pos"]

    for key in required:
        if key not in settings:
            raise AimsWrapperError(key + "' missing from: '" + filename + "'")

    for key in optional:
        if key not in settings:
            if key.lower() in ["extra_input_files", "remove", "compress", "backup"]:
                settings[key] = []
            elif key.lower() in ["move"]:
                settings[key] = DEFAULT_AIMS_MOVE_LIST
            elif key.lower() in ["copy"]:
                settings[key] = DEFAULT_AIMS_COPY_LIST
            else:
                settings[key] = None

    if settings["basis"] != "light" and settings["basis"] != "tight":
        raise AimsWrapperError("Basis setting for FHI-aims must be >light< or >tight<")

    if settings["is_slab"] == "True":
        if 'fix_pos' not in settings:
            err_str = "Error: When running a slab you need to define a 'fix_pos' key"
            err_str += "which defines the position below which atoms in the slab will be fixed (all coordinates)."
            err_str += "fix_pos: Units in absolute Angstrom coordinates (NOT fractional!)"
            raise AimsWrapperError(err_str)
        settings["is_slab"] = True
    else:
        settings["is_slab"] = False

    if 'strict_kpoints' not in settings:
        settings['strict_kpoints'] = False
    else:
        settings['strict_kpoints'] = True

    return settings


def write_settings(settings, filename):
    """ Write 'settings' as json file, 'filename' """
    with open(filename, 'w') as f:
        json.dump(settings, f, indent=4)


def aims_input_file_names(workdir, configname, clex):
    """
    Collect casm.aimswrapper input files from the CASM project hierarchy

    Looks for:

      control.in:   
        The base input file used for calculations. Found via:
          DirectoryStructure.settings_path_crawl

      POS: 
        The CASM-generated POS file giving the initial structure to be calculated.

      SPECIES:
        The SPECIES file specifying FHI-aims basis settings for each species in the structure.

    Arguments
    ---------
      
      workdir: casm.project.DirectoryStructure instance
        CASM project directory hierarchy
      
      configname: str
        The name of the configuration to be calculated
      
      clex: casm.project.ClexDescription instance
        The cluster expansion being worked on. Used for the 'calctype' settings.
    
    
    Returns
    -------
      
      filepaths: tuple(control.in, POS, SPECIES)
        A tuple containing the paths to the aimswrapper input files
    
    
    Raises
    ------
      If any required file is not found.
    
    """
    # Find required input files in CASM project directory tree
    controlfile = workdir.settings_path_crawl("control.skel", configname, clex)
    prim_geometryfile = workdir.settings_path_crawl("geometry.skel", configname, clex)
    super_poscarfile = workdir.POS(configname)
    basisfile = workdir.settings_path_crawl("basis", configname, clex)

    # Verify that required input files exist
    if controlfile is None:
        raise AimsError("aims_input_file_names failed. No control.skel file found in CASM project.")
    if prim_geometryfile is None:
        raise AimsError("No reference geometry.skel found in CASM project.")
    if super_poscarfile is None:
        raise AimsError("aims_input_file_names failed. No POS file found for this configuration.")
    if basisfile is None:
        raise AimsError("aims_input_file_names failed. No SPECIES file found in CASM project.")

    return controlfile, prim_geometryfile, super_poscarfile, basisfile
