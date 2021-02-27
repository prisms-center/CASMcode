from __future__ import absolute_import, division, print_function, unicode_literals

import json
import six
from casm.misc import noindent


class ProjectIOError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def write_eci(proj, eci, fit_details=None, clex=None, verbose=False):
    """
    Write eci.json
    
    Arguments
    ---------
    
      proj: casm.project.Project instance
        The CASM project
      
      eci: List[(index, value)]
        index (int): linear index of basis function
        value (float): ECI value

      fit_details: Dict
        Description of the fitting method used to generate the ECI,
        usually as output by casm.learn.to_json
      
      clex: ClexDescription instance, optional, default=proj.settings.default_clex
        Specifies where to write the ECI

      verbose: be verbose?
    
    """
    wdir = proj.dir
    if clex is None:
      clex = proj.settings.default_clex
    
    # read basis.json
    filename = wdir.basis(clex)
    with open(filename, 'rb') as f:
        j = json.loads(f.read().decode('utf-8'))
    #print(json.dumps(j, indent=2))
    
    # edit to add fitting settings
    j["fit"] = fit_details
    
    # edit to add eci
    for index, value in eci:
      j["cluster_functions"][index]["eci"] = value
    
    # pretty printing
    for entry in j["site_functions"]:
      if entry["basis"] is not None:
        basis = entry["basis"]
        for key, val in basis.items():
          basis[key] = noindent.NoIndent(val)
    for entry in j["cluster_functions"]:
      entry["orbit"] = noindent.NoIndent(entry["orbit"])
      sites = entry["prototype"]["sites"]
      for i in range(len(sites)):
        sites[i] = noindent.NoIndent(sites[i])
    
    # write eci.json
    filename = wdir.eci(clex)
    
    if verbose:
      print("Writing:", filename, "\n")
    with open(filename, 'wb') as f:
      f.write(six.u(json.dumps(j, indent=2, cls=noindent.NoIndentEncoder)).encode('utf-8'))
    
    # refresh proj to reflect new eci
    proj.refresh(clear_clex=True)


def read_project_settings(filename):
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
        "dft_cmd": DFT execution command (NO default)
        "ncpus": number of cpus (cores) to run on (default $PBS_NP)
        "run_limit": number of quantum espresso runs until "not_converging" (default 10)
        "nrg_convergence": converged if last two runs complete and differ in energy by less
        "move": files to move at the end of a run (ex. [ ".wfc"], default [])
        "copy": files to copy from run to run ( default [infilename])
        "remove": files to remove at the end of a run ( default [".wfc",".igk",".save"]
        "compress": files to compress at the end of a run (ex. [outfilename], default [])
        "backup": files to compress to backups at the end of a run, used in conjunction with move (ex. [".wfc"])
        "encut": [START, STOP, STEP] values for converging ecutwfc to within nrg_convergence (ex. ["450", "Auto", "10"],
                 default ["Auto", "Auto", "10"] where "Auto" is either the largest ENMAX in all
                 UPFS called in SPECIES for START, or 2.0 * largest ENMAX for STOP)
        "kpoints": [start, stop, step] values for converging KPOINTS to within nrg_convergence (ex. ["5", "50", "1"],
                 default ["5", "Auto", "1"] <---- Needs to be adjusted for grid convergence
        "extra_input_files": extra input files to be copied from the settings directory, e.g., OCCUPATIONS file.
        "initial" : location of infile with tags for the initial run, if desired
        "final" : location of infile with tags for the final run, if desired
        "err_types" : list of errors to check for. Allowed entries are "IbzkptError"
                      and "SubSpaceMatrixError". Default: ["SubSpaceMatrixError"] <---- STILL NEED TO IMPLEMENT
    """
    try:
        with open(filename, 'rb') as file:
            settings = json.loads(file.read().decode('utf-8'))
    except IOError as e:
        print("Error reading settings file:", filename)
        raise e

#    if settings['software'] != 'vasp':
#        raise IOError('ONLY VASP in experimental status...')

    required = ["queue", "ppn", "walltime", "software", "run_cmd"]

    optional = ["account", "pmem", "priority", "message", "email", "qos", "npar", "ncore", "kpar",
                "ncpus", "run_limit", "nrg_convergence", "prerun", "postrun",
                "encut", "kpoints", "extra_input_files", "move", "copy", "remove", "compress",
                "backup", "initial", "final", "strict_kpoints", "err_types", "preamble",
                "infilename", "outfilename", "atom_per_proc", "nodes", "is_slab", "fix_pos", "basis"]

    for key in required:
        if key not in settings:
            raise ProjectIOError(key + "' missing from: '" + filename + "'")

    for key in optional:
        if key not in settings:
            if key.lower() in ["extra_input_files", "remove", "compress", "backup", "move", "copy"]:
                settings[key] = []
            else:
                settings[key] = None

    # TODO: Define a COMMON DEFAULT_XXXX_LIST for all codes and append / handle these, don't use separates for codes

    if settings['software'] == 'vasp':
        from casm.vasp.io.io import DEFAULT_VASP_COPY_LIST as DEFAULT_COPY_LIST
        from casm.vasp.io.io import DEFAULT_VASP_MOVE_LIST as DEFAULT_MOVE_LIST
        from casm.vasp.io.io import DEFAULT_VASP_REMOVE_LIST as DEFAULT_REMOVE_LIST
        DEFAULT_GZIP_LIST = []
    elif settings['software'] == 'qe':
        from casm.quantumespresso.qeio import DEFAULT_QE_COPY_LIST as DEFAULT_COPY_LIST
        from casm.quantumespresso.qeio import DEFAULT_QE_MOVE_LIST as DEFAULT_MOVE_LIST
        from casm.quantumespresso.qeio import DEFAULT_QE_REMOVE_LIST as DEFAULT_REMOVE_LIST
        DEFAULT_GZIP_LIST = []
    elif settings['software'] == 'aims':
        from casm.aims.io.io import DEFAULT_AIMS_COPY_LIST as DEFAULT_COPY_LIST
        from casm.aims.io.io import DEFAULT_AIMS_MOVE_LIST as DEFAULT_MOVE_LIST
        from casm.aims.io.io import DEFAULT_AIMS_REMOVE_LIST as DEFAULT_REMOVE_LIST
        from casm.aims.io.io import DEFAULT_AIMS_GZIP_LIST as DEFAULT_GZIP_LIST
    else:
        raise ProjectIOError('Software needs to be defined in relax.json.')

    if type(settings["remove"]) == list:
        if 'default' in settings["remove"]:
            settings["remove"] += DEFAULT_REMOVE_LIST
    elif type(settings["remove"]) == str:
        if settings["remove"].lower() == 'default':
            settings["remove"] = DEFAULT_REMOVE_LIST
        else:
            settings["remove"] = [settings["remove"]]

    if "priority" not in settings:
        settings["priority"] = 0
    if "extra_input_files" not in settings:
        settings["extra_input_files"] = []
    if "strict_kpoints" not in settings:
        settings["strict_kpoints"] = False

    for k in settings.keys():
        if k not in required:
            if k not in optional:
                print("unknown key '" + k + "' found in: '" + filename + "'")
                raise ProjectIOError("unknown key '" + k + "' found in: '" + filename + "'")

    settings['DEF_CP'] = DEFAULT_COPY_LIST
    settings['DEF_MV'] = DEFAULT_MOVE_LIST
    settings['DEF_RM'] = DEFAULT_REMOVE_LIST
    settings['DEF_GZ'] = DEFAULT_GZIP_LIST

    return settings


def read_band_settings(filename):
    """Returns a JSON object reading JSON files containing settings for band structure plots.

    Returns:
         settings = a JSON object containing the settings file contents
                      This can be accessed like a dict: settings["account"], etc.
                      ** All values are expected to be 'str' type. **
    """
    try:
        with open(filename, 'rb') as file:
            settings = json.loads(file.read().decode('utf-8'))
    except IOError as e:
        print("Error reading settings file:", filename)
        raise e

    band_args = ["band_subdiv", "band_bs_projection", "band_dos_projection", "band_vb_energy_range",
                 "band_cb_energy_range", "band_fixed_cb_energy", "band_egrid_interval", "band_font",
                 "band_axis_fontsize", "band_tick_fontsize", "band_legend_fontsize", "band_bs_legend",
                 "band_dos_legend", "band_rgb_legend", "band_fig_size", "band_plot_name", "band_plot_dpi",
                 "band_kpt_dens"]

    for key in band_args:
        if key not in settings:
            if key == 'band_kpt_dens':
                settings['band_kpt_dens'] = 1000
            if key == 'band_subdiv':
                settings['band_subdiv'] = 100
            if key == 'band_bs_projection':
                settings['band_bs_projection'] = 'elements'
            if key == 'band_dos_projection':
                settings['band_dos_projection'] = 'elements'
            if key == 'band_vb_energy_range':
                settings['band_vb_energy_range'] = 4
            if key == 'band_cb_energy_range':
                settings['band_cb_energy_range'] = 4
            if key == 'band_fixed_cb_energy':
                settings['band_fixed_cb_energy'] = False
            if key == 'band_egrid_interval':
                settings['band_egrid_interval'] = 1
            if key == 'band_font':
                settings['band_font'] = 'Times New Roman'
            if key == 'band_axis_fontsize':
                settings['band_axis_fontsize'] = 20
            if key == 'band_tick_fontsize':
                settings['band_tick_fontsize'] = 15
            if key == 'band_legend_fontsize':
                settings['band_legend_fontsize'] = 14
            if key == 'band_bs_legend':
                settings['band_bs_legend'] = 'best'
            if key == 'band_dos_legend':
                settings['band_dos_legend'] = 'best'
            if key == 'band_rgb_legend':
                settings['band_rgb_legend'] = True
            if key == 'band_fig_size':
                settings['band_fig_size'] = (11, 8.5)
            if key == 'band_plot_name':
                settings['band_plot_name'] = 'BSDOS.png'
            if key == 'band_plot_dpi':
                settings['band_plot_dpi'] = 170

    return settings


def read_feff_settings(filename):
    """Returns a JSON object reading JSON files containing settings for XAS computations.

    Returns:
         settings = a JSON object containing the settings file contents
                      This can be accessed like a dict: settings["account"], etc.
                      ** All values are expected to be 'str' type. **
    """
    try:
        with open(filename, 'rb') as file:
            settings = json.loads(file.read().decode('utf-8'))
    except IOError as e:
        print("Error reading settings file:", filename)
        raise e

    feff_args = ["feff_cmd", "feff_nkpts", "feff_radius", "feff_user_tags",
                 "feff_plot_sigma", "feff_plot_use_omega", "feff_base_dir",
                 "feff_base_mpi"]

    for key in feff_args:
        if key not in settings:
            if key == 'feff_cmd':
                settings['feff_cmd'] = 'feff'
            if key == 'feff_nkpts':
                settings['feff_nkpts'] = 1000
            if key == 'feff_radius':
                settings['feff_radius'] = 12
            if key == 'feff_user_tags':
                settings['feff_user_tags'] = {}
            if key == 'feff_plot_sigma':
                settings['feff_plot_sigma'] = 1.0
            if key == 'feff_plot_use_omega':
                settings['feff_plot_use_omega'] = False

    return settings
