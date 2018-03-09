"""See class info"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

### External ###
import os
import shutil
import six
from math import ceil, sqrt
import sys
import json
# import re
# import warnings
try:
  from prisms_jobs import Job, JobDB, error_job, complete_job, JobsError, JobDBError, EligibilityError
except ImportError:
  # use of the pbs module is deprecated after CASM v0.2.1
  from pbs import Job, JobDB, error_job, complete_job, JobDBError, EligibilityError
  from pbs import PBSError as JobsError


### Local ###
from casm import vasp, wrapper
from casm.misc import noindent
from casm.project import DirectoryStructure, ProjectSettings
from casm.vaspwrapper import VaspWrapperError, read_settings, write_settings, \
  vasp_input_file_names

### Globals ###
VALID_PROP_TYPES = ["KPOINTS", "ENCUT", "NBANDS", "SIGMA"]
VALID_TOL_TYPES = ["relaxed_energy", "relaxed_volume", "relaxed_lattice"]

#pylint: disable=line-too-long, too-many-instance-attributes, too-many-arguments, too-many-branches, too-many-statements, too-many-locals, too-many-nested-blocks, too-many-lines


class ConvergeError(Exception):
    """ Exception sub-class for errors related to convergence """
    pass

class Converge(object):
    """The Converge class contains functions for setting up, executing, and parsing a VASP convergence test.

        The convergence creates the following directory structure:
        config/
            calctype.name/
                property/
                    value.0/
                    value.1/
                    ...
                    value.n/

        'value.i' directories are created for each discrete step in 'value'
          being explored for 'property', e.g., every 10 eV for ENCUT, or every
          1 for KPOINTS.
        If runs are static runs, 'value.i' contains INCAR, POSCAR, etc.
        If runs are relaxation runs, 'value.i' contains 'run.0', 'run.1', ...,
          'run.final', where 'run.final' contains a static run.

        This automatically looks for VASP settings files using:
          casm.project.DirectoryStructure.settings_path_crawl

    Attributes
    ----------

      casm_settings: casm.project.ProjectSettings instance
        CASM project settings

      casm_directories: casm.project.DirectoryStructure instance
        CASM project directory hierarchy

      settings: dict
        Settings for job submission and the relaxation, see vaspwrapper.read_settings

      configdir: str
        Directory where configuration results are stored. The result of:
          casm.project.DirectoryStructure.configuration_dir(self.configname)

      configname: str
        The name of the configuration to be calculated

      auto: boolean
        True if using prisms_jobs module's JobDB to manage jobs

      sort: boolean
        True if sorting atoms in POSCAR by type

      clex: casm.project.ClexDescription instance
        The cluster expansion being worked on. Used for the 'calctype' settings.
        Currently, fixed to self.casm_settings.default_clex.

    """
    def __init__(self, configdir=None, auto=True, sort=True, propdir=None, prop=None):
        """
        Construct a VASP convergence job object.

        Arguments
        ----------

            configdir: str, optional, default=None
              Path to configuration directory. If None, uses the current working
              directory
            auto: boolean, optional, default=True,
              Use True to use the prisms_jobs module's JobDB to manage jobs

            sort: boolean, optional, default=True,
              Use True to sort atoms in POSCAR by type

            propdir: str, optional, default=None
              Name of the directory within the configdir for converging the
                current property-of-interest

            prop: str, optional, default=None
              Name of an INCAR property that is being converged

        """
        print("Construct a casm.vaspwrapper.Converge instance:")

        if configdir is None:
            configdir = os.getcwd()
        print("  Input directory:", configdir)

        self.propdir = propdir
        self.prop = prop

        # get the configname from the configdir path
        _res = os.path.split(configdir)
        self.configname = os.path.split(_res[0])[1] + "/" + _res[1]
        print("  Configuration:", self.configname)

        print("Reading CASM settings")
        self.casm_directories = DirectoryStructure(configdir)
        self.casm_settings = ProjectSettings(configdir)
        if self.casm_settings is None:
            raise VaspWrapperError("Not in a CASM project. The file '.casm' directory was not found.")

        if os.path.abspath(configdir) != self.configdir:
            print("")
            print("input configdir:", configdir)
            print("determined configname:", self.configname)
            print("expected configdir given configname:", self.configdir)
            raise VaspWrapperError("Mismatch between configname and configdir")

        # fixed to default_clex for now
        self.clex = self.casm_settings.default_clex

        # store path to .../config/calctype.name, and create if not existing
        self.calcdir = self.casm_directories.calctype_dir(self.configname, self.clex)
        try:
            os.mkdir(self.calcdir)
        except:
            pass
        print("  Calculations directory:", self.calcdir)

        # read the settings json file
        print("  Reading converge.json settings file")
        sys.stdout.flush()
        setfile = self.casm_directories.settings_path_crawl("converge.json",
                                                            self.configname,
                                                            self.clex)

        if setfile is None:
            raise VaspWrapperError("Could not find \"converge.json\" in an appropriate \"settings\" directory")

        else:
            print("  Read settings from:", setfile)
        self.settings = read_settings(setfile)

        # add required keys to settings if not present
        if not "ncore" in self.settings:
            self.settings["ncore"] = None
        if not "npar" in self.settings:
            self.settings["npar"] = None
        if not "kpar" in self.settings:
            self.settings["kpar"] = None
        if not "vasp_cmd" in self.settings:
            self.settings["vasp_cmd"] = None
        if not "ncpus" in self.settings:
            self.settings["ncpus"] = None
        if not "run_limit" in self.settings:
            self.settings["run_limit"] = None
        if not "prerun" in self.settings:
            self.settings["prerun"] = None
        if not "postrun" in self.settings:
            self.settings["postrun"] = None
        if not "prop" in self.settings:
            self.settings["prop"] = None
        if not "prop_start" in self.settings:
            self.settings["prop_start"] = None
        if not "prop_step" in self.settings:
            self.settings["prop_step"] = None
        if not "prop_stop" in self.settings:
            self.settings["prop_stop"] = None
        if not "tol" in self.settings:
            self.settings["tol"] = "relaxed_energy"
        if not "tol_amount" in self.settings:
            self.settings["tol_amount"] = 0.001

        self.auto = auto
        self.sort = sort

        # Invoke the settings checker
        self._settings_checker()

        # Check if a name was specified
        if self.settings["name"] is None:
            self.name = self.settings["prop"]+"_converge"
        else:
            self.name = self.settings["name"]

        # If we have a list-type prop_start (e.g. a k-point string), we need to be careful with making our list
        if isinstance(self.settings["prop_start"], list):
            # For constructing directory names, we'll want a 3-position sub
            self.prop_name = "%i_%i_%i"
            self.prop_list = []
            i, j, k = self.settings["prop_start"]
            while [i, j, k] <= self.settings["prop_stop"]:
                self.prop_list += [[i, j, k]]
                i += self.settings["prop_step"][0]
                j += self.settings["prop_step"][1]
                k += self.settings["prop_step"][2]
        # Anything integer or float, however, is easily range'd
        elif float(self.settings["prop_start"]) != int(self.settings["prop_start"]) or float(self.settings["prop_stop"]) != int(self.settings["prop_stop"]) or float(self.settings["prop_step"]) != int(self.settings["prop_step"]):
            # For constructing directory names, we'll want a 1-position sub
            self.prop_name = "%g"
            self.prop_list = list(self.frange(self.settings["prop_start"], self.settings["prop_stop"] + self.settings["prop_step"], self.settings["prop_step"]))
        else:
            # For constructing directory names, we'll want a 1-position sub
            self.prop_name = "%i"
            self.prop_list = range(self.settings["prop_start"], self.settings["prop_stop"] + self.settings["prop_step"], self.settings["prop_step"])

        # Making all the directories
        self.prop_dir_list = []
        for prop in self.prop_list:
            if isinstance(prop, list):
                self.prop_dir_list += [os.path.join(self.calcdir, self.prop_name % tuple(prop))]
            else:
                self.prop_dir_list += [os.path.join(self.calcdir, self.prop_name % prop)]
            try:
                if isinstance(prop, list):
                    os.mkdir(os.path.join(self.calcdir, self.prop_name % tuple(prop)))
                else:
                    os.mkdir(os.path.join(self.calcdir, self.prop_name % prop))
            except OSError:
                pass

        print("  DONE\n")
        sys.stdout.flush()

    def collect(self):
        """ Collect the results of a convergence run """

        nrg_data = []
        vol_data = []
        lat_a_data = []
        lat_b_data = []
        lat_c_data = []
        # Iterate over each individual value of the property, since parallelism!
        for prop, propdir in zip(self.prop_list, self.prop_dir_list):
            # construct the Converge object
            convergence = vasp.Converge(propdir, self.run_settings(), prop)

            # check the current status
            (status, task) = convergence.status()   #pylint: disable=unused-variable

            # If any job isn't done, we should abort
            if status == "incomplete":
                print("Directory " + propdir + " Has Status:" + status + "  Please wait for the job to finish.")
                return
            elif status != "complete":
                print("Directory " + propdir + " Has Status:" + status + "  Please inspect for errors, then re-submit.")
                return
            else:
                # ensure results report written
                if not os.path.isfile(os.path.join(propdir, "properties.calc.json")):
                    ###TODO: Why did I write it like this? Why not call convergence.finalize()?
                    store_prop = self.prop
                    store_propdir = self.propdir

                    self.prop = prop
                    self.propdir = propdir

                    self.finalize()

                    self.prop = store_prop
                    self.propdir = store_propdir

                # Collect job data
                my_properties = read_properties(os.path.join(propdir, "properties.calc.json"))

                nrg_data += [float(my_properties["relaxed_energy"])/sum(my_properties["atoms_per_type"])] #Energy per atom!!!
                vol_data += [self.volume(my_properties["relaxed_lattice"])]
                lat_a, lat_b, lat_c = self.lengths(my_properties["relaxed_lattice"])    #pylint: disable=invalid-name
                lat_a_data += [lat_a]
                lat_b_data += [lat_b]
                lat_c_data += [lat_c]

        # Calculate the diffs
        nrg_diff = []
        vol_diff = []
        lat_a_diff = []
        lat_b_diff = []
        lat_c_diff = []
        prior_nrg = 0.
        prior_vol = 0.
        prior_lat_a = 0.
        prior_lat_b = 0.
        prior_lat_c = 0.
        for nrg, vol, lat_a, lat_b, lat_c in zip(nrg_data, vol_data, lat_a_data, lat_b_data, lat_c_data):
            nrg_diff += [nrg - prior_nrg]
            vol_diff += [vol - prior_vol]
            lat_a_diff += [lat_a - prior_lat_a]
            lat_b_diff += [lat_b - prior_lat_b]
            lat_c_diff += [lat_c - prior_lat_c]
            prior_nrg = nrg
            prior_vol = vol
            prior_lat_a = lat_a
            prior_lat_b = lat_b
            prior_lat_c = lat_c
        abs_nrg_diff = [abs(x) for x in nrg_diff]
        abs_vol_diff = [abs(x) for x in vol_diff]
        abs_lat_a_diff = [abs(x) for x in lat_a_diff]
        abs_lat_b_diff = [abs(x) for x in lat_b_diff]
        abs_lat_c_diff = [abs(x) for x in lat_c_diff]

        abs_nrg_diff_per = [100.*abs(x/y) for x, y in zip(abs_nrg_diff, nrg_data)]
        abs_vol_diff_per = [100.*abs(x/y) for x, y in zip(abs_vol_diff, vol_data)]
        abs_lat_a_diff_per = [100.*abs(x/y) for x, y in zip(abs_lat_a_diff, lat_a_data)]
        abs_lat_b_diff_per = [100.*abs(x/y) for x, y in zip(abs_lat_b_diff, lat_b_data)]
        abs_lat_c_diff_per = [100.*abs(x/y) for x, y in zip(abs_lat_c_diff, lat_c_data)]

        # Determine where, if anywhere, we have reached convergence
        prop_conv = None
        prop_idx = -1
        if self.settings["tol"] is None:
            pass
        else:
            # Enumerate over our specified tolerances
            for tol, tol_amount in zip(self.settings["tol"], self.settings["tol_amount"]):
                if tol.lower() == "relaxed_energy":

                    # If even the last step doesn't hit tol, we have NOT converged, and should abort this loop now!!!
                    if abs_nrg_diff[-1] > tol_amount:
                        prop_idx = -1
                        break

                    # Walk forward in steps and see when we've converged
                    for idx, nrg_diff_i in enumerate(abs_nrg_diff):

                        # Check if THIS value fulfills the tol (the value may oscillate, so this is NOT a break condition)
                        if nrg_diff_i <= tol_amount:
                            continue

                        # Increment our "global" prop index, as appropriate
                        else:
                            if idx > prop_idx:
                                prop_idx = idx

                    continue

                elif tol.lower() == "relaxed_volume":

                    # If even the last step doesn't hit tol, we have NOT converged, and should abort this loop now!!!
                    if abs_vol_diff[-1] > tol_amount:
                        prop_idx = -1
                        break

                    # Walk forward in steps and see when we've converged
                    for idx, vol_diff_i in enumerate(abs_vol_diff):

                        # Check if THIS value fulfills the tol (the value may oscillate, so this is NOT a break condition)
                        if vol_diff_i <= tol_amount:
                            continue

                        # Increment our "global" prop index, as appropriate
                        else:
                            if idx > prop_idx:
                                prop_idx = idx

                    continue

                elif tol.lower() == "relaxed_lattice":
                    # If even the last step doesn't hit tol, we have NOT converged, and should abort this loop now!!!
                    if abs_lat_a_diff[-1] > tol_amount or abs_lat_b_diff[-1] > tol_amount or abs_lat_c_diff[-1] > tol_amount:
                        prop_idx = -1
                        break

                    # Walk forward in steps and see when we've converged
                    for idx, (lat_a_diff_i, lat_b_diff_i, lat_c_diff_i) in enumerate(zip(abs_lat_a_diff, abs_lat_b_diff, abs_lat_c_diff)):

                        # Check if THIS value fulfills the tol (the value may oscillate, so this is NOT a break condition)
                        if lat_a_diff_i <= tol_amount and lat_b_diff_i <= tol_amount and lat_c_diff_i <= tol_amount:
                            continue

                        # Increment our "global" prop index, as appropriate
                        else:
                            if idx > prop_idx:
                                prop_idx = idx

                    continue

                else:
                    raise ConvergeError("Invalid tolerance type %s!" % tol)

            if prop_idx >= 0:
                # We want to do prop_idx + 1 because we stop incrementing the idx counter on the step BEFORE
                # the tol, and we want the step AT the tol
                prop_conv = self.prop_list[prop_idx + 1]

        # Print our results to convergence.calc.json
        conv_dict = {}
        conv_dict["property"] = {"property_name" : self.settings["prop"], "property_converged_value" : prop_conv}
        # conv_dict["property"] = self.settings["prop"]
        # conv_dict["property_converged_value"] = prop_conv
        conv_dict["tolerance"] = {"tolerance_type" : self.settings["tol"], "tolerance_value" : self.settings["tol_amount"]}
        # conv_dict["tolerances"] = self.settings["tol"]
        # conv_dict["tolerance_amounts"] = self.settings["tol_amount"]
        conv_dict["data"] = {
            "property_values" : self.prop_list,
            "energy" : {"relaxed" : nrg_data, "abs_relaxed_diff" : abs_nrg_diff, "percent_relaxed_diff" : abs_nrg_diff_per},
            "volume" : {"relaxed" : vol_data, "abs_relaxed_diff" : abs_vol_diff, "percent_relaxed_diff" : abs_vol_diff_per},
            "lattice_vector" : {
                "relaxed" : [[a, b, c] for a, b, c in zip(lat_a_data, lat_b_data, lat_c_data)],
                "abs_relaxed_diff" : [[a, b, c] for a, b, c in zip(abs_lat_a_diff, abs_lat_b_diff, abs_lat_c_diff)],
                "percent_relaxed_diff" : [[a, b, c] for a, b, c in zip(abs_lat_a_diff_per, abs_lat_b_diff_per, abs_lat_c_diff_per)]
                }
            }

        # conv_dict["values"] = self.prop_list
        # conv_dict["relaxed_energies"] = nrg_data
        # conv_dict["relaxed_energy_diffs"] = nrg_diff
        # conv_dict["abs_relaxed_energy_diffs"] = abs_nrg_diff
        # conv_dict["relaxed_volumes"] = vol_data
        # conv_dict["relaxed_volume_diffs"] = vol_diff
        # conv_dict["abs_relaxed_volume_diffs"] = abs_vol_diff
        # conv_dict["relaxed_lattice_vector_lengths"] = [[a, b, c] for a, b, c in zip(lat_a_data, lat_b_data, lat_c_data)]
        # conv_dict["relaxed_lattice_vector_length_diffs"] = [[a, b, c] for a, b, c in zip(lat_a_diff, lat_b_diff, lat_c_diff)]
        # conv_dict["abs_relaxed_lattice_vector_length_diffs"] = [[a, b, c] for a, b, c in zip(abs_lat_a_diff, abs_lat_b_diff, abs_lat_c_diff)]
        # write properties.calc.json

        conv_file_name = "convergence.calc.json"
        outputfile = os.path.join(self.calcdir, conv_file_name)
        with open(outputfile, 'wb') as my_file:
            my_file.write(six.u(json.dumps(conv_dict, cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)).encode('utf-8'))
        print("Wrote " + outputfile)
        sys.stdout.flush()

        # Copy our converge.json so we can revisit it later
        write_settings(self.settings, os.path.join(self.calcdir, "converge.json"))
        print("Wrote " + os.path.join(self.calcdir, "converge.json"))
        sys.stdout.flush()

    def setup(self):
        """ Setup initial convergence run

            Uses the following files from the most local .../settings/calctype.name directory:
                INCAR: VASP input settings
                KPOINTS: VASP kpoints settings
                POSCAR: reference for KPOINTS if KPOINTS mode is not A/AUTO/Automatic
                SPECIES: info for each species such as which POTCAR files to use, MAGMOM, GGA+U, etc.

            Uses the following files from the .../config directory:
                POS: structure of the configuration to be converged

        """
        # Find required input files in CASM project directory tree
        vaspfiles = vasp_input_file_names(self.casm_directories, self.configname, self.clex)
        incarfile, kpointsfile, _, poscarfile, speciesfile = vaspfiles

        # Verify that required input files exist
        if incarfile is None:
            raise vasp.VaspError("Converge.setup failed. No INCAR file found in CASM project.")
        if kpointsfile is None:
            raise vasp.VaspError("Converge.setup failed. No KPOINTS file found in CASM project.")
        if poscarfile is None:
            raise vasp.VaspError("Converge.setup failed. No POS file found for this configuration.")
        if speciesfile is None:
            raise vasp.VaspError("Converge.setup failed. No SPECIES file found in CASM project.")

        # Find optional input files
        extra_input_files = []
        for s in self.settings["extra_input_files"]:
            extra_input_files.append(self.casm_directories.settings_path_crawl(s, self.configname, self.clex))
            if extra_input_files[-1] is None:
                raise vasp.VaspError("Converge.setup failed. Extra input file " + s + " not found in CASM project.")
        if self.settings["initial"]:
            extra_input_files += [ self.casm_directories.settings_path_crawl(self.settings["initial"], self.configname, self.clex) ]
            if extra_input_files[-1] is None:
                raise vasp.VaspError("Converge.setup failed. No initial INCAR file " + self.settings["initial"] + " found in CASM project.")
        if self.settings["final"]:
            extra_input_files += [ self.casm_directories.settings_path_crawl(self.settings["final"], self.configname, self.clex) ]
            if extra_input_files[-1] is None:
                raise vasp.VaspError("Converge.setup failed. No final INCAR file " + self.settings["final"] + " found in CASM project.")
        sys.stdout.flush()

        print("Setting up VASP input files:", self.propdir)

        # read prim and prim kpoints
        print("  Reading KPOINTS:", kpointsfile)
        kpoints = vasp.io.Kpoints(kpointsfile)

        # read species, super poscar, incar, and generate super kpoints
        print("  Reading SPECIES:", speciesfile)
        species_settings = vasp.io.species_settings(speciesfile)
        print("  Reading supercell POS:", poscarfile)
        poscar = vasp.io.Poscar(poscarfile, species_settings)
        print("  Reading INCAR:", incarfile)
        incar = vasp.io.Incar(incarfile, species_settings, poscar, True)


        # Modify as appropriate
        if self.settings["prop"].upper() == "ENCUT":
            try:
                incar.tags["ENCUT"] = int(self.prop)
            except:
                raise ConvergeError("Error in Cconverge.collect: something has gone wrong and the run-specific property %s could not be cast as int!" % self.prop)
        elif self.settings["prop"].upper() == "NBANDS":
            try:
                incar.tags["NBANDS"] = int(self.prop)
            except:
                raise ConvergeError("Error in Converge.collect: something has gone wrong and the run-specific property %s could not be cast as int!" % self.prop)
        elif self.settings["prop"].upper() == "SIGMA":
            try:
                incar.tags["SIGMA"] = float(self.prop)
            except:
                raise ConvergeError("Error in Converge.collect: something has gone wrong and the run-specific property %s could not be cast as int!" % self.prop)
        elif self.settings["prop"].upper() == "KPOINTS":
            if isinstance(self.prop, list):
                if kpoints.automode[0].lower() == "a":
                    raise ConvergeError("Error in Converge.collect: An \"Automatic\" k-point generation scheme was supplied in %s, but you gave me 3 k-points values: %i, %i, %i!" % (kpointsfile, self.prop[0], self.prop[1], self.prop[2]))
                kpoints.subdivisions = self.prop
            else:
                if kpoints.automode[0].lower() == "a":
                    kpoints.subdivisions = self.prop
                else:
                    kpoints.subdivisions = [self.prop, self.prop, self.prop]

        else:
            raise ConvergeError("Error in Converge.collect: \"prop: %s\" not a valid convergence prop type!\nCurrently supported convergence prop types are %s" % (self.settings["prop"], VALID_PROP_TYPES))

        # write main input files
        print("  Writing supercell POSCAR:", os.path.join(self.propdir, 'POSCAR'))
        poscar.write(os.path.join(self.propdir, 'POSCAR'), True)
        print("  Writing INCAR:", os.path.join(self.propdir, 'INCAR'))
        incar.write(os.path.join(self.propdir, 'INCAR'))
        print("  Writing supercell KPOINTS:", os.path.join(self.propdir, 'KPOINTS'))
        kpoints.write(os.path.join(self.propdir, 'KPOINTS'))
        print("  Writing POTCAR:", os.path.join(self.propdir, 'POTCAR'))
        vasp.io.write_potcar(os.path.join(self.propdir, 'POTCAR'), poscar, species_settings, True)

        # copy extra input files
        print("  Copying extra input files", end=' ')
        for my_input_file in extra_input_files:
            print(my_input_file, end=' ')
            shutil.copy(my_input_file, self.propdir)
        print("")

        print("VASP input files complete\n")
        sys.stdout.flush()

    def submit(self):
        """Submit jobs for this VASP convergence"""

        print("Submitting...")
        print("Configuration:", self.configname)
        print("Property:", self.settings["prop"])
        # Iterate over each individual value of the property, since parallelism!
        for prop, propdir in zip(self.prop_list, self.prop_dir_list):
            try:
                # Store our propdir before editing it
                store_propdir = self.propdir
                store_prop = self.prop
                self.propdir = propdir
                self.prop = prop
                print("Value:", self.prop)
                # first, check if the job has already been submitted and is not completed
                db = JobDB()    #pylint: disable=invalid-name
                print("Calculation directory:", propdir)
                jobid = db.select_regex_id("rundir", propdir)
                print("JobID:", jobid)
                sys.stdout.flush()
                if jobid != []:
                    for j in jobid:
                        job = db.select_job(j)
                        # taskstatus = ["Incomplete", "Complete", "Continued", "Check", "Error:.*", "Aborted"]
                        # jobstatus = ["C", "Q", "R", "E", "W", "H", "M"]
                        if job["jobstatus"] != "C":
                            print("JobID:", job["jobid"], "  Jobstatus:", job["jobstatus"], "  Not submitting.")
                            sys.stdout.flush()
                            raise UserWarning()
                        #elif job["taskstatus"] in ["Complete", "Check"] or re.match("Error:.*", job["taskstatus"]):
                        #    print "JobID:", job["jobid"], "  Taskstatus:", job["taskstatus"], "  Not submitting."
                        #    sys.stdout.flush()
                        #    return


                # second, only submit a job if convergence status is "incomplete"

                # construct the Converge object
                convergence = vasp.Converge(propdir, self.run_settings(), prop)
                # convergence = vasp.Converge(self.calcdir, self.run_settings())

                # check the current status
                (status, task) = convergence.status()

                if status == "complete":
                    print("Status:", status, "  Not submitting.")
                    sys.stdout.flush()

                    # ensure job marked as complete in db
                    if self.auto:
                        for j in jobid:
                            job = db.select_job(j)
                            if job["taskstatus"] == "Incomplete":
                                try:
                                    complete_job(jobid=j)
                                except (JobsError, JobDBError, EligibilityError) as my_error:
                                    print(str(my_error))
                                    sys.stdout.flush()

                    # ensure results report written
                    if not os.path.isfile(os.path.join(propdir, "properties.calc.json")):
                        self.finalize()

                    continue

                elif status == "not_converging":
                    print("Status:", status, "  Not submitting.")
                    sys.stdout.flush()
                    continue

                elif status != "incomplete":
                    raise VaspWrapperError("unexpected convergence status: '" + status + "' and task: '" + task + "'")
                    # sys.stdout.flush()        ### Liz: Why is this here? Its unreachable.
                    # continue


                print("Preparing to submit a VASP Convergence job")
                sys.stdout.flush()

                # cd to configdir, submit jobs from configdir, then cd back to currdir
                currdir = os.getcwd()
                os.chdir(propdir)

                # determine the number of atoms in the configuration
                print("  Counting atoms in the POSCAR")
                sys.stdout.flush()
                pos = vasp.io.Poscar(os.path.join(self.configdir, "POS"))
                N = len(pos.basis) #pylint: disable=invalid-name

                # Construct the run command
                cmd = ""
                if self.settings["preamble"] is not None:
                # Append any instructions given in the 'preamble' file, if given
                    preamble = self.casm_directories.settings_path_crawl(self.settings["preamble"], self.configname, self.clex)
                    with open(preamble, 'rb') as my_preamble:
                        cmd += "".join(my_preamble.read().decode('utf-8'))
                # Or just execute a single prerun line, if given
                if self.settings["prerun"] is not None:
                    cmd += self.settings["prerun"] + "\n"
                cmd += "python -c \"import casm.vaspwrapper; casm.vaspwrapper.Converge(configdir='" + self.configdir + "', propdir='" + propdir + "', prop=" + str(prop) + ").run()\"\n"

                if self.settings["postrun"] is not None:
                    cmd += self.settings["postrun"] + "\n"

                print("  Constructing a job")
                sys.stdout.flush()
                
                # construct a Job
                job = Job(name=wrapper.jobname(self.configname) + "_" + '.'.join(propdir.split(os.sep)[-2:]),\
                              account=self.settings["account"],\
                              nodes=int(ceil(float(N)/float(self.settings["atom_per_proc"])/float(self.settings["ppn"]))),\
                              ppn=int(self.settings["ppn"]),\
                              walltime=self.settings["walltime"],\
                              pmem=self.settings["pmem"],\
                              qos=self.settings["qos"],\
                              queue=self.settings["queue"],\
                              message=self.settings["message"],\
                              email=self.settings["email"],\
                              priority=self.settings["priority"],\
                              command=cmd,\
                              auto=self.auto)

                print("  Submitting")
                sys.stdout.flush()
                # submit the job
                job.submit()
                self.report_status("submitted")

                # return to current directory and restore settings
                os.chdir(currdir)
                self.propdir = store_propdir
                self.prop = store_prop

            except UserWarning:
                continue

        print("CASM VASPWrapper Convergence job submission complete\n")
        sys.stdout.flush()


    def run_settings(self):
        """ Set default values based on runtime environment"""
        settings = dict(self.settings)

        # set default values

        if settings["npar"] == "CASM_DEFAULT":
            if "PBS_NUM_NODES" in os.environ:
                settings["npar"] = int(os.environ["PBS_NUM_NODES"])
            elif "SLURM_JOB_NUM_NODES" in os.environ:
                settings["npar"] = int(os.environ["SLURM_JOB_NUM_NODES"])
            else:
                settings["npar"] = None
        elif settings["npar"] == "VASP_DEFAULT":
            settings["npar"] = None

        if settings["npar"] is None:
            if settings["ncore"] == "CASM_DEFAULT":
                if "PBS_NUM_PPN" in os.environ:
                    settings["ncore"] = int(os.environ["PBS_NUM_PPN"])
                elif "SLURM_CPUS_ON_NODE" in os.environ:
                    settings["ncore"] = int(os.environ["SLURM_CPUS_ON_NODE"])
                else:
                    settings["ncore"] = None
            elif settings["ncore"] == "VASP_DEFAULT":
                settings["ncore"] = 1
        else:
            settings["ncore"] = None

        if settings["ncpus"] is None or settings["ncpus"] == "CASM_DEFAULT":
            if "PBS_NP" in os.environ:
                settings["ncpus"] = int(os.environ["PBS_NP"])
            elif "SLURM_NTASKS" in os.environ:
                settings["ncpus"] = int(os.environ["SLURM_NTASKS"])
            else:
                settings["ncpus"] = None

        if settings["run_limit"] is None or settings["run_limit"] == "CASM_DEFAULT":
            settings["run_limit"] = 10

        return settings


    def run(self):
        """ Setup input files, run a vasp convergence, and report results """

        # construct the Converge object
        convergence = vasp.Converge(self.propdir, self.run_settings(), self.prop)

        # check the current status
        (status, task) = convergence.status()


        if status == "complete":
            print("Status:", status)
            sys.stdout.flush()

            # mark job as complete in db
            if self.auto:
                try:
                    complete_job()
                except (JobsError, JobDBError, EligibilityError) as my_error:
                    print(str(my_error))
                    sys.stdout.flush()

            # write results to properties.calc.json
            self.finalize()
            return

        elif status == "not_converging":
            print("Status:", status)
            self.report_status("failed", "run_limit")
            print("Returning")
            sys.stdout.flush()
            return

        elif status == "incomplete":

            if task == "setup":
                self.setup()

            self.report_status("started")
            (status, task) = convergence.run()

        else:
            self.report_status("failed", "unknown")
            raise VaspWrapperError("unexpected convergence status: '" + status + "' and task: '" + task + "'")
            # sys.stdout.flush()


        # once the run is done, update database records accordingly

        if status == "not_converging":

            # mark error
            if self.auto:
                try:
                    error_job("Not converging")
                except (JobsError, JobDBError) as my_error:
                    print(str(my_error))
                    sys.stdout.flush()

            print("Not Converging!")
            sys.stdout.flush()
            self.report_status("failed", "run_limit")

            # print a local settings file, so that the run_limit can be extended if the
            #   convergence problems are fixed
            try:
                os.makedirs(os.path.join(self.configdir, "settings", self.casm_settings["curr_calctype"], self.settings["prop"]+"_converge"))
            except OSError:
                pass
            settingsfile = os.path.join(self.configdir, "settings", self.casm_settings["curr_calctype"], self.settings["prop"] + "_converge", "converge.json")
            write_settings(self.settings, settingsfile)

            print("Writing:", settingsfile)
            print("Edit the 'run_limit' property if you wish to continue.")
            sys.stdout.flush()
            return

        elif status == "complete":

            # mark job as complete in db
            if self.auto:
                try:
                    complete_job()
                except (JobsError, JobDBError, EligibilityError) as my_error:
                    print(str(my_error))
                    sys.stdout.flush()

            # write results to properties.calc.json
            self.finalize()

        else:
            self.report_status("failed", "unknown")
            raise VaspWrapperError("vasp convergence complete with unexpected status: '" + status + "' and task: '" + task + "'")
            #sys.stdout.flush()

    def report_status(self, status, failure_type=None):
        """Report calculation status to status.json file in configuration directory.

        Args:
            status: string describing calculation status. Currently used values are
                 not_submitted
                 submitted
                 complete
                 failed
             failure_type: optional string describing reason for failure. Currently used values are
                 unknown
                 electronic_convergence
                 run_limit"""

        output = dict()
        output["status"] = status
        if failure_type is not None:
            output["failure_type"] = failure_type

        outputfile = os.path.join(self.propdir, "status.json")
        with open(outputfile, 'wb') as my_file:
            my_file.write(six.u(json.dumps(output, cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)).encode('utf-8'))
        print("Wrote " + outputfile)
        sys.stdout.flush()

    def finalize(self):
        """ Checks if a series of runs is converged, and writes properties.calc.json as appropriate """
        if self.is_converged():
            # write properties.calc.json
            vaspdir = os.path.join(self.propdir, "run.final")
            super_poscarfile = os.path.join(self.configdir, "POS")
            speciesfile = self.casm_directories.settings_path_crawl("SPECIES", self.configname, self.clex)
            output = self.properties(vaspdir, super_poscarfile, speciesfile)
            prop_file_name = "properties.calc.json"
            outputfile = os.path.join(self.propdir, prop_file_name)
            with open(outputfile, 'wb') as my_file:
                my_file.write(six.u(json.dumps(output, cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)).encode('utf-8'))
            print("Wrote " + outputfile)
            sys.stdout.flush()
            self.report_status('complete')

    def is_converged(self):
        """Check for electronic convergence in completed calculations. Returns True or False."""

        # Verify that the last relaxation reached electronic convergence
        convergence = vasp.Converge(self.propdir, self.run_settings())
        for i in range(len(convergence.rundir)):
            try:
                print(self.propdir)
                vrun = vasp.io.Vasprun(os.path.join(self.propdir, convergence.rundir[-i-1], "vasprun.xml"))
                if len(vrun.all_e_0[-1]) >= vrun.nelm:
                    print(('The last convergence run (' +
                          os.path.basename(convergence.rundir[-i-1]) +
                          ') failed to achieve electronic convergence; properties.calc.json will not be written.\n'))
                    self.report_status('failed', 'electronic_convergence')
                    return False
                break
            ### WHY IS THIS HERE
            except: #pylint: disable=bare-except
                pass

        # Verify that the final static run reached electronic convergence
        vrun = vasp.io.Vasprun(os.path.join(self.propdir, "run.final", "vasprun.xml"))
        if len(vrun.all_e_0[0]) >= vrun.nelm:
            print('The final run failed to achieve electronic convergence; properties.calc.json will not be written.\n')
            self.report_status('failed', 'electronic_convergence')
            return False

        return True

    def _settings_checker(self):
        """ Checks that settings given in converge.json are valid """
        # Check for required settings having sane values
        try:
            if not self.settings["prop"].upper() in VALID_PROP_TYPES:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"prop: %s\" not a valid convergence prop type!\nCurrently supported convergence prop types are %s" % (self.settings["prop"], VALID_PROP_TYPES))
        except:
            raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"prop: %s\" missing, or could not be converted to a string!\nCurrently supported convergence prop types are %s" % (self.settings["prop"], VALID_PROP_TYPES))
        # ENCUT requires an int for start, stop, and step
        if self.settings["prop"].upper() == "ENCUT":
            try:
                self.settings["prop_start"] = int(self.settings["prop_start"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain integer prop_start for prop ENCUT. I found: %s" % self.settings["prop_start"])
            try:
                self.settings["prop_step"] = int(self.settings["prop_step"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain integer prop_step for prop ENCUT. I found: %s" % self.settings["prop_step"])
            try:
                self.settings["prop_stop"] = int(self.settings["prop_stop"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain integer prop_stop for prop ENCUT. I found: %s" % self.settings["prop_step"])
        # NBANDS requires an int for start, stop, and step
        elif self.settings["prop"].upper() == "NBANDS":
            try:
                self.settings["prop_start"] = int(self.settings["prop_start"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain integer prop_start for prop NBANDS. I found: %s" % self.settings["prop_start"])
            try:
                self.settings["prop_step"] = int(self.settings["prop_step"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain integer prop_step for prop NBANDS. I found: %s" % self.settings["prop_step"])
            try:
                self.settings["prop_stop"] = int(self.settings["prop_stop"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain integer prop_stop for prop NBANDS. I found: %s" % self.settings["prop_step"])
        # SIGMA requires a float for start, stop, and step
        elif self.settings["prop"].upper() == "SIGMA":
            try:
                self.settings["prop_start"] = float(self.settings["prop_start"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain float prop_start for prop SIGMA. I found: %s" % self.settings["prop_start"])
            try:
                self.settings["prop_step"] = float(self.settings["prop_step"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain float prop_step for prop SIGMA. I found: %s" % self.settings["prop_step"])
            try:
                self.settings["prop_stop"] = float(self.settings["prop_stop"])
            except:
                raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain float prop_stop for prop SIGMA. I found: %s" % self.settings["prop_step"])
        # KPOINTS requires either an int or a length-3 list of ints for start, stop, and step
        elif self.settings["prop"].upper() == "KPOINTS":
            try:
                self.settings["prop_start"] = int(self.settings["prop_start"])
            except TypeError:
                try:
                    self.settings["prop_start"] = [int(k) for k in self.settings["prop_start"]]
                except:
                    raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain integer or 3-list of integer prop_start for prop KPOINTS. I found: %s" % self.settings["prop_start"])
            if isinstance(self.settings["prop_start"], list):
                try:
                    self.settings["prop_step"] = [int(k) for k in self.settings["prop_step"]]
                except:
                    raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain prop_step of the same type as prop_start for prop KPOINTS. I found: %s" % self.settings["prop_step"])
                try:
                    self.settings["prop_stop"] = [int(k) for k in self.settings["prop_stop"]]
                except:
                    raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain prop_stop of the same type as prop_start for prop KPOINTS. I found: %s" % self.settings["prop_stop"])
            else:
                try:
                    self.settings["prop_step"] = int(self.settings["prop_step"])
                except:
                    raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain prop_step of the same type as prop_start for prop KPOINTS. I found: %s" % self.settings["prop_step"])
                try:
                    self.settings["prop_stop"] = int(self.settings["prop_stop"])
                except:
                    raise ConvergeError("Error in casm.vaspwrapper.Converge(): converge.json must contain prop_stop of the same type as prop_start for prop KPOINTS. I found: %s" % self.settings["prop_stop"])

        # If "tol" is present, check for a valid tol type and "tol_amount" value
        if self.settings["tol"] is not None:
            if isinstance(self.settings["tol"], list):
                if len(self.settings["tol"]) != len(self.settings["tol_amount"]) or not isinstance(self.settings["tol_amount"], list):
                    raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"tol\" and \"tol_amount\" must have the same number of entries! However, you have %s and %s" % (str(self.settings["tol"]), str(self.settings["tol_amount"])))
                for my_tol, my_tol_amount in zip(self.settings["tol"], self.settings["tol_amount"]):
                    try:
                        if not my_tol.lower() in VALID_TOL_TYPES:
                            raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"tol: %s\" not a valid convergence tolerance type!\nCurrently supported convergence tolerance types are %s" % (self.settings["tol"], VALID_TOL_TYPES))
                    except:
                        raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"tol: %s\" not a valid convergence tolerance type!\nCurrently supported convergence tolerance types are %s" % (self.settings["tol"], VALID_TOL_TYPES))
                    try:
                        my_tol_amount = abs(float(my_tol_amount))
                    except:
                        raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"tol_amount: %s\" cannot be converted to float, but a float is needed!" % self.settings["tol_amount"])
            else:
                try:
                    if not self.settings["tol"].lower() in VALID_TOL_TYPES:
                        raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"tol: %s\" not a valid convergence tolerance type!\nCurrently supported convergence tolerance types are %s" % (self.settings["tol"], VALID_TOL_TYPES))
                    else:
                        self.settings["tol"] = [self.settings["tol"]]
                except:
                    raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"tol: %s\" not a valid convergence tolerance type!\nCurrently supported convergence tolerance types are %s" % (self.settings["tol"], VALID_TOL_TYPES))

                try:
                    self.settings["tol_amount"] = [abs(float(self.settings["tol_amount"]))]
                except:
                    raise ConvergeError("Error in casm.vaspwrapper.Converge(): \"tol_amount: %s\" cannot be converted to float, but a float is needed!" % self.settings["tol_amount"])

    @property
    def configdir(self):
        """ Produces a configname from the configdir """
        return self.casm_directories.configuration_dir(self.configname)

    @staticmethod
    def properties(vaspdir, super_poscarfile=None, speciesfile=None):
        """Report results to properties.calc.json file in configuration directory, after checking for electronic convergence."""

        output = dict()
        vrun = vasp.io.Vasprun(os.path.join(vaspdir, "vasprun.xml"))

        # the calculation is run on the 'sorted' POSCAR, need to report results 'unsorted'

        if (super_poscarfile is not None) and (speciesfile is not None):
            species_settings = vasp.io.species_settings(speciesfile)
            super_poscar = vasp.io.Poscar(super_poscarfile, species_settings)
            unsort_dict = super_poscar.unsort_dict()
        else:
            # fake unsort_dict (unsort_dict[i] == i)
            unsort_dict = dict(zip(range(0, len(vrun.basis)), range(0, len(vrun.basis))))
        super_poscar = vasp.io.Poscar(os.path.join(vaspdir, "POSCAR"))

        # unsort_dict:
        #   Returns 'unsort_dict', for which: unsorted_dict[orig_index] == sorted_index;
        #   unsorted_dict[sorted_index] == orig_index
        #   For example:
        #     'unsort_dict[0]' returns the index into the unsorted POSCAR of the first atom in the sorted POSCAR


        output["atom_type"] = super_poscar.type_atoms
        output["atoms_per_type"] = super_poscar.num_atoms
        output["coord_mode"] = vrun.coord_mode

        # as lists
        output["relaxed_forces"] = [None for i in range(len(vrun.forces))]
        for i, v in enumerate(vrun.forces): #pylint: disable=invalid-name
            output["relaxed_forces"][unsort_dict[i]] = noindent.NoIndent(vrun.forces[i])

        output["relaxed_lattice"] = [noindent.NoIndent(v) for v in vrun.lattice]

        output["relaxed_basis"] = [None for i in range(len(vrun.basis))]
        for i, v in enumerate(vrun.basis):  #pylint: disable=invalid-name
            output["relaxed_basis"][unsort_dict[i]] = noindent.NoIndent(vrun.basis[i])

        output["relaxed_energy"] = vrun.total_energy

        return output

    @staticmethod
    def volume(lat):
        """ Computes the volume of a parallelpiped given a 3-list of 3-lists (i.e., 3 3-vectors)"""
        try:
            return (
                lat[0][0]*(lat[1][1]*lat[2][2] - lat[1][2]*lat[2][1])
                - lat[0][1]*(lat[1][0]*lat[2][2] - lat[1][2]*lat[2][0])
                + lat[0][2]*(lat[1][0]*lat[2][1] - lat[1][1]*lat[2][0]))
        except:
            raise RuntimeError("The given lattice %s could not be parsed as a 3-list of 3-vectors." % str(lat))

    @staticmethod
    def lengths(lat):
        """ Computes the lengths of the 3 basis vectors in a lattice """
        try:
            return [
                sqrt(lat[0][0]**2 + lat[0][1]**2 + lat[0][2]**2),
                sqrt(lat[1][0]**2 + lat[1][1]**2 + lat[1][2]**2),
                sqrt(lat[2][0]**2 + lat[2][1]**2 + lat[2][2]**2)
                ]
        except:
            raise RuntimeError("The given lattice %s could not be parsed as a 3-list of 3-vectors." % str(lat))

    @staticmethod
    def frange(x, y, jump=1.0): #pylint: disable=invalid-name
        """Robust Range for floats, without using numpy."""
        i = 0.0
        x = float(x)  # Prevent yielding integers.
        x0 = x  #pylint: disable=invalid-name
        epsilon = jump / 2.0
        yield x  # yield always first value
        while x + epsilon < y:
            i += 1.0
            x = x0 + i * jump
            yield x
