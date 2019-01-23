from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import math
import os
import re
import six
import sys
import warnings

try:
  from prisms_jobs import Job, JobDB, error_job, complete_job, JobsError, JobDBError, EligibilityError
except ImportError:
  # use of the pbs module is deprecated after CASM v0.2.1
  from pbs import Job, JobDB, error_job, complete_job, JobDBError, EligibilityError
  from pbs import PBSError as JobsError

from casm import quantumespresso, wrapper
from casm.misc import noindent
from casm.project import DirectoryStructure, ProjectSettings
from casm.qewrapper import QEWrapperError, qe_input_file_names, read_settings, write_settings

class Relax(object):
    """The Relax class contains functions for setting up, executing, and parsing a Quantum Espresso relaxation.

        The relaxation creates the following directory structure:
        config/
          calctype.name/
              run.0/
              run.1/
              ...
              run.final/

        'run.i' directories are only created when ready.
        'run.final' is a final constant volume run {"calculation":'relax'}

        This automatically looks for Quantum Espresso settings files in .../settings/calctype.name,
        where '...' is the nearest parent directory of 'self.configdir' in the CASM project repository

        Contains:
            self.configdir (.../config)
            self.calcdir   (.../config/calctype.name)

            self.settings = dictionary of settings for job submission and the relaxation, see qewrapper.read_settings

            self.auto = True if using prisms_jobs module's JobDB to manage jobs
            self.sort = True if sorting atoms in POSCAR by type
    """
    def __init__(self, configdir=None, auto = True, sort = True):
        """
        Construct a Quantum Espresso relaxation job object.

        Args:
            configdir: path to configuration
            auto: True if using prisms_jobs module's JobDB to manage jobs

        """
        if configdir == None:
            configdir = os.getcwd()

        print("Working on directory "+str(configdir))

        # get the configname from the configdir path
        _res = os.path.split(configdir)
        self.configname = os.path.split(_res[0])[1] + "/" + _res[1]
        print("  Configuration:", self.configname)

        print("Reading CASM settings")
        self.casm_settings = ProjectSettings(configdir)
        if self.casm_settings == None:
            raise QEWrapperError("Not in a CASM project. The file '.casm' directory was not found.")

        self.casm_directories=DirectoryStructure(configdir)

        print("Constructing a CASM QEWrapper Relax object")
        sys.stdout.flush()

        print("  Setting up directories")
        sys.stdout.flush()

        # store path to .../config, if not existing raise
        self.configdir = os.path.abspath(configdir)
        if not os.path.isdir(self.configdir):
            raise quantumespresso.QuantumEspressoError("Error in casm.qewrapper.Relax: Did not find directory: " + self.configdir)
            sys.stdout.flush()

        # store path to .../config/calctype.name, and create if not existing
        self.calcdir = self.casm_directories.calctype_dir(self.configname,self.casm_settings.default_clex)
        try:
            os.mkdir(self.calcdir)
        except:
            pass


        # read the settings json file
        print("  Reading relax.json settings file")
        sys.stdout.flush()
        setfile = self.casm_directories.settings_path_crawl("relax.json",self.configname,self.casm_settings.default_clex)

        if setfile == None:
            raise QEWrapperError("Could not find \"relax.json\" in an appropriate \"settings\" directory")
            sys.stdout.flush()

        else:
            print("Using "+str(setfile)+" as settings...")

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
        if not "strict_kpoints" in self.settings:
            self.settings["strict_kpoints"]=False
        if not "infilename" in self.settings:
            print("WARNING: No input file specified in relax.json using default infilename of std.in")
            self.settings["infilename"]="std.in"
        if not "outfilename" in self.settings:
            print("WARNING: No output file specified in relax.json using default outfilename of std.out")
            self.settings["outfilename"]="std.out"


        self.auto = auto
        self.sort = sort
        print("Quantum Espresso Relax object constructed\n")
        sys.stdout.flush()


    def setup(self):
        """ Setup initial relaxation run

            Uses the following files from the most local .../settings/calctype.name directory:
                infile: Quantum Espresso input settings
                SPECIES: info for each species such as which UPF files to use, atomic specific tags etc GGA+U, etc.

            Uses the following files from the .../config directory:
                POS: structure of the configuration to be relaxed

        """
        # Find required input files in CASM project directory tree
        infilename=self.settings["infilename"]
        qefiles=qe_input_file_names(self.casm_directories,self.configname,self.casm_settings.default_clex,infilename)
        infilename,super_poscarfile,speciesfile=qefiles


        # Find optional input files
        extra_input_files = []
        for s in self.settings["extra_input_files"]:
            extra_input_files.append(self.casm_directories.settings_path_crawl(s,self.configname,self.casm_settings.default_clex))
            if extra_input_files[-1] is None:
                raise quantumespresso.QuantumEspressoError("Relax.setup failed. Extra input file " + s + " not found in CASM project.")
        if self.settings["initial"]:
            extra_input_files += [ self.casm_directories.settings_path_crawl(self.settings["initial"],self.configname,self.casm_settings.default_clex) ]
            if extra_input_files[-1] is None:
                raise quantumespresso.QuantumEspressoError("Relax.setup failed. No initial  Infile " + self.settings["initial"] + " found in CASM project.")
        if self.settings["final"]:
            extra_input_files += [ self.casm_directories.settings_path_crawl(self.settings["final"],self.configname,self.casm_settings.default_clex) ]
            if extra_input_files[-1] is None:
                raise quantumespresso.QuantumEspressoError("Relax.setup failed. No final Infile " + self.settings["final"] + " found in CASM project.")


        sys.stdout.flush()

        quantumespresso.qeio.write_quantum_espresso_input(self.calcdir, infilename, super_poscarfile, speciesfile, self.sort, extra_input_files,self.settings["strict_kpoints"])


    def submit(self):
        """Submit a PBS job for this Quantum Espresso relaxation"""

        # first, check if the job has already been submitted and is not completed
        db = JobDB()
        print("rundir", self.calcdir)
        id = db.select_regex_id("rundir", self.calcdir)
        print("id:", id)
        sys.stdout.flush()
        if id != []:
            db.update()
            for j in id:
                job = db.select_job(j)
                # taskstatus = ["Incomplete","Complete","Continued","Check","Error:.*","Aborted"]
                # jobstatus = ["C","Q","R","E","W","H","M"]
                if job["jobstatus"] != "C":
                    print("JobID:", job["jobid"], "  Jobstatus:", job["jobstatus"], "  Not submitting.")
                    sys.stdout.flush()
                    return
                #elif job["taskstatus"] in ["Complete", "Check"] or re.match( "Error:.*", job["taskstatus"]):
                #    print "JobID:", job["jobid"], "  Taskstatus:", job["taskstatus"], "  Not submitting."
                #    sys.stdout.flush()
                #    return


        # second, only submit a job if relaxation status is "incomplete"

        # construct the Relax object
        relaxation = quantumespresso.Relax(self.calcdir, self.run_settings())

        # load infilename and outfilename
        infilename=self.settings["infilename"]
        outfilename=self.settings["outfilename"]

        # check the current status
        (status, task) = relaxation.status()

        if status == "complete":
            print("Status:", status, "  Not submitting.")
            sys.stdout.flush()

            # ensure job marked as complete in db
            if self.auto:
                for j in id:
                  job = db.select_job(j)
                  if job["taskstatus"] == "Incomplete":
                      try:
                          complete_job(jobid=j)
                      except (JobsError, JobDBError, EligibilityError) as e:
                          print(str(e))
                          sys.stdout.flush()

            # ensure results report written
            if not os.path.isfile(os.path.join(self.calcdir, "properties.calc.json")):
                self.finalize()

            return

        elif status == "not_converging":
            print("Status:", status, "  Not submitting.")
            sys.stdout.flush()
            return

        elif status != "incomplete":
            raise QEWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()
            return


        print("Preparing to submit a Quantum Espresso relaxation PBS job")
        sys.stdout.flush()

        # cd to configdir, submit jobs from configdir, then cd back to currdir
        currdir = os.getcwd()
        os.chdir(self.calcdir)

        # determine the number of atoms in the configuration
        print("  Counting atoms in the POSCAR")
        sys.stdout.flush()
        pos = quantumespresso.qeio.Poscar(os.path.join(self.configdir,"POS"))
        N = len(pos.basis)

        print("  Constructing a PBS job")
        sys.stdout.flush()
        # construct a Job
        job = Job(name=wrapper.jobname(self.configname),\
                      account=self.settings["account"],\
                      nodes=int(math.ceil(float(N)/float(self.settings["atom_per_proc"])/float(self.settings["ppn"]))),\
                      ppn=int(self.settings["ppn"]),\
                      walltime=self.settings["walltime"],\
                      pmem=self.settings["pmem"],\
                      qos=self.settings["qos"],\
                      queue=self.settings["queue"],\
                      message=self.settings["message"],\
                      email=self.settings["email"],\
                      priority=self.settings["priority"],\
                      command="python -c \"import casm.qewrapper; casm.qewrapper.Relax('" + self.configdir + "').run()\"",\
                      auto=self.auto)

        print("  Submitting")
        sys.stdout.flush()
        # submit the job
        job.submit()
        self.report_status("submitted")

        # return to current directory
        os.chdir(currdir)

        print("CASM QEWrapper relaxation PBS job submission complete\n")
        sys.stdout.flush()


    def run_settings(self):
        """ Set default values based on runtime environment"""
        settings = dict(self.settings)

        # set default values

        if settings["npar"] == "CASM_DEFAULT":
            if "PBS_NUM_NODES" in os.environ:
                settings["npar"] = int(os.environ["PBS_NUM_NODES"])
            else:
                settings["npar"] = None

        if settings["npar"] == None:
            if settings["ncore"] == "CASM_DEFAULT":
                if "PBS_NUM_PPN" in os.environ:
                    settings["ncore"] = int(os.environ["PBS_NUM_PPN"])
                else:
                    settings["ncore"] = None
        else:
            settings["ncore"] = None

        if settings["ncpus"] == None or settings["ncpus"] == "CASM_DEFAULT":
            if "PBS_NP" in os.environ:
                settings["ncpus"] = int(os.environ["PBS_NP"])
            else:
                settings["ncpus"] = None

        if settings["run_limit"] == None or settings["run_limit"] == "CASM_DEFAULT":
            settings["run_limit"] = 10

        return settings


    def run(self):
        """ Setup input files, run a qe relaxation, and report results """

        # load infilename and outfilename
        infilename=self.settings["infilename"]
        outfilename=self.settings["outfilename"]
        
        # construct the Relax object
        relaxation = quantumespresso.Relax(self.calcdir, self.run_settings())

        # check the current status
        (status, task) = relaxation.status()


        if status == "complete":
            print("Status:", status)
            sys.stdout.flush()

            # mark job as complete in db
            if self.auto:
                try:
                    complete_job()
                except (JobsError, JobDBError, EligibilityError) as e:
                    print(str(e))
                    sys.stdout.flush()

            # write results to properties.calc.json
            self.finalize()
            return

        elif status == "not_converging":
            print("Status:", status)
            self.report_status("failed","run_limit")
            print("Returning")
            sys.stdout.flush()
            return

        elif status == "incomplete":

            if task == "setup":
                self.setup()

            self.report_status("started")
            (status, task) = relaxation.run()

        else:
            self.report_status("failed","unknown")
            raise QEWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()


        # once the run is done, update database records accordingly

        if status == "not_converging":

            # mark error
            if self.auto:
                try:
                    error_job("Not converging")
                except (JobsError, JobDBError) as e:
                    print(str(e))
                    sys.stdout.flush()

            print("Not Converging!")
            sys.stdout.flush()
            self.report_status("failed","run_limit")

            # print a local settings file, so that the run_limit can be extended if the
            #   convergence problems are fixed
            try:
                os.makedirs(self.casm_directories.configuration_calc_settings_dir(self.casm_settings.default_clex))
            except:
                pass
            settingsfile = os.path.join(self.casm_directories.configuration_calc_settings_dir(self.casm_settings.default_clex), "relax.json")
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
                except (JobsError, JobDBError, EligibilityError) as e:
                    print(str(e))
                    sys.stdout.flush()

            # write results to properties.calc.json
            self.finalize()

        else:
            self.report_status("failed","unknown")
            raise QEWrapperError("Quantum Espresso relaxation complete with unexpected status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()

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

        outputfile = os.path.join(self.calcdir, "status.json")
        with open(outputfile, 'wb') as file:
            file.write(six.u(json.dumps(output, cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)).encode('utf-8'))
        print("Wrote " + outputfile)
        sys.stdout.flush()

    def finalize(self):
        outfilename=self.settings["outfilename"]
        if self.is_converged():
            # write properties.calc.json
            qedir = os.path.join(self.calcdir, "run.final")
            super_poscarfile = os.path.join(self.configdir,"POS")
            speciesfile = self.casm_directories.settings_path_crawl("SPECIES",self.configname,self.casm_settings.default_clex)
            output = self.properties(qedir, outfilename)
            outputfile = os.path.join(self.calcdir, "properties.calc.json")
            with open(outputfile, 'wb') as file:
                file.write(six.u(json.dumps(output, cls=noindent.NoIndentEncoder, indent=4, sort_keys=True)).encode('utf-8'))
            print("Wrote " + outputfile)
            sys.stdout.flush()
            self.report_status('complete')

    def is_converged(self):
      # Check for electronic convergence in completed calculations. Returns True or False.
      outfilename=self.settings["outfilename"]
      # Verify that the last relaxation reached electronic convergence
      relaxation = quantumespresso.Relax(self.calcdir, self.run_settings())
      for i in range(len(relaxation.rundir)):
        try:
          qrun = quantumespresso.qeio.QErun( os.path.join(self.calcdir, relaxation.rundir[-i-1], outfilename))
          if not qrun.elec_conv:
            print(('The last relaxation run (' +
                os.path.basename(relaxation.rundir[-i-1]) +
                ') failed to achieve electronic convergence; properties.calc.json will not be written.\n'))
            self.report_status('failed','electronic_convergence')
            return False
          break
        except:
          pass

      # Verify that the final static run reached electronic convergence
      qrun = quantumespresso.qeio.QErun( os.path.join(self.calcdir, "run.final", outfilename) )
      if not qrun.elec_conv:
          print('The final run failed to achieve electronic convergence; properties.calc.json will not be written.\n')
          self.report_status('failed','electronic_convergence')
          return False

      return True

    @staticmethod
    def properties(qedir, outfilename):
        """Report results to properties.calc.json file in configuration directory, after checking for electronic convergence."""

        output = dict()
        qrun = quantumespresso.qeio.QErun( os.path.join(qedir, outfilename) )


        output["atom_type"] = qrun.atom_type
        output["atoms_per_type"] = qrun.atoms_per_type
        if qrun.coord_mode == "crystal":
            output["coord_mode"] = "direct"
        else:
            output["coord_mode"] = "cartesian" + qrun.coord_mode

        output["relaxed_forces"] = [noindent.NoIndent(v) for v in list(map(lambda y: list(map(lambda x: x*13.605698066/0.52918,y)), qrun.forces ))] #convert Ry/bohr to eV/angst

        output["relaxed_lattice"] = [noindent.NoIndent(v) for v in list(map(lambda y: list(map(lambda x: x* 0.52918,y)),qrun.lattice) )] #convert bohr to angst

        output["relaxed_basis"] = [noindent.NoIndent(v) for v in qrun.basis]

        output["relaxed_energy"] = qrun.total_energy * 13.605698066 #convert Ry to eV

        return output

