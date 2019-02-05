import os
import sys
import json

from prisms_jobs import Job, JobDB, error_job, complete_job, JobsError, JobDBError, EligibilityError

from casm.wrapper.misc import confname_as_jobname

from casm import wrapper
from casm.project import DirectoryStructure, ProjectSettings
from casm.misc.noindent import NoIndent, NoIndentEncoder

from casm.aims.relax import AimsRelax
from casm.aims.io.io import write_aims_input
from casm.aimswrapper import AimsWrapperError, read_settings, write_settings, aims_input_file_names
from casm.aims.io.aimsrun import AimsRun
from casm.aims.io.geometry import Geometry


class Relax(object):
    """The Relax class contains functions for setting up, executing, and parsing a relaxation.

        The relaxation creates the following directory structure:
        config/
          calctype.name/
              run.0/
              run.1/
              ...
              run.final/

        'run.i' directories are only created when ready.
        'run.final' is a final run

        This automatically looks for the settings files using casm.project.DirectoryStructure.settings_path_crawl

    Attributes
    ----------
      
      casm_settings: casm.project.ProjectSettings instance
        CASM project settings
      
      casm_directories: casm.project.DirectoryStructure instance
        CASM project directory hierarchy
      
      settings: dict
        Settings for pbs and the relaxation, see casm.project.io.read_project_settings
      
      configdir: str
        Directory where configuration results are stored. The result of:
          casm.project.DirectoryStructure.configuration_dir(self.configname)
      
      configname: str
        The name of the configuration to be calculated
      
      sort: boolean
        True if sorting atoms in POSCAR by type
      
      clex: casm.project.ClexDescription instance
        The cluster expansion being worked on. Used for the 'calctype' settings.
        Currently, fixed to self.casm_settings.default_clex.
    
    """
    def __init__(self, configdir=None, auto=True, sort=True):
        """
        Construct a relaxation job object.

        Arguments
        ----------
    
            configdir: str, optional, default=None
              Path to configuration directory. If None, uses the current working directory
            
            sort: boolean, optional, default=True,
              Use True to sort atoms in POSCAR by type

        """
        print("Construct a casm.aimswrapper.Relax instance:")
        
        if configdir is None:
            configdir = os.getcwd()
        print("  Input directory:" + configdir)
        
        # get the configname from the configdir path
        _res = os.path.split(configdir)
        self.configname = os.path.split(_res[0])[1] + "/" + _res[1]
        print("  Configuration:" + self.configname)
        
        print("  Reading CASM settings")
        self.casm_directories = DirectoryStructure(configdir)
        self.casm_settings = ProjectSettings(configdir)
        if self.casm_settings is None:
            raise AimsWrapperError("Not in a CASM project. The '.casm' directory was not found.")
        
        if os.path.abspath(configdir) != self.configdir:
            print("")
            print("input configdir:" + configdir)
            print("determined configname:" + self.configname)
            print("expected configdir given configname:" + self.configdir)
            raise AimsWrapperError("Mismatch between configname and configdir")
        
        # fixed to default_clex for now
        self.clex = self.casm_settings.default_clex
        
        # store path to .../config/calctype.name, and create if not existing
        self.calcdir = self.casm_directories.calctype_dir(self.configname, self.clex)

        print("  Calculations directory:" + self.calcdir)
        if not os.path.isdir(self.calcdir):
            os.mkdir(self.calcdir)

        # read the settings json file
        print("  Reading relax.json...")
        sys.stdout.flush()
        setfile = self.casm_directories.settings_path_crawl("relax.json", self.configname, self.clex)

        if setfile is None:
            raise AimsWrapperError("Could not find relax.json in settings directory")
        else:
            print("  Read settings from:" + setfile)
        self.settings = read_settings(setfile)

        # set default settings if not present
        if "aims_cmd" not in self.settings:
            self.settings["aims_cmd"] = None
        if "ppn" not in self.settings:
            self.settings["ppn"] = None
        if "nodes" not in self.settings:
            self.settings["nodes"] = None
        if "run_limit" not in self.settings:
            self.settings["run_limit"] = None
        if "prerun" not in self.settings:
            self.settings["prerun"] = None
        if "postrun" not in self.settings:
            self.settings["postrun"] = None

        self.auto = auto
        self.sort = sort
        self.strict_kpoints = self.settings['strict_kpoints']

        print("  DONE\n")
        sys.stdout.flush()

    @property
    def configdir(self):
        return self.casm_directories.configuration_dir(self.configname)

    def setup(self):
        """ Setup initial relaxation run
            Uses the following files from the most local .../settings/calctype.name directory:
                control.skel: input settings
                geometry.skel: reference for KPOINTS object
                basis: species info
            Uses the following files from the .../config directory:
                POS: structure of the configuration to be relaxed
        """
        # Find required input files in CASM project directory tree
        aimsfiles = aims_input_file_names(self.casm_directories, self.configname, self.clex)
        controlfile, prim_posfile, super_posfile, basisfile = aimsfiles

        write_aims_input(self.calcdir, controlfile, prim_posfile, super_posfile, basisfile, self.strict_kpoints)

    def submit(self):
        """Submit a job for this relaxation"""
        print("Submitting configuration: " + self.configname)
        # first, check if the job has already been submitted and is not completed
        db = JobDB()
        print("Calculation directory: ", self.calcdir)
        sub_id = db.select_regex_id("rundir", self.calcdir)

        if sub_id is not []:
            for j in sub_id:
                job = db.select_job(j)
                if job["jobstatus"] != "?":
                    print("JobID: " + job["jobid"],
                          "  Jobstatus:" + job["jobstatus"] + "  Not submitting.")
                    sys.stdout.flush()
                    return

        # second, only submit a job if relaxation status is "incomplete"
        # construct the Relax object
        relaxation = AimsRelax(self.calcdir, self.run_settings())
        # check the current status
        status, task = relaxation.status()

        if status == "complete":
            print("Status:", status, "  Not submitting.")
            sys.stdout.flush()

            # ensure job marked as complete in db
            if self.auto:
                for j in sub_id:
                    job = db.select_job(j)
                    if job["taskstatus"] == "Incomplete":
                        try:
                            complete_job(jobid=j)
                        except IOError as e:
                            raise AimsWrapperError(e)

            # ensure results report written
            if not os.path.isfile(os.path.join(self.calcdir, "properties.calc.json")):
                self.finalize()

            return

        elif status == "not_converging":
            print("Status:" + status + "  Not submitting.")
            sys.stdout.flush()
            return

        elif status != "incomplete":
            raise AimsWrapperError("unexpected relaxation status: '"
                                   + status + "' and task: '" + task + "'")

        print("Preparing to submit a FHI-aims relaxation job")
        sys.stdout.flush()

        # cd to configdir, submit jobs from configdir, then cd back to currdir
        currdir = os.getcwd()
        os.chdir(self.calcdir)

        # determine the number of atoms in the configuration
        print("Counting atoms in the Supercell")
        sys.stdout.flush()
#        pos = vasp.io.Poscar(os.path.join(self.configdir, "POS"))
#        N = len(pos.basis)
        
        # construct command to be run
        cmd = ""
        if self.settings["prerun"] is not None:
            cmd += self.settings["prerun"] + "\n"
        cmd += "python -c \"import casm.aimswrapper; casm.aimswrapper.Relax('" + self.configdir + "').run()\"\n"
        if self.settings["postrun"] is not None:
            cmd += self.settings["postrun"] + "\n"

        ncpus = int(self.settings['nodes']) * int(self.settings['ppn'])

        print("Constructing the job")
        sys.stdout.flush()

        # construct a pbs.Job
        job = Job(name=confname_as_jobname(self.configdir),
                  account=self.settings["account"],
                  nodes=int(self.settings["nodes"]),
                  ppn=int(ncpus),
                  walltime=self.settings["walltime"],
                  pmem=self.settings["pmem"],
                  qos=self.settings["qos"],
                  queue=self.settings["queue"],
                  message=self.settings["message"],
                  email=self.settings["email"],
                  priority=self.settings["priority"],
                  command=cmd,
                  auto=self.auto)

        print("Submitting")
        sys.stdout.flush()
        # submit the job
        job.submit()
        self.report_status("submitted")

        # return to current directory
        os.chdir(currdir)

        print("CASM AimsWrapper relaxation job submission complete\n")
        sys.stdout.flush()

    def run_settings(self):
        """ Set default values based on runtime environment"""
        settings = dict(self.settings)

        # set default values
        if settings["ncpus"] is None:
                settings["ncpus"] = int(settings["nodes"]) * int(settings["ppn"])

        if settings["run_limit"] is None or settings["run_limit"] == "CASM_DEFAULT":
            settings["run_limit"] = 10

        return settings

    def run(self):
        """ Setup input files, run a FHI-aims relaxation, and report results """
        # construct the Relax object
        relaxation = AimsRelax(self.calcdir, self.run_settings())

        # check the current status
        status, task = relaxation.status()

        if status == "complete":
            print("Status: " + status)
            sys.stdout.flush()
            # mark job as complete in db
            if self.auto:
                try:
                    complete_job()
                except IOError as e:
                    raise AimsWrapperError(e)

            # write results to properties.calc.json
            self.finalize()
            return

        elif status == "not_converging":
            print("Status:" + status)
            self.report_status("failed", "run_limit")
            print("Returning")
            sys.stdout.flush()
            return

        elif status == "incomplete":
            if task == "setup":
                self.setup()
            self.report_status("started")
            status, task = relaxation.run()

        else:
            self.report_status("failed", "unknown")
            raise AimsWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")

        # once the run is done, update database records accordingly
        if status == "not_converging":
            # mark error
            if self.auto:
                try:
                    error_job("Not converging")
                except (JobsError, JobDBError) as e:
                    raise AimsWrapperError(e)

            print("Not Converging!")
            sys.stdout.flush()
            self.report_status("failed", "run_limit")

            # print a local settings file, so that the run_limit can be extended if the
            # convergence problems are fixed
            
            config_set_dir = self.casm_directories.configuration_calc_settings_dir(self.configname, self.clex)
            
            try:
                os.makedirs(config_set_dir)
            except IOError:
                pass
            settingsfile = os.path.join(config_set_dir, "relax.json")
            write_settings(self.settings, settingsfile)

            print("Writing:" + settingsfile)
            print("Edit the 'run_limit' property if you wish to continue.")
            sys.stdout.flush()
            return

        elif status == "complete":
            # mark job as complete in db
            try:
                complete_job()
            except (JobsError, JobDBError, EligibilityError) as e:
                raise AimsWrapperError(e)
            self.finalize()

        else:
            self.report_status("failed", "unknown")
            raise AimsWrapperError("Relaxation complete with unexpected status: '" +
                                   status + "' and task: '" + task + "'")

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
#        output["basis_type"] = basis_type
        if failure_type is not None:
            output["failure_type"] = failure_type

        outputfile = os.path.join(self.calcdir, "status.json")
        with open(outputfile, 'w') as file:
            file.write(json.dumps(output, cls=NoIndentEncoder, indent=4, sort_keys=True))
        print("Wrote " + outputfile)
        sys.stdout.flush()

    def finalize(self):
        if self.is_converged():
            # write properties.calc.json
            finaldir = "run." + str(self.settings["basis"])
            aimsdir = os.path.join(self.calcdir, finaldir)
            super_posfile = os.path.join(self.configdir, "POS")
            speciesfile = self.casm_directories.settings_path_crawl("basis", self.configname, self.clex)
            output = self.properties(aimsdir, super_posfile, speciesfile)
            outputfile = os.path.join(self.calcdir, "properties.calc.json")
            with open(outputfile, 'w') as file:
                file.write(json.dumps(output, cls=NoIndentEncoder, indent=4, sort_keys=True))
            print("Wrote " + outputfile)
            sys.stdout.flush()
            self.report_status(status='complete')

    def is_converged(self):
        # Check for electronic convergence in completed calculations. Returns True or False.
        relaxation = AimsRelax(self.calcdir, self.run_settings())
        rundirname = 'run.' + str(self.settings["basis"])
        dname = os.path.join(self.configdir, 'calctype.default', rundirname)
        relaxation.rundir.append(dname)
        for i in range(len(relaxation.rundir)):
            try:
                with open(os.path.join(self.calcdir, relaxation.rundir[-i-1], "std.out")) as runfile:
                    for line in runfile:
                        if "Have a nice day." in line:
                            return True
            except IOError:
                pass

        return False

    @staticmethod
    def properties(aimsdir, super_posfile=None, basisfile=None):
        """Report results to properties.calc.json file in configuration directory,
           after checking for electronic convergence."""
        output = dict()
        arun = AimsRun(os.path.join(aimsdir, "std.out"))

        if super_posfile is not None and basisfile is not None:
            # basis_settings_loc = basis_settings(basisfile)
            super_cell = Geometry(super_posfile)
            unsort_dict = super_cell.unsort_dict()
        else:
            # fake unsort_dict (unsort_dict[i] == i)
            unsort_dict = dict(zip(range(0, len(arun.basis)), range(0, len(arun.basis))))
            super_cell = Geometry(os.path.join(aimsdir, "geometry.in"))

        output["atom_type"] = super_cell.type_atoms
        output["atoms_per_type"] = super_cell.num_atoms
        output["coord_mode"] = arun.coord_mode

        # as lists
        output["relaxed_forces"] = [None in range(len(arun.forces))]
        for i, v in enumerate(arun.forces):
            output["relaxed_forces"][unsort_dict[i]] = NoIndent(arun.forces[i])

        output["relaxed_lattice"] = [NoIndent(v) for v in arun.lattice]
        output["relaxed_basis"] = [None in range(len(arun.basis))]

        for i, v in enumerate(arun.basis):
            output["relaxed_basis"][unsort_dict[i]] = NoIndent(arun.basis[i])

        output["relaxed_energy"] = arun.total_energy
        return output
