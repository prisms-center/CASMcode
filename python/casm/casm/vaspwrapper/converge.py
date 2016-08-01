import os, math, sys, json, re, warnings
import pbs
import vasp
import casm
import casm.project
import vaspwrapper

class Converge(object):
    """The Converge class contains functions for setting up, executing, and parsing a VASP convergence.

        The convergence creates the following directory structure:
        config/
          calctype.name/
              converge.convtype/
                  value_0/
                  value_1/
                  ...
                  value_n/
        'value_i' directories represent different values of convtype (e.g. k-point lengths for KPOINTS,
        or energy cutoffs for ENCUT) and are populated, from lower accuracy to higher accuracy, using
        settings in .casmroot (if provided) or by the default behavior, which depends on the convtype.

        This automatically looks for VASP settings files in .../settings/calctype.name,
        where '...' is the nearest parent directory of 'self.configdir' in the CASM project repository

        Contains:
            self.configdir (.../config)
            self.calcdir   (.../config/calctype.name)

            self.convtype = What is being converged (e.g. KPOINTS, ENCUT, etc)
            self.settings = dictionary of settings for pbs and the convergence, see vaspwrapper.read_settings

            self.auto = True if using pbs module's JobDB to manage pbs jobs
            self.sort = True if sorting atoms in POSCAR by type
    """
    def __init__(self, configdir = None, convtype = None, auto = True, sort = True):
        """
        Construct a VASP convergence job object.

        Args:
            configdir: path to configuration
            convtype: type of convergence (e.g. 'kpoints' or 'encut')
            auto: True if using pbs module's JobDB to manage pbs jobs

        """
        if configdir == None:
            configdir = os.getcwd()

        print "Reading CASM settings"
        self.casm_settings = casm.project.ProjectSettings()
        if self.casm_settings == None:
            raise vaspwrapper.VaspWrapperError("Not in a CASM project. The file '.casm/project_settings.json' was not found.")

        self.casm_directories=casm.project.DirectoryStructure()

        print "Constructing a CASM VASPWrapper Converge object"
        sys.stdout.flush()

        print "  Setting up directories"
        sys.stdout.flush()

        # store path to .../config, if not existing raise
        self.configdir = os.path.abspath(configdir)
        if not os.path.isdir(self.configdir):
            raise vasp.VaspError("Error in casm.vasp.relax: Did not find directory: " + self.configdir)
            sys.stdout.flush()

        # store path to .../config/calctype.name, and create if not existing
        self.calcdir = self.casm_directories.calctype_dir(configdir,self.casm_settings.default_clex)
        try:
            os.mkdir(self.calcdir)
        except:
            pass

        # read the settings json file
        print "  Reading converge.json settings file"
        sys.stdout.flush()
        setfile = self.casm_directories.settings_path_crawl("converge.json",self.casm_settings.default_clex,self.configdir)

        if setfile == None:
            raise vaspwrapper.VaspWrapperError("Could not find \"converge.json\" in an appropriate \"settings\" directory")
            sys.stdout.flush()
        self.settings = vaspwrapper.read_settings(setfile)

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
        if not "encut" in self.settings:
            self.settings["encut"] = None
        if not "kpoints" in self.settings:
            self.settings["kpoints"] = None
        if not "nrg_convergence" in self.settings:
            self.settings["nrg_convergence"] = None
        if not "run_limit" in self.settings:
            self.settings["run_limit"] = None

        self.auto = auto
        self.sort = sort
        print "VASP Converge object constructed\n"
        sys.stdout.flush()


    def setup(self):
        """ Setup initial relaxation run

            Uses the following files from the most local .../settings/calctype.name directory:
                INCAR: VASP input settings
                KPOINTS: VASP kpoints settings
                POSCAR: reference for KPOINTS if KPOINTS mode is not A/AUTO/Automatic
                SPECIES: info for each species such as which POTCAR files to use, MAGMOM, GGA+U, etc.

            Uses the following files from the .../config directory:
                POS: structure of the configuration to be relaxed

        """
        vaspfiles=casm.vaspwrapper.vasp_input_file_names(self.casm_directories,self.casm_settings.default_clex,self.configdir)
        incarfile,prim_kpointsfile,prim_poscarfile,super_poscarfile,speciesfile=vaspfiles

        vasp.io.write_vasp_input(self.calcdir, incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, self.sort)


    def submit(self):
        """Submit a PBS job for this VASP convergence"""

        # first, check if the job has already been submitted and is not completed
        db = pbs.JobDB()
        print "rundir", self.calcdir
        id = db.select_regex_id("rundir", self.calcdir)
        print "id:", id
        sys.stdout.flush()
        if id != []:
            for j in id:
                job = db.select_job(j)
                # taskstatus = ["Incomplete","Complete","Continued","Check","Error:.*","Aborted"]
                # jobstatus = ["C","Q","R","E","W","H","M"]
                if job["jobstatus"] != "C":
                    print "JobID:", job["jobid"], "  Jobstatus:", job["jobstatus"], "  Not submitting."
                    sys.stdout.flush()
                    return
                #elif job["taskstatus"] in ["Complete", "Check"] or re.match( "Error:.*", job["taskstatus"]):
                #    print "JobID:", job["jobid"], "  Taskstatus:", job["taskstatus"], "  Not submitting."
                #    sys.stdout.flush()
                #    return


        # second, only submit a job if relaxation status is "incomplete"

        # construct the Relax object
        convergence = vasp.Converge(self.calcdir, self.run_settings(), convtype=self.convtype)

        # check the current status
        (status, task) = convergence.status()

        if status == "complete":
            print "Status:", status, "  Not submitting."
            sys.stdout.flush()
            return

        elif status == "not_converging":
            print "Status:", status, "  Not submitting."
            sys.stdout.flush()
            return

        elif status != "incomplete":
            raise vaspwrapper.VaspWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()
            return


        print "Preparing to submit a VASP convergence PBS job"
        sys.stdout.flush()

        # cd to configdir, submit jobs from configdir, then cd back to currdir
        currdir = os.getcwd()
        os.chdir(self.calcdir)

        # determine the number of atoms in the configuration
        print "  Counting atoms in the POSCAR"
        sys.stdout.flush()
        pos = vasp.io.Poscar(os.path.join(self.configdir,"POS"))
        N = len(pos.basis)

        print "  Constructing a PBS job"
        sys.stdout.flush()
        # construct a pbs.Job
        job = pbs.Job(name=casm.jobname(self.configdir),\
                      account=self.settings["account"],\
                      nodes=int(math.ceil(float(N)/float(self.settings["atom_per_proc"])/float(self.settings["ppn"]))),\
                      ppn=int(self.settings["ppn"]),\
                      walltime=self.settings["walltime"],\
                      pmem=self.settings["pmem"],\
                      queue=self.settings["queue"],\
                      message=self.settings["message"],\
                      email=self.settings["email"],\
                      priority=self.settings["priority"],\
                      command="python -c \"import casm; casm.vaspwrapper.Converge('" + self.configdir + "', '" + self.convtype + "').run()\"",\
                      auto=self.auto)

        print "  Submitting"
        sys.stdout.flush()
        # submit the job
        job.submit()

        # return to current directory
        os.chdir(currdir)

        print "CASM VASPWrapper relaxation PBS job submission complete\n"
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
        elif settings["npar"] == "VASP_DEFAULT":
            settings["npar"] = None

        if settings["npar"] == None:
            if settings["ncore"] == "CASM_DEFAULT":
                if "PBS_NUM_PPN" in os.environ:
                    settings["ncore"] = int(os.environ["PBS_NUM_PPN"])
                else:
                    settings["ncore"] = None
            elif settings["ncore"] == "VASP_DEFAULT":
                settings["ncore"] = 1
        else:
            settings["ncore"] = None

        if settings["ncpus"] == None or settings["ncpus"] == "CASM_DEFAULT":
            if "PBS_NP" in os.environ:
                settings["ncpus"] = int(os.environ["PBS_NP"])
            else:
                settings["ncpus"] = None

        if settings["run_limit"] == None or settings["run_limit"] == "CASM_DEFAULT":
            settings["run_limit"] = 10

        if settings["encut"] == None:
            settings["encut"] = ["Auto", "Auto", 10]

        if settings["kpoints"] == None:
            settings["kpoints"] = [5, "Auto", 1]

        if settings["nrg_convergence"] == None:
            settings["nrg_convergence"] = 1E-3

        return settings


    def run(self):
        """ Setup input files, run a vasp convergence, and report results """

        # construct the Relax object
        convergence = vasp.Converge(self.calcdir, self.run_settings(), convtype=self.convtype)

        # check the current status
        (status, task) = convergence.status()


        if status == "complete":
            print "Status:", status
            sys.stdout.flush()

            # mark job as complete in db
            if self.auto:
                try:
                    pbs.complete_job()
                except (pbs.PBSError, pbs.JobDBError, pbs.EligibilityError) as e:
                    print str(e)
                    sys.stdout.flush()

###<HERE>###

            # write results to output.VASP
            self.report()
            return

        elif status == "not_converging":
            print "Status:", status
            print "Returning"
            sys.stdout.flush()
            return

        elif status == "incomplete":

            if task == "setup":
                self.setup()

            (status, task) = relaxation.run()

        else:
            raise vaspwrapper.VaspWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()


        # once the run is done, update database records accordingly

        if status == "not_converging":

            # mark error
            if self.auto:
                try:
                    pbs.error_job("Not converging")
                except (pbs.PBSError, pbs.JobDBError) as e:
                    print str(e)
                    sys.stdout.flush()

            print "Not Converging!"
            sys.stdout.flush()

            # print a local settings file, so that the run_limit can be extended if the
            #   convergence problems are fixed
            try:
                os.makedirs(self.casm_directories.configuration_calc_settings_dir(self.casm_settings.default_clex))
            except:
                pass
            settingsfile = os.path.join(self.casm_directories.configuration_calc_settings_dir(self.casm_settings.default_clex), "relax.json")
            vaspwrapper.write_settings(self.settings, settingsfile)

            print "Writing:", settingsfile
            print "Edit the 'run_limit' property if you wish to continue."
            sys.stdout.flush()
            return

        elif status == "complete":

            # mark job as complete in db
            if self.auto:
                try:
                    pbs.complete_job()
                except (pbs.PBSError, pbs.JobDBError, pbs.EligibilityError) as e:
                    print str(e)
                    sys.stdout.flush()

            # write results to output.VASP
            self.report()

        else:
            raise vaspwrapper.VaspWrapperError("vasp relaxation complete with unexpected status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()



    def report(self):
        """Report final results to properties.calc.json file in configuration directory."""

        outputfile = os.path.join(self.calcdir, "properties.calc.json")

        if os.path.isfile(outputfile):
            with open(outputfile, 'r') as file:
                output = json.load(file)
        else:
            output = dict()

        vrun = vasp.io.Vasprun( os.path.join(self.calcdir, "run.final", "vasprun.xml") )

        # the calculation is run on the 'sorted' POSCAR, need to report results 'unsorted'

        if self.sort:
            super_poscarfile = os.path.join(self.configdir,"POS")
            speciesfile = casm.settings_path_crawl("SPECIES",self.casm_settings.default_clex,self.configdir)
            species_settings = vasp.io.species_settings(speciesfile)
            super = vasp.io.Poscar(super_poscarfile, species_settings)
            unsort_dict = super.unsort_dict()
        else:
            # fake unsort_dict (unsort_dict[i] == i)
            unsort_dict = dict(zip(range(0,len(vrun.basis)),range(0,len(vrun.basis))))

        print unsort_dict

        # unsort_dict:
        #   Returns 'unsort_dict', for which: unsorted_dict[orig_index] == sorted_index;
        #   unsorted_dict[sorted_index] == orig_index
        #   For example:
        #     'unsort_dict[0]' returns the index into the unsorted POSCAR of the first atom in the sorted POSCAR


        results = dict()
        results["is_complete"] = vrun.is_complete
        results["atom_type"] = super.type_atoms
        results["atoms_per_type"] = super.num_atoms
        results["coord_mode"] = vrun.coord_mode


        # as lists
        results["relaxed_forces"] = [ None for i in range(len(vrun.forces))]
        for i, v in enumerate(vrun.forces):
            results["relaxed_forces"][unsort_dict[i] ] = vrun.forces[i]

        results["relaxed_lattice"] = vrun.lattice

        results["relaxed_basis"] = [ None for i in range(len(vrun.basis))]
        for i, v in enumerate(vrun.basis):
            results["relaxed_basis"][unsort_dict[i] ] = vrun.basis[i]

        results["relaxed_energy"] = vrun.total_energy

        with open(outputfile, 'w') as file:
            json.dump(results, file, indent=4, sort_keys=True)


        print "Wrote " + outputfile
        sys.stdout.flush()


