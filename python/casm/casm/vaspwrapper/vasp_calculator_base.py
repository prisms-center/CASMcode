"""implements the parent class for vasp calculations"""

import os, math, sys, json
import pandas
import numpy as np
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
import pbs

class VaspCalculatorBase(object):
    """
    Base class containing all the basic functions that method classes can inherit
    
    Attributes
    ----------
    selection : casm.project.Selection
        selection of configuration
    calctype : string
        calctype to setup and run the neb calculations
    auto : bool
    sort : bool

    """
    def __init__(self, selection, calctype=None, auto=True, sort=True):
        """set up attributes for the base class"""
        self.selection = selection
        self.calctype = calctype
        self.auto = auto
        self.sort = sort
        self.casm_directories = self.selection.proj.dir
        self.casm_settings = self.selection.proj.settings
        if self.casm_settings is None:
            raise vaspwrapper.VaspWrapperError("Not in a CASM project. The file '.casm' directory was not found.")

        self.clex = self.casm_settings.default_clex
        if calctype:
            self.clex.calctype = calctype
        self.calc_subdir = ""
        self.results_subdir = '' #everything between $(calcdir)/run.*/ and OSZICAR and OUTCAR files
        self.calculator = None
        self.append_selection_data()

    @classmethod
    def from_configuration_dir(cls, configuration_dir, calctype, auto=True, sort=True):
        """returns a instance of the Neb class instantited with a single configuration"""
        # change config_dir to configuration_dir all over
        proj = Project(configuration_dir)
        sel = Selection(proj, "EMPTY", "config", False)
        split_path = configuration_dir.split(os.path.sep)
        index = split_path.index("training_data")
        name = '/'.join(split_path[index+1:])
        sel.data = pandas.DataFrame({"name":name, "selected":1}, index=range(1))
        # try:
        #     os.mkdir(os.path.join(proj.path, ".casm/tmp"))
        # except:
        #     pass
        # sel_config = sel.saveas(os.path.join(proj.path, ".casm/tmp", configname.replace('/', '.')), True)
        obj = cls(sel, calctype, auto, sort)
        return obj

    def config_properties(self, config_data):
        """read properties directories of a specific configuration"""
        config_dict = dict(config_data)
        config_dict["configdir"] =  self.casm_directories.configuration_dir(config_data["name"], self.calc_subdir)
        config_dict["calcdir"] = self.casm_directories.calctype_dir(config_data["name"], self.clex, self.calc_subdir)
        config_dict["setfile"] = self.casm_directories.settings_path_crawl("calc.json", config_data["name"], self.clex, self.calc_subdir)
        return config_dict

    def append_selection_data(self):
        """append configproperties to selection.data"""
        config_dicts = []
        for index,config_data in self.selection.data.iterrows():
            config_dicts.append(self.config_properties(config_data))
        properites_dict = {}
        for key in config_dicts[0].keys():
            if key in self.selection.data.columns:
                continue
            properites_dict[key] = []
            for config_dict in config_dicts:
                properites_dict[key].append([config_dict[key]])
            self.selection.add_data(key, properites_dict[key])

    def pre_setup(self):
        """has to be overloaded in the child class"""
        pass

    def setup(self):
        """Setup initial relaxation run for the selection"""
        self.pre_setup()
        for index, config_data in self.selection.data.iterrows():
            self.config_setup(config_data)

    def config_setup(self, config_data):
        """ Setup initial relaxation run

            Uses the following files from the most local .../settings/calctype.name directory:
                INCAR: VASP input settings
                KPOINTS: VASP kpoints settings
                POSCAR: reference for KPOINTS if KPOINTS mode is not A/AUTO/Automatic
                SPECIES: info for each species such as which POTCAR files to use, MAGMOM, GGA+U, etc.

            Uses the following files from the .../config_calcdir/00 directory:
                POSCAR: sample structure of the configuration to be relaxed

        """
        settings = self.read_settings(config_data["setfile"])
        vaspfiles = self.get_vasp_input_files(config_data, settings)
        incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, extra_input_files = vaspfiles
        vasp.io.write_vasp_input(config_data["calcdir"], incarfile,
                                 prim_kpointsfile, prim_poscarfile,
                                 super_poscarfile, speciesfile,
                                 self.sort, extra_input_files,
                                 settings["strict_kpoints"])
        if settings["initial_deformation"] != None:
            deformation = self.get_deformation(settings)
            self.apply_deformation(deformation, config_data["calcdir"])

    def get_deformation(self, settings):
        """ Either reads or queries for deformation matrix from settings dict"""
        if settings["initial_deformation"]["method"] == 'manual':
            deformation = np.array(settings["initial_deformation"]["deformation"])
        elif settings["initial_deformation"]["method"] == 'auto':
            configname = settings["initial_deformation"]["configname"]
            calctype = settings["initial_deformation"]["calctype"]
            sel_tmp = Selection(self.selection.proj, "EMPTY", "config", False)
            sel_tmp.data = pandas.DataFrame({"configname":configname, "selected":1},
                                            index=range(1))
            try:
                os.mkdir(os.path.join(proj.path, ".casm/tmp"))
            except:
                pass
            sel_config = sel_tmp.saveas(os.path.join(self.selection.proj.path, ".casm/tmp",
                                                     configname.replace('/', '.')), True)
            sel_config.query(["relaxation_strain(U,0:5,{})".format(calctype)])
            deformation = np.array([float(sel_config.data["relaxation_strain(U,{},{})".format(i, calctype)].loc[0]) for i in range(6)])
            os.remove(os.path.join(self.selection.proj.path, ".casm/tmp",
                                   configname.replace('/', '.')))
        else:
            raise vaspwrapper.VaspWrapperError("use manual or auto mode to set initial deformation. see casm format --vasp for settings")

        if deformation.ndim == 1:
            deformation = np.array([[deformation[0], deformation[5], deformation[4]],
                                    [deformation[5], deformation[1], deformation[3]],
                                    [deformation[4], deformation[3], deformation[2]]])
        return deformation

    def apply_deformation(self, deformation, calcdir):
        """applies the deformation and write the poscar to calcdir"""
        poscarfile = os.path.join(calcdir, "POSCAR")
        poscar_obj = vasp.io.poscar.Poscar(poscarfile)
        poscar_obj.apply_deformation(deformation)
        poscar_obj.write(poscarfile)

    @staticmethod
    def read_settings(setfile):
        """ Read settings from a settings calc.json file"""

        settings = vaspwrapper.read_settings(setfile)
        # set default settings if not present
        if not "ncore" in settings:
            settings["ncore"] = None
        if not "npar" in settings:
            settings["npar"] = None
        if not "kpar" in settings:
            settings["kpar"] = None
        if not "vasp_cmd" in settings:
            settings["vasp_cmd"] = None
        if not "ncpus" in settings:
            settings["ncpus"] = None
        if not "run_limit" in settings:
            settings["run_limit"] = None
        if not "prerun" in settings:
            settings["prerun"] = None
        if not "postrun" in settings:
            settings["postrun"] = None

        return settings

    def get_vasp_input_files(self, config_data, settings):
        # Find required input files in CASM project directory tree
        vaspfiles = casm.vaspwrapper.vasp_input_file_names(self.casm_directories,
                                                           config_data["name"],
                                                           self.clex,
                                                           self.calc_subdir)
        incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile = vaspfiles
        # Find optional input files
        extra_input_files = []
        for s in settings["extra_input_files"]:
            extra_input_files.append(self.casm_directories.settings_path_crawl(s, config_data["name"],
                                                                               self.clex, self.calc_subdir))
            if extra_input_files[-1] is None:
                raise vasp.VaspError("Neb.setup failed. Extra input file " + s + " not found in CASM project.")
        if settings["initial"]:
            extra_input_files += [self.casm_directories.settings_path_crawl(settings["initial"],
                                                                            config_data["name"],
                                                                            self.clex, self.calc_subdir)]
            if extra_input_files[-1] is None:
                raise vasp.VaspError("Neb.setup failed. No initial INCAR file " + settings["initial"] + " found in CASM project.")
        if settings["final"]:
            extra_input_files += [self.casm_directories.settings_path_crawl(settings["final"],
                                                                            config_data["name"],
                                                                            self.clex, self.calc_subdir)]
            if extra_input_files[-1] is None:
                raise vasp.VaspError("Neb.setup failed. No final INCAR file " + settings["final"] + " found in CASM project.")
        return incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, extra_input_files

    def submit(self):
        """ submit jobs for a selection"""
        self.pre_setup()
        db = pbs.JobDB()
        for index,config_data in self.selection.data.iterrows():
            print "Submitting..."
            print "Configuration:", config_data["name"]
            #first, check if the job has already been submitted and is not completed
            print "Calculation directory:", config_data["calcdir"]
            id = db.select_regex_id("rundir", config_data["calcdir"])
            print "JobID:", id
            sys.stdout.flush()
            try:
                if id != []:
                    db.update()
                    for j in id:
                        job = db.select_job(j)
                        if job["jobstatus"] != "C":
                            print "JobID:", job["jobid"], "  Jobstatus:", job["jobstatus"], "  Not submitting."
                            sys.stdout.flush()
                            raise BreakException
            except BreakException:
                continue
            settings = self.read_settings(config_data["setfile"])
            # construct the Relax object
            calculation = self.calculator(config_data["calcdir"], self.run_settings(settings))
            # check the current status
            (status, task) = calculation.status()

            if status == "complete":
                print "Status:", status, "  Not submitting."
                sys.stdout.flush()

                # ensure job marked as complete in db
                if self.auto:
                    for j in id:
                        job = db.select_job(j)
                        if job["taskstatus"] == "Incomplete":
                            try:
                                pbs.complete_job(jobid=j)
                            except (pbs.PBSError, pbs.JobDBError, pbs.EligibilityError) as e:
                                print str(e)
                                sys.stdout.flush()

                # ensure results report written
                if not os.path.isfile(os.path.join(config_data["calcdir"], "properties.calc.json")):
                    self.finalize(config_data)

                continue

            elif status == "not_converging":
                print "Status:", status, "  Not submitting."
                sys.stdout.flush()
                continue

            elif status != "incomplete":
                raise vaspwrapper.VaspWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")
                sys.stdout.flush()
                continue

            print "Preparing to submit a VASP relaxation PBS job"
            sys.stdout.flush()

            # cd to configdir, submit jobs from configdir, then cd back to currdir
            currdir = os.getcwd()
            os.chdir(config_data["calcdir"])

            nodes, ppn = self._calc_submit_node_info(settings, config_data)

            # construct command to be run
            cmd = ""
            if settings["preamble"] is not None:
                # Append any instructions given in the 'preamble' file, if given
                preamble = self.casm_directories.settings_path_crawl(settings["preamble"],
                                                                     config_data["name"],
                                                                     self.clex,
                                                                     self.calc_subdir)
                with open(preamble) as my_preamble:
                    cmd += "".join(my_preamble)
            # Or just execute a single prerun line, if given
            if settings["prerun"] is not None:
                cmd += settings["prerun"] + "\n"
            cmd += self.run_cmd(config_data["configdir"], self.calctype)
            if settings["postrun"] is not None:
                cmd += settings["postrun"] + "\n"

            print "Constructing a PBS job"
            sys.stdout.flush()
            # construct a pbs.Job
            job = pbs.Job(name=casm.jobname(config_data["configdir"]),\
                          account=settings["account"],\
                          nodes=nodes, ppn=ppn,\
                          walltime=settings["walltime"],\
                          pmem=settings["pmem"],\
                          qos=settings["qos"],\
                          queue=settings["queue"],\
                          message=settings["message"],\
                          email=settings["email"],\
                          priority=settings["priority"],\
                          command=cmd,\
                          auto=self.auto,
			  software=db.config["software"])

            print "Submitting"
            sys.stdout.flush()
            # submit the job
            job.submit()
            self.report_status(config_data["calcdir"], "submitted")

            # return to current directory
            os.chdir(currdir)

            print "CASM VASPWrapper relaxation PBS job submission complete\n"
            sys.stdout.flush()

    def run_cmd(self, configdir, calctype):
        """has to be overloaded in the method class"""
        return None

    @staticmethod
    def _calc_submit_node_info(settings, config_data):
        """return nodes, ppn from settings of a configuration"""
        if settings["nodes"] != None and settings["ppn"] != None:
            return int(settings["nodes"]), int(settings["ppn"])
        elif settings["atom_per_proc"] != None and settings["ppn"] != None:
            pos = vasp.io.Poscar(os.path.join(config_data["calcdir"], "POSCAR"))
            num = len(pos.basis)
            nodes = int(math.ceil(float(num)/float(settings["atom_per_proc"])/float(settings["ppn"])))
            return nodes, int(settings["ppn"])
        elif settings["nodes_per_image"] != None and settings["ppn"] != None:
            nodes = int(config_data["n_images"]) * float(settings["nodes_per_image"])
            return nodes, int(settings["ppn"])
        else:
            raise vaspwrapper.VaspWrapperError("Not enough information to determine nodes and ppn information")

    @staticmethod
    def run_settings(settings):
        """ Set default values based on runtime environment"""

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
        """run the job of a selection"""
        self.pre_setup()
        for index,config_data in self.selection.data.iterrows():
            settings = self.read_settings(config_data["setfile"])
            calculation = self.calculator(config_data["calcdir"], self.run_settings(settings))

            # check the current status
            (status, task) = calculation.status()

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

                # write results to properties.calc.json
                self.finalize(config_data)
                continue

            elif status == "not_converging":
                print "Status:", status
                self.report_status(config_data["calcdir"], "failed", "run_limit")
                print "Returning"
                sys.stdout.flush()
                continue

            elif status == "incomplete":

                if task == "setup":
                    self.config_setup(config_data)

                self.report_status(config_data["calcdir"], "started")
                (status, task) = calculation.run()

            else:
                self.report_status(config_data["calcdir"], "failed", "unknown")
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
                self.report_status(config_data["calcdir"], "failed", "run_limit")

                # print a local settings file, so that the run_limit can be extended if the
                #   convergence problems are fixed

                config_set_dir = self.casm_directories.configuration_calc_settings_dir(config_data["name"],
                                                                                       self.clex,
                                                                                       self.calc_subdir)

                try:
                    os.makedirs(config_set_dir)
                except:
                    pass
                settingsfile = os.path.join(config_set_dir, "calc.json")
                vaspwrapper.write_settings(settings, settingsfile)

                print "Writing:", settingsfile
                print "Edit the 'run_limit' property if you wish to continue."
                sys.stdout.flush()
                continue

            elif status == "complete":

                # mark job as complete in db
                if self.auto:
                    try:
                        pbs.complete_job()
                    except (pbs.PBSError, pbs.JobDBError, pbs.EligibilityError) as e:
                        print str(e)
                        sys.stdout.flush()

                # write results to properties.calc.json
                self.finalize(config_data)

            else:
                self.report_status(config_data["calcdir"], "failed", "unknown")
                raise vaspwrapper.VaspWrapperError("vasp relaxation complete with unexpected status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()


    def report_status(self, calcdir, status, failure_type=None):
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

        outputfile = os.path.join(calcdir, "status.json")
        with open(outputfile, 'w') as file:
            file.write(json.dumps(output, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))
        print "Wrote " + outputfile
        sys.stdout.flush()

    def report(self):
        """
        report status for the selection
        Notes
        -----
        checks for convergence
        calls the finalize function to write the approprite properties files.
        """
        for index,config_data in self.selection.data.iterrows():
            try:
                settings = self.read_settings(config_data["setfile"])
                calculation = self.calculator(config_data["calcdir"], self.run_settings(settings))
                if self.is_converged(calculation):
                    self.finalize(config_data)
            except:
                print("Unable to report properties for directory {}.\n"
                      "Please verify that it contains a completed VASP calculation.".format(config_data["configdir"]))
                raise

    def finalize(self, config_data, super_poscarfile=None):
        # write properties.calc.json
        vaspdir = os.path.join(config_data["calcdir"], "run.final")
        speciesfile = self.casm_directories.settings_path_crawl("SPECIES",
                                                                config_data["name"],
                                                                self.clex,
                                                                self.calc_subdir)
        output = self.properties(vaspdir, super_poscarfile, speciesfile)
        outputfile = os.path.join(config_data["calcdir"], "properties.calc.json")
        with open(outputfile, 'w') as file:
            file.write(json.dumps(output, file, cls=casm.NoIndentEncoder,
                                  indent=4, sort_keys=True))
        print "Wrote " + outputfile
        sys.stdout.flush()
        self.report_status(config_data["calcdir"], 'complete')

    def is_converged(self, calculation):
        # Check for electronic convergence in completed calculations. Returns True or False.
        # Verify that the last relaxation reached electronic convergence
        for i in range(len(calculation.rundir)):
            try:
                vrun_oszicar = vasp.io.Oszicar(os.path.join(calculation.calcdir, calculation.rundir[-i-1],
                                                            self.results_subdir, "OSZICAR"))
                vrun_nelm = vasp.io.get_incar_tag("NELM", os.path.join(calculation.calcdir, calculation.rundir[-i-1]))
                if vrun_nelm == None: vrun_nelm = 60 ##pushing the default. may be write an addon to get it from outcar
                if len(vrun_oszicar.num_elm[-1]) >= vrun_nelm:
                    print('The last relaxation run (' +
                          os.path.basename(relaxation.rundir[-i-1]) +
                          ') failed to achieve electronic convergence; properties.calc.json will not be written.\n')
                    self.report_status(calculation.calcdir, 'failed', 'electronic_convergence')
                    return False
                break
            except:
                pass

        # Verify that the final static run reached electronic convergence
        vrun_oszicar = vasp.io.Oszicar(os.path.join(calculation.calcdir, "run.final",
                                                    self.results_subdir, "OSZICAR"))
        vrun_nelm = vasp.io.get_incar_tag("NELM", os.path.join(calculation.calcdir, "run.final"))
        if vrun_nelm == None: vrun_nelm = 60 ##pushing the default. may be write an addon to get it from outcar
        if vrun_oszicar.num_elm[-1] >= vrun_nelm:
            print('The final run failed to achieve electronic convergence; properties.calc.json will not be written.\n')
            self.report_status(calculation.calcdir, 'failed', 'electronic_convergence')
            return False

        return True

    @staticmethod
    def properties(vaspdir, super_poscarfile=None, speciesfile=None):
        """ return a dict of output form a vasp directory"""

        output = dict()
        # load the OSZICAR and OUTCAR
        zcar = vasp.io.Oszicar(os.path.join(vaspdir, "OSZICAR"))
        ocar = vasp.io.Outcar(os.path.join(vaspdir, "OUTCAR"))

        # the calculation is run on the 'sorted' POSCAR, need to report results 'unsorted'

        if (super_poscarfile is not None) and (speciesfile is not None):
            species_settings = vasp.io.species_settings(speciesfile)
            super_poscar = vasp.io.Poscar(super_poscarfile, species_settings)
            unsort_dict = super_poscar.unsort_dict()
        else:
            # fake unsort_dict (unsort_dict[i] == i)
            super_poscar = vasp.io.Poscar(os.path.join(vaspdir, "POSCAR"))
            unsort_dict = dict(zip(range(0, len(super_poscar.basis)),
                                   range(0, len(super_poscar.basis))))
        super_contcar = vasp.io.Poscar(os.path.join(vaspdir, "CONTCAR"))

        # unsort_dict:
        #   Returns 'unsort_dict', for which: unsorted_dict[orig_index] == sorted_index;
        #   unsorted_dict[sorted_index] == orig_index
        #   For example:
        #     'unsort_dict[0]' returns the index into the unsorted POSCAR of the first atom in the sorted POSCAR


        output["atom_type"] = super_poscar.type_atoms
        output["atoms_per_type"] = super_poscar.num_atoms
        output["coord_mode"] = super_poscar.coord_mode

        # as lists
        output["relaxed_forces"] = [None for i in range(len(ocar.forces))]
        for i, force in enumerate(ocar.forces):
            output["relaxed_forces"][unsort_dict[i]] = casm.NoIndent(force)

        output["relaxed_lattice"] = [casm.NoIndent(list(v)) for v in super_contcar.lattice()]
        output["relaxed_basis"] = [None for i in range(len(super_contcar.basis))]
        for i, ba in enumerate(super_contcar.basis):
            output["relaxed_basis"][unsort_dict[i]] = casm.NoIndent(list(ba.position))

        output["relaxed_energy"] = zcar.E[-1]

        if ocar.ispin == 2:
            output["relaxed_magmom"] = zcar.mag[-1]
            if ocar.lorbit in [1, 2, 11, 12]:
                output["relaxed_mag_basis"] = [None for i in range(len(super_contcar.basis))]
                for i, v in enumerate(super_contcar.basis):
                    output["relaxed_mag_basis"][unsort_dict[i]] = casm.NoIndent(ocar.mag[i])

        return output

class BreakException(Exception):
    """use this exception to break an outer loop"""
    pass
