import os, math, sys, json, re, warnings, shutil
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
import pbs
from casm.vaspwrapper import Submit, Run, Report

class VaspCalculatorBase(Submit, Run, Report):
    """ Base class containing all the basic functions that methods can inherit
    """
    def __init__(self, selection, calctype=None, auto=True, sort=True):
        """ set up attributes for the base class
        """
        self.selection = selection
        self.calctype = calctype
        self.auto = auto
        self.sort = sort
        self.casm_directories = self.selection.proj.dir
        self.casm_settings = self.selection.settings
        if self.casm_settings is None:
            raise vaspwrapper.VaspWrapperError("Not in a CASM project. The file '.casm' directory was not found.")

        self.clex = self.casm_settings.default_clex
        if calctype:
            self.clex.calctype = calctype
        self.calc_subdir = ""
        self.results_subdir = '' #everything between $(calcdir)/run.*/ and OSZICAR and OUTCAR files

    def config_properties(self, config_data):
        """read properties directories of a specific configuration"""
        config_dict = dict(config_data)
        config_dict["configdir"] =  self.casm_directories.configuration_dir(configname, self.calc_subdir)
        config_dict["calcdir"] = self.casm_directories.calctype_dir(configname, self.clex, self.calc_subdir)
        config_dict["setfile"] = self.casm_directories.settings_path_crawl("calc.json", configname, self.clex, self.calc_subdir)

        return config_dict

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

    def submit(self, calculator):
        """ submit a job
        """

        db = pbs.JobDB()
        db.update()
        for config_data in self.selection.data:
            print "Submitting..."
            print "Configuration:", config_data["configname"]
            #first, check if the job has already been submitted and is not completed
            print "Calculation directory:", config_data["calcdir"]
            id = db.select_regex_id("rundir", config_data["calcdir"])
            print "JobID:", id
            sys.stdout.flush()
            if id != []:
                for j in id:
                    job = db.select_job(j)
                    if job["jobstatus"] != "C":
                        print "JobID:", job["jobid"], "  Jobstatus:", job["jobstatus"], "  Not submitting."
                        sys.stdout.flush()
                        continue
            settings = self.read_settings(config_data["setfile"])
            # construct the Relax object
            calculation = calculator(config_data["calcdir"], run_settings(settings))
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
                    self.finalize() #TODO

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
                                                                     config_data["configname"],
                                                                     self.clex,
                                                                     self.calc_subdir)
                with open(preamble) as my_preamble:
                    cmd += "".join(my_preamble)
            # Or just execute a single prerun line, if given
            if settings["prerun"] is not None:
                cmd += settings["prerun"] + "\n"
            cmd += "python -c \"import casm.vaspwrapper; casm.vaspwrapper.Run('" + config_obj.configdir + "')\"\n" #TODO
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
                          auto=self.auto)

            print "Submitting"
            sys.stdout.flush()
            # submit the job
            job.submit()
            self.report_status("submitted")

            # return to current directory
            os.chdir(currdir)

            print "CASM VASPWrapper relaxation PBS job submission complete\n"
            sys.stdout.flush()

    @staticmethod
    def _calc_submit_node_info(config_data, settings):
        if "nodes" in settings and "ppn" in settings:
            return int(settings["nodes"]), int(settings["ppn"])
        elif "atoms_per_proc" in settings and "ppn" in settings:
            pos = vasp.io.Poscar(os.path.join(config_data["calcdir"], "POSCAR"))
            N = len(pos.basis)
            nodes = int(math.ceil(float(N)/float(settings["atom_per_proc"])/float(settings["ppn"])))
            return nodes, int(settings["ppn"])
        elif "nodes_per_image" in settings and "ppn" in settings:
            nodes = int(config_data["n_images"]) * float(settings["nodes_per_image"])
            return nodes, int(settings["ppn"])
        else:
            raise vaspwrapper.VaspWrapperError("Not enough information to determine nodes and ppn information")


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
