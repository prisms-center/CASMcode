import os, math, sys
import pbs
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
import vasp.Neb,vasp.Relax

class Submit(object):
    """Submit a PBS job for this VASP relaxation"""

    def __init__(self, selection=None, method='relax', auto=True):
        """
        Construct a VASP submit job object.

        Arguments
        ----------

            selection: casm.project.Selection object, default= yet to be implemented #todo
              Selection of all DiffTransConfigurations to submit calculations.
              default should be MASTER selection and yet to be implemented

            method: name of the method to perform the calculation
              Options include 'relax' and 'Neb' etc.

            auto: boolean, optional, default=True,
              Use True to use the pbs module's JobDB to manage pbs jobs

        """

        self.selection = selection
        self.auto = auto
        available_methods = {'relax': vasp.Relax, 'neb': vasp.Neb}
        from available_methods[method] import ConfigProperties #TODO
        db = pbs.JobDB()
        db.update()
        for config_data in self.selection.data:
            config_obj = ConfigProperties(config_data)
            print "Submitting..."
            print "Configuration:", config_obj.configname
            #first, check if the job has already been submitted and is not completed
            print "Calculation directory:", config_obj.calcdir
            id = db.select_regex_id("rundir", config_obj.calcdir)
            print "JobID:", id
            sys.stdout.flush()
            if id != []:
                for j in id:
                    job = db.select_job(j)
                    if job["jobstatus"] != "C":
                        print "JobID:", job["jobid"], "  Jobstatus:", job["jobstatus"], "  Not submitting."
                        sys.stdout.flush()
                        continue

            # construct the Relax object
            calculation = available_methods[method](config_obj.calcdir, run_settings(config_obj.settings))
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
                if not os.path.isfile(os.path.join(config_obj.calcdir, "properties.calc.json")):
                    report.finalize() #TODO

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
            os.chdir(config_obj.calcdir)

            nodes, ppn = self._calc_submit_node_info(config_obj)

            # construct command to be run
            cmd = ""
            if config_obj.settings["preamble"] is not None:
                # Append any instructions given in the 'preamble' file, if given
                preamble = config_obj.casm_directories.settings_path_crawl(config_obj.settings["preamble"],
                                                                           config_obj.configname,
                                                                           config_obj.clex,
                                                                           config_obj.calc_subdir)
                with open(preamble) as my_preamble:
                    cmd += "".join(my_preamble)
            # Or just execute a single prerun line, if given
            if config_obj.settings["prerun"] is not None:
                cmd += config_obj.settings["prerun"] + "\n"
            cmd += "python -c \"import casm.vaspwrapper; casm.vaspwrapper.Run('" + config_obj.configdir + "')\"\n"
            if config_obj.settings["postrun"] is not None:
                cmd += config_obj.settings["postrun"] + "\n"

            print "Constructing a PBS job"
            sys.stdout.flush()
            # construct a pbs.Job
            job = pbs.Job(name=casm.jobname(config_obj.configdir),\
                          account=config_obj.settings["account"],\
                          nodes=nodes, ppn=ppn,\
                          walltime=config_obj.settings["walltime"],\
                          pmem=config_obj.settings["pmem"],\
                          qos=config_obj.settings["qos"],\
                          queue=config_obj.settings["queue"],\
                          message=config_obj.settings["message"],\
                          email=config_obj.settings["email"],\
                          priority=config_obj.settings["priority"],\
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


    def _calc_submit_node_info(self, config_obj):
        settings = config_obj.settings
        if "nodes" in settings and "ppn" in settings:
            return int(settings["nodes"]), int(settings["ppn"])
        elif "atoms_per_proc" in settings and "ppn" in settings:
            pos = vasp.io.Poscar(os.path.join(settings.configdir, "POS"))
            N = len(pos.basis)
            nodes = int(math.ceil(float(N)/float(settings["atom_per_proc"])/float(settings["ppn"])))
            return nodes, int(settings["ppn"])
        elif "nodes_per_image" in settings and "ppn" in settings:
            nodes = int(config_obj.n_images) * float(settings["nodes_per_image"])
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
