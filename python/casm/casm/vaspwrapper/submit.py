import os, math, sys, json, re, warnings, shutil
import pbs
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper

class Submit(object):
    """Submit a PBS job for this VASP relaxation"""

    def __init__(self, selection = None, method = 'relax'):
        self.selection = selection
        for config_data in self.selection.data:
            print "Submitting..."
            print "Configuration:", config_data.configname
            #first, check if the job has already been submitted and is not completed
            db = pbs.JobDB()
            print "Calculation directory:", config_data.calcdir
            id = db.select_regex_id("rundir", config_data.calcdir)
            print "JobID:", id
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

            available_methods = ['relax':vasp.Relax ,'neb':vasp.Neb]
            # construct the Relax object
            calculation = available_methods[method](self.calcdir, self.run_settings())

            # check the current status
            (status, task) = calculation.status()
##edit from here
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
            if not os.path.isfile(os.path.join(self.calcdir, "properties.calc.json")):
                self.finalize()

            return

        elif status == "not_converging":
            print "Status:", status, "  Not submitting."
            sys.stdout.flush()
            return

        elif status != "incomplete":
            raise vaspwrapper.VaspWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()
            return


        print "Preparing to submit a VASP relaxation PBS job"
        sys.stdout.flush()

        # cd to configdir, submit jobs from configdir, then cd back to currdir
        currdir = os.getcwd()
        os.chdir(self.calcdir)

        # determine the number of atoms in the configuration
        #print "Counting atoms in the POSCAR"
        #sys.stdout.flush()
        #pos = vasp.io.Poscar(os.path.join(self.configdir,"POS"))
        #N = len(pos.basis)

        # determine the number of nodes from settings #TODO
        # still has to implement a option to use either atoms_per_proc or nodes_per_image or images_per_node
        nodes = int(self.n_images) * float(self.settings["nodes_per_image"])
        # construct command to be run
        cmd = ""
        if self.settings["preamble"] is not None:
        # Append any instructions given in the 'preamble' file, if given
            preamble = self.casm_directories.settings_path_crawl(self.settings["preamble"],
                                                                 self.configname, self.clex,
                                                                 self.calc_subdir)
            with open(preamble) as my_preamble:
                cmd += "".join(my_preamble)
        # Or just execute a single prerun line, if given
        if self.settings["prerun"] is not None:
            cmd += self.settings["prerun"] + "\n"
        cmd += "python -c \"import casm.vaspwrapper; casm.vaspwrapper.Neb('" + self.configdir + "').run()\"\n"
        if self.settings["postrun"] is not None:
            cmd += self.settings["postrun"] + "\n"

        print "Constructing a PBS job"
        sys.stdout.flush()
        # construct a pbs.Job
        job = pbs.Job(name=casm.jobname(self.configdir),\
                      account=self.settings["account"],\
                      nodes=int(nodes),
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

        print "Submitting"
        sys.stdout.flush()
        # submit the job
        job.submit()
        self.report_status("submitted")

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
