import os, math, sys, json, re, warnings, shutil
import pbs
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
from casm.vaspwrapper.submit import run_settings

class Run(object):
    """ Setup input files, run a vasp relaxation, and report results 
    """
    def __init__(self, configname, method='relax', calctype=None, auto=True):

        available_methods = {'relax': vasp.Relax, 'neb': vasp.Neb}
        from available_methods[method] import Container ##TODO
        self.config_obj, self.finalize_obj , self.properties_obj = Container(configname, calctype)
        self.calculation = available_methods[method](self.config_obj.calcdir, run_settings(self.config_obj.settings))

    def run(self):
        """run the job"""
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
            self.finalize_obj.finalize()
            return

        elif status == "not_converging":
            print "Status:", status
            self.finalize_obj.report_status("failed","run_limit")
            print "Returning"
            sys.stdout.flush()
            return

        elif status == "incomplete":

            if task == "setup":
                self.setup() #TODO

            self.finalize_obj.report_status("started")
            (status, task) = calculation.run()

        else:
            self.finalize_obj.report_status("failed", "unknown")
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
            self.finalize_obj.report_status("failed", "run_limit")

            # print a local settings file, so that the run_limit can be extended if the
            #   convergence problems are fixed

            config_set_dir = self.config_obj.casm_directories.configuration_calc_settings_dir(self.config_obj.configname,
                                                                                              self.config_obj.clex,
                                                                                              self.config_obj.calc_subdir)

            try:
                os.makedirs(config_set_dir)
            except:
                pass
            settingsfile = os.path.join(config_set_dir, "calc.json")
            vaspwrapper.write_settings(self.config_obj.settings, settingsfile)

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

            # write results to properties.calc.json
            self.finalize_obj.finalize()

        else:
            self.finalize_obj.report_status("failed", "unknown")
            raise vaspwrapper.VaspWrapperError("vasp relaxation complete with unexpected status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()
