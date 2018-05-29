"""Defines the neb module methods"""
import pbs
import numpy as np
import os
import sys
import json
import shutil
import vasp
import casm
import pandas
from casm.project import Project, Selection
import vaspwrapper
from casm.vaspwrapper.vasp_calculator_base import VaspCalculatorBase
from vasp import Neb as calculator
from vasp import Relax as fake_calculator
    
def shuffle_endpoint_props(myjson,pos):
    positions=map(lambda site: site.position.tolist(),pos.basis)
    sort_inds=[ p[1] for p in sorted(zip(positions,range(len(positions))),key=lambda pair: pair[0])]
    unsort_inds=np.argsort(sort_inds)
    prev_basis=myjson["relaxed_basis"]
    prev_basis_sorted=sorted(prev_basis)
    myjson["relaxed_basis"]=[prev_basis_sorted[i] for i in unsort_inds]
    if "relaxed_forces" in myjson.keys():
        prev_forces=myjson["relaxed_forces"]
        prev_forces_sorted=sorted(prev_forces)
        myjson["relaxed_forces"]=[prev_forces_sorted[i] for i in unsort_inds]
    return myjson


class Neb(VaspCalculatorBase):
    """
    The Neb class contains functions for setting up, executing, and parsing a VASP neb calculation.

    Attributes
    ----------
    selection : casm.project.Selection
        selection of configuration
    calctype : string
        calctype to setup and run the neb calculations
    auto : bool
    sort : bool

    Methods
    -------
    from_configuration_dir(configdir='string', calctype='string', bool, bool)
        returns a instance of the Neb class instantited with a single configuration
    config_properties(config_data=dict/Pandas.DataFrame)
        return a dict of the properties required to setup a configuration
    pre_setup
        creates folder and makes POS files for each image
    setup
        sets up the input vasp files for the selection
    config_setup
        sets up the input vasp files for a single configuration
    get_vasp_input_files(config_data=dict/Pandas.DataFrame, settings=dict)
        returns filenames of a vasp neb calculation
    submit
        submit a job for each configuration
    run
        runs the neb calcutation on the selection
    report
        reports results for the selection
    run_cmd(configdir='string', calctype='string')
        return a string of command to run a single configuration
    finalize(config_data=dict/pandas_data, super_poscarfile='string')
        checks convergnce and write a properties file for the selection
    properties(calcdir='string', super_poscarfile='string', speciesfile='string')
        return a dict containing all the relaxed properties for a configuration

    Notes
    -----
    The calculation creates the following directory structure for each configuration:
    config/
        calctype.name/
            N_images_`n_images`\
                run.0/
                ....

    'run.i' directories are only created when ready.

    This automatically looks for VASP settings files using:
    casm.project.DirectoryStructure.settings_path_crawl

    The class inhertes from VaspCalculatorbase and methods overload the functionality in the parent

    """
    def __init__(self, selection, calctype=None, auto=True, sort=True):
        """Construct a VASP neb job object."""
        print "Construct a casm.vaspwrapper.Neb instance:"
        VaspCalculatorBase.__init__(self, selection, calctype, auto, sort)
        self.results_subdir = '01'
        self.calculator = calculator
        self.fake_calculator = fake_calculator

    def config_properties(self, config_data):
        """return configuration properties as a dict"""
        config_dict = super(Neb, self).config_properties(config_data)
        try:
            n_images = json.load(open(config_dict["setfile"]))["n_images"]
            endstate_calctype = json.load(open(config_dict["setfile"]))["endstate_calctype"]
        except:
            ## Error message if "n_images" not present in settings
            raise vaspwrapper.VaspWrapperError("Could not find \"n_images\" in \"calc.json\" in an appropriate \"settings\" directory")
        config_dict["n_images"] = n_images
        config_dict["endstate_calctype"] = endstate_calctype
        config_dict["calcdir"] = os.path.join(config_dict["calcdir"], "N_images_{}".format(n_images))
        return config_dict

    def pre_setup(self):
        """Makes settings file for interpolation and POS folder/files for each image"""
        proj = self.selection.proj
        try:
            os.mkdir(os.path.join(proj.path, ".casm/tmp"))
        except:
            pass
        sel_tmp = self.selection.saveas(os.path.join(proj.path, ".casm/tmp", "neb_interpolation_selection_tmp"), force=True)
        dict = {}
        dict["selection"] = os.path.join(proj.path, ".casm/tmp", "neb_interpolation_selection_tmp")
        for index, config_data in self.selection.data.iterrows():
            conf_dict = {"n_images" : config_data["n_images"],
                         "endstate_calctype" : config_data["endstate_calctype"]}
            dict[config_data["name"]] = conf_dict
            #try:
            #    os.makedirs(config_data["calcdir"])
            #except:
            #    pass
            #try:
            #    for i in range(config_data["n_images"] + 2):
                   #os.makedirs(os.path.join(config_data["calcdir"], 'poscars', str(i).zfill(2)))
            #except:
            #    pass
        dict["n_images"] = conf_dict["n_images"]
        dict["endstate_calctype"] = conf_dict["endstate_calctype"]
        dict["calctype"]=self.calctype
        tmp_folder = os.path.join(proj.path, '.casm/tmp')
        filename = "neb_interpolation_settings.json"
        with open(os.path.join(tmp_folder, filename), 'w') as file:
            file.write(json.dumps(dict, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))

        ## write the selection Interpolation command
        args = "enum --method DiffTransConfigInterpolation -s {}".format(os.path.join(tmp_folder, filename))
        output = proj.command(args)
        print "Interpolation complete"
        if os.path.isfile(os.path.join(tmp_folder, filename)):
            os.remove(os.path.join(tmp_folder, filename))

    def setup(self):
        """ Setup initial relaxation run for the selection"""
        super(Neb, self).setup()

    def config_setup(self, config_data):
        """ Setup initial relaxation run for a configuration

            Uses the following files from the most local .../settings/calctype.name directory:
                INCAR: VASP input settings
                KPOINTS: VASP kpoints settings
                POSCAR: reference for KPOINTS if KPOINTS mode is not A/AUTO/Automatic
                SPECIES: info for each species such as which POTCAR files to use, MAGMOM, GGA+U, etc.

            Uses the following files from the .../config_calcdir/00 directory:
                POSCAR: sample structure of the configuration to be relaxed

        """
        settings = self.read_settings(config_data["setfile"])
        difftransconfigname=config_data["name"]
        sel_tmp = Selection(self.selection.proj, "EMPTY", "diff_trans_config", False)
        sel_tmp.data = pandas.DataFrame({"name":difftransconfigname, "selected":1},index=range(1))
        sel_tmp = sel_tmp.saveas(os.path.join(self.selection.proj.path, ".casm/tmp","mirrors"),True)
        sel_config = Selection(self.selection.proj,os.path.join(self.selection.proj.path, ".casm/tmp","mirrors"), "diff_trans_config", False)
        sel_config.query(["and(rs(from_configname,to_configname),rs(to_configname,from_configname))"])
        symmetric=bool(sel_config.data["and(rs(from_configname,to_configname),rs(to_configname,from_configname))"].loc[0])
        override_mirrors=False
        if "override_mirrors" in settings.keys():
            override_mirrors=settings["override_mirrors"]
        if (not override_mirrors and symmetric):
            config_data["calcdir"] = config_data["calcdir"][:-1] + "1"
        super(Neb, self).config_setup(config_data)
        ## settings the images tag in incar file
        
        tmp_dict = {"images": config_data["n_images"]}
        vasp.io.set_incar_tag(tmp_dict, config_data["calcdir"])

    def get_vasp_input_files(self, config_data, settings):
        """returns filenames of a vasp neb calculation"""
        vaspfiles = super(Neb, self).get_vasp_input_files(config_data, settings)
        incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, extra_input_files = vaspfiles
        super_poscarfile = os.path.join(config_data["calcdir"], "poscars", "01", "POSCAR")
        return incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile, extra_input_files

    def submit(self):
        """submit a job for each configuration"""
        self.pre_setup()
        db = pbs.JobDB()
        for index,config_data in self.selection.data.iterrows():
            settings = self.read_settings(config_data["setfile"])
            difftransconfigname=config_data["name"]
            sel_tmp = Selection(self.selection.proj, "EMPTY", "diff_trans_config", False)
            sel_tmp.data = pandas.DataFrame({"name":difftransconfigname, "selected":1},index=range(1))
            sel_tmp = sel_tmp.saveas(os.path.join(self.selection.proj.path, ".casm/tmp","mirrors"),True)
            sel_config = Selection(self.selection.proj,os.path.join(self.selection.proj.path, ".casm/tmp","mirrors"), "diff_trans_config", False)
            sel_config.query(["and(rs(from_configname,to_configname),rs(to_configname,from_configname))"])
            symmetric=bool(sel_config.data["and(rs(from_configname,to_configname),rs(to_configname,from_configname))"].loc[0])
            override_mirrors=False
            if "override_mirrors" in settings.keys():
                override_mirrors=settings["override_mirrors"]
            if (not override_mirrors and symmetric):
                config_data["calcdir"] = config_data["calcdir"][:-1] + "1"
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
            if (not override_mirrors and symmetric):
                settings["subdir"] ="01"
                calculation = self.fake_calculator(config_data["calcdir"], self.run_settings(settings))
            else:
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
                    if (is_converged(calculation)):
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

            self.config_setup(config_data)
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
            job = pbs.Job(name=config_data["name"].replace("/",".")),\
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
               

    def run(self):
        """runs the neb calculation on the selection"""
        for index,config_data in self.selection.data.iterrows():
            settings = self.read_settings(config_data["setfile"])
            difftransconfigname=config_data["name"]
            sel_tmp = Selection(self.selection.proj, "EMPTY", "diff_trans_config", False)
            sel_tmp.data = pandas.DataFrame({"name":difftransconfigname, "selected":1},index=range(1))
            sel_tmp = sel_tmp.saveas(os.path.join(self.selection.proj.path, ".casm/tmp","mirrors"),True)
            sel_config = Selection(self.selection.proj,os.path.join(self.selection.proj.path, ".casm/tmp","mirrors"), "diff_trans_config", False)
            sel_config.query(["and(rs(from_configname,to_configname),rs(to_configname,from_configname))"])
            symmetric=bool(sel_config.data["and(rs(from_configname,to_configname),rs(to_configname,from_configname))"].loc[0])
            override_mirrors=False
            if "override_mirrors" in settings.keys():
                override_mirrors=settings["override_mirrors"]
            if (not override_mirrors and symmetric):
                config_data["calcdir"] = config_data["calcdir"][:-1] + "1"
            if (not override_mirrors and symmetric):
                settings["subdir"] ="01"
                calculation = self.fake_calculator(config_data["calcdir"], self.run_settings(settings))
            else:
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
                if (is_converged(calculation)):
                    self.finalize(config_data)
                continue

            elif status == "not_converging":
                print "Status:", status
                self.report_status(config_data["calcdir"], "failed", "run_limit")
                print "Returning"
                sys.stdout.flush()
                continue

            elif status == "incomplete":


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
                if is_converged(calculation):
                    self.finalize(config_data)

            else:
                self.report_status(config_data["calcdir"], "failed", "unknown")
                raise vaspwrapper.VaspWrapperError("vasp relaxation complete with unexpected status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()

    def report(self):
        """reports results for the selection"""
        super(Neb, self).report()

    @staticmethod
    def run_cmd(configdir, calctype):
        """return a string of command to run a single configuration"""
        return "python -c \"import casm.vaspwrapper; obj = casm.vaspwrapper.Neb.from_configuration_dir('{0}', '{1}'); obj.run()\"\n".format(configdir,
                                                                                                                                            calctype)

    def finalize(self, config_data, super_poscarfile=None):
        """checks convergnce and write a properties file for the selection"""
        if super_poscarfile is None:
            super_poscarfile = os.path.join(config_data["calcdir"], "run.final/00/POSCAR")
        super(Neb, self).finalize(config_data, super_poscarfile)
        images_calculated = [int(i.strip().split('_')[-1]) for i in os.listdir(os.path.split(config_data["calcdir"])[0]) if "N_images" in i]
        num_images = config_data["n_images"]
        if num_images == max(images_calculated):
            shutil.copy(os.path.join(config_data["calcdir"], "properties.calc.json"),
                        os.path.join(os.path.split(config_data["calcdir"])[0],
                                     "properties.calc.json"))
            print "As the present run has highest number of images copied {0} to {1}".format(os.path.join(config_data["calcdir"],
                                                                                                          "properties.calc.json"),
                                                                                             os.path.join(os.path.split(config_data["calcdir"])[0],
                                                                                                          "properties.calc.json"))
        sys.stdout.flush()

    def properties(self, calcdir, super_poscarfile=None, speciesfile=None):
        """Make properties output as a list of dict of each image properties"""
        final_output = {}
        endpts = json.load(open(os.path.join(calcdir,"../endpoint_props.json")))
        from_pos=os.path.join(calcdir,"00/POSCAR")
        final_output["00"]=shuffle_endpoint_props(endpts["0"],vasp.io.poscar.Poscar(from_pos))
        num_images = vasp.io.get_incar_tag("IMAGES", calcdir)
        to_pos=os.path.join(calcdir,str(num_images+1).zfill(2)+"/POSCAR")
        for img in [str(j).zfill(2) for j in range(1, num_images+1)]:
            vaspdir = calcdir + "/{}".format(img)
            output = super(Neb, self).properties(vaspdir, super_poscarfile, speciesfile)
            final_output[img] = output
        final_output[str(num_images+1).zfill(2)]=shuffle_endpoint_props(endpts[str(num_images+1)],vasp.io.poscar.Poscar(to_pos))
        return final_output
