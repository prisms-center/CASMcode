""" FIXME """

import os
import math
import sys
import json
import pbs
import seqquest
import casm
from . import questwrapper

class Relax(object):
    """The Relax class contains functions for setting up, executing, and parsing a SeqQuest relaxation.

        The relaxation creates the following directory structure:
        config/
          calctype.name/
              run.0/
              run.1/
              ...
              run.final/

        'run.i' directories are only created when ready.
        'run.final' is a final constant volume run

        This automatically looks for VASP settings files in .../settings/calctype.name,
        where '...' is the nearest parent directory of 'self.configdir' in the CASM project repository

        Contains:
            self.configdir (.../config)
            self.calcdir   (.../config/calctype.name)

            self.settings = dictionary of settings for pbs and the relaxation, see questwrapper.read_settings

            self.auto = True if using pbs module's JobDB to manage pbs jobs
            self.sort = True if sorting atoms in POSCAR by type
    """
    def __init__(self, configdir=None, auto=True, sort=True):
        """
        Construct a SeqQuest relaxation job object.

        Args:
            configdir: path to configuration
            auto: True if using pbs module's JobDB to manage pbs jobs

        """
        if configdir is None:
            configdir = os.getcwd()

        print "Reading CASM settings"
        self.casm_settings = casm.casm_settings(configdir)
        if self.casm_settings is None:
            raise questwrapper.QuestWrapperError("Not in a CASM project.\
                                                  The file '.casm' directory was not found.")

        print "Constructing a CASM QuestWrapper Relax object"
        sys.stdout.flush()

        print "  Setting up directories"
        sys.stdout.flush()

        # store path to .../config, if not existing raise
        self.configdir = os.path.abspath(configdir)
        if not os.path.isdir(self.configdir):
            raise seqquest.SeqQuestError("Error in casm.quest.relax: Did not find directory: "
                                         + self.configdir)

        # store path to .../config/calctype.name, and create if not existing
        self.calcdir = os.path.join(self.configdir, self.casm_settings["curr_calctype"])
        try:
            os.mkdir(self.calcdir)
        except OSError:
            pass

        # read the settings json file
        print "  Reading relax.json settings file"
        sys.stdout.flush()
        setfile = casm.settings_path("relax.json", self.casm_settings["curr_calctype"],
                                     self.configdir)
        if setfile is None:
            raise questwrapper.QuestWrapperError("Could not find .../settings/"
                                                 + self.casm_settings["curr_calctype"]
                                                 + "/relax.json file.")
        self.settings = questwrapper.read_settings(setfile)

        # add required keys to settings if not present
        if not "run_cmd" in self.settings:
            self.settings["run_cmd"] = None
        if not "ncpus" in self.settings:
            self.settings["ncpus"] = None
        if not "run_limit" in self.settings:
            self.settings["run_limit"] = None

        self.auto = auto
        self.sort = sort
        print "SeqQuest Relax object constructed\n"
        sys.stdout.flush()

    def setup(self):
        """ Setup initial relaxation run

            Uses the following files from the most local .../settings/calctype.name directory:
                lcao.in: SeqQuest input settings
                SPECIES: info for each species such as which *.atm files to use, MAGMOM, GGA+U, etc.

            Uses the following files from the .../config directory:
                POS: structure of the configuration to be relaxed

        """
        # Find required input files in CASM project directory tree
        # self.species = species_settings(speciesfile)
        # self.lcao_in = LcaoIN(lcao_in_file, speciesfile=speciesfile, POS=super_poscarfile)
        lcao_in = casm.settings_path("lcao.in", self.casm_settings["curr_calctype"], self.configdir)
        super_poscarfile = os.path.join(self.configdir, "POS")
        speciesfile = casm.settings_path("SPECIES", self.casm_settings["curr_calctype"],
                                         self.configdir)

        # Verify that required input files exist
        if lcao_in is None:
            raise seqquest.SeqQuestError("Relax.setup failed. No lcao.in file found in CASM\
					 project.")
        if super_poscarfile is None:
            raise seqquest.SeqQuestError("Relax.setup failed. No POS file found for this\
					 configuration.")
        if speciesfile is None:
            raise seqquest.SeqQuestError("Relax.setup failed. No SPECIES file found in CASM\
 					 project.")

        # Find optional input files
        extra_input_files = []
        for s in self.settings["extra_input_files"]:
            extra_input_files.append(casm.settings_path(s, self.casm_settings["curr_calctype"],
                                                        self.configdir))
            if extra_input_files[-1] is None:
                raise seqquest.SeqQuestError("Relax.setup failed. Extra input file "
                                             + s + " not found in CASM project.")
        if self.settings["initial"]:
            extra_input_files += [casm.settings_path(self.settings["initial"],
                                                     self.casm_settings["curr_calctype"],
                                                     self.configdir)]
            if extra_input_files[-1] is None:
                raise seqquest.SeqQuestError("Relax.setup failed. No initial lcao.in file "
                                             + self.settings["initial"] + " found in CASM project.")
        if self.settings["final"]:
            extra_input_files += [casm.settings_path(self.settings["final"],
                                                     self.casm_settings["curr_calctype"],
                                                     self.configdir)]
            if extra_input_files[-1] is None:
                raise seqquest.SeqQuestError("Relax.setup failed. No final lcao.in file "
                                             + self.settings["final"] + " found in CASM project.")

        if self.settings["cont_relax"] is not None:
            if os.path.isfile(os.path.join("..", "calctype"+self.settings["cont_relax"],
                                           "properties.calc.json")):
                with open(os.path.isfile(os.path.join("..", "calctype"+self.settings["cont_relax"],
                                                      "properties.calc.json"))) as stream:
                    final = json.load(stream)
                    with open(os.path.join(self.calcdir, self.settings["cont_relax"] + "_POS",
                                           "w")) as fake_pos:
                        fake_pos.write("Cont from " + self.settings["cont_relax"] + "\n")
                        fake_pos.write("1.0000000000 \n")
                        fake_pos.write("\n".join(["     " + "  ".join(map(str, x))
                                                  for x in final['relaxed_lattice']]) + "\n")
                        fake_pos.write(" ".join(final['atom_type'])+"\n")
                        fake_pos.write(" ".join(map(str, final['atoms_per_type']))+"\n")
                        fake_pos.write("Cartesian\n")
                        basis_idx = 0
                        for t, n in zip(final['atom_type'], final['atoms_per_type']):   #pylint: disable=invalid-name
                            for _ in range(n):
                                fake_pos.write("     "
                                               + "  ".join(map(str,
                                                               final['relaxed_basis'][basis_idx]))
                                               + "  " +  t + "\n")
                                basis_idx += 1

            if os.path.isfile(os.path.join(self.calcdir, self.settings["cont_relax"] + "_POS")):
                super_poscarfile = os.path.join(self.calcdir, self.settings["cont_relax"] + "_POS")

        sys.stdout.flush()

        seqquest.seqquest_io.SeqquestIO(lcao_in, super_poscarfile, speciesfile, extra_input_files).write(self.calcdir)

    def submit(self):   #pylint: disable=too-many-statements
        """Submit a PBS job for this SeqQuest relaxation"""

        # first, check if the job has already been submitted and is not completed
        db = pbs.JobDB()
        print "rundir", self.calcdir
        id = db.select_regex_id("rundir", self.calcdir)
        print "id:", id
        sys.stdout.flush()
        if id != []:
            for j in id:
                job = db.select_job(j)
                if job["jobstatus"] != "C":
                    print "JobID:", job["jobid"], "  Jobstatus:", job["jobstatus"], "  Not submitting."
                    sys.stdout.flush()
                    return

        # second, only submit a job if relaxation status is "incomplete"

        # construct the Relax object
        relaxation = seqquest.Relax(self.calcdir, self.run_settings())

        # check the current status
        (status, task) = relaxation.status()

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
            raise questwrapper.QuestWrapperError("unexpected relaxation status: '" + status
                                                 + "' and task: '" + task + "'")

        print "Preparing to submit a VASP relaxation PBS job"
        sys.stdout.flush()

        # cd to configdir, submit jobs from configdir, then cd back to currdir
        currdir = os.getcwd()
        os.chdir(self.calcdir)

        # determine the number of atoms in the configuration
        print "  Counting atoms in the POSCAR"
        sys.stdout.flush()
        geom = seqquest.seqquest_io.Geom.POS(os.path.join(self.configdir, "POS"))
        N = len(geom.basis)

        # Construct the run command
        command = "python -c \"import casm.questwrapper; casm.questwrapper.Relax('" + self.configdir + "').run()\""
        if self.settings["preamble"] is not None:
        # Append any instructions given in the 'preamble' file, if given
            preamble = casm.settings_path(self.settings["preamble"],
                                          self.casm_settings["curr_calctype"],
                                          self.configdir)
            with open(preamble) as my_preamble:
                command = "".join(my_preamble) + command
                # for line in my_preamble:
                #     command = line.rstrip() + "; " + command

        print "  Constructing a PBS job"
        sys.stdout.flush()
        # construct a pbs.Job
        job = pbs.Job(name=casm.jobname(self.configdir),\
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
                      command=command,\
                      auto=self.auto)

        print "  Submitting"
        sys.stdout.flush()
        # submit the job
        job.submit()
        self.report_status("submitted")

        # return to current directory
        os.chdir(currdir)

        print "CASM questwrapper relaxation PBS job submission complete\n"
        sys.stdout.flush()

    def run_settings(self):
        """ Set default values based on runtime environment"""
        settings = dict(self.settings)

        # set default values

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
        """ Setup input files, run a vasp relaxation, and report results """

        # construct the Relax object
        relaxation = seqquest.Relax(self.calcdir, self.run_settings())

        # check the current status
        (status, task) = relaxation.status()


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
            self.finalize()
            return

        elif status == "not_converging":
            print "Status:", status
            self.report_status("failed", "run_limit")
            print "Returning"
            sys.stdout.flush()
            return

        elif status == "incomplete":

            if task == "setup":
                self.setup()

            self.report_status("started")
            (status, task) = relaxation.run()

        else:
            self.report_status("failed", "unknown")
            raise questwrapper.QuestWrapperError("unexpected relaxation status: '"
                                                 + status + "' and task: '" + task + "'")


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
            self.report_status("failed", "run_limit")

            # print a local settings file, so that the run_limit can be extended if the
            #   convergence problems are fixed
            try:
                os.makedirs(os.path.join(self.configdir, "settings",
                                         self.casm_settings["curr_calctype"]))
            except OSError:
                pass
            settingsfile = os.path.join(self.configdir, "settings",
                                        self.casm_settings["curr_calctype"], "relax.json")
            questwrapper.write_settings(self.settings, settingsfile)

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
            self.finalize()

        else:
            self.report_status("failed", "unknown")
            raise questwrapper.QuestWrapperError("vasp relaxation complete with unexpected status: '"
                                                 + status + "' and task: '" + task + "'")

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
        with open(outputfile, 'w') as stream:
            stream.write(json.dumps(output, stream, cls=casm.NoIndentEncoder, indent=4,
                                    sort_keys=True))
        print "Wrote " + outputfile
        sys.stdout.flush()

    def finalize(self):
        """ Collect final run things """
        if self.is_converged():
            # write properties.calc.json
            rundir = os.path.join(self.calcdir, "run.final")
            output = self.properties(rundir)
            outputfile = os.path.join(self.calcdir, "properties.calc.json")
            with open(outputfile, 'w') as stream:
                stream.write(json.dumps(output, stream, cls=casm.NoIndentEncoder, indent=4,
                                        sort_keys=True))
            print "Wrote " + outputfile
            sys.stdout.flush()
            self.report_status('complete')

    def is_converged(self):
        """Check for electronic convergence in completed calculations. Returns True or False."""
        # Not currently implemented for SeqQuest

        # # Verify that the last relaxation reached electronic convergence
        # relaxation = vasp.Relax(self.calcdir, self.run_settings())
        # for i in range(len(relaxation.rundir)):
        #   try:
        #     vrun = vasp.io.Vasprun( os.path.join(self.calcdir, relaxation.rundir[-i-1], "vasprun.xml"))
        #     if len(vrun.all_e_0[-1]) >= vrun.nelm:
        #       print('The last relaxation run (' +
        #           os.path.basename(relaxation.rundir[-i-1]) +
        #           ') failed to achieve electronic convergence; properties.calc.json will not be written.\n')
        #       self.report_status('failed','electronic_convergence')
        #       return False
        #     break
        #   except:
        #     pass

        # # Verify that the final static run reached electronic convergence
        # vrun = vasp.io.Vasprun( os.path.join(self.calcdir, "run.final", "vasprun.xml") )
        # if len(vrun.all_e_0[0]) >= vrun.nelm:
        #     print('The final run failed to achieve electronic convergence; properties.calc.json will not be written.\n')
        #     self.report_status('failed','electronic_convergence')
        #     return False

        return True

    @staticmethod
    def properties(rundir):
        """Report results to properties.calc.json file in configuration directory, after checking for electronic convergence."""

        output = dict()
        ofile = seqquest.seqquest_io.LcaoOUT(os.path.join(rundir, "lcao.out"))
        ogeom = seqquest.seqquest_io.Geom.geom(os.path.join(rundir, "lcao.geom"))

        # the calculation is run on the 'sorted' POSCAR, need to report results 'unsorted'

        output["atom_type"] = ogeom.type_atoms
        output["atoms_per_type"] = ogeom.num_atoms
        output["coord_mode"] = ofile.coord_mode

        # fake unsort_dict (unsort_dict[i] == i)
        unsort_dict = dict(zip(range(0, len(ogeom.basis)), range(0, len(ogeom.basis))))

        # as lists
        output["relaxed_forces"] = [None for i in range(len(ofile.forces))]

        print casm.NoIndent(ofile.forces[0])
        for i, v in enumerate(ofile.forces):
            output["relaxed_forces"][unsort_dict[i]] = casm.NoIndent(ofile.forces[i])

        output["relaxed_lattice"] = [casm.NoIndent(v) for v in ofile.cell.lattice]

        output["relaxed_basis"] = [None for i in range(len(ogeom.basis))]
        for i, v in enumerate(ogeom.basis):
            output["relaxed_basis"][unsort_dict[i]] = casm.NoIndent(ogeom.basis[i].position)

        output["relaxed_energy"] = ofile.total_energy

        return output


