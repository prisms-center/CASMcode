import os
import shutil
import re
import subprocess
import sys
import time
import gzip
import warnings

import casm.project
import casm.aimswrapper
import casm.aims
import casm.aims.io.geometry

from casm.aims.io.io import DEFAULT_AIMS_GZIP_LIST, DEFAULT_AIMS_COPY_LIST, \
                            DEFAULT_AIMS_REMOVE_LIST, DEFAULT_AIMS_MOVE_LIST


class AimsError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class AimsWarning(Warning):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


def continue_job(jobdir, contdir, settings):
    """Use the files in job directory 'jobdir', to setup a job in directory 'contdir'.

       Args:
         jobdir: path to current job directory
         contdir: path to new directory to continue the job
         settings:
             copy: Files copied from 'jobdir' to 'contdir'
                It also copies either geometry.in.next_step or if that does not exist geometry.in
             move: Files moved from 'jobdir' to 'contdir'
             keep: Files (along with those in 'copy') to keep in 'jobdir'. The rest are removed.
    """

    configdir = os.getcwd()
    _res = os.path.split(configdir)
    cfgname = os.path.split(_res[0])[1] + "/" + _res[1]
    casm_dirs = casm.project.DirectoryStructure(configdir)
    casm_sets = casm.project.ProjectSettings(configdir)
    clex = casm_sets.default_clex
    aimsfiles = casm.aimswrapper.aims_input_file_names(casm_dirs, cfgname, clex)
    controlfile, prim_posfile, super_posfile, basisfile = aimsfiles

    print("Continue FHI-aims job:\n  Original: " + jobdir + "\n  Continuation: " + contdir)
    sys.stdout.flush()

    # remove duplicates
    move = list(set(settings['move'] + DEFAULT_AIMS_MOVE_LIST))
    copy = list(set(settings['copy'] + DEFAULT_AIMS_COPY_LIST))
    remove = list(set(settings['remove'] + DEFAULT_AIMS_REMOVE_LIST))
    compress = list(set(settings['compress']))

    # Check that the user isn't being contradictory or silly
    for f in move:
        if f == "geometry.in" or f == "geometry.in.next_step":
            raise AimsError("Error in casm.vasp.general.continue_job(). "
                            "Do not include geo.in and geo.in.next_step 'move'; use 'backup' if you want a backup")
        if f in remove:
            if f in DEFAULT_AIMS_MOVE_LIST:
                raise AimsError("Error in casm.aims.general.continue_job(). "
                                "%s cannot be removed, FHI-aims will not run!!!" % f)
            else:
                warnings.warn("Warning: %s found in both 'move' and 'remove'. "
                              "The file will not be removed." % f, AimsWarning)
                remove = list(set(remove) - set(f))
        if f in copy:
            if f in DEFAULT_AIMS_MOVE_LIST:
                warnings.warn("Warning: %s found in both 'move' and 'copy'. "
                              "The fill will be moved only." % f, AimsWarning)
                copy = list(set(copy) - set(f))
            else:
                warnings.warn("Warning: %s found in both 'move' and 'copy'. "
                              "The file will be copied only." % f, AimsWarning)
                move = list(set(move) - set(f))
        if f in compress:
            if f in DEFAULT_AIMS_GZIP_LIST:
                raise AimsError("Error in casm.aims.general.continue_job(). "
                                "%s cannot be compressed, FHI-aims will not run!!!" % f)
            else:
                raise AimsError("Error in casm.aims.general.continue_job(). "
                                "%s found in both 'move' and 'compress', but these options contradict. "
                                "Did you mean 'backup'?" % f)

    # make the new contdir
    if not os.path.isdir(contdir):
        os.mkdir(contdir)

    # copy/move files
    if os.path.isfile(os.path.join(jobdir, "geometry.in.next_step")) and \
            os.path.getsize(os.path.join(jobdir, "geometry.in.next_step")) > 0:
        shutil.move(os.path.join(jobdir, "control.in"), os.path.join(contdir, "control.in"))
        my_basis = casm.aims.io.basis.basis_settings(basisfile)
        newgeo = casm.aims.io.geometry.Geometry(os.path.join(jobdir, "geometry.in.next_step"), my_basis)
        newgeo.write(os.path.join(contdir, "geometry.in"), my_basis)
        print("  cp geometry.in.next_step -> geometry.in (and added default moments)")
    else:
        shutil.copyfile(os.path.join(jobdir, "geometry.in"), os.path.join(contdir, "geometry.in"))
        print("  no geometry.in.next_step, took old geometry.in (this means no relaxation step was performed)")
        shutil.move(os.path.join(jobdir, "control.in"), os.path.join(contdir, "control.in"))


def complete_job(jobdir, settings):
    """Clean up and compress output

       Args:
         jobdir: path to current job directory
         settings: casm settings
    """
    print("Completing FHI-aims job: " + jobdir)
    sys.stdout.flush()

    # compress files
    print(" gzip: ", end='')
    for file in settings["compress"]:
        if os.path.isfile(os.path.join(jobdir, file)):
            print(file, end='')
            # Open target file, target file.gz
            f_in = open(os.path.join(jobdir, file), 'rb')
            f_out = gzip.open(os.path.join(jobdir, file)+'.gz', 'wb')
            # Compress, close files
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
            # Remove original target file
            os.remove(os.path.join(jobdir, file))
    print(" gzipping DONE")
    print("\n")
    sys.stdout.flush()


class FreezeError(object):
    """VASP appears frozen"""
    def __init__(self):
        self.pattern = None

    def __str__(self):
        return "FHI-aims appears to be frozen"

    @staticmethod
    def error(jobdir=None):
        """ Check if aims is frozen

        Args:
            jobdir: job directory
        Returns:
            True if:
                1) no file has been modified for 5 minutes
                2) 'LOOP+' exists in OUTCAR and no output file has been modified
                   in 5x the time for the slowest loop
        """

        # Check if any files modified in last 300 s
        most_recent = None
        most_recent_file = None
        for f in os.listdir(jobdir):
            t = time.time() - os.path.getmtime(os.path.join(jobdir, f))
            if most_recent is None:
                most_recent = t
                most_recent_file = f
            elif t < most_recent:
                most_recent = t
                most_recent_file = f

        print("Most recent file output (" + most_recent_file + "):", most_recent, " seconds ago.")
        sys.stdout.flush()
        if most_recent < 300:
            return False

    @staticmethod
    def fix(err_jobdir, new_jobdir, settings):
        """ Fix by killing the job and resubmitting."""
        print("  Kill job and continue...")
        continue_job(err_jobdir, new_jobdir, settings)


def error_check(jobdir, stdoutfile):
    """ Check vasp stdout for errors """
    err = dict()

    # Error to check line by line, only look for first of each type
    sout = open(stdoutfile, 'r')

    # Error to check for once
    possible = [FreezeError()]
    for p in possible:
        if p.error(jobdir=jobdir):
            err[p.__class__.__name__] = p

    sout.close()
    if len(err) == 0:
        return None
    else:
        return err


def run(jobdir=None, stdout="std.out", stderr="std.err", command=None, ncpus=None,
        poll_check_time=5.0, err_check_time=60.0):
    """ Run FHI-aims using subprocess.

        The 'command' is executed in the directory 'jobdir'.

        Args:
            jobdir:     directory to run.  If jobdir is None, the current directory is used.
            stdout:     filename to write to.  If stdout is None, "std.out" is used.
            stderr:     filename to write to.  If stderr is None, "std.err" is used.
            command:    (str or None) FHI-aims execution command
                        If command != None: then 'command' is run in a subprocess
                        Else, if ncpus == 1, then command = "aims"
                        Else, command = "mpirun -np {NCPUS} aims"
            ncpus:      (int) if '{NCPUS}' is in 'command' string, then 'ncpus' is substituted in the command.
                        if ncpus==None, $PBS_NP is used if it exists, else 1
            poll_check_time: how frequently to check if the vasp job is completed
            err_check_time: how frequently to parse vasp output to check for errors
    """
    print("Begin FHI-aims run:")
    sys.stdout.flush()

    if jobdir is None:
        jobdir = os.getcwd()

    currdir = os.getcwd()
    os.chdir(jobdir)

    if ncpus is None:
        if "PBS_NP" in os.environ:
            ncpus = os.environ["PBS_NP"]
        else:
            ncpus = 1

    if command is None:
        if ncpus == 1:
            command = "aims"
        else:
            command = "mpirun -np {NCPUS} aims"

    if re.search("NCPUS", command):
        command = command.format(NCPUS=str(ncpus))

    print("  jobdir:", jobdir)
    print("  exec:", command)
    sys.stdout.flush()

    sout = open(os.path.join(jobdir, stdout), 'w')
    serr = open(os.path.join(jobdir, stderr), 'w')
    err = None
    p = subprocess.Popen(command.split(), stdout=sout, stderr=serr)

    # wait for process to end, and periodically check for errors
    last_check = time.time()
    while p.poll() is None:
        time.sleep(poll_check_time)
        if time.time() - last_check > err_check_time:
            last_check = time.time()
            err = error_check(jobdir, os.path.join(jobdir, stdout))
            if err is not None:
                # FreezeErrors are fatal and usually not helped with abort_scf
                if "FreezeError" in err.keys():
                    print("  FHI-aims seems frozen, killing job")
                    sys.stdout.flush()
                    p.kill()

    # close output files
    sout.close()
    serr.close()

    os.chdir(currdir)

    print("Run complete")
    sys.stdout.flush()

    # check finished job for errors
    if err is None:
        err = error_check(jobdir, os.path.join(jobdir, stdout))
        if err is not None:
            print("  Found errors:", end='')
            for e in err:
                print(e, end='')
        print("\n")
    return err
