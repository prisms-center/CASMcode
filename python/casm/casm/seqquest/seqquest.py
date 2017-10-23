""" FIXME """
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os
import shutil
import re
import subprocess
import sys
import time
import gzip
import warnings
import signal
from . import seqquest_io
from .seqquest_io import (QUEST_INPUT_FILE_LIST, DEFAULT_QUEST_MOVE_LIST,
                          DEFAULT_QUEST_COPY_LIST, DEFAULT_QUEST_REMOVE_LIST)
class SeqQuestError(Exception):
    """ Class for SeqQuest errors """
    pass

class SeqQuestWarning(Warning):
    """ Class for SeqQuest warnings """
    pass

def continue_job(jobdir, contdir, settings):
    """Use the files in seqquest job directory 'jobdir', to setup a seqquest job in directory
        'contdir'.

       Args:
         jobdir: path to current job directory
         contdir: path to new directory to continue the job
       Settings:
         copy: Files copied from 'jobdir' to 'contdir'
                    It also copies either CONTCAR, or if that does not exist POSCAR.
         move: Files moved from 'jobdir' to 'contdir'
         keep: Files (along with those in 'copy') to keep in 'jobdir'. The rest are removed.

         Do not include "POSCAR" or "CONTCAR" in 'move' or 'copy'. jobdir/CONTCAR is copied
         to contdir/POSCAR, unless it does not exist, in which case jobdir/POSCAR is copied
         to contdir/POSCAR. jobdir/POSCAR and jobdir/CONTCAR are always kept.
    """

    print("Continue SeqQuest job:\n  Original: " + jobdir + "\n  Continuation: " + contdir)
    sys.stdout.flush()

    # remove duplicates
    move = list(set(settings['move']+DEFAULT_QUEST_MOVE_LIST+settings['extra_input_files']))
    copy = list(set(settings['copy']+DEFAULT_QUEST_COPY_LIST))
    remove = list(set(settings['remove']))
    compress = list(set(settings['compress']))
    backup = list(set(settings['backup']))

    # Check that necessary files are being moved/copied: lcao.in, lcao.geom_in
    if not "lcao.in" in move + copy:
        warnings.warn("Warning: lcao.in not found in either 'move' or 'copy'.\
                       Moving POTCAR by default", SeqQuestWarning)
        move += ["lcao.in"]

    # make the new contdir
    try:
        os.mkdir(contdir)
    except:
        pass

    # make compressed backups of files sensitive to corruption (e.g. WAVECAR)
    print(" backup:")
    for file in backup:
        if os.path.isfile(os.path.join(jobdir, file)):
            print(file, end=' ')
            # Open target file, target file.gz
            f_in = open(os.path.join(jobdir, file), 'rb')
            f_out = gzip.open(os.path.join(jobdir, file)+'_BACKUP.gz', 'wb')
            # Compress, close files
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
    print("")


    # copy geom/lcao.in/etc
    # if os.path.isfile(os.path.join(jobdir, "lcao.geom")) and os.path.getsize(os.path.join(jobdir, "lcao.geom")) > 0:
    if "lcao.geom" in os.listdir(jobdir):
        if os.path.getsize(os.path.join(jobdir, "lcao.geom")) > 0:
            shutil.copyfile(os.path.join(jobdir, "lcao.geom"), os.path.join(contdir, "lcao.geom_in"))
            print("  cp lcao.geom -> lcao.geom_in")
        else:
            shutil.copyfile(os.path.join(jobdir, "lcao.geom_in"), os.path.join(contdir, "lcao.geom_in"))
            print("  no lcao.geom : cp lcao.geom_in -> lcao.geom_in")
    else:
        shutil.copyfile(os.path.join(jobdir, "lcao.geom_in"), os.path.join(contdir, "lcao.geom_in"))
        print("  no lcao.geom : cp lcao.geom_in -> lcao.geom_in")

    # copy .atm files
    for atm in os.listdir(jobdir):
        if re.search(r".*\.atm", atm):
            if os.path.islink(os.path.join(jobdir, atm)):
                os.symlink(os.readlink(os.path.join(jobdir, atm)), os.path.join(contdir, atm))
            else:
                shutil.copy(os.path.join(jobdir, atm), os.path.join(contdir, atm))

    # move files
    print("  mv:", end=' ')
    # for file in move:
    #     print file,
    #     os.rename(os.path.join(jobdir, file), os.path.join(contdir, file))
    for mname in move:
        for mfile in os.listdir(jobdir):
            if re.match(mname, mfile):
                print(mfile, end=' ')
                os.rename(os.path.join(jobdir, mfile), os.path.join(contdir, mfile))
    print("")

    # copy files
    print("  cp:", end=' ')
    # for file in copy:
    #     print file,
    #     shutil.copyfile(os.path.join(jobdir, file), os.path.join(contdir, file))
    for cname in copy:
        for cfile in os.listdir(jobdir):
            if re.match(cname, cfile):
                print(cfile, end=' ')
                shutil.copyfile(os.path.join(jobdir, cfile), os.path.join(contdir, cfile))
    print("")

    # remove files
    print("  rm:", end=' ')
    # for file in remove:
    #     if os.path.isfile(os.path.join(jobdir, file)):
    #         print file,
    #         os.remove(os.path.join(jobdir, file))
    for rname in remove:
        for rfile in os.listdir(jobdir):
            if re.match(rname, rfile):
                print(rfile, end=' ')
                os.remove(os.path.join(jobdir, rfile))
    print("")

    # compress files
    print(" gzip:", end=' ')
    # for file in compress:
    #     if os.path.isfile(os.path.join(jobdir, file)):
    #         print file,
    #         # Open target file, target file.gz
    #         f_in = open(os.path.join(jobdir, file), 'rb')
    #         f_out = gzip.open(os.path.join(jobdir, file)+'.gz', 'wb')
    #         # Compress, close files
    #         f_out.writelines(f_in)
    #         f_out.close()
    #         f_in.close()
    #         # Remove original target file
    #         os.remove(os.path.join(jobdir, file))
    for cname in compress:
        for cfile in os.listdir(jobdir):
            if re.match(cname, cfile):
                print(cfile, end=' ')
                # Open target file, target file.gz
                f_in = open(os.path.join(jobdir, cfile), 'rb')
                f_out = gzip.open(os.path.join(jobdir, cfile)+'.gz', 'wb')
                # Compress, close files
                f_out.writelines(f_in)
                f_out.close()
                f_in.close()
                # Remove original target file
                os.remove(os.path.join(jobdir, cfile))
    print("")
    print("")
    sys.stdout.flush()

def complete_job(jobdir, settings):
# def complete_job(jobdir,\
#                  remove = io.DEFAULT_VASP_REMOVE_LIST + ["POTCAR"],\
#                  compress = []):
    """Remove files from a quest job directory

       Args:
         jobdir: path to current job directory
       Settings:
         copy: Does nothing
         move: Does nothing
         compress: Compresses listed files
         backup: Does nothing
         remove: Deletes listed files
    """

    print("Complete SeqQuest job: " + jobdir)
    sys.stdout.flush()

    # remove files
    print("  rm:", end=' ')
    for f in settings["remove"]:
        if not f in (settings["copy"] + settings["move"] + settings["compress"] + settings["backup"]):
            if os.path.isfile(os.path.join(jobdir,f)):
                print(f, end=' ')
                os.remove(os.path.join(jobdir,f))
    for f in settings["extra_input_files"]:
        if os.path.isfile(os.path.join(jobdir, f)):
            print(f, end=' ')
            os.remove(os.path.join(jobdir,f))
    print("")

    # compress files
    print(" gzip:", end=' ')
    for file in settings["compress"]:
        if os.path.isfile(os.path.join(jobdir,file)):
            print(file, end=' ')
            # Open target file, target file.gz
            f_in = open(os.path.join(jobdir, file), 'rb')
            f_out = gzip.open(os.path.join(jobdir, file)+'.gz', 'wb')
            # Compress, close files
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
            # Remove original target file
            os.remove(os.path.join(jobdir,file))
    print("")
    print("")
    sys.stdout.flush()

def run(jobdir=None, stdout="std.out", stderr="std.err", command=None, ncpus=None, poll_check_time=5.0, err_check_time=60.0, err_types=None):
    """ Run seqquest using subprocess.

        The 'command' is executed in the directory 'jobdir'.

        Args:
            jobdir:     directory to run vasp.  If jobdir == None, the current directory is used.
            stdout:     filename to write to.  If stdout == None, "std.out" is used.
            stderr:     filename to write to.  If stderr == None, "std.err" is used.
            command:    (str or None) vasp execution command
                        If command != None: then 'command' is run in a subprocess
                        Else, if ncpus == 1, then command = "quest"
                        Else, command = "mpirun -np {NCPUS} quest"
            ncpus:      (int) if '{NCPUS}' is in 'command' string, then 'ncpus' is substituted in the command.
                        if ncpus==None, $PBS_NP is used if it exists, else 1
            poll_check_time: how frequently to check if the vasp job is completed
            err_check_time: how frequently to parse vasp output to check for errors
            err_types:  List of error types to check for. Supported errors: 'IbzkptError', 'SubSpaceMatrixError'. Default: None, in which case only SubSpaceMatrixErrors are checked.

    """
    print("Begin quest run:")
    sys.stdout.flush()

    if jobdir is None:
        jobdir = os.getcwd()

    currdir = os.getcwd()
    os.chdir(jobdir)

    if ncpus is None:
        if "PBS_NP" in os.environ:
            ncpus = os.environ["PBS_NP"]
        elif "SLURM_NTASKS" in os.environ:
            ncpus = os.environ["SLURM_NTASKS"]
        else:
            ncpus = 1

    if command is None:
        if ncpus == 1:
            command = "quest"
        else:
            command = "mpirun -np {NCPUS} quest"

    if re.search(r"\{NCPUS\}", command):
        command = command.format(NCPUS=str(ncpus))

    ### Expand remaining environment variables
    command = os.path.expandvars(command)

    print("  jobdir:", jobdir)
    print("  exec:", command)
    sys.stdout.flush()

    sout = open(os.path.join(jobdir, stdout), 'w')
    serr = open(os.path.join(jobdir, stderr), 'w')
    err = None

    proc = subprocess.Popen(command.split(), stdout=sout, stderr=serr)

    # wait for process to end, and periodically check for errors
    poll = proc.poll()
    last_check = time.time()
    # stopcar_time = None
    while poll is None:
        time.sleep(poll_check_time)

        if time.time() - last_check > err_check_time:
            last_check = time.time()
            # Error handling not yet implemented
            # err = error_check(jobdir, os.path.join(jobdir,stdout), err_types)
            # if err != None:
            #     # FreezeErrors are fatal and usually not helped with STOPCAR
            #     if "FreezeError" in err.keys():
            #         print "  VASP is frozen, killing job"
            #         sys.stdout.flush()
            #         # Sometimes p.kill doesn't work if the process is on multiple nodes
            #         os.kill(p.pid, signal.SIGKILL)
            #         p.kill()
            #         # If the job is re-invoked (e.g. via mpirun or srun) too quickly
            #         #   after the previous job ended, infinitiband clusters can have
            #         #   some issues with resource allocation. A 30s sleep solves this.
            #         time.sleep(60)
            #     # Other errors can be killed with STOPCAR, which is safer
            #     elif stopcar_time == None:
            #         print "  Found errors:",
            #         for e in err:
            #             print e,
            #         print "\n  Killing job with STOPCAR"
            #         sys.stdout.flush()
            #         io.write_stopcar('e', jobdir)
            #         stopcar_time = time.time()
            #         time.sleep(60)
            #     # If the STOPCAR exists, wait 5 min before manually killing the job
            #     elif time.time() - stopcar_time > 300:
            #         print "  VASP is non-responsive, killing job"
            #         sys.stdout.flush()
            #         os.kill(p.pid, signal.SIGKILL)
            #         p.kill()
            #         # If the job is re-invoked (e.g. via mpirun or srun) too quickly
            #         #   after the previous job ended, infinitiband clusters can have
            #         #   some issues with resource allocation. A 30s sleep solves this.
            #         time.sleep(60)

        poll = proc.poll()

    # close output files
    sout.close()
    serr.close()

    os.chdir(currdir)

    print("Run complete")
    sys.stdout.flush()

    # check finished job for errors
    # if err == None:
    #     err = error_check(jobdir, os.path.join(jobdir,stdout), err_types)
    #     if err != None:
    #         print "  Found errors:",
    #         for e in err:
    #             print e,
    #     print "\n"

    # return err
    return None


