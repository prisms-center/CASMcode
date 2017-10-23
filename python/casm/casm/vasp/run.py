""" Job manipuation routines for VASP"""
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

from casm.vasp.error import VaspError, VaspWarning, error_check, crash_check
from casm.vasp import io

def complete_job(jobdir, settings):
    """Remove files from a vasp job directory

       Args:
         jobdir: path to current job directory
       Settings:
         copy: Does nothing
         move: Does nothing
         compress: Compresses listed files
         backup: Does nothing
         remove: Deletes listed files
    """

    print("Complete VASP job: " + jobdir)
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

def run(jobdir = None, stdout = "std.out", stderr = "std.err", npar=None, ncore=None, command=None, ncpus=None, kpar=None, poll_check_time = 5.0, err_check_time = 60.0, err_types=None):
    """ Run vasp using subprocess.

        The 'command' is executed in the directory 'jobdir'.

        Args:
            jobdir:     directory to run vasp.  If jobdir is None, the current directory is used.
            stdout:     filename to write to.  If stdout is None, "std.out" is used.
            stderr:     filename to write to.  If stderr is None, "std.err" is used.
            npar:       (int or None) VASP INCAR NPAR setting. If npar is None, then NPAR is removed from INCAR
            kpar:       (int or None) VASP INCAR KPAR setting. If kpar is None, then KPAR is removed from INCAR
            ncore:      (int or None) VASP INCAR NCORE setting. If not npar is None or ncore is None, then NCORE is removed from INCAR
            command:    (str or None) vasp execution command
                        If command != None: then 'command' is run in a subprocess
                        Else, if ncpus == 1, then command = "vasp"
                        Else, command = "mpirun -np {NCPUS} vasp"
            ncpus:      (int) if '{NCPUS}' is in 'command' string, then 'ncpus' is substituted in the command.
                        if ncpus==None, $PBS_NP is used if it exists, else 1
            poll_check_time: how frequently to check if the vasp job is completed
            err_check_time: how frequently to parse vasp output to check for errors
            err_types:  List of error types to check for. Supported errors: 'IbzkptError', 'SubSpaceMatrixError', 'NbandsError'. Default: None, in which case only SubSpaceMatrixErrors are checked.

    """
    print("Begin vasp run:")
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
            command = "vasp"
        else:
            command = "mpirun -np {NCPUS} vasp"

    if re.search("\{NCPUS\}",command):
        command = command.format(NCPUS=str(ncpus))

    ### Expand remaining environment variables
    command = os.path.expandvars(command)

    if npar is not None:
        ncore = None

    if npar is not None or ncore is not None:
        io.set_incar_tag({"NPAR":npar, "NCORE":ncore}, jobdir)

    if kpar is not None:
        io.set_incar_tag({"KPAR":kpar}, jobdir)

    print("  jobdir:", jobdir)
    print("  exec:", command)
    sys.stdout.flush()

    sout = open(os.path.join(jobdir,stdout),'w')
    serr = open(os.path.join(jobdir,stderr),'w')
    err = None
    p = subprocess.Popen(command.split(),stdout=sout, stderr=serr)

    # wait for process to end, and periodically check for errors
    poll = p.poll()
    last_check = time.time()
    stopcar_time = None
    while poll  is None:
        time.sleep(poll_check_time)

        if time.time() - last_check > err_check_time:
            last_check = time.time()
            err = error_check(jobdir, os.path.join(jobdir, stdout), err_types)
            if err != None:
                # FreezeErrors are fatal and usually not helped with STOPCAR
                if "FreezeError" in err.keys():
                    print("  VASP is frozen, killing job")
                    sys.stdout.flush()
                    # Sometimes p.kill doesn't work if the process is on multiple nodes
                    os.kill(p.pid, signal.SIGKILL)
                    p.kill()
                    # If the job is re-invoked (e.g. via mpirun or srun) too quickly
                    #   after the previous job ended, infinitiband clusters can have
                    #   some issues with resource allocation. A 30s sleep solves this.
                    time.sleep(30)
                # Other errors can be killed with STOPCAR, which is safer
                elif stopcar_time is None:
                    print("  Found errors:", end=' ')
                    for e in err:
                        print(e, end=' ')
                    print("\n  Killing job with STOPCAR")
                    sys.stdout.flush()
                    io.write_stopcar('e', jobdir)
                    stopcar_time = time.time()
                    time.sleep(30)
                # If the STOPCAR exists, wait 5 min before manually killing the job
                elif time.time() - stopcar_time > 300:
                    print("  VASP is non-responsive, killing job")
                    sys.stdout.flush()
                    os.kill(p.pid, signal.SIGKILL)
                    p.kill()
                    # If the job is re-invoked (e.g. via mpirun or srun) too quickly
                    #   after the previous job ended, infinitiband clusters can have
                    #   some issues with resource allocation. A 30s sleep solves this.
                    time.sleep(30)

        poll = p.poll()

    # close output files
    sout.close()
    serr.close()

    os.chdir(currdir)

    print("Run complete")
    sys.stdout.flush()

    # check finished job for errors
    if err is None:
        # Crash-type errors take priority over any other error that may show up
        err = crash_check(jobdir, os.path.join(jobdir, stdout), err_types)
        if err is None:
            err = error_check(jobdir, os.path.join(jobdir,stdout), err_types)
    if err != None:
        print("  Found errors:", end=' ')
        for e in err:
            print(e, end=' ')
    print("\n")

    return err
