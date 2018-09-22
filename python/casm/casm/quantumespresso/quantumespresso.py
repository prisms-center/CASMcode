from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, shutil, re, subprocess, sys, time, gzip, warnings
from casm.quantumespresso import qeio


class QuantumEspressoError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class QuantumEspressoWarning(Warning):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def continue_job(jobdir, contdir, settings):
# def continue_job(jobdir,\
#                  contdir,\
#                  move = io.DEFAULT_VASP_MOVE_LIST,\
#                  copy = io.DEFAULT_VASP_COPY_LIST,\
#                  remove = io.DEFAULT_VASP_REMOVE_LIST,\
#                  compress = [], backup = []):
    """Use the files in vasp job directory 'jobdir', to setup a vasp job in directory 'contdir'.

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

    print("Continue Quantum Espresso job:\n  Original: " + jobdir + "\n  Continuation: " + contdir)
    sys.stdout.flush()

    infilename=settings["infilename"]
    outfilename=settings["outfilename"]
    # remove duplicates
    move = list(set(settings['move']+settings['extra_input_files']))
    copy = list(set(settings['copy']+[infilename]))
    remove = list(set(settings['remove']))
    compress = list(set(settings['compress']))
    backup = list(set(settings['backup']))

    # Check that necessary files are being moved/copied 
    if not infilename in (move + copy):
        warnings.warn("Warning: Infile not found in either 'move' or 'copy'. Copying " + infilename + " by default", QuantumEspressoWarning)
        copy += [infilename]

    # Check that the user isn't being contradictory or silly
    for f in move:
        if f in remove:
            warnings.warn("Warning: %s found in both 'move' and 'remove'. The file will not be removed." % f, QuantumEspressoWarning)
            remove = list(set(remove) - set([f]))
        if f in copy:
            warnings.warn("Warning: %s found in both 'move' and 'copy'. The file will be copied only." % f, QuantumEspressoWarning)
            move = list(set(move) - set([f]))
        if f in compress:
            raise QuantumEspressoError("Error in casm.quantumespresso.continue_job(). %s found in both 'move' and 'compress', but these options contradict. Did you mean 'backup'?" % f)
    for f in copy:
        if re.match(".*.wfc.*",f) or re.match(".*.igk.*",f):
            warnings.warn("Warning: %s can be rather huge. It is suggested to include %s in 'move', rather than 'copy'. If you're worried about corruption due to interrupted runs, add %s to 'backup' as well." % (f, f, f), QuantumEspressoWarning)
        if f in remove:
            warnings.warn("Warning: %sfound in 'copy' and 'remove', which is the same as 'move'." % f, QuantumEspressoWarning)
            move += f
            copy = list(set(copy) - set([f]))
            remove = list(set(copy) - set([f]))
        if f in compress:
            warnings.warn("Warning: %s found in 'copy' and 'compress'. It is suggested to include %s in 'move' and 'backup' instead." % (f, f), QuantumEspressoWarning)
        if f in backup:
            warnings.warn("Warning: %s found in 'copy' and 'backup'. This will result in three copies of your file being made (old, old_BACKUP.gz, and copy)!" % f, QuantumEspressoWarning)
    for f in remove:
        if f in compress:
            warnings.warn("Warning: %s found in both 'compress' and 'remove'. Defaulting to 'compress' only!" % f, QuantumEspressoWarning)
            remove = list(set(remove) - set([f]))
    for f in compress:
        if f in backup:
            warnings.warn("Warning: %s found in 'compress' and 'backup'. Defaulting to ompressing only!" % f, QuantumEspressoWarning)
            backup = list(set(backup) - set([f]))
    for f in backup:
        if (f not in move) and (f in remove):
            warnings.warn("Warning: 'backup' is meant for files being moved, but %s is found under 'remove' only. Did you maybe mean 'compress'?", QuantumEspressoWarning)
        elif (f not in move):
            warnings.warn("Warning: %s not found in 'move', so you'll just end up with 2 copies of the file..." % f, QuantumEspressoWarning)

    # make the new contdir
    try:
        os.mkdir(contdir)
    except:
        pass

    # make compressed backups of files sensitive to corruption
    print(" backup:")
    for file in backup:
        if os.path.isfile(os.path.join(jobdir,file)):
            print(file, end=' ')
            # Open target file, target file.gz
            f_in = open(os.path.join(jobdir, file), 'rb')
            f_out = gzip.open(os.path.join(jobdir, file)+'_BACKUP.gz', 'wb')
            # Compress, close files
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
    print("")


    # move files
    print("  mv:", end=' ')
    for file in move:
        print(file, end=' ')
        os.rename(os.path.join(jobdir,file),os.path.join(contdir,file))
    print("")

    # copy files
    print("  cp:", end=' ')
    for file in copy:
        print(file, end=' ')
        shutil.copyfile(os.path.join(jobdir,file),os.path.join(contdir,file))
    print("")

    # Rewrite new positions into infile if necessary 
    if os.path.isfile(os.path.join(jobdir,outfilename)) and os.path.isfile(os.path.join(contdir,infilename)):
        postpos=qeio.Poscar(os.path.join(jobdir,outfilename)) # see poscar.py regarding run transfer errors
        newinfile=qeio.Infile(os.path.join(contdir,infilename))
        newinfile.rewrite_poscar_info(postpos) # see infile.py regarding run transfer errors
        newinfile.write(os.path.join(contdir,infilename))
        print("Overwrote " + infilename + " positions and lattice with positions and lattice from " + outfilename)
    else:
        print("  no calculated positions to overwrite")

    # remove files
    print("  rm:", end=' ')
    for file in remove:
        if os.path.isfile(os.path.join(jobdir,file)):
            print(file, end=' ')
            os.remove(os.path.join(jobdir,file))
    print("")

    # compress files
    print(" gzip:", end=' ')
    for file in compress:
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


def complete_job(jobdir, settings):
# def complete_job(jobdir,\
#                  remove = io.DEFAULT_VASP_REMOVE_LIST + ["POTCAR"],\
#                  compress = []):
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

    print("Complete Quantum Espresso job: " + jobdir)
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


class FreezeError(object):
    """Quantum Espresso appears frozen"""
    def __init__(self):
        self.pattern = None


    def __str__(self):
        return "Quantum Espresso appears to have frozen"


    def error(self, outfilename ,line=None, jobdir=None):
        """ Check if QE appears frozen

            Returns true if:
            1) no file has been modified for 5 minutes
            2) 'total cpu time+' exists in Outfile and no output file has been modified
            in 5x the time for the slowest loop
        """

        # Check if any files modified in last 300 s
        most_recent = None
        most_recent_file = None
        for f in os.listdir(jobdir):
            t = time.time() - os.path.getmtime(os.path.join(jobdir,f))
            if most_recent == None:
                most_recent = t
                most_recent_file = f
            elif t < most_recent:
                most_recent = t
                most_recent_file = f

        print("Most recent file output (" + f + "):", most_recent, " seconds ago.")
        sys.stdout.flush()
        if t < 300:
            return False

        myoutfile = qeio.Outfile(os.path.join(jobdir, outfilename))
        if myoutfile.complete:
            print("outfile.complete:", myoutfile.complete)
            sys.stdout.flush()
            return False
        elif myoutfile.slowest_loop != None and most_recent > 5.0*myoutfile.slowest_loop:
            print("slowest_loop:", myoutfile.slowest_loop)
            print("5.0*slowest_loop:", 5.0*myoutfile.slowest_loop)
            print("most_recent:", most_recent)
            sys.stdout.flush()
            return True
        return False


    def fix(self, err_jobdir, new_jobdir, settings):
        """ Fix by killing the job and resubmitting.
        """
        continue_job(err_jobdir, new_jobdir, settings)
        print("  Kill job and try to continue")


class NbandsError(object):
    """
    Your highest band is occupied at some k-points! Unless you are
    """
    def __init__(self):
        self.pattern = "Your highest band is occupied at some k-points! Unless you are"

    def __str__(self):
        return self.pattern

    def error(self, line=None, jobdir=None):
        """ Check if pattern found in line """
        if re.search(self.pattern,line):
            return True
        return False

    def fix(self, err_jobdir, new_jobdir, settings):
        """ Try to fix the error by increasing the number of bands"""
        infilename=settings["infilename"]
        outfilename=settings["outfilename"]
        continue_job(err_jobdir, new_jobdir, settings)
        with open(os.path.join(err_jobdir,outfilename)) as f:
            err_outcar = f.read().splitlines()
        nbands_line = []
        for line in err_outcar:
            if "NBANDS" in line:
                nbands_line.append(line)
        if len(nbands_line)<1:
            print("SERIOUS WARNING :  ")
            print("        Couldn't find any reference to nbands in the Outfile. Continuing without fixing")
        else:
            for l in nbands_line:
                if 'k-points' in l.strip().split():
                    print("  Set NBANDS = "+str(int(1.1 * float(l.strip().split()[-1]))))
                    sys.stdout.flush()
                    qeio.set_infile_tag({"NBANDS":int(1.1 * float(l.strip().split()[-1]))}, infilename,jobdir=new_jobdir)
                    break
        
def error_check(jobdir, stdoutfile, err_types):
    """ Check quantum espresso stdout for errors """
    err = dict()
    if err_types is None:
        possible = []
    else:
        err_objs = {'NbandsError' : NbandsError()}
        for s in err_types:
            if s not in err_objs.keys():
                raise QuantumEspressoError('Invalid err_type: %s'%s)
        possible = [err_objs[s] for s in err_types]

    # Error to check line by line, only look for first of each type
    sout = open(stdoutfile, 'r')
    for line in sout:
        for p in possible:
            if not p.__class__.__name__ in err:
                if p.error(line=line, jobdir=jobdir):
                    err[ p.__class__.__name__] = p

    # Error to check for once
    possible = [FreezeError()]
    for p in possible:
        if p.error(line=None, jobdir=jobdir):
            err[ p.__class__.__name__] = p

    sout.close()
    if len(err) == 0:
        return None
    else:
        return err

def run(infilename, outfilename, jobdir = None, stdout = "std.out", stderr = "std.err", command=None, ncpus=None, poll_check_time = 5.0, err_check_time = 60.0, err_types=None):
    """ Run quantum espresso using subprocess.

        The 'command' is executed in the directory 'jobdir'.

        Args:
            jobdir:     directory to run quantumespresso.  If jobdir == None, the current directory is used.
            stdout:     filename to write to.  If stdout == None, "std.out" is used.
            stderr:     filename to write to.  If stderr == None, "std.err" is used.
            command:    (str or None) qe execution command
                        If command != None: then 'command' is run in a subprocess
                        Else, if ncpus == 1, then command = "pw.x < infilename > outfilename"
                        Else, command = "mpirun -np {NCPUS} pw.x < infilename > outfilename"
            ncpus:      (int) if '{NCPUS}' is in 'command' string, then 'ncpus' is substituted in the command.
                        if ncpus==None, $PBS_NP is used if it exists, else 1
            poll_check_time: how frequently to check if the quantum espresso job is completed
            err_check_time: how frequently to parse quantum espresso output to check for errors
            

    """
    print("Begin quantum espresso run:")
    sys.stdout.flush()

    if jobdir == None:
        jobdir = os.getcwd()

    currdir = os.getcwd()
    os.chdir(jobdir)

    if ncpus == None:
        if "PBS_NP" in os.environ:
            ncpus = os.environ["PBS_NP"]
        else:
            ncpus = 1

    if command == None:
        if ncpus == 1:
            command = "pw.x < {INFILE} > {OUTFILE}"
        else:
            command = "mpirun -np {NCPUS} pw.x < {INFILE} > {OUTFILE}"

    if re.search("\{NCPUS\}",command):
        command = command.format(NCPUS=str(ncpus),INFILE=infilename,OUTFILE=outfilename)

    command = command.format(INFILE=infilename,OUTFILE=outfilename)


    print("  jobdir:", jobdir)
    print("  exec:", command)
    sys.stdout.flush()

    sout = open(os.path.join(jobdir,stdout),'w')
    serr = open(os.path.join(jobdir,stderr),'w')
    err = None

    p = subprocess.Popen(command,stdout=sout, stderr=serr,shell=True)
    # wait for process to end, and periodically check for errors
    poll = p.poll()
    last_check = time.time()
    stopcar_time = None
    while poll  == None:
        time.sleep(poll_check_time)

        if time.time() - last_check > err_check_time:
            last_check = time.time()
            #err = error_check(jobdir, os.path.join(jobdir,stdout), err_typevims)
            if err != None:
                # FreezeErrors are fatal and usually not helped with STOPCAR
                if "FreezeError" in err.keys():
                    print("  Quantum Espresso is frozen, killing job")
                    sys.stdout.flush()
                    p.kill()
        poll = p.poll()

    # close output files
    sout.close()
    serr.close()

    os.chdir(currdir)

    print("Run complete")
    sys.stdout.flush()

    # check finished job for errors
    if err == None:
        #err = error_check(jobdir, os.path.join(jobdir,stdout), err_types)
        if err != None:
            print("  Found errors:", end=' ')
            for e in err:
                print(e, end=' ')
        print("\n")

    return err

