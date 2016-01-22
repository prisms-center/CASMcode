import os, shutil, re, subprocess, sys, time, gzip, warnings
import io


class VaspError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class VaspWarning(Warning):
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

    print "Continue VASP job:\n  Original: " + jobdir + "\n  Continuation: " + contdir
    sys.stdout.flush()

    # remove duplicates
    move = list(set(settings['move']+io.DEFAULT_VASP_MOVE_LIST+settings['extra_input_files']))
    copy = list(set(settings['copy']+io.DEFAULT_VASP_COPY_LIST))
    remove = list(set(settings['remove']))
    compress = list(set(settings['compress']))
    backup = list(set(settings['backup']))

    # Check that necessary files are being moved/copied: INCAR, POTCAR, KPOINTS
    if not "POTCAR" in (move + copy):
        warnings.warn("Warning: POTCAR not found in either 'move' or 'copy'. Moving POTCAR by default", VaspWarning)
        move += ["POTCAR"]
    if not "INCAR" in (move + copy):
        warnings.warn("Warning: INCAR not found in either 'move' or 'copy'. Copying INCAR by default", VaspWarning)
        copy += ["INCAR"]
    if not "KPOINTS" in (move + copy):
        warnings.warn("Warning: KPOINTS not found in either 'move' or 'copy'. Copying KPOINTS by default", VaspWarning)
        copy += ["KPOINTS"]

    # Check that the user isn't being contradictory or silly
    for f in move:
        if f == "POSCAR" or f == "CONTCAR":
            raise VaspError("Error in casm.vasp.general.continue_job().  Do not include POSCAR or CONTCAR in 'move'; use 'backup' if you want a backup")
        if f in remove:
            if f in io.DEFAULT_VASP_MOVE_LIST:
                raise VaspError("Error in casm.vasp.general.continue_job(). %s cannot be removed, VASP will not run!!!" % f)
            else:
                warnings.warn("Warning: %s found in both 'move' and 'remove'. The file will not be removed." % f, VaspWarning)
                remove = list(set(remove) - set([f]))
        if f in copy:
            if f in io.DEFAULT_VASP_MOVE_LIST:
                warnings.warn("Warning: %s found in both 'move' and 'copy'. The fill will be moved only." % f, VaspWarning)
                copy = list(set(copy) - set([f]))
            else:
                warnings.warn("Warning: %s found in both 'move' and 'copy'. The file will be copied only." % f, VaspWarning)
                move = list(set(move) - set([f]))
        if f in compress:
            if f in io.DEFAULT_VASP_MOVE_LIST:
                raise VaspError("Error in casm.vasp.general.continue_job(). %s cannot be compressed, VASP will not run!!!" % f)
            else:
                raise VaspError("Error in casm.vasp.general.continue_job(). %s found in both 'move' and 'compress', but these options contradict. Did you mean 'backup'?" % f)
    for f in copy:
        if f == "POSCAR" or f == "CONTCAR":
            raise VaspError("Error in casm.vasp.general.continue_job().  Do not include POSCAR or CONTCAR in 'copy'; use 'backup' if you want a backup")
        if f == "WAVECAR" or f == "CHGCAR":
            warnings.warn("Warning: %s can be rather huge. It is suggested to include %s in 'move', rather than 'copy'. If you're worried about corruption due to interrupted runs, add %s to 'backup' as well." % (f, f, f), VaspWarning)
        if f in remove:
            warnings.warn("Warning: %sfound in 'copy' and 'remove', which is the same as 'move'." % f, VaspWarning)
            move += f
            copy = list(set(copy) - set([f]))
            remove = list(set(copy) - set([f]))
        if f in compress:
            warnings.warn("Warning: %s found in 'copy' and 'compress'. It is suggested to include %s in 'move' and 'backup' instead." % (f, f), VaspWarning)
        if f in backup:
            warnings.warn("Warning: %s found in 'copy' and 'backup'. This will result in three copies of your file being made (old, old_BACKUP.gz, and copy)!" % f, VaspWarning)
    for f in remove:
        if f in compress:
            warnings.warn("Warning: %s found in both 'compress' and 'remove'. Defaulting to 'compress' only!" % f, VaspWarning)
            remove = list(set(remove) - set([f]))
    for f in compress:
        if f in backup:
            warnings.warn("Warning: %s found in 'compress' and 'backup'. Defaulting to ompressing only!" % f, VaspWarning)
            backup = list(set(backup) - set([f]))
    for f in backup:
        if (f not in move) and (f in remove):
            warnings.warn("Warning: 'backup' is meant for files being moved, but %s is found under 'remove' only. Did you maybe mean 'compress'?", VaspWarning)
        elif (f not in move):
            warnings.warn("Warning: %s not found in 'move', so you'll just end up with 2 copies of the file..." % f, VaspWarning)

    # make the new contdir
    try:
        os.mkdir(contdir)
    except:
        pass

    # make compressed backups of files sensitive to corruption (e.g. WAVECAR)
    print " backup:"
    for file in backup:
        if os.path.isfile(os.path.join(jobdir,file)):
            print file,
            # Open target file, target file.gz
            f_in = open(os.path.join(jobdir, file), 'rb')
            f_out = gzip.open(os.path.join(jobdir, file)+'_BACKUP.gz', 'wb')
            # Compress, close files
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
    print ""

    # copy CONTCAR/POSCAR/etc
    if os.path.isfile(os.path.join(jobdir,"CONTCAR")) and os.path.getsize(os.path.join(jobdir,"CONTCAR")) > 0:
        shutil.copyfile(os.path.join(jobdir,"CONTCAR"),os.path.join(contdir,"POSCAR"))
        print "  cp CONTCAR -> POSCAR"
    else:
        shutil.copyfile(os.path.join(jobdir,"POSCAR"),os.path.join(contdir,"POSCAR"))
        print "  no CONTCAR: cp POSCAR -> POSCAR"

    # move files
    print "  mv:",
    for file in move:
        print file,
        os.rename(os.path.join(jobdir,file),os.path.join(contdir,file))
    print ""

    # copy files
    print "  cp:",
    for file in copy:
        print file,
        shutil.copyfile(os.path.join(jobdir,file),os.path.join(contdir,file))
    print ""

    # remove files
    print "  rm:",
    for file in remove:
        if os.path.isfile(os.path.join(jobdir,file)):
            print file,
            os.remove(os.path.join(jobdir,file))
    print ""

    # compress files
    print " gzip:",
    for file in compress:
        if os.path.isfile(os.path.join(jobdir,file)):
            print file,
            # Open target file, target file.gz
            f_in = open(os.path.join(jobdir, file), 'rb')
            f_out = gzip.open(os.path.join(jobdir, file)+'.gz', 'wb')
            # Compress, close files
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
            # Remove original target file
            os.remove(os.path.join(jobdir,file))
    print ""
    print ""
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

    print "Complete VASP job: " + jobdir
    sys.stdout.flush()

    # remove files
    print "  rm:",
    for f in settings["remove"]:
        if not f in (settings["copy"] + settings["move"] + settings["compress"] + settings["backup"]):
            if os.path.isfile(os.path.join(jobdir,f)):
                print f,
                os.remove(os.path.join(jobdir,f))
    for f in settings["extra_input_files"]:
        if os.path.isfile(os.path.join(jobdir, f)):
            print f,
            os.remove(os.path.join(jobdir,f))
    print ""

    # compress files
    print " gzip:",
    for file in settings["compress"]:
        if os.path.isfile(os.path.join(jobdir,file)):
            print file,
            # Open target file, target file.gz
            f_in = open(os.path.join(jobdir, file), 'rb')
            f_out = gzip.open(os.path.join(jobdir, file)+'.gz', 'wb')
            # Compress, close files
            f_out.writelines(f_in)
            f_out.close()
            f_in.close()
            # Remove original target file
            os.remove(os.path.join(jobdir,file))
    print ""
    print ""
    sys.stdout.flush()


class IbzkptError(object):
    """ VERY BAD NEWS! internal error in subroutine IBZKPT

        This is not necessarily a fatal error.  Also, these fixes are not working.
    """
    def __init__(self):
        self.pattern = "VERY BAD NEWS! internal error in subroutine IBZKPT"
#        self.pattern = "DONOTMATCH"


    def __str__(self):
        return self.pattern


    def error(self, line=None, jobdir=None):
        """ Check if pattern found in line, and KPOINTS subdivisions are odd or not Gamma-centered """
        if re.search(self.pattern, line):
            kpt = io.Kpoints(os.path.join(jobdir,"KPOINTS"))
            odd = True
            for i in range(len(kpt.subdivisions)):
                if kpt.subdivisions[i] % 2 == 0:
                    odd = False
            if odd == False or kpt.automode[0].lower() != "g":
                return True
        return False


    def fix(self, err_jobdir, new_jobdir, settings):
        """ Try to fix by making kpoints subdivisions odd or changing kpoint length"""
        continue_job(err_jobdir, new_jobdir, settings)
        kpt = io.Kpoints(os.path.join(new_jobdir,"KPOINTS"))

        if kpt.automode[0].lower() == "a":
            old_kpt = kpt.subdivisions[0]
            ocr = io.Outcar(os.path.join(err_jobdir,"OUTCAR"))
            max_k = max(ocr.kpts)
            kpt.automode = "GAMMA"
            kpt.subdivisions = [max_k, max_k, max_k]
            kpt.shift = [0.0, 0.0, 0.0]
            print "  Changed KPOINTS from AUTO (%s) to GAMMA (%s)" % (old_kpt, kpt.subdivisions)
#           kpt.subdivisions[0]+=1
#           print "  Increased KPOINTS length from ", kpt.subdivisions[0]-1, " to ", kpt.subdivisions[0]
            kpt.write(os.path.join(new_jobdir,"KPOINTS"))

        else:
            odd = True
            for i in range(len(kpt.subdivisions)):
                if kpt.subdivisions[i] % 2 == 0:
                    odd = False
                    kpt.subdivisions[i] += 1

            if odd == False:
                print "  Set KPOINTS subdivision to be odd:", kpt.subdivisions, "and mode to Gamma"
                kpt.automode = "Gamma"
                kpt.write(os.path.join(new_jobdir,"KPOINTS"))
#        else:
#            symprec = io.get_incar_tag("SYMPREC", jobdir = new_jobdir)
#            if symprec == None or symprec < 1.1e-8:
#                print "  Set SYMPREC = 1e-8"
#                io.set_incar_tag({"SYMPREC": 1e-8}, jobdir = new_jobdir)
#            elif io.get_incar_tag("ISYM", jobdir = new_jobdir) != 0:
#                print "  Set ISYM = 0"
#                io.set_incar_tag({"ISYM": 0}, jobdir = new_jobdir)


class SubSpaceMatrixError(object):
    """WARNING: Sub-Space-Matrix is not hermitian"""
    def __init__(self):
        self.pattern = "WARNING: Sub-Space-Matrix is not hermitian"


    def __str__(self):
        return self.pattern


    def error(self, line=None, jobdir=None):
        """ Check if pattern found in line """
        if re.search(self.pattern, line):
            return True
        return False


    def fix(self, err_jobdir, new_jobdir, settings):
        """ First attempt:
                Set ALGO = VeryFast
                Unset IALGO
            Second attempt:
                Set IBRION = 1 and POTIM = 0.1
            Final attempt:
                Set LREAD = .FALSE.
        """
        continue_job(err_jobdir, new_jobdir, settings)
        if io.get_incar_tag("ALGO", jobdir=new_jobdir) != "VeryFast":
            print "  Set Algo = VeryFast, and Unset IALGO"
            io.set_incar_tag({"IALGO":None, "ALGO":"VeryFast"}, jobdir=new_jobdir)
        elif io.get_incar_tag("IBRION",jobdir=new_jobdir) != 1:
            print "  Set IBRION = 1 and POTIM = 0.1"
            sys.stdout.flush()
            io.set_incar_tag({"IBRION":1, "POTIM":0.1}, jobdir=new_jobdir)
        elif io.get_incar_tag("LREAL", jobdir = new_jobdir) != False:
            print "  Set LREAL = .FALSE."
            sys.stdout.flush()
            io.set_incar_tag({"LREAL":False}, jobdir=new_jobdir)


class FreezeError(object):
    """VASP appears frozen"""
    def __init__(self):
        self.pattern = None


    def __str__(self):
        return "VASP appears to have frozen"


    def error(self, line=None, jobdir=None):
        """ Check if VASP appears frozen

            Returns true if:
            1) no file has been modified for 5 minutes
            2) 'LOOP+' exists in OUTCAR and no output file has been modified
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

        print "Most recent file output (" + f + "):", most_recent, " seconds ago."
        sys.stdout.flush()
        if t < 300:
            return False

        outcar = io.Outcar(os.path.join(jobdir, "OUTCAR"))
        if outcar.complete:
            print "outcar.complete:", outcar.complete
            sys.stdout.flush()
            return False
        elif outcar.slowest_loop != None and most_recent > 5.0*outcar.slowest_loop:
            print "slowest_loop:", outcar.slowest_loop
            print "5.0*slowest_loop:", 5.0*outcar.slowest_loop
            print "most_recent:", most_recent
            sys.stdout.flush()
            return True
        return False


    def fix(self, err_jobdir, new_jobdir, settings):
        """ Fix by killing the job and resubmitting.
        """
        continue_job(err_jobdir, new_jobdir, settings)
        print "  Kill job and try to continue"


def error_check(jobdir, stdoutfile, err_types):
    """ Check vasp stdout for errors """
    err = dict()
    if err_types is None:
        possible = [SubSpaceMatrixError()]
    else:
        err_objs = {'IbzkptError' : IbzkptError(), 'SubSpaceMatrixError' : SubSpaceMatrixError()}
        for s in err_types:
            if s not in err_objs.keys():
                raise VaspError('Invalid err_type: %s'%s)
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


def run(jobdir = None, stdout = "std.out", stderr = "std.err", npar=None, ncore=None, command=None, ncpus=None, kpar=None, poll_check_time = 5.0, err_check_time = 60.0, err_types=None):
    """ Run vasp using subprocess.

        The 'command' is executed in the directory 'jobdir'.

        Args:
            jobdir:     directory to run vasp.  If jobdir == None, the current directory is used.
            stdout:     filename to write to.  If stdout == None, "std.out" is used.
            stderr:     filename to write to.  If stderr == None, "std.err" is used.
            npar:       (int or None) VASP INCAR NPAR setting. If npar == None, then NPAR is removed from INCAR
            kpar:       (int or None) VASP INCAR KPAR setting. If kpar == None, then KPAR is removed from INCAR
            ncore:      (int or None) VASP INCAR NCORE setting. If not npar == None or ncore == None, then NCORE is removed from INCAR
            command:    (str or None) vasp execution command
                        If command != None: then 'command' is run in a subprocess
                        Else, if ncpus == 1, then command = "vasp"
                        Else, command = "mpirun -np {NCPUS} vasp"
            ncpus:      (int) if '{NCPUS}' is in 'command' string, then 'ncpus' is substituted in the command.
                        if ncpus==None, $PBS_NP is used if it exists, else 1
            poll_check_time: how frequently to check if the vasp job is completed
            err_check_time: how frequently to parse vasp output to check for errors
            err_types:  List of error types to check for. Supported errors: 'IbzkptError', 'SubSpaceMatrixError'. Default: None, in which case only SubSpaceMatrixErrors are checked.

    """
    print "Begin vasp run:"
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
            command = "vasp"
        else:
            command = "mpirun -np {NCPUS} vasp"

    if re.search("\{NCPUS\}",command):
        command = command.format(NCPUS=str(ncpus))

    if not npar == None:
        ncore = None

    io.set_incar_tag({"NPAR":npar, "NCORE":ncore, "KPAR":kpar}, jobdir)

    print "  jobdir:", jobdir
    print "  exec:", command
    sys.stdout.flush()

    sout = open(os.path.join(jobdir,stdout),'w')
    serr = open(os.path.join(jobdir,stderr),'w')
    err = None
    p = subprocess.Popen(command.split(),stdout=sout, stderr=serr)

    # wait for process to end, and periodically check for errors
    poll = p.poll()
    last_check = time.time()
    stopcar_time = None
    while poll  == None:
        time.sleep(poll_check_time)

        if time.time() - last_check > err_check_time:
            last_check = time.time()
            err = error_check(jobdir, os.path.join(jobdir,stdout), err_types)
            if err != None:
                # FreezeErrors are fatal and usually not helped with STOPCAR
                if "FreezeError" in err.keys():
                    print "  VASP is frozen, killing job"
                    sys.stdout.flush()
                    p.kill()
                # Other errors can be killed with STOPCAR, which is safer
                elif stopcar_time == None:
                    print "  Found errors:",
                    for e in err:
                        print e,
                    print "\n  Killing job with STOPCAR"
                    sys.stdout.flush()
                    io.write_stopcar('e', jobdir)
                    stopcar_time = time.time()
                # If the STOPCAR exists, wait 5 min before manually killing the job
                elif time.time() - stopcar_time > 300:
                    print "  VASP is non-responsive, killing job"
                    sys.stdout.flush()
                    p.kill()

        poll = p.poll()

    # close output files
    sout.close()
    serr.close()

    os.chdir(currdir)

    print "Run complete"
    sys.stdout.flush()

    # check finished job for errors
    if err == None:
        err = error_check(jobdir, os.path.join(jobdir,stdout), err_types)
        if err != None:
            print "  Found errors:",
            for e in err:
                print e,
        print "\n"

    return err

