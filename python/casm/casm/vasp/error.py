""" Job manipuation routines for VASP"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import abc
import os
import shutil
import re
import sys
import time

from casm.vasp import io

class VaspError(Exception):
    """ VASP related errors """
    pass

class VaspWarning(Warning):
    """ VASP related warnings"""
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg

def continue_job(jobdir, contdir, settings):
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

    print("Continue VASP job:\n  Original: " + jobdir + "\n  Continuation: " + contdir)
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

    # copy CONTCAR/POSCAR/etc
    if os.path.isfile(os.path.join(jobdir, "CONTCAR")) and os.path.getsize(os.path.join(jobdir, "CONTCAR")) > 0:
        shutil.copyfile(os.path.join(jobdir, "CONTCAR"), os.path.join(contdir, "POSCAR"))
        print("  cp CONTCAR -> POSCAR")
    else:
        shutil.copyfile(os.path.join(jobdir, "POSCAR"), os.path.join(contdir, "POSCAR"))
        print("  no CONTCAR: cp POSCAR -> POSCAR")

    # move files
    print("  mv:", end=' ')
    for file in move:
        print(file, end=' ')
        # This prevents a missing file from crashing the vasprun (e.g. a missing WAVECAR)
        if os.path.isfile(os.path.join(jobdir, file)):
            os.rename(os.path.join(jobdir, file), os.path.join(contdir, file))
        else:
            print("Could not find file %s, skipping!" % file)
    print("")

    # copy files
    print("  cp:", end=' ')
    for file in copy:
        print(file, end=' ')
        # This prevents a missing file from crashing the vasprun (e.g. a missing WAVECAR)
        if os.path.isfile(os.path.join(jobdir, file)):
            shutil.copyfile(os.path.join(jobdir, file), os.path.join(contdir, file))
        else:
            print("Could not find file %s, skipping!" % file)
    print("")

    # remove files
    print("  rm:", end=' ')
    for file in remove:
        if os.path.isfile(os.path.join(jobdir, file)):
            print(file, end=' ')
            os.remove(os.path.join(jobdir, file))
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


class _RunError(object):
    """ Template class for all run errors to be caught """
    ___metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def error(self, line=None, jobdir=None):
        """ Check for the error in the OUTCAR """

    @abc.abstractmethod
    def fix(self, err_jobdir, new_jobdir, settings):
        """ Fix the error """

    @abc.abstractproperty
    def pattern(self):
        """ The pattern the Class matches to """

    def __str__(self):
        return self.pattern

class _CrashError(object):
    """ Template class for all crash-related errors to be caught """
    ___metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def error(self, line=None, jobdir=None):
        """ Check for the error in the OUTCAR """

    @abc.abstractmethod
    def fix(self, err_jobdir, new_jobdir, settings):
        """ Fix the error """

    @abc.abstractproperty
    def pattern(self):
        """ The pattern the Class matches to """

    def __str__(self):
        return self.pattern

class _FreezeError(object):
    """ Template class for all crash-related errors to be caught """
    ___metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def error(self, line=None, jobdir=None):
        """ Check for the error in the OUTCAR """

    @abc.abstractmethod
    def fix(self, err_jobdir, new_jobdir, settings):
        """ Fix the error """

    @abc.abstractproperty
    def pattern(self):
        """ The pattern the Class matches to """

    def __str__(self):
        return self.pattern

class IbzkptError(_RunError):
    """ VERY BAD NEWS! internal error in subroutine IBZKPT

        This is not necessarily a fatal error.  Also, these fixes are not working.
    """
    # This should be a class member: all IbzkptErrors have it
    pattern = "VERY BAD NEWS! internal error in subroutine IBZKPT"

    def __str__(self):
        return self.pattern



    def error(self, line=None, jobdir=None):
        """ Check if pattern found in line, and KPOINTS subdivisions are odd or not Gamma"""
        if re.search(self.pattern, line):
            kpt = io.Kpoints(os.path.join(jobdir, "KPOINTS"))
            odd = True
            for i in range(len(kpt.subdivisions)):
                if kpt.subdivisions[i] % 2 == 0:
                    odd = False
            if odd is False or kpt.automode[0].lower() != "g":
                return True
        return False

    def fix(self, err_jobdir, new_jobdir, settings):
        """ Try to fix by making kpoints subdivisions odd or changing kpoint length"""
        continue_job(err_jobdir, new_jobdir, settings)
        kpt = io.Kpoints(os.path.join(new_jobdir, "KPOINTS"))

        if kpt.automode[0].lower() == "a":
            old_kpt = kpt.subdivisions[0]
            ocr = io.Outcar(os.path.join(err_jobdir, "OUTCAR"))
            max_k = max(ocr.kpts)
            kpt.automode = "GAMMA"
            kpt.subdivisions = [max_k, max_k, max_k]
            kpt.shift = [0.0, 0.0, 0.0]
            print("  Changed KPOINTS from AUTO (%s) to GAMMA (%s)" % (old_kpt, kpt.subdivisions))
            kpt.write(os.path.join(new_jobdir, "KPOINTS"))

        else:
            odd = True
            for i in range(len(kpt.subdivisions)):
                if kpt.subdivisions[i] % 2 == 0:
                    odd = False
                    kpt.subdivisions[i] += 1

            if odd == False:
                print("  Set KPOINTS subdivision to be odd:", kpt.subdivisions, "and mode to Gamma")
                kpt.automode = "Gamma"
                kpt.write(os.path.join(new_jobdir, "KPOINTS"))
#        else:
#            symprec = io.get_incar_tag("SYMPREC", jobdir = new_jobdir)
#            if symprec is None or symprec < 1.1e-8:
#                print "  Set SYMPREC = 1e-8"
#                io.set_incar_tag({"SYMPREC": 1e-8}, jobdir = new_jobdir)
#            elif io.get_incar_tag("ISYM", jobdir = new_jobdir) != 0:
#                print "  Set ISYM = 0"
#                io.set_incar_tag({"ISYM": 0}, jobdir = new_jobdir)

class FEXCFError(_RunError):
    """ERROR FEXCF: supplied exchange-correlation table"""

    pattern = "ERROR FEXCF: supplied exchange-correlation table"

    def __str__(self):
        return self.pattern

    def error(self, line=None, jobdir=None):
        """ Check if pattern found in line """
        if re.search(self.pattern, line):
            return True
        return False

    def fix(self, err_jobdir, new_jobdir, settings):
        """ First attempt:
                Set IBRION = 2
            Second attempt:
                Reduce POTIM to 0.1
            Final attempt:
                Reduce POTIM to 0.01
        """
        continue_job(err_jobdir, new_jobdir, settings)
        if io.get_incar_tag("IBRION", jobdir=new_jobdir) != 2:
            print("  Set IBRION = 2")
            io.set_incar_tag({"IBRION":2}, jobdir=new_jobdir)
        elif io.get_incar_tag("POTIM", jobdir=new_jobdir) > 0.1:
            print("  Set POTIM = 0.1")
            sys.stdout.flush()
            io.set_incar_tag({"POTIM":0.1}, jobdir=new_jobdir)
        elif io.get_incar_tag("POTIM", jobdir=new_jobdir) > 0.01:
            print("  Set POTIM = 0.01")
            sys.stdout.flush()
            io.set_incar_tag({"POTIM":0.01}, jobdir=new_jobdir)


class SubSpaceMatrixError(_RunError):
    """WARNING: Sub-Space-Matrix is not hermitian"""

    pattern = "WARNING: Sub-Space-Matrix is not hermitian"

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
                Set LREAL = .FALSE.
        """
        continue_job(err_jobdir, new_jobdir, settings)
        if io.get_incar_tag("ALGO", jobdir=new_jobdir) != "VeryFast":
            print("  Set Algo = VeryFast, and Unset IALGO")
            io.set_incar_tag({"IALGO":None, "ALGO":"VeryFast"}, jobdir=new_jobdir)
        elif io.get_incar_tag("IBRION", jobdir=new_jobdir) != 1:
            print("  Set IBRION = 1 and POTIM = 0.1")
            sys.stdout.flush()
            io.set_incar_tag({"IBRION":1, "POTIM":0.1}, jobdir=new_jobdir)
        elif io.get_incar_tag("LREAL", jobdir=new_jobdir) != "False":
            print("  Set LREAL = .FALSE.")
            sys.stdout.flush()
            io.set_incar_tag({"LREAL":False}, jobdir=new_jobdir)
        else:
            print("  Set Algo = Normal, and Unset IALGO")
            sys.stdout.flush()
            io.set_incar_tag({"IALGO":None, "ALGO":"Normal"}, jobdir=new_jobdir)

class InisymError(_CrashError):
    """ VASP can't figure out the (magnetic) symmetry """
    pattern = " INISYM: ERROR: Unable to resolve symmetry "

    def fix(self, err_jobdir, new_jobdir, settings):
        """ Up symprec, or turn off symmetry"""
        continue_job(err_jobdir, new_jobdir, settings)
        symprec = io.get_incar_tag("SYMPREC", jobdir = new_jobdir)
        if symprec is None or symprec > 1.1e-8:
            print("  Set SYMPREC = 1e-8")
            io.set_incar_tag({"SYMPREC": 1e-8}, jobdir = new_jobdir)
        elif io.get_incar_tag("ISYM", jobdir = new_jobdir) != 0:
            print("  Set ISYM = 0")
            io.set_incar_tag({"ISYM": 0}, jobdir = new_jobdir)

class SgrconError(_CrashError):
    """ VASP is having yet another symmetry problem """
    pattern = "VERY BAD NEWS! internal error in subroutine SGRCON:"

    def fix(self, err_jobdir, new_jobdir, settings):
        """ Up symprec, or turn off symmetry"""
        continue_job(err_jobdir, new_jobdir, settings)
        symprec = io.get_incar_tag("SYMPREC", jobdir = new_jobdir)
        if symprec is None or symprec > 1.1e-8:
            print("  Set SYMPREC = 1e-8")
            io.set_incar_tag({"SYMPREC": 1e-8}, jobdir = new_jobdir)
        elif io.get_incar_tag("ISYM", jobdir = new_jobdir) != 0:
            print("  Set ISYM = 0")
            io.set_incar_tag({"ISYM": 0}, jobdir = new_jobdir)

class WavecarError(_CrashError):
    """ A bad WAVECAR is causing the job to crash/abort """
    pattern = "  ERROR: while reading WAVECAR, plane wave coefficients changed"

    def error(self, line=None, jobdir=None):
        """ Check if pattern found in line """
        if re.search(self.pattern, line):
            return True
        return False

    def fix(self, err_jobdir, new_jobdir, settings):
        """ Delete WAVECAR and retry """
        continue_job(err_jobdir, new_jobdir, settings)
        if os.path.isfile(os.path.join(new_jobdir, "WAVECAR")):
            os.remove(os.path.join(new_jobdir, "WAECAR"))
        with open(os.path.join(new_jobdir, "WAVECAR"), 'a') as f:
            pass

class NbandsError(_RunError):
    """Your highest band is occupied at some k-points! Unless you are"""
    pattern = "Your highest band is occupied at some k-points! Unless you are"

    def __str__(self):
        return self.pattern

    def error(self, line=None, jobdir=None):
        """ Check if pattern found in line """
        if re.search(self.pattern,line):
            return True
        return False

    def fix(self, err_jobdir, new_jobdir, settings):
        """ Try to fix the error by increasing the number of bands"""
        continue_job(err_jobdir, new_jobdir, settings)
        with open(os.path.join(err_jobdir,'OUTCAR')) as f:
            err_outcar = f.read().splitlines()
        nbands_line = []
        for line in err_outcar:
            if "NBANDS" in line:
                nbands_line.append(line)
        if len(nbands_line)<1:
            print("SERIOUS WARNING :  ")
            print("        Couldn't find any reference to nbands in the OUTCAR. Continuing without fixing")
        else:
            for l in nbands_line:
                if 'k-points' in l.strip().split():
                    print("  Set NBANDS = "+str(int(1.1 * float(l.strip().split()[-1]))))
                    sys.stdout.flush()
                    io.set_incar_tag({"NBANDS":int(1.1 * float(l.strip().split()[-1]))}, jobdir=new_jobdir)
                    break

class NoConvergeError(_RunError):
    """An ionic step was performed without a fully-converged density"""
    pattern = r"[A-Za-z]+:\s+%i\s+"

    def __str__(self):
        return "VASP ran out of electronic steps before convergence was achieved!"

    def error(self, line=None, jobdir=None):
        """ Check if pattern found in line """
        # I don't like having to open the INCAR every time...
        nelm = io.get_incar_tag("NELM", jobdir=jobdir)
        if nelm is None:
            nelm = 40
        ediff = io.get_incar_tag("EDIFF", jobdir=jobdir)
        if ediff is None:
            ediff = 1e-4
        # We know the SCF ended at exactly NELM steps
        if re.search(self.pattern % nelm, line):
            # This may not be a SCF line at all, so pass to-float errors
            try:
                # Determining if the SCF came *close* to converging
                if ediff*10. <= float(line.split()[3]):
                    return True
                return False
            except:
                return False

    def fix(self, err_jobdir, new_jobdir, settings):
        """ Try to fix the error by changing the algo"""
        continue_job(err_jobdir, new_jobdir, settings)
        # Replace the potentially bad POSCAR
        shutil.copyfile(os.path.join(err_jobdir, "POSCAR"), os.path.join(new_jobdir, "POSCAR"))
        # First, see if a change of ALGO helps
        curr_algo = io.get_incar_tag("ALGO", new_jobdir).upper()
        if curr_algo == 'FAST':
            io.set_incar_tag({"ALGO":"Normal", "IALGO":None}, jobdir=new_jobdir)
            print("  Set ALGO = Normal")
        elif curr_algo == 'NORMAL':
            io.set_incar_tag({"ALGO":"All", "IALGO":None}, jobdir=new_jobdir)
            print("  Set ALGO = All")
        elif curr_algo == 'ALL':
            io.set_incar_tag({"ALGO":"Damped", "IALGO":None}, jobdir=new_jobdir)
            io.set_incar_tag({"TIME":"0.4"}, jobdir=new_jobdir)
            print("  Set ALGO = Damped, TIME = 0.4")

class FreezeError(_FreezeError):
    """VASP appears frozen"""
    pattern = None

    def __str__(self):
        return "VASP appears to have frozen"

    @staticmethod
    def error(line=None, jobdir=None):  #pylint: disable=unused-argument
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
            t = time.time() - os.path.getmtime(os.path.join(jobdir, f))
            if most_recent is None:
                most_recent = t
                most_recent_file = f
            elif t < most_recent:
                most_recent = t
                most_recent_file = f

        print("Most recent file output (" + most_recent_file + "):", most_recent, " seconds ago.")
        sys.stdout.flush()
        if t < 300:
            return False

        outcar = io.Outcar(os.path.join(jobdir, "OUTCAR"))
        if outcar.complete:
            print("outcar.complete:", outcar.complete)
            sys.stdout.flush()
            return False
        elif outcar.slowest_loop != None and most_recent > 5.0*outcar.slowest_loop:
            print("slowest_loop:", outcar.slowest_loop)
            print("5.0*slowest_loop:", 5.0*outcar.slowest_loop)
            print("most_recent:", most_recent)
            sys.stdout.flush()
            return True
        return False

    @staticmethod
    def fix(err_jobdir, new_jobdir, settings):
        """ Fix by killing the job and resubmitting.
        """
        continue_job(err_jobdir, new_jobdir, settings)
        print("  Kill job and try to continue")

def error_check(jobdir, stdoutfile, err_types):
    """ Check vasp stdout for errors """
    err = dict()
    err_objs = {}
    for i_err in _RunError.__subclasses__():
        err_objs[i_err.__name__] = i_err()
    if err_types is None:
        possible = [SubSpaceMatrixError()]
    else:
        # err_objs = {'IbzkptError' : IbzkptError(), 'SubSpaceMatrixError' : SubSpaceMatrixError(), 'NbandsError' : NbandsError()}
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
    possible = [i_err() for i_err in _FreezeError.__subclasses__()]
    for p in possible:
        if p.error(line=None, jobdir=jobdir):
            err[ p.__class__.__name__] = p

    sout.close()
    if len(err) == 0:
        return None
    else:
        return err

def crash_check(jobdir, stdoutfile, crash_types):
    """ Check vasp stdout for evidence of a crash """
    err = None
    possible = [i_err() for i_err in _CrashError.__subclasses__()]
    ocar = io.Outcar(os.path.join(jobdir, "OUTCAR"))
    if ocar.complete:
        return None
    # Error to check line by line, only track most-recent error
    sout = open(stdoutfile, 'r')
    for line in sout:
        for p in possible:
            if p.error(line=line, jobdir=jobdir):
                err = p
    sout.close()
    if err is None:
        return None
    else:
        return {err.__class__.__name__: err}

