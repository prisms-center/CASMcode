import unittest
import fixtures
import os
from os.path import join,exists
import casm
import shutil
import json
import pbs
import time

class TestRun(unittest.TestCase):
  
  def setUp(self):
    """
    Read test case data
    """
    self.cases = fixtures.read_vasp_test_cases(".")["run"]
    print ""
  
  def test_run_single(self):
    """
    Test vasp.Relax.run()
    """
    curr=os.getcwd()
    cases = self.cases["run_single"]
    for case in cases:
      # do tests
      
      # load Relax directory
      filesdir = join(fixtures.projects_dir, case["proj"])
      if os.path.isdir(join(filesdir,"relaxdir")):
        shutil.rmtree(join(filesdir,"relaxdir"))
      shutil.copytree(filesdir,join(filesdir,"relaxdir"))
      relaxdir = join(filesdir,"relaxdir")
      ##CREATE POTCAR##
      if "VASP_POTENTIAL_DIR" not in os.environ.keys():
        print "NO SPECIFIED ENVIRONMENT VARIABLE VASP_POTENTIAL_DIR\n"
        print "USING DEFAULT OF $HOME/vasp_potentials\n"
        vasp_pot="$HOME/vasp_potentials"
      else:
        vasp_pot=os.environ["VASP_POTENTIAL_DIR"].encode('ascii','ignore')
      with open(join(relaxdir,"POTCAR"),'wb') as wfd:
        for f in case["potcars"]:
          with open(join(join(vasp_pot,f),"POTCAR"),'rb') as fd:
            shutil.copyfileobj(fd, wfd)
          fd.close()
      wfd.close()

      print "Relax dir is",relaxdir

      db=pbs.JobDB()

      #initialize Relax object and begin relaxation
      cmd = "python -c \"import vasp; relaxation=vasp.Relax('" + relaxdir +"');relaxation.run();\"\n"

      job = pbs.Job(name="casm_unit_test",\
                    account=None,\
                    nodes=1,\
                    ppn=2,\
                    walltime="1:00:00",\
                    pmem=None,\
                    qos=None,\
                    queue="batch",\
                    message=None,\
                    email=None,\
                    priority=0,\
                    command=cmd,\
                    auto=True)
      os.chdir(relaxdir)
      job.submit()
      os.chdir(curr)
      db.update();
      id = db.select_regex_id("rundir",relaxdir)
      j=db.select_job(id[-1])
      while j["jobstatus"]=="Q":
        time.sleep(10)
        db.update();
        j=db.select_job(id[-1])
      while j["jobstatus"]=="R":
        time.sleep(10)
        db.update();
        j=db.select_job(id[-1])
      self.assertTrue(exists(join(relaxdir,"run.final/OUTCAR")))
      if relaxdir==join(filesdir,"relaxdir"):
        shutil.rmtree(relaxdir)
      print "done!"

  def test_setup_many(self):
    """
    Test casm.vaspwrapper.setup() using casm-calc --setup
    """
    cases = self.cases["setup_many"]
    from casm.project import Project,Selection
    import subprocess,shlex,sys,fileinput
    curr_dir=os.getcwd()
    for case in cases:
      # do tests
      # create casm project
      filesdir = join(fixtures.projects_dir, case["proj"])
      if os.path.isdir(join(filesdir,"casmproj_vasp")):
        shutil.rmtree(join(filesdir,"casmproj_vasp"))
      os.mkdir(join(filesdir,"casmproj_vasp"))
      for file in os.listdir(filesdir):
        if os.path.isfile(join(filesdir,file)):
            shutil.copy(join(filesdir,file),join(filesdir,"casmproj_vasp"))
      casmproj = join(filesdir,"casmproj_vasp")
      if "VASP_POTENTIAL_DIR" not in os.environ.keys():
        print "NO SPECIFIED ENVIRONMENT VARIABLE VASP_POTENTIAL_DIR\n"
        print "USING DEFAULT OF $HOME/vasp_potentials\n"
        vasp_pot="$HOME/vasp_potentials"
      else:
        vasp_pot=os.environ["VASP_POTENTIAL_DIR"].encode('ascii','ignore')
      #proj = casm.project.Project(path=join(filesdir, "casmproj"))
      print "casm project dir is", casmproj
      os.chdir(casmproj)
      pinit=subprocess.Popen(shlex.split("casm init"))
      pinit.wait()
      pcomp=subprocess.Popen(shlex.split("casm composition -c"))
      pcomp.wait()
      pcompsel=subprocess.Popen(shlex.split("casm composition -s 0"))
      pcompsel.wait()
      penums=subprocess.Popen(shlex.split("casm enum -s --max " + str(case["max_vol"])))
      penums.wait()
      penumc=subprocess.Popen(shlex.split("casm enum -c -a"))
      penumc.wait()
      psel=subprocess.Popen(shlex.split("casm select --set-on"))
      psel.wait()
      settings_dir = join(casmproj,"training_data/settings/calctype.default")

      #move files to settings directory
      shutil.copyfile("INCAR",join(settings_dir,"INCAR"))
      shutil.copyfile("KPOINTS",join(settings_dir,"KPOINTS"))
      shutil.copyfile("POSCAR",join(settings_dir,"POSCAR"))
      for line in fileinput.input("VSPECIES",inplace=1):
        if "$!%&" in line:
          line = line.replace("$!%&",vasp_pot)
        sys.stdout.write(line)
      shutil.copyfile("VSPECIES",join(settings_dir,"SPECIES"))
      shutil.copyfile("vrelax.json",join(settings_dir,"relax.json"))

      psetup=subprocess.Popen(shlex.split("casm-calc --setup"))
      psetup.wait()
      train_data=join(casmproj,"training_data")
      for dirName in os.listdir(train_data):
        if "SCEL" in dirName and os.path.isdir(join(train_data,dirName)):
          print "Found:",dirName
          for config in os.listdir(join(train_data,dirName)):
            if os.path.isdir(join(join(train_data,dirName),config)):
              print "Checking", config
              self.assertTrue(exists(join(join(join(join(train_data,dirName),config)),"POS")))
              copiedfiles=os.listdir(join(join(join(join(train_data,dirName),config)),"calctype.default"))
              self.assertTrue("INCAR" in copiedfiles)
              self.assertTrue("POSCAR" in copiedfiles)
              self.assertTrue("POTCAR" in copiedfiles)
              self.assertTrue("KPOINTS" in copiedfiles)
      os.chdir(curr_dir)
      if os.path.isdir(join(filesdir,"casmproj_vasp")):
        shutil.rmtree(join(filesdir,"casmproj_vasp"))
      print "done!"   

  def test_run_many(self):
    """
    Test casm.vaspwrapper.run() using casm-calc --submit
    """
    cases = self.cases["run_many"]
    from casm.project import Project,Selection
    import subprocess,shlex,sys,fileinput
    curr_dir=os.getcwd()
    for case in cases:
      # do tests
      # create casm project
      filesdir = join(fixtures.projects_dir, case["proj"])
      if os.path.isdir(join(filesdir,"casmproj_vasp")):
        shutil.rmtree(join(filesdir,"casmproj_vasp"))
      os.mkdir(join(filesdir,"casmproj_vasp"))
      for file in os.listdir(filesdir):
        if os.path.isfile(join(filesdir,file)):
            shutil.copy(join(filesdir,file),join(filesdir,"casmproj_vasp"))
      casmproj = join(filesdir,"casmproj_vasp")
      if "VASP_POTENTIAL_DIR" not in os.environ.keys():
        print "NO SPECIFIED ENVIRONMENT VARIABLE VASP_POTENTIAL_DIR\n"
        print "USING DEFAULT OF $HOME/vasp_potentials\n"
        vasp_pot="$HOME/vasp_potentials"
      else:
        vasp_pot=os.environ["VASP_POTENTIAL_DIR"].encode('ascii','ignore')
      #proj = casm.project.Project(path=join(filesdir, "casmproj"))
      print "casm project dir is", casmproj
      os.chdir(casmproj)
      pinit=subprocess.Popen(shlex.split("casm init"))
      pinit.wait()
      pcomp=subprocess.Popen(shlex.split("casm composition -c"))
      pcomp.wait()
      pcompsel=subprocess.Popen(shlex.split("casm composition -s 0"))
      pcompsel.wait()
      penums=subprocess.Popen(shlex.split("casm enum -s --max " + str(case["max_vol"])))
      penums.wait()
      penumc=subprocess.Popen(shlex.split("casm enum -c -a"))
      penumc.wait()
      psel=subprocess.Popen(shlex.split("casm select --set-on"))
      psel.wait()
      settings_dir = join(casmproj,"training_data/settings/calctype.default")
      #Move files to setting directory
      shutil.copyfile("INCAR",join(settings_dir,"INCAR"))
      shutil.copyfile("KPOINTS",join(settings_dir,"KPOINTS"))
      shutil.copyfile("POSCAR",join(settings_dir,"POSCAR"))
      for line in fileinput.input("VSPECIES",inplace=1):
        if "$!%&" in line:
          line = line.replace("$!%&",vasp_pot)
        sys.stdout.write(line)
      shutil.copyfile("VSPECIES",join(settings_dir,"SPECIES"))
      shutil.copyfile("vrelax.json",join(settings_dir,"relax.json"))
      psubmit=subprocess.Popen(shlex.split("casm-calc --submit"))
      psubmit.wait()
      db=pbs.JobDB()
      db.update();

      train_data=join(casmproj,"training_data")

      for dirName in os.listdir(train_data):
        if "SCEL" in dirName and os.path.isdir(join(train_data,dirName)):
          print "Found:",dirName
          for config in os.listdir(join(train_data,dirName)):
            if os.path.isdir(join(join(train_data,dirName),config)):
              print "Checking", config
              self.assertTrue(exists(join(join(join(join(train_data,dirName),config)),"POS")))
              id = db.select_regex_id("rundir",join(join(join(join(join(train_data,dirName),config)),"calctype.default")))
              j=db.select_job(id[-1])
              while (j["jobstatus"]=="Q" or j["jobstatus"]=="R"):
                time.sleep(10)
                db.update();
                j=db.select_job(id[-1])
              filename=join(join(join(join(join(join(train_data,dirName),config)),"calctype.default")),"status.json")
              file=open(filename)
              status=json.load(file)
              file.close()
              if status["status"]=="complete":
                self.assertTrue(os.path.isdir(join(join(join(join(join(train_data,dirName),config)),"calctype.default"),"run.final")))
                self.assertTrue(exists(join(join(join(join(join(train_data,dirName),config)),"calctype.default"),"properties.calc.json")))
                print "verified complete calculation!"
              else:
                self.assertFalse(exists(join(join(join(join(join(train_data,dirName),config)),"calctype.default"),"properties.calc.json")))
                print dirName + config + " failed and did not produce properties.calc.json as expected"
      os.chdir(curr_dir)
      if os.path.isdir(join(filesdir,"casmproj_vasp")):
        shutil.rmtree(join(filesdir,"casmproj_vasp"))
      print "done!" 

if __name__ == '__main__':
  unittest.main()