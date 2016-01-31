import project
import os, subprocess

class Selection(object):
    """
    A Selection object contains information about a CASM project
    
    Attributes:
      proj: the CASM Project for this selection
      path: absolute path to selection file, or "MASTER" for master list
    
    """
    def __init__(self, proj=None, path="MASTER"):
        """
        Construct a CASM Project representation.

        Args:
            proj: a Project
            path: path to selection file (Default="MASTER"). 

        """
        if proj == None:
          proj = project.Project()
        elif not isinstance(proj, project.Project):
          raise Exception("Error constructing Selection: proj argument is not a CASM project")
        self.proj = proj
        
        if path == "MASTER":
          self.path = "MASTER"
        elif os.path.isfile(path):
          self.path = os.path.abspath(path)
        else:
          raise Exception("No file exists named: " + path)
    
    
    def df(self):
        """
        Read Selection file as a pandas.DataFrame
        """
        if self.path[-5:].lower() == ".json":
          return pandas.read_json(self.path, orient='records')
        else:
          f = open(self.path, 'r')
          f.read(1)
          return pandas.read_csv(f, sep=' *', engine='python')
    
    
    def set_on(self, criteria="", output=None, force=False):
        """
        Perform 'casm select --set-on' using this Selection as input
        
        Args:
            criteria: selection criteria string (Default="").
            output: Name of output file to write result to
        
        """
        self._generic_set("set-on", criteria, output, force)
    
    
    def set_off(self, criteria="", output=None, force=False):
        """
        Perform 'casm select --set-off' using this Selection as input
        
        Args:
            criteria: selection criteria string (Default="").
            output: Name of output file to write result to
        
        """
        self._generic_set("set-off", criteria, output, force)
    
    
    def set(self, criteria="", output=None, force=False):
        """
        Perform 'casm select --set' using this Selection as input
        
        Args:
            criteria: selection criteria string (Default="").
            output: Name of output file to write result to
        
        """
        self._generic_set("set", criteria, output, force)

    
    def _generic_set(self, command, criteria="", output=None, force=False):
        """
        Perform 'casm select --command' using this Selection as input, where
        command='set-on','set-off', or 'set'
        
        Args:
            criteria: selection criteria string (Default="").
            output: Name of output file to write result to
        
        """
        args = "select --" + command + " " + criteria + " -c " + self.path
        if output != None:
          args += " -o " + output
        if force:
          args += " -f"
        self.proj.command(args)
        


        