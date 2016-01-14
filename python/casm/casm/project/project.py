import casm
import casm.project
import os, subprocess, pandas
import StringIO


class Project(object):
    """The Project class contains information about a CASM project
    """
    def __init__(self, path=None):
        """
        Construct a CASM Project representation.

        Args:
            path: path to CASM project (Default=None). 
            auto: True if using pbs module's JobDB to manage pbs jobs

        """
        if path == None:
          if casm.project_path(path) == None:
            if path == None:
              raise Exception("No CASM project found using " + os.getcwd())
            else:
              raise Exception("No CASM project found using " + path)
        self.path = casm.project_path(path) 
    
    
    def command(self, args):
        """
        Execute a command via the command line interface. 
        
        Args:
          args: A string containing the command to be executed. Ex: "select --set-on -o /abspath/to/my_selection"
        
        Returns:
          (stdout, stderr): The result of running the command via the command line iterface
        """
        cwd = os.getcwd()
        os.chdir(self.path)
        result = subprocess.Popen(args.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
        os.chdir(cwd)
        return result

    def query(self, columns, selection=None, verbatim=True):
        """Return a pandas DataFrame object containing the output of a 
           'casm query' command.
           
           Args:
             columns: iterable of strings corresponding to 'casm query -k' args
             proj: Project to query (default is CASM project containing the current working directory)
             verbatim: if True, use 'casm query --verbatim' option (default is True)
             selection: a Selection to query (default is "MASTER" selection)
          
           Returns:
             data: a pandas DataFrame containing the query results
        """
        if selection == None:
          selection = casm.project.Selection(self)
        elif not isinstance(selection, casm.project.Selection):
          raise Exception("Error, argument 'selection' must be None or a Selection")
        
        args = "casm query -k "
        for k in columns:
          args += k + " "
        if selection.path != "MASTER":
          args += " -c " + selection.path
        if verbatim == True:
          args += " -v"
        args += " -o STDOUT"
        
        (stdout, stderr) = self.command(args)
        
        return pandas.read_csv(StringIO.StringIO(stdout[1:]), sep=' *') 
    
