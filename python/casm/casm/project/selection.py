import project
import query
import os, subprocess, json
import pandas

class Selection(object):
    """
    A Selection object contains information about a CASM project
    
    Attributes:
      proj: the CASM Project for this selection
      path: absolute path to selection file, or "MASTER" for master list
    
    Properties:
      data: A pandas.DataFrame describing the selected configurations. Has at least
        'configname' and 'selected' (as bool) columns.
      
    """
    def __init__(self, proj=None, path="MASTER"):
        """
        Construct a CASM Project representation.

        Args:
            proj: a Project
            path: path to selection file, or "MASTER" (Default="MASTER") 

        """
        if proj == None:
          proj = project.Project()
        elif not isinstance(proj, project.Project):
          raise Exception("Error constructing Selection: proj argument is not a CASM project")
        self.proj = proj
        
        self.path = path
        if os.path.isfile(path):
          self.path = os.path.abspath(path)

#        if path == "MASTER":
#          self.path = "MASTER"
#        elif os.path.isfile(path):
#          self.path = os.path.abspath(path)
#        else:
#          raise Exception("No file exists named: " + path)
        
        self._data = None
    
    
    @property
    def data(self):
        """
        Get Selection data as a pandas.DataFrame
        
        If the data is modified, 'save' must be called for CASM to use the modified selection.
        """
        if self._data is None:
          if self.path in ["MASTER", "ALL", "CALCULATED"]:
            self._data = query.query(self.proj, ['configname', 'selected'], self)
          elif self._is_json():
            self._data = pandas.read_json(self.path, orient='records')
          else:
            f = open(self.path, 'r')
            f.read(1)
            self._data = pandas.read_csv(f, sep=' *', engine='python')
          self._clean_data()
        return self._data
    
    
    def save(self, data=None, force=False):
        """
        Save the current selection. Also allows completely replacing the 'data'
        describing the selected configurations.
        
        Args:
          data: None (default), or pandas.DataFrame describing the Selection with 
            'configname' and 'selected' columns. If path=="MASTER", Configurations 
            not included in 'data' will be set to not selected.
          force: Boolean, force overwrite existing files
        """
        if self.path == "MASTER":
          
          if data is not None:
            self._data = data
            self._clean_data()
        
          if self._data is None:
            return
          
          clist = self.proj.dir.config_list()
          backup = clist + ".tmp"
          if os.path.exists(backup):
            raise Exception("File: " + backup + " already exists")
          
          # read
          j = json.load(open(clist, 'r'))
          
          for sk, sv in j["supercells"].iteritems():
            for ck, cv in sv.iteritems():
              sv[ck]["selected"] = False
          
          # set selection
          for index, row in self._data.iterrows():
            scelname, configid = row["configname"].split('/')
            j["supercells"][scelname][configid]["selected"] = row["selected"]
            
          # write
          f = open(backup, 'w')
          json.dump(j, f)
          os.rename(backup, clist)
          
        elif self.path in ["ALL", "CALCULATED"]:
          raise Exception("Cannot save the '" + self.path + "' Selection")
        
        else:
          
          if data is not None:
            self._data = data
            self._clean_data()
        
          if os.path.exists(self.path) and not force:
            raise Exception("File: " + self.path + " already exists")
        
          backup = self.path + ".tmp"
          if os.path.exists(backup):
            raise Exception("File: " + backup + " already exists")
          
          if self._is_json():
            self._data.to_json(path_or_buf=backup, orient='records')
          else:
            f = open(backup, 'w')
            f.write('#')
            self._data["selected"] = self._data["selected"].astype(int)
            self._data.to_csv(path_or_buf=f, sep=' ', index=False)
            self._clean_data()
          os.rename(backup, self.path)
        
        
    
    def saveas(self, path,force=False):
        """
        Create a new Selection from this one, save and return it
        
        Args:
          path: path to selection file (Default="MASTER")
          force: Boolean, force overwrite existing files
        
        Returns:
          sel: the new Selection created from this one
        """
        sel = copy.deepcopy(self)
        sel.path = path
        sel.save(force=force)
        return sel
    
    
    def _is_json(self):
        return self.path[-5:].lower() == ".json"
    
    def _clean_data(self):
        self._data['selected'] = self._data['selected'].astype(bool)
        
    
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
        self._data = None


        