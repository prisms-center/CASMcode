import warnings
from casm import project_path, API
import os, subprocess, json
from os.path import join


class ClexDescription(object):
    """
    Settings for a cluster expansion
    
    Attributes
    ----------
      
      name: str
        Cluster expansion name
      
      property: str
        Name of the property being cluster expanded
      
      calctype: str
        Calctype name
      
      ref: str
        Reference state name
      
      bset: str
        Basis set
      
      eci: str
        ECI set name
      
    """
    def __init__(self, name, property, calctype, ref, bset, eci):
      self.name = name
      self.property = property
      self.calctype = calctype
      self.ref = ref
      self.bset = bset
      self.eci = eci


class ProjectSettings(object):
    """
    Settings for a CASM project
    
    Attributes:
      path: Path to CASM project
      data: Dict storing contents of project_settings.json
    
    """
    def __init__(self, path=None):
        """
        Construct a CASM ProjectSettings representation.

        Args:
            path: path to CASM project (Default=None, uses project containing current directory). 

        """
        if project_path(path) is None:
          if path is None:
            raise Exception("No CASM project found using " + os.getcwd())
          else:
            raise Exception("No CASM project found using " + path)
        self.path = project_path(path)
        dir = DirectoryStructure(self.path)
        self.data = json.load(open(dir.project_settings()))

        d = self.data["cluster_expansions"][self.data["default_clex"]]
        self._default_clex = ClexDescription(d["name"], d["property"], d["calctype"], d["ref"], d["bset"], d["eci"])

        self._clex = [
          ClexDescription(d[1]["name"], d[1]["property"], d[1]["calctype"], d[1]["ref"], d[1]["bset"], d[1]["eci"]) 
          for d in self.data["cluster_expansions"].iteritems()
        ]
        
        d = self.data["cluster_expansions"].get("formation_energy", None)
        self._formation_energy_clex = None
        if d is not None:
          self._formation_energy_clex = ClexDescription(d["name"], d["property"], d["calctype"], d["ref"], d["bset"], d["eci"])
    
    # -- Accessors --
    
    @property
    def cluster_expansions(self):
        return self._clex
    
    @property
    def default_clex(self):
        return self._default_clex
    
    @property
    def formation_energy_clex(self):
        return self._formation_energy_clex
    

class DirectoryStructure(object):
    """Standard file and directory locations for a CASM project"""
    def __init__(self, path=None):
        """
        Construct a CASM Project DirectoryStructure representation.

        Args:
            path: path to CASM project (Default=None, uses project containing current directory). 

        """
        if project_path(path) is None:
          if path is None:
            raise Exception("No CASM project found using " + os.getcwd())
          else:
            raise Exception("No CASM project found using " + path)
        self.path = project_path(path)
        self.__casm_dir = ".casm"
        self.__bset_dir = "basis_sets"
        self.__calc_dir = "training_data"
        self.__set_dir = "settings"
        self.__sym_dir = "symmetry"
        self.__clex_dir = "cluster_expansions"
    
    
    # ** Query filesystem **

    def all_bset(self):
      """Check filesystem directory structure and return list of all basis set names"""
      return self.__all_settings("bset", join(self.path, self.__bset_dir))

    def all_calctype(self):
      """Check filesystem directory structure and return list of all calctype names"""
      return self.__all_settings("calctype", join(self.path, self.__calc_dir, self.__set_dir))
    
    def all_ref(self, calctype):
      """Check filesystem directory structure and return list of all ref names for a given calctype"""
      return self.__all_settings("ref", self.calc_settings_dir(calctype))

    def all_clex_name(self):
      """Check filesystem directory structure and return list of all cluster expansion names"""
      return self.__all_settings("clex", join(self.path, self.__clex_dir))

    def all_eci(self, property, calctype, ref, bset):
      """Check filesystem directory structure and return list of all eci names"""
      return self.__all_settings("eci", join(self.path, self.__clex_dir, self.__clex_name(property), self.__calctype(calctype), self.__ref(ref), self.__bset(bset)))


    # ** File and Directory paths **


    # -- Project directory --------

    def root_dir(self):
      """Return casm project directory path"""
      return self.path

    def prim(self):
      """Return prim.json path"""
      return join(self.path, "prim.json")

    # -- Hidden .casm directory --------

    def casm_dir(self):
      """Return hidden .casm dir path"""
      return join(self.path, self.__casm_dir)

    def project_settings(self):
      """Return project_settings.json path"""
      return join(self.casm_dir(), "project_settings.json")

    def scel_list(self, scelname):
      """Return master scel_list.json path"""
      return join(self.casm_dir(), "scel_list.json")

    def config_list(self):
      """Return master config_list.json file path"""
      return join(self.casm_dir(), "config_list.json")


    # -- Symmetry --------

    def symmetry_dir(self):
      """Return symmetry directory path"""
      return join(self.path, self.sym_dir)

    def lattice_point_group(self):
      """Return lattice_point_group.json path"""
      return join(self.symmetry_dir(), "lattice_point_group.json")

    def factor_group(self):
      """Return factor_group.json path"""
      return join(self.symmetry_dir(), "factor_group.json")

    def crystal_point_group(self):
      """Return crystal_point_group.json path"""
      return join(self.symmetry_dir(), "crystal_point_group.json")


    # -- Basis sets --------

    def bset_dir(self, clex):
      """Return path to directory contain basis set info"""
      return join(self.path, self.__bset_dir, self.__bset(clex.bset))

    def bspecs(self, clex):
      """Return basis function specs (bspecs.json) file path"""
      return join(self.bset_dir(clex), "bspecs.json")

    def clust(self, clex):
      """Returns path to the clust.json file"""
      return join(self.bset_dir(clex), "clust.json")
    
    def basis(self, clex):
      """Returns path to the basis.json file"""
      return join(self.bset_dir(clex), "basis.json")

    def clexulator_dir(self, clex):
      """Returns path to directory containing global clexulator"""
      return join(self.bset_dir(clex))

    def clexulator_src(self, project, clex):
      """Returns path to global clexulator source file"""
      return join(self.bset_dir(clex), (project + "_Clexulator.cc"))

    def clexulator_o(self, project, clex):
      """Returns path to global clexulator.o file"""
      return join(self.bset_dir(clex), (project + "_Clexulator.o"))
    
    def clexulator_so(self, project, clex):
      """Returns path to global clexulator so file"""
      return join(self.bset_dir(clex), (project + "_Clexulator.so"))


    # -- Calculations and reference --------

    def settings_path_crawl(self, filename, configname, clex):
        """
        Returns the path to the first file named 'filename' found in the settings 
        directories.
        
        Searches:
          1) self.configuration_calc_settings_dir(configname, clex)
          2) self.supercell_calc_settings_dir(scelname, clex)
          3) self.calc_settings_dir(clex)
          DirectoryStructure.configuration_calc_settings_dir(configname, clex)Crawl casm directory structure starting at configdir and moving upwards 
        
        Returns None if file named 'filename' not found in any of the three directories.
        
        
        Arguments
        ---------
          filename: str
            The name of the file being searched for
          
          configname: str
            The name of the configuration
            
          clex: a casm.project.ClexDescription instance
            Used to specify the calctype to find settings for
          
        
        Returns
        ---------
          filepath: str or None
            The path to the first file named 'filename' found in the settings
            directories, or None if not found.
      
        """
        filepath = join(self.configuration_calc_settings_dir(configname, clex), filename)
        if os.path.exists(filepath):
          return filepath
        
        scelname = configname.split('/')[0]
        filepath = join(self.supercell_calc_settings_dir(scelname, clex), filename)
        if os.path.exists(filepath):
          return filepath
        
        filepath = join(self.calc_settings_dir(clex), filename)
        if os.path.exists(filepath):
          return filepath
        
        return None

    def supercell_dir(self, scelname):
      """Return supercell directory path (scelname has format SCELV_A_B_C_D_E_F)"""
      return join(self.path, self.__calc_dir, scelname)

    def configuration_dir(self, configname):
      """Return configuration directory path (configname has format SCELV_A_B_C_D_E_F/I)"""
      return join(self.path, self.__calc_dir, configname)
    
    def POS(self, configname):
      """Return path to POS file"""
      return join(self.configuration_dir(configname), "POS")
    
    def calctype_dir(self, configname, clex):
      """Return calctype directory path (e.g. training_data/SCEL_...../0/calctype.default"""
      return join(self.configuration_dir(configname),self.__calctype(clex.calctype))

    def calc_settings_dir(self, clex):
      """Return calculation settings directory path, for global settings"""
      return join(self.path, self.__calc_dir, self.__set_dir, self.__calctype(clex.calctype))

    def supercell_calc_settings_dir(self, scelname, clex):
      """Return calculation settings directory path, for supercell specific settings"""
      return join(self.supercell_dir(scelname), self.__set_dir, self.__calctype(clex.calctype))

    def configuration_calc_settings_dir(self, configname, clex):
      """Return calculation settings directory path, for configuration specific settings"""
      return join(self.configuration_dir(configname), self.__set_dir, self.__calctype(clex.calctype))

    def calculated_properties(self, configname, clex):
      """Return calculated properties file path"""
      return join(self.configuration_dir(configname), self.__calctype(clex.calctype), "properties.calc.json")


    def ref_dir(self, clex):
      """Return calculation reference settings directory path, for global settings"""
      return join(self.calc_settings_dir(clex.calctype), self.__ref(clex.ref))

    def composition_axes(self):
      """Return composition axes file path"""
      return join(self.casm_dir(), "composition_axes.json")
    
    def chemical_reference(self, clex):
      """Return chemical reference file path"""
      return join(self.ref_dir(clex), "chemical_reference.json")
    

    # -- Cluster expansions --------

    def property_dir(self, clex):
      """Returns path to eci directory"""
      return join(self.path, self.__clex_dir, self.__clex_name(clex.property))

    def eci_dir(self, clex):
      """
      Returns path to eci directory
      
      Arguments
      ---------
        clex: a casm.project.ClexDescription instance
          Specifies the cluster expansion to get the eci directory for
        
      Returns
      -------
        p: str
          Path to the eci directory
      """
      return join(self.property_dir(clex), self.__calctype(clex.calctype), self.__ref(clex.ref), self.__bset(clex.bset), self.__eci(clex.eci))

    def eci(self, clex):
      """
      Returns path to eci.json
      
      Arguments
      ---------
        clex: a casm.project.ClexDescription instance
          Specifies the cluster expansion to get the eci.json for
        
      Returns
      -------
        p: str
          Path to the eci directory
      """
      return join(self.eci_dir(clex), "eci.json")


    # private:

    def __bset(self, bset):
      return "bset." + bset

    def __calctype(self, calctype):
      return "calctype." + calctype

    def __ref(self, ref):
      return "ref." + ref

    def __clex_name(self, clex_name):
      return "clex." + clex_name

    def __eci(self, eci):
      return "eci." + eci

    
    def __all_settings(self, pattern, location):
      """
      Find all directories at 'location' that match 'pattern.something'
      and return a std::vector of the 'something'
      """
      
      all = [];
      pattern += ".";

      # get all
      if not os.path.exists(location):
        return all
      
      for item in os.listdir(location):
        if os.path.isdir(os.path.join(location,item)) and item[:len(pattern)] == pattern:
          all.append(item[len(pattern):])
      return sorted(all)

    

class Project(object):
    """The Project class contains information about a CASM project
    
    Attributes
    ----------
      
      path: str
        Path to project root directory
      
      settings: casm.project.ProjectSettings instance
        Contains project settings
      
      dir: casm.project.DirectoryStructure instance
        Provides file and directory locations within the project
      
      casm_exe: str
        The casm CLI executable to use when necessary
      
      verbose: bool
        How much to print to stdout
      
    """
    def __init__(self, path=None, casm_exe=None, verbose=True):
      """
      Construct a CASM Project representation.

      Arguments
      ----------
        
        path: str, optional, default=None 
          Path to project root directory. Default=None uses project containing 
          current working directory
        
        casm_exe: str, optional, default=None
          CASM executable to use for command line interface. Default
          uses $CASM if it exists in the environment, else "casm".
        
        verbose: bool, optional, default=True
          How much to print to stdout
      
      """
      
      # will hold a ctypes.c_void_p when loading CASM project into memory
      self._ptr = None
      
      # will keep a casm.API instance
      self._api = None
      
      # set path to this CASM project
      if project_path(path) is None:
        if path is None:
          raise Exception("No CASM project found using " + os.getcwd())
        else:
          raise Exception("No CASM project found using " + path)
      
      # set executable name
      if casm_exe is None:
        if "CASM" in os.environ:
          casm_exe = os.environ["CASM"]
        else:
          casm_exe = "casm"
      
      self.path = project_path(path)
      self.__refresh()
      self.casm_exe = casm_exe
      self.verbose = verbose
      
      
    
    def __del__(self):
      self.__unload()
    
    
    def __load(self):
      """
      Explicitly load CASM project into memory.
      """
      if self._ptr is None:
        self._api = API()
        if self.verbose:
          streamptr = self._api.stdout()
        else:
          streamptr = self._api.nullstream()
        
        if self.verbose:
          errstreamptr = self._api.stderr()
        else:
          errstreamptr = self._api.nullstream()
        
        self._ptr = self._api.primclex_new(self.path, streamptr, streamptr, errstreamptr)
        
    
    def __unload(self):
      """
      Explicitly unload CASM project from memory.
      """
      if self._ptr is not None:
        self._api.primclex_delete(self._ptr)
        self._ptr = None
    
    
    def __refresh(self):
      """
      Reload self.settings and self.dir
      
      Use this after adding or modifying files in the CASM project but no
      special call to refresh PrimClex properties is required
      """
      self.dir = DirectoryStructure(self.path)
      self.settings = ProjectSettings(self.path)
    
    def refresh(self, read_settings=False, read_composition=False, read_chem_ref=False, read_configs=False, clear_clex=False):
      """
      Refresh PrimClex properties to reflect changes to CASM project files.
      """
      if read_settings:
        self.dir = DirectoryStructure(self.path)
        self.settings = ProjectSettings(self.path)
      if self._ptr is not None:
        self._api.primclex_refresh(
          self.data(), 
          read_settings,
          read_composition,
          read_chem_ref,
          read_configs,
          clear_clex)
      
    
    def data(self):
      """
      Returns a 'ctypes.c_void_p' that points to a CASM project. (PrimClex)
      """
      self.__load()
      return self._ptr
    
    
    def command(self, args):
      """
      Execute a command via the c api. 
      
      Args:
        args: A string containing the command to be executed. Ex: "select --set-on -o /abspath/to/my_selection"
      
      Returns:
        (stdout, stderr, returncode): The result of running the command via the command line iterface
      """
      return self.command_via_capi(args)
    
        
    def command_via_cli(self, args):
      """
      Execute a command via the command line interface. 
      
      Args:
        args: A string containing the command to be executed. Ex: "select --set-on -o /abspath/to/my_selection"
      
      Returns:
        (stdout, stderr, returncode): The result of running the command via the command line iterface
      """
      child = subprocess.Popen([self.casm_exe] + args.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE, cwd=self.path)
      result = child.communicate()
      self.__refresh()
      return (result[0], result[1], child.returncode)
    
    
    def command_via_capi(self, args):
      """
      Execute a command via the c api. 
      
      Args:
        args: A string containing the command to be executed. Ex: "select --set-on -o /abspath/to/my_selection"
      
      Returns:
        (stdout, stderr, returncode): The result of running the command via the command line iterface
      """
      # this also ensures self._api is not None
      data = self.data()
      
      # construct stringstream objects to capture stdout, debug, stderr
      ss = self._api.ostringstream_new()
      ss_debug = self._api.ostringstream_new()
      ss_err = self._api.ostringstream_new()
      
      prev = os.getcwd()
      os.chdir(self.path)
      try:
        res = self._api(args, self.data(), ss, ss_debug, ss_err)
      finally:
        os.chdir(prev)
      
      # copy strings and delete stringstreams
      stdout = self._api.ostringstream_to_str(ss)
      self._api.ostringstream_delete(ss)
      
      self._api.ostringstream_delete(ss_debug)
      
      stderr = self._api.ostringstream_to_str(ss_err)
      self._api.ostringstream_delete(ss_err)
      
      self.__refresh()
      return (stdout, stderr, res)
      
