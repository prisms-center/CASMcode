import casm
import os, subprocess, pandas, json
import StringIO


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
        if path == None:
          if casm.project_path(path) == None:
            if path == None:
              raise Exception("No CASM project found using " + os.getcwd())
            else:
              raise Exception("No CASM project found using " + path)
        self.path = casm.project_path(path)
        dir = DirectoryStructure(self.path)
        self.data = json.load(open(dir.project_settings()))
        
    
    # -- Accessors --
    
    def bset(self):
        return self.data["curr_bset"]
    
    def calctype(self):
        return self.data["curr_calctype"]
    
    def eci(self):
        return self.data["curr_eci"]
    
    def clex(self):
        return self.data["curr_clex"]
    
    def ref(self):
        return self.data["curr_ref"]
    
    def properties(self):
        return self.data["curr_properties"]
    
    def name(self):
        return self.data["name"]
    
    
    # -- Mutators --
    
    def set_bset(self, bset):
        self.data["curr_bset"] = bset
    
    def set_calctype(self, calctype):
        self.data["curr_calctype"] = calctype
    
    def set_eci(self, eci):
        self.data["curr_eci"] = eci
    
    def set_clex(self, clex):
        self.data["curr_clex"] = clex
    
    def set_ref(self, ref):
        self.data["curr_ref"] = ref
    
    def set_properties(self, properties):
        self.data["curr_properties"] = list(properties)
    
    
    def commit(self):
        dir = DirectoryStructure(self.path)
        f = open(dir.project_settings(),'w')
        json.dump(self.data, f)
        f.close()
    

class DirectoryStructure(object):
    """Standard file and directory locations for a CASM project"""
    def __init__(self, path=None):
        """
        Construct a CASM Project DirectoryStructure representation.

        Args:
            path: path to CASM project (Default=None, uses project containing current directory). 

        """
        if path == None:
          if casm.project_path(path) == None:
            if path == None:
              raise Exception("No CASM project found using " + os.getcwd())
            else:
              raise Exception("No CASM project found using " + path)
        self.path = casm.project_path(path)
        self.__casm_dir = ".casm"
        self.__bset_dir = "basis_sets"
        self.__calc_dir = "training_data"
        self.__set_dir = "settings"
        self.__sym_dir = "symmetry"
        self.__clex_dir = "cluster_expansions"
    
    
    # ** Query filesystem **

    def all_bset(self):
      """Check filesystem directory structure and return list of all basis set names"""
      return __all_settings("bset", os.path.join(self.path, self.__bset_dir))

    def all_calctype(self):
      """Check filesystem directory structure and return list of all calctype names"""
      return __all_settings("calctype", os.path.join(self.path, self.__calc_dir, self.__set_dir))
    
    def all_ref(self, calctype):
      """Check filesystem directory structure and return list of all ref names for a given calctype"""
      return __all_settings("ref", calc_settings_dir(calctype))

    def all_clex(self):
      """Check filesystem directory structure and return list of all cluster expansion names"""
      return __all_settings("clex", os.path.join(self.path, self.__clex_dir))

    def all_eci(self, clex, calctype, ref, bset):
      """Check filesystem directory structure and return list of all eci names"""
      return __all_settings("eci", os.path.join(self.path, self.__clex_dir, __clex(clex), __calctype(calctype), __ref(ref), __bset(bset)))


    # ** File and Directory paths **


    # -- Project directory --------

    def root_dir(self):
      """Return casm project directory path"""
      return self.path

    def prim(self):
      """Return prim.json path"""
      return os.path.join(self.path, "prim.json")

    # -- Hidden .casm directory --------

    def casm_dir(self):
      """Return hidden .casm dir path"""
      return os.path.join(self.path, self.__casm_dir)

    def project_settings(self):
      """Return project_settings.json path"""
      return os.path.join(self.casm_dir(), "project_settings.json")

    def scel_list(self, scelname):
      """Return master scel_list.json path"""
      return os.path.join(self.casm_dir(), "scel_list.json")

    def config_list(self):
      """Return master config_list.json file path"""
      return os.path.join(self.casm_dir(), "config_list.json")


    # -- Symmetry --------

    def symmetry_dir(self):
      """Return symmetry directory path"""
      return os.path.join(self.path, self.sym_dir)

    def lattice_point_group(self):
      """Return lattice_point_group.json path"""
      return os.path.join(self.symmetry_dir(), "lattice_point_group.json")

    def factor_group(self):
      """Return factor_group.json path"""
      return os.path.join(self.symmetry_dir(), "factor_group.json")

    def crystal_point_group(self):
      """Return crystal_point_group.json path"""
      return os.path.join(self.symmetry_dir(), "crystal_point_group.json")


    # -- Basis sets --------

    def bset_dir(self, bset):
      """Return path to directory contain basis set info"""
      return os.path.join(self.path, self.__bset_dir, self.__bset(bset))

    def bspecs(self, bset):
      """Return basis function specs (bspecs.json) file path"""
      return os.path.join(self.bset_dir(bset), "bspecs.json")

    def clust(self, bset):
      """Returns path to the clust.json file"""
      return os.path.join(self.bset_dir(bset), "clust.json")
    
    def basis(self, bset):
      """Returns path to the basis.json file"""
      return os.path.join(self.bset_dir(bset), "basis.json")

    def clexulator_dir(self, bset):
      """Returns path to directory containing global clexulator"""
      return os.path.join(self.bset_dir(bset))

    def clexulator_src(self, project, bset):
      """Returns path to global clexulator source file"""
      return os.path.join(self.bset_dir(bset), (project + "_Clexulator.cc"))

    def clexulator_o(self, project, bset):
      """Returns path to global clexulator.o file"""
      return os.path.join(self.bset_dir(bset), (project + "_Clexulator.o"))
    
    def clexulator_so(self, project, bset):
      """Returns path to global clexulator so file"""
      return os.path.join(self.bset_dir(bset), (project + "_Clexulator.so"))


    # -- Calculations and reference --------

    def supercell_dir(self, scelname):
      """Return supercell directory path (scelname has format SCELV_A_B_C_D_E_F)"""
      return os.path.join(self.path, self.__calc_dir, scelname)

    def configuration_dir(self, configname):
      """Return configuration directory path (configname has format SCELV_A_B_C_D_E_F/I)"""
      return os.path.join(self.path, self.__calc_dir, configname)


    def calc_settings_dir(self, calctype):
      """Return calculation settings directory path, for global settings"""
      return os.path.join(self.path, self.__calc_dir, self.__set_dir, self.__calctype(calctype))

    def supercell_calc_settings_dir(self, scelname, calctype):
      """Return calculation settings directory path, for supercell specific settings"""
      return os.path.join(self.supercell_dir(scelname), self.__set_dir, self.__calctype(calctype))

    def configuration_calc_settings_dir(self, configname, calctype):
      """Return calculation settings directory path, for configuration specific settings"""
      return os.path.join(self.configuration_dir(configname), self.__set_dir, self.__calctype(calctype))

    def calculated_properties(self, configname, calctype):
      """Return calculated properties file path"""
      return os.path.join(self.configuration_dir(configname), self.__calctype(calctype), "properties.calc.json")


    def ref_dir(self, calctype, ref):
      """Return calculation reference settings directory path, for global settings"""
      return os.path.join(self.calc_settings_dir(calctype), self.__ref(ref))

    def composition_axes(self, calctype, ref):
      """Return composition axes file path"""
      return os.path.join(self.ref_dir(calctype, ref), "composition_axes.json")
    
    def chemical_reference(self, calctype, ref):
      """Return chemical reference file path"""
      return os.path.join(self.ref_dir(calctype, ref), "chemical_reference.json")
    

    # -- Cluster expansions --------

    def clex_dir(self, clex):
      """Returns path to eci directory"""
      return os.path.join(self.path, self.__clex_dir, self.__clex(clex))

    def eci_dir(self, clex, calctype, ref, bset, eci):
      """Returns path to eci directory"""
      return os.path.join(self.clex_dir(clex), self.__calctype(calctype), self.__ref(ref), self.__bset(bset), self.__eci(eci))

    def eci(self, clex, calctype, ref, bset, eci):
      """Returns path to eci.json"""
      return os.path.join(self.eci_dir(clex, calctype, ref, bset, eci), "eci.out")


    # private:

    def __bset(self, bset):
      return "bset." + bset

    def __calctype(self, calctype):
      return "calctype." + calctype

    def __ref(self, ref):
      return "ref." + ref

    def __clex(self, clex):
      return "clex." + clex

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
        if os.path.isdir(item) and item[:len(pattern)+1] == pattern:
          all.append(item[len(pattern)+1:])
          
      return sorted(all)

    

class Project(object):
    """The Project class contains information about a CASM project
    """
    def __init__(self, path=None):
        """
        Construct a CASM Project representation.

        Args:
            path: path to CASM project (Default=None, uses project containing current directory). 

        """
        if path == None:
          if casm.project_path(path) == None:
            if path == None:
              raise Exception("No CASM project found using " + os.getcwd())
            else:
              raise Exception("No CASM project found using " + path)
        self.path = casm.project_path(path)
        self.dir = DirectoryStructure(path)
        self.settings = ProjectSettings(path)
        
    def command(self, args, execcasm = "casm"):
        """
        Execute a command via the command line interface. 
        
        Args:
          args: A string containing the command to be executed. Ex: "select --set-on -o /abspath/to/my_selection"
        
        Returns:
          (stdout, stderr): The result of running the command via the command line iterface
        """
        cwd = os.getcwd()
        os.chdir(self.path)
        result = subprocess.Popen([execcasm] + args.split(' '),stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
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
        
        args = "query -k "
        for k in columns:
          args += k + " "
        if selection.path != "MASTER":
          args += " -c " + selection.path
        if verbatim == True:
          args += " -v"
        args += " -o STDOUT"
        
        (stdout, stderr) = self.command(args)
        
        return pandas.read_csv(StringIO.StringIO(stdout[1:]), sep=' *') 

