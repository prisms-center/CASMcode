import warnings
import casm
import vasp
import os, subprocess, json
from os.path import join
import ctypes, glob

if 'LIBCASM' in os.environ:
  libname = os.environ['LIBCASM']
elif 'CASMPREFIX' in os.environ:
  libname = glob.glob(join(os.environ['CASMPREFIX'], 'lib', 'libcasm.*'))[0]
else:
  libname = glob.glob(join('/usr', 'local', 'lib', 'libcasm.*'))[0]
lib_casm = ctypes.CDLL(libname, mode=ctypes.RTLD_GLOBAL)

if 'LIBCCASM' in os.environ:
  libname = os.environ['LIBCCASM']
elif 'CASMPREFIX' in os.environ:
  libname = glob.glob(join(os.environ['CASMPREFIX'], 'lib', 'libccasm.*'))[0]
else:
  libname = glob.glob(join('/usr', 'local', 'lib', 'libccasm.*'))[0]
lib_ccasm = ctypes.CDLL(libname, mode=ctypes.RTLD_GLOBAL)

#### Argument types

lib_ccasm.casm_STDOUT.restype = ctypes.c_void_p

lib_ccasm.casm_STDERR.restype = ctypes.c_void_p


lib_ccasm.casm_nullstream_new.restype = ctypes.c_void_p

lib_ccasm.casm_nullstream_delete.argtypes = [ctypes.c_void_p]
lib_ccasm.casm_nullstream_delete.restype = None


lib_ccasm.casm_ostringstream_new.restype = ctypes.c_void_p

lib_ccasm.casm_ostringstream_delete.argtypes = [ctypes.c_void_p]
lib_ccasm.casm_ostringstream_delete.restype = None

lib_ccasm.casm_ostringstream_size.argtypes = [ctypes.c_void_p]
lib_ccasm.casm_ostringstream_size.restype = ctypes.c_ulong

lib_ccasm.casm_ostringstream_strcpy.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_char)]
lib_ccasm.casm_ostringstream_strcpy.restype = ctypes.POINTER(ctypes.c_char)


lib_ccasm.casm_primclex_new.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
lib_ccasm.casm_primclex_new.restype = ctypes.c_void_p

lib_ccasm.casm_primclex_delete.argtypes = [ctypes.c_void_p]
lib_ccasm.casm_primclex_delete.restype = None

lib_ccasm.casm_capi.argtypes = [ctypes.c_char_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
lib_ccasm.casm_capi.restype = ctypes.c_int


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

    def calctype_dir(self):
      return "calctype."+self.calctype

    def ref_dir(self):
      return "ref."+self.ref

    def bset_dir(self):
      return "bset."+self.bset

    def eci_dir(self):
      return "eci."+self.eci
      

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
        if path is None:
          if casm.project_path(path) is None:
            if path is None:
              raise Exception("No CASM project found using " + os.getcwd())
            else:
              raise Exception("No CASM project found using " + path)
        self.path = casm.project_path(path)
        dir = DirectoryStructure(self.path)
        print dir.project_settings()
        self.data = json.load(open(dir.project_settings()))

        d = self.data["cluster_expansions"][self.data["default_clex"]]
        self._default_clex = ClexDescription(d["name"], d["property"], d["calctype"], d["ref"], d["bset"], d["eci"])

        self._clex = [
          ClexDescription(d[1]["name"], d[1]["property"], d[1]["calctype"], d[1]["ref"], d[1]["bset"], d[1]["eci"]) 
          for d in self.data["cluster_expansions"].iteritems()
        ]

    
    # -- Accessors --
    
    @property
    def cluster_expansions(self):
        return self._clex
    
    @property
    def default_clex(self):
        return self._default_clex
    

class DirectoryStructure(object):
    """Standard file and directory locations for a CASM project"""
    def __init__(self, path=None):
        """
        Construct a CASM Project DirectoryStructure representation.

        Args:
            path: path to CASM project (Default=None, uses project containing current directory). 

        """
        if path is None:
          if casm.project_path(path) is None:
            if path is None:
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
      return __all_settings("bset", join(self.path, self.__bset_dir))

    def all_calctype(self):
      """Check filesystem directory structure and return list of all calctype names"""
      return __all_settings("calctype", join(self.path, self.__calc_dir, self.__set_dir))
    
    def all_ref(self, calctype):
      """Check filesystem directory structure and return list of all ref names for a given calctype"""
      return __all_settings("ref", calc_settings_dir(calctype))

    def all_clex_name(self):
      """Check filesystem directory structure and return list of all cluster expansion names"""
      return __all_settings("clex", join(self.path, self.__clex_dir))

    def all_eci(self, property, calctype, ref, bset):
      """Check filesystem directory structure and return list of all eci names"""
      return __all_settings("eci", join(self.path, self.__clex_dir, __clex_name(property), __calctype(calctype), __ref(ref), __bset(bset)))


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

    def supercell_dir(self, scelname):
      """Return supercell directory path (scelname has format SCELV_A_B_C_D_E_F)"""
      return join(self.path, self.__calc_dir, scelname)

    def configuration_dir(self, configname):
      """Return configuration directory path (configname has format SCELV_A_B_C_D_E_F/I)"""
      return join(self.path, self.__calc_dir, configname)


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
        if os.path.isdir(item) and item[:len(pattern)+1] == pattern:
          all.append(item[len(pattern)+1:])
          
      return sorted(all)

    

class Project(object):
    """The Project class contains information about a CASM project
    """
    def __init__(self, path=None, casm_exe=None, verbose=True):
      """
      Construct a CASM Project representation.

      Args:
          path: path to CASM project (Default=None, uses project containing 
            current directory). 
          case_exe: CASM executable to use for command line interface. (Default
            uses $CASM if it exists in the environment, else "casm")
      """
      # set path to this CASM project
      if path is None:
        if casm.project_path(path) is None:
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
      
      self.path = casm.project_path(path)
      self.dir = DirectoryStructure(path)
      self.settings = ProjectSettings(path)
      self.casm_exe = casm_exe
      self.verbose = verbose
      
      # will hold a ctypes.c_void_p when loading CASM project into memory
      self._ptr = None
    
    
    def __del__(self):
      self.__unload()
    
    
    def __load(self):
      """
      Explicitly load CASM project into memory.
      """
      if self._ptr is None:
        if self.verbose:
          streamptr = lib_ccasm.casm_STDOUT()
        else:
          streamptr = lib_ccasm.casm_ostringstream_new()
        
        self._ptr = lib_ccasm.casm_primclex_new(self.path, streamptr)
        
        if not self.verbose:
          #qstr = ctypes.create_string_buffer(lib_ccasm.casm_ostringstream_size(ss))
          lib_ccasm.casm_ostringstream_delete(streamptr)
        
    
    def __unload(self):
      """
      Explicitly unload CASM project from memory.
      """
      if self._ptr is not None:
        lib_ccasm.casm_primclex_delete(self._ptr)
        self._ptr = None
        
    
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
      cwd = os.getcwd()
      os.chdir(self.path)
      child = subprocess.Popen([self.casm_exe] + args.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      result = child.communicate()
      os.chdir(cwd)
      return (result[0], result[1], child.returncode)
    
    
    def command_via_capi(self, args):
      """
      Execute a command via the c api. 
      
      Args:
        args: A string containing the command to be executed. Ex: "select --set-on -o /abspath/to/my_selection"
      
      Returns:
        (stdout, stderr, returncode): The result of running the command via the command line iterface
      """
      cwd = os.getcwd()
      os.chdir(self.path)
      
      # construct stringstream objects to capture stdout, stderr
      ss = lib_ccasm.casm_ostringstream_new()
      ss_err = lib_ccasm.casm_ostringstream_new()
      
      res = lib_ccasm.casm_capi(args, self.data(), ss, ss_err)
      
      # copy string and delete stringstream
      qstr = ctypes.create_string_buffer(lib_ccasm.casm_ostringstream_size(ss))
      lib_ccasm.casm_ostringstream_strcpy(ss, qstr)
      lib_ccasm.casm_ostringstream_delete(ss)
      
      # copy string and delete stringstream
      qstr_err = ctypes.create_string_buffer(lib_ccasm.casm_ostringstream_size(ss_err))
      lib_ccasm.casm_ostringstream_strcpy(ss_err, qstr_err)
      lib_ccasm.casm_ostringstream_delete(ss_err)
      
      return (qstr, qstr_err, res)
      
def vasp_input_file_names(settings, configdir):
    # Find required input files in CASM project directory tree
    curr_calctype=settings.default_clex.calctype_dir()

    incarfile = casm.settings_path("INCAR",curr_calctype,configdir)
    prim_kpointsfile = casm.settings_path("KPOINTS",curr_calctype,configdir)
    prim_poscarfile = casm.settings_path("POSCAR",curr_calctype,configdir)
    super_poscarfile = os.path.join(configdir,"POS")
    speciesfile = casm.settings_path("SPECIES",curr_calctype,configdir)

    # Verify that required input files exist
    if incarfile is None:
        raise vasp.VaspError("vasp_input_file_names failed. No INCAR file found in CASM project.")
    if prim_kpointsfile is None:
        raise vasp.VaspError("vasp_input_file_names failed. No KPOINTS file found in CASM project.")
    if prim_poscarfile is None:
        warnings.warn("No reference POSCAR file found in CASM project. I hope your KPOINTS mode is A/AUTO/Automatic or this will fail!", vasp.VaspWarning)
    if super_poscarfile is None:
        raise vasp.VaspError("vasp_input_file_names failed. No POS file found for this configuration.")
    if speciesfile is None:
        raise vasp.VaspError("vasp_input_file_names failed. No SPECIES file found in CASM project.")

    return (incarfile, prim_kpointsfile, prim_poscarfile, super_poscarfile, speciesfile)

