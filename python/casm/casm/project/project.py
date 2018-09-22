from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import os
import math
import warnings
from os.path import join
from string import ascii_lowercase

import numpy as np
import six

from casm.project import syminfo
from casm.api import API, casm_command, casm_capture

def project_path(dir=None):
    """
    Crawl up from dir to find '.casm'. If found returns the directory containing the '.casm' directory.
    If not found, return None.

    Args:
    If dir == None, set to os.getcwd()
    """
    if dir == None:
      dir = os.getcwd()
    else:
      dir = os.path.abspath(dir)
    if not os.path.isdir(dir):
      raise Exception("Error, no directory named: " + dir)
    curr = dir
    cont = True
    while cont == True:
        test_path = os.path.join(curr,".casm")
        if os.path.isdir(test_path):
            return curr
        elif curr == os.path.dirname(curr):
            return None
        else:
            curr = os.path.dirname(curr)
    return None

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
          for d in six.iteritems(self.data["cluster_expansions"])
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

      name: str
        Project name

      settings: casm.project.ProjectSettings instance
        Contains project settings

      dir: casm.project.DirectoryStructure instance
        Provides file and directory locations within the project

      prim: casm.project.Prim instance
        Represents the primitive crystal structure

      composition_axes: casm.project.CompositionAxes or None
        Currently selected composition axes, or None

      all_composition_axes: dict(str:casm.project.CompositionAxes)
        Dict containing name:CompositionAxes pairs, including both standard and
        custom composition axes

      verbose: bool
        How much to print to stdout

    """
    def __init__(self, path=None, verbose=True):
      """
      Construct a CASM Project representation.

      Arguments
      ----------

        path: str, optional, default=None
          Path to project root directory. Default=None uses project containing
          current working directory

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

      self.path = project_path(path)
      self.__refresh()
      self.verbose = verbose

      self.all_composition_axes = {}
      if os.path.exists(self.dir.composition_axes()):
          with open(self.dir.composition_axes(), 'r') as f:
              data = json.load(f)
              if "standard_axes" in data:
                  for key, val in six.iteritems(data["standard_axes"]):
                      self.all_composition_axes[key] = CompositionAxes(key, val)
              if "custom_axes" in data:
                  for key, val in six.iteritems(data["custom_axes"]):
                      self.all_composition_axes[key] = CompositionAxes(key, val)
              self.composition_axes = None
              if "current_axes" in data:
                  self.composition_axes = self.all_composition_axes[data["current_axes"]]


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
      self._prim = None

    @property
    def prim(self):
        if self._prim is None:
            self._prim = Prim(self)
        return self._prim

    @property
    def name(self):
        return self.settings.data['name']

    def refresh(self, read_settings=False, read_composition=False, read_chem_ref=False, read_configs=False, clear_clex=False):
      """
      Refresh PrimClex properties to reflect changes to CASM project files.
      """
      if read_settings:
        self.__refresh()
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
      Execute a command via the c api, writing output to stdout/stderr.

      Args:
        args: A string containing the command to be executed. Ex: "select --set-on -o
            /abspath/to/my_selection"

      Returns:
        returncode: The returncode of the command via the
            CASM C API.

      """
      # this also ensures self._api is not None
      data = self.data()
      returncode = self._api.capi_call(args, self.data())
      self.__refresh()
      return returncode

    def capture(self, args, combine_output=False):
      """
      Execute a command via the c api and store stdout/stderr result as str.

      Args:
        args: A string containing the command to be executed. Ex: "select --set-on -o
        /abspath/to/my_selection"

      Returns
      -------
        (stdout, stderr, returncode): The result of running the command via the
            command line iterface. 'stdout' and 'stderr' are in text type ('unicode'/'str'). If
            'combine_output' is True, then returns (combined_output, returncode).

      """
      # this also ensures self._api is not None
      data = self.data()

      # construct stringstream objects to capture stdout, debug, stderr
      ss = self._api.ostringstream_new()
      if combine_output:
          ss_debug = ss
          ss_err = ss
      else:
          ss_debug = self._api.ostringstream_new()
          ss_err = self._api.ostringstream_new()

      self._api.primclex_set_logging(self.data(), ss, ss_debug, ss_err)
      returncode = self._api.capi_call(args, self.data())

      # copy strings and delete stringstreams
      stdout = self._api.ostringstream_to_str(ss)
      self._api.ostringstream_delete(ss)

      if combine_output:
          res = (stdout.decode('utf-8'), returncode)
      else:
          stderr = self._api.ostringstream_to_str(ss_err)
          self._api.ostringstream_delete(ss_err)

          res = (stdout.decode('utf-8'), stderr.decode('utf-8'), returncode)

      # reset logging to write to stdout/stderr
      self._api.primclex_set_logging(self.data(), self._api.stdout(), self._api.stdout(), self._api.stderr())

      self.__refresh()
      return res

    @classmethod
    def init(cls, root, verbose=True):
        """ Calls `casm init` to create a new CASM project in the given directory

        Arguments
        ---------

          root: str (optional, default=os.getcwd())
            A string giving the path to the root directory of the new CASM project. A `prim.json` file must be present in the directory.

          verbose: bool (optional, default=True)
            Passed to casm.project.Project constructor. How much to print to stdout.

        Returns
        -------
          proj: A casm.project.Project instance for the new CASM project.

        Raises
        ------
          An exception is raised if a new project could not be initialized. This could be due to an already existing project, bad or missing input file, or other cause.

        """
        output, returncode = casm_capture("init", root=root, combine_output=True)
        if returncode != 0:
            print(output)
            raise Exception("Could not initialize the project")
        output, returncode = casm_capture("composition ", root=root, combine_output=True)

        return Project(root, verbose=verbose)


class Prim(object):
    """The Primitive Crystal Structure

    Attributes
    ----------

        proj: casm.Project
          the CASM project the Prim belongs to

        lattice_matrix: np.array of shape (3, 3)
          lattice vectors as matrix columns

        lattice_parameters: dict
          Lattice parameters and angles (in degrees), as:
            {'a':a, 'b':b, 'c':c, 'alpha':alpha, 'beta':beta, 'gamma':gamma}

        basis: List(dict)
          crystal basis, as read directly from prim.json (format may change)

        coordinate_mode: str
          crystal basis coordinate_mode, as read directly from prim.json (format
          may change)

        lattice_symmetry_s: str
          lattice point group, in Schoenflies notation

        lattice_symmetry_hm: str
          lattice point group, in Hermann-Mauguin notation

        lattice_system: str
          lattice system name, ('cubic', 'hexagonal', 'rhombohedral', etc.)

        crystal_symmetry_s: str
          crystal point group, in Schoenflies notation

        crystal_symmetry_hm: str
          crystal point group, in Hermann-Mauguin notation

        crystal_system: str
          crystal system name, ('cubic', 'hexagonal', 'trigonal', etc.)

        crystal_family: str
          crystal family name, ('cubic', 'hexagonal', etc.)

        space_group_number: str
          range of possible space group number

        components: List[str]
          occupational components

        elements: List[str]
          all allowed elements

        n_independent_compositions: int
          number of independent composition axes

        degrees_of_freedom: List[str]
          allowed degrees of freedom, from:
            'occupation'

    """

    def __init__(self, proj):
        """
        Construct a CASM Prim

        Arguments
        ---------

          proj: casm.Project, optional, default=Project containing the current working directory
            the CASM project the Prim belongs to

        """
        if proj == None:
            proj = Project()
        elif not isinstance(proj, Project):
            raise Exception("Error constructing Prim: proj argument is not a CASM Project")
        self.proj = proj

        # raw prim.json (for some properties not yet supported in the API)
        with open(self.proj.dir.prim()) as f:
            raw_prim = json.load(f)
        self.lattice_matrix = np.array(raw_prim['lattice_vectors']).transpose()
        self.basis = raw_prim['basis']
        self.coordinate_mode = raw_prim['coordinate_mode']

        def _angle(a, b):
            return math.degrees(math.acos(
                np.dot(a,b) / (np.linalg.norm(a) * np.linalg.norm(b))
            ))

        def _lattice_parameters(L):
            a = np.linalg.norm(L[:,0])
            b = np.linalg.norm(L[:,1])
            c = np.linalg.norm(L[:,2])
            alpha = _angle(L[:,1], L[:,2])
            beta = _angle(L[:,0], L[:,2])
            gamma = _angle(L[:,0], L[:,1])
            return {'a':a, 'b':b, 'c':c, 'alpha':alpha, 'beta':beta, 'gamma':gamma}
        self.lattice_parameters = _lattice_parameters(self.lattice_matrix)

        (stdout, stderr, returncode) = proj.capture("sym")

        # lattice symmetry
        self.lattice_symmetry_s = syminfo.lattice_symmetry(stdout)
        self.lattice_symmetry_hm = syminfo.hm_symmetry(self.lattice_symmetry_s)
        self.lattice_system = syminfo.lattice_system(self.lattice_symmetry_s)

        # crystal symmetry
        self.crystal_symmetry_s = syminfo.crystal_symmetry(stdout)
        self.crystal_symmetry_hm = syminfo.hm_symmetry(self.crystal_symmetry_s)
        self.crystal_system = syminfo.crystal_system(self.crystal_symmetry_s)
        self.crystal_family = syminfo.crystal_family(self.crystal_symmetry_s)
        self.space_group_number = syminfo.space_group_number_map[self.crystal_symmetry_s]

        # composition (v0.2.X: elements and components are identical, only 'occupation' allowed)
        with open(self.proj.dir.composition_axes()) as f:
            raw_composition_axes = json.load(f)

        self.components = raw_composition_axes['standard_axes']['0']['components']
        self.elements = self.components
        self.n_independent_compositions = raw_composition_axes['standard_axes']['0']['independent_compositions']
        self.degrees_of_freedom = ['occupation']


class CompositionAxes(object):
    """A composition axes object

    Attributes
    ----------

        name: str
          composition axes name

        components: List[str]
          occupational components

        n_independent_compositions: int
          number of independent composition axes

        mol_formula: str
          number of each component in terms of the parametric compositions

        param_formula: str
          parametric compositions in terms of the number of components

        end_members: dict of np.array of shape=(n_components,)
          the number of components per unit cell in each end member state, in form:
          {'origin':np.array, 'a':np.array, 'b', np.array, ...}. Order matches
          that given by self.components.


    """
    def __init__(self, name, data):
        self._name = name
        self._data = data

        self._end_members = {}
        for c in ascii_lowercase:
            if c in self._data:
                self.end_members[c] = np.array(self._data[c])[:,0]
            else:
                break
        self._end_members['origin'] = np.array(self._data['origin'])[:,0]

    @property
    def name(self):
        return self._name

    @property
    def components(self):
        return self._data['components']

    @property
    def n_independent_compositions(self):
        return self._data['independent_compositions']

    @property
    def mol_formula(self):
        return self._data['mol_formula']

    @property
    def param_formula(self):
        return self._data['param_formula']

    @property
    def origin(self):
        return self._origin

    @property
    def end_members(self):
        return self._end_members
