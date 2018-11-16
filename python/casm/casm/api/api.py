from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import ctypes
import glob
import json
import os
import six
from distutils.spawn import find_executable
from os.path import dirname, join
from sys import platform

import sh

class API(object):
  """
  Class to provide access to the libccasm C API.

  Acts like a singleton by loading libcasm and libccasm into a class member the
  first time (and only the first time) that a new API instance is constructed.

  Each API instance uses the same shared class member to make calls.
  """
  class __API(object):
    """
    Hidden class of which there will be only one instance

    Attributes
    ----------

      lib_casm: ctypes.CDLL
        Handle to libcasm.*

      lib_ccasm: ctypes.CDLL
        Handle to libccasm.*

    """
    def __init__(self):
      """Loads dynamic libraries"""

      try:
          if 'LIBCASM' in os.environ:
              libcasm_path = os.environ['LIBCASM']
          else:
              casm_path = find_executable('ccasm')
              if platform == 'darwin':
                  libcasm_path = sh.grep(sh.otool('-L', casm_path), 'libcasm').split()[0]
                  libcasm_path = libcasm_path.replace('@loader_path', dirname(casm_path))
              else:
                  libcasm_path = sh.grep(sh.ldd(casm_path), 'libcasm').split()[2]
          libccasm_path = libcasm_path.replace('libcasm', 'libccasm')

          self.lib_casm = ctypes.CDLL(libcasm_path, mode=ctypes.RTLD_GLOBAL)
          self.lib_ccasm = ctypes.CDLL(libccasm_path, mode=ctypes.RTLD_GLOBAL)
      except Exception as e:
          print("\n~~~ Error loading casm libraries ~~~")
          if 'LIBCASM' in os.environ:
              libcasm_path = os.environ['LIBCASM']
              print("Looking for libcasm at LIBCASM:", libcasm_path)
              if not os.path.exists(libcasm_path):
                  print("File does not exist")
                  print("Install CASM if it is not installed, or update your PATH, or set LIBCASM to the location of libcasm.")
              else:
                  print("File exists, but for unknown reason could not be loaded.")
          else:
              casm_path = find_executable('ccasm')
              print("find_executable('ccasm'):", casm_path)
              if casm_path is None:
                  print("Could not find 'ccasm' executable. CASM is not installed on your PATH.")
                  print("Install CASM if it is not installed, or update your PATH, or set LIBCASM to the location of libcasm.")
              elif basename(dirname(casm_path)) == ".libs":
                  print("Found 'ccasm' executable in a '.libs' directory. Are you running tests?")
                  if platform == 'darwin':
                      check = join(dirname(casm_path), "libcasm.dylib")
                  else:
                      check = join(dirname(casm_path), "libcasm.so")
                  if os.path.exists(check):
                      print("You probably need to set LIBCASM="+check)
                  else:
                      print("You probably need to re-make or update your PATH")
              else:
                  print("Found 'ccasm', but for unknown reason could not determine libcasm location.")
                  if platform == 'darwin':
                      print("otool -L:")
                      res = sh.otool('-L', casm_path)
                      for val in res:
                          print(val.strip())
                  else:
                      print("ldd:")
                      res = sh.ldd(casm_path)
                      for val in res:
                          print(val.strip())
          print("")
          raise e

      #### Argument types

      self.lib_ccasm.casm_STDOUT.restype = ctypes.c_void_p

      self.lib_ccasm.casm_STDERR.restype = ctypes.c_void_p

      self.lib_ccasm.casm_nullstream.restype = ctypes.c_void_p

      self.lib_ccasm.casm_ostringstream_new.restype = ctypes.c_void_p

      self.lib_ccasm.casm_ostringstream_delete.argtypes = [ctypes.c_void_p]
      self.lib_ccasm.casm_ostringstream_delete.restype = None

      self.lib_ccasm.casm_ostringstream_size.argtypes = [ctypes.c_void_p]
      self.lib_ccasm.casm_ostringstream_size.restype = ctypes.c_ulong

      self.lib_ccasm.casm_ostringstream_strcpy.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_char)]
      self.lib_ccasm.casm_ostringstream_strcpy.restype = ctypes.POINTER(ctypes.c_char)


      self.lib_ccasm.casm_primclex_null.argtypes = None
      self.lib_ccasm.casm_primclex_null.restype = ctypes.c_void_p

      self.lib_ccasm.casm_primclex_new.argtypes = [ctypes.c_char_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
      self.lib_ccasm.casm_primclex_new.restype = ctypes.c_void_p

      self.lib_ccasm.casm_primclex_delete.argtypes = [ctypes.c_void_p]
      self.lib_ccasm.casm_primclex_delete.restype = None

      self.lib_ccasm.casm_primclex_refresh.argtypes = [ctypes.c_void_p, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool, ctypes.c_bool]
      self.lib_ccasm.casm_primclex_refresh.restype = None

      self.lib_ccasm.casm_primclex_set_logging.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
      self.lib_ccasm.casm_primclex_set_logging.restype = None

      self.lib_ccasm.casm_command_list.argtypes = [ctypes.c_void_p]
      self.lib_ccasm.casm_command_list.restype = None

      self.lib_ccasm.casm_capi.argtypes = [ctypes.c_char_p, ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
      self.lib_ccasm.casm_capi.restype = ctypes.c_int

      self.lib_ccasm.casm_capi_call.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
      self.lib_ccasm.casm_capi_call.restype = ctypes.c_int

  __api = None

  def __init__(self):
    """
    Acts like a singleton by loading libcasm and libccasm into a class member the
    first time (and only the first time) that a new API instance is constructed.

    Each API instance uses the same shared class member to make calls.
    """
    if API.__api is None:
      API.__api = API.__API()
    return

  def stdout(self):
    """
    Get a pointer to a CASM::Log that writes to std::cout.

    This does not need to be deleted manually.
    """
    return API.__api.lib_ccasm.casm_STDOUT()

  def stderr(self):
    """
    Get a pointer to a CASM::Log that writes to std::cerr

    This does not need to be deleted manually.
    """
    return API.__api.lib_ccasm.casm_STDERR()

  def nullstream(self):
    """
    Get a pointer to a CASM::Log that writes to null.

    This does not need to be deleted manually.
    """
    return API.__api.lib_ccasm.casm_nullstream()

  def ostringstream_new(self):
    """
    Construct a CASM::OStringStreamLog that writes to std::ostringstream.

    Returns
    -------
      ptr: CASM::OStringStreamLog pointer
        Used for capturing CASM output in a string.
        This ptr needs to be deleted manually by using API.ostringstream_delete(ptr)
    """
    return API.__api.lib_ccasm.casm_ostringstream_new()

  def ostringstream_to_str(self, ptr):
    """
    Copy the value of a CASM::OStringStreamLog to a Python string.

    Arguments
    ---------

      ptr: CASM::OStringStreamLog pointer


    Returns
    -------
      s: str
        The contents of a CASM::OStringStreamLog

    """
    c_str = ctypes.create_string_buffer(API.__api.lib_ccasm.casm_ostringstream_size(ptr))
    API.__api.lib_ccasm.casm_ostringstream_strcpy(ptr, c_str)
    return c_str.value

  def ostringstream_delete(self, ptr):
    """
    Delete a CASM::OStringStreamLog

    Arguments
    ---------

      ptr: CASM::OStringStreamLog pointer

    """
    API.__api.lib_ccasm.casm_ostringstream_delete(ptr)
    return

  def primclex_null(self):
    """
    Construct a CASM::PrimClex nullptr

    Returns
    -------
      ptr: CASM::PrimClex nullptr

    """
    return API.__api.lib_ccasm.casm_primclex_null()

  def primclex_new(self, path, log, debug_log, err_log):
    """
    Construct a new CASM::PrimClex

    Arguments
    ---------

      path: str
        The path to a CASM project

      log: CASM::Log pointer
        A pointer to a CASM::Log to write standard output

      debug_log: CASM::Log pointer
        A pointer to a CASM::Log to write debug output

      err_log: CASM::Log pointer
        A pointer to a CASM::Log to write error output


    Notes
    -------

    CASM::Log pointers can be obtained from:
      1) API.stdout()
      2) API.stderr()
      3) API.nullstream()
      4) API.ostringstream()


    Returns
    -------
      ptr: CASM::PrimClex pointer

    """
    return API.__api.lib_ccasm.casm_primclex_new(six.b(path), log, debug_log, err_log)

  def primclex_refresh(self, ptr, read_settings=False, read_composition=False, read_chem_ref=False, read_configs=False, clear_clex=False):
    """
    Call CASM::PrimClex::refresh to reload PrimClex data from settings files

    Arguments
    ---------

      ptr: CASM::PrimClex pointer

      read_settings: bool, optional, default=False
        If True, read project_settings.json

      read_composition: bool, optional, default=False
        If True, read composition_axes.json

      read_chem_ref: bool, optional, default=False
        If True, read chemical_reference.json

      read_configs: bool, optional, default=False
        If True, read SCEL and config_list.json

      clear_clex: bool, optional, default=False
        If True, clear stored orbitrees, clexulators, and eci


    Note
    ----

      This does not check if what you request will cause problems.

    """
    API.__api.lib_ccasm.casm_primclex_refresh(ptr, read_settings, read_composition, read_chem_ref, read_configs, clear_clex)
    return

  def primclex_delete(self, ptr):
    """
    Delete a CASM::PrimClex

    Arguments
    ---------

      ptr: CASM::PrimClex pointer

    """
    API.__api.lib_ccasm.casm_primclex_delete(ptr)
    return

  def primclex_set_logging(self, primclex, log, debug_log, err_log):
    """
    Set the PrimClex Logging
    """
    API.__api.lib_ccasm.casm_primclex_set_logging(primclex, log, debug_log, err_log)

  def command_list(self):
    """
    Get list of recognized casm commands implemented at the libcasm level

    Returns
    -------
      s: str
        JSON array containing the list of recognized casm commands

    """
    ptr = self.ostringstream_new()
    API.__api.lib_ccasm.casm_command_list(ptr)
    s = self.ostringstream_to_str(ptr)
    self.ostringstream_delete(ptr)
    return s

  def capi(self, args, primclex, root, log, debug_log, err_log):
    """
    Make an API call

    Arguments
    ---------

      args: str
        A string containing the arguments for the casm command to be executed.

          Ex: "select --set-on -o /abspath/to/my_selection"
          Ex: "query -k 'configname selected' -v -o STDOUT"

      primclex: CASM::PrimClex pointer
        A pointer to a CASM::PrimClex, as obtained from API.primclex_new()

      root: str
        A string giving the path to a root directory of a CASM project, typically
        casm.project.Project.path

      log: CASM::Log pointer
        A pointer to a CASM::Log to write standard output

      debug_log: CASM::Log pointer
        A pointer to a CASM::Log to write debug output

      err_log: CASM::Log pointer
        A pointer to a CASM::Log to write error output

    Returns
    -------
      returncode: int

        Possible values:

          0: No error

          1: ERR_INVALID_ARG
            Command line input is not valid

          2: ERR_UNKNOWN
            Misc. or Unknown error

          3: ERR_NO_PROJ
            No CASM project can be found in expected location

          4: ERR_INVALID_INPUT_FILE
            An expected input file is invalid

          5: ERR_MISSING_INPUT_FILE
            An expected input file can not be found

          6: ERR_EXISTING_FILE
            A file might be overwritten

          7: ERR_MISSING_DEPENDS
            Requested command can not be performed because some dependency needs
            to be done first (i.e. no basis set, so can't use clexulator)

          8: ERR_OTHER_PROJ
            Unknown attempting to overwrite another CASM project

    """
    return API.__api.lib_ccasm.casm_capi(six.b(args), primclex, six.b(root), log, debug_log, err_log)

  def capi_call(self, args, primclex):
    """
    Make an API call, using existing PrimClex's path and logging

    Arguments
    ---------

      args: str
        A string containing the arguments for the casm command to be executed.

          Ex: "select --set-on -o /abspath/to/my_selection"
          Ex: "query -k 'configname selected' -v -o STDOUT"

      primclex: CASM::PrimClex pointer
        A pointer to a CASM::PrimClex, as obtained from API.primclex_new()

    Returns
    -------
      returncode: int

        Possible values:

          0: No error

          1: ERR_INVALID_ARG
            Command line input is not valid

          2: ERR_UNKNOWN
            Misc. or Unknown error

          3: ERR_NO_PROJ
            No CASM project can be found in expected location

          4: ERR_INVALID_INPUT_FILE
            An expected input file is invalid

          5: ERR_MISSING_INPUT_FILE
            An expected input file can not be found

          6: ERR_EXISTING_FILE
            A file might be overwritten

          7: ERR_MISSING_DEPENDS
            Requested command can not be performed because some dependency needs
            to be done first (i.e. no basis set, so can't use clexulator)

          8: ERR_OTHER_PROJ
            Unknown attempting to overwrite another CASM project

    """
    return API.__api.lib_ccasm.casm_capi_call(six.b(args), primclex)


def command_list():
    """
    Get list of recognized casm commands implemented at the libcasm level

    Returns
    -------
      s: str
        JSON array containing the list of recognized casm commands

    """
    _api = API()
    return _api.command_list()

def casm_command(args, root=None, combine_output=False):
    """
    Execute a command via the c api, writing output to stdout/stderr. If required by the specific
    command, a temporary PrimClex instance will be constrcuted. Use the Project.command or Project.capture to avoid re-initializations.

    Arguments
    ---------

      args: str
        A string containing the arguments for the casm command to be executed.

          Ex: "select --set-on -o /abspath/to/my_selection"
          Ex: "query -k 'configname selected' -v -o STDOUT"

      root: str (optional, default=os.getcwd())
        A string giving the path to a root directory of a CASM project, typically
        casm.project.Project.path

      combine_output: bool (optional, default=False)
        If True, print stdout and stderr to same str and only ret

    Returns
    -------
      returncode: The result of running the command via the
          command line iterface. 'stdout' and 'stderr' are in text type ('unicode'/'str'). If
          'combine_output' is True, then returns (combined_output, returncode).
    """
    _api = API()

    # set default root
    if root is None:
        root = os.getcwd()

    # construct stringstream objects to capture stdout, debug, stderr
    ss = _api.stdout()
    ss_debug = ss
    if combine_output:
        ss_err = ss
    else:
        ss_err = _api.stderr()

    returncode = _api.capi(args, _api.primclex_null(), root, ss, ss_debug, ss_err)

    return returncode

def casm_capture(args, root=None, combine_output=False):
    """
    Execute a command via the c api and store stdout/stderr result as str. If required by the
    specific command, a temporary PrimClex instance will be constrcuted. Use the Project.command or
    Project.capture to avoid re-initializations.

    Arguments
    ---------

      args: str
        A string containing the arguments for the casm command to be executed.

          Ex: "select --set-on -o /abspath/to/my_selection"
          Ex: "query -k 'configname selected' -v -o STDOUT"

      root: str (optional, default=os.getcwd())
        A string giving the path to a root directory of a CASM project, typically
        casm.project.Project.path

      combine_output: bool (optional, default=False)
        If True, print stdout and stderr to same str and only ret

    Notes
    -------
      Prefer not to use this function if you have an existing casm.project.Project instance that
      the command will operate on. Instead use Project.command() which will execute the command and
      then update the Project instance's ProjectSettings and DirectoryStructure to reflect any
      changes that have occurred.

    Returns
    -------
      (stdout, stderr, returncode): The result of running the command via the
          command line iterface. 'stdout' and 'stderr' are in text type ('unicode'/'str'). If
          'combine_output' is True, then returns (combined_output, returncode).

    """
    _api = API()

    # set default root
    if root is None:
        root = os.getcwd()

    # construct stringstream objects to capture stdout, debug, stderr
    ss = _api.ostringstream_new()
    if combine_output:
        ss_debug = ss
        ss_err = ss
    else:
        ss_debug = ss
        ss_err = _api.ostringstream_new()

    print("_api.capi begin")
    returncode = _api.capi(args, _api.primclex_null(), root, ss, ss_debug, ss_err)
    print("_api.capi done")

    # copy strings and delete stringstreams
    stdout = _api.ostringstream_to_str(ss)
    _api.ostringstream_delete(ss)

    if combine_output:
        res = (stdout.decode('utf-8'), returncode)
    else:
        stderr = self._api.ostringstream_to_str(ss_err)
        _api.ostringstream_delete(ss_err)

        res = (stdout.decode('utf-8'), stderr.decode('utf-8'), returncode)

    return res
