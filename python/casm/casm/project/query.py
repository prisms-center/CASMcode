from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

# conda's current version of pandas raises these warnings, but they are safe
# see: https://stackoverflow.com/questions/40845304
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

from io import StringIO
import pandas
import six
import casm
from casm.misc import compat

def query(proj, columns, selection=None, verbatim=True, all=False):
  """Return a pandas DataFrame object containing the output of a
     'casm query' command.

     Args:
       proj: Project to query (default is CASM project containing the current working directory)
       columns: iterable of strings corresponding to 'casm query -k' args
       selection: a Selection to query (default is "MASTER" selection)
       verbatim: if True, use 'casm query --verbatim' option (default is True)
       all: if True, use 'casm query --all' option (default is False)

     Returns:
       data: a pandas DataFrame containing the query results
  """
  args = _query_args(proj, columns, selection, verbatim, all, api=True)

  stdout, stderr, returncode = proj.capture(args)

  try:
    sout = StringIO(stdout)
    if compat.peek(sout) == '#':
        sout.read(1)
    return pandas.read_csv(sout, sep=compat.str(' +'), engine='python')
  except:
    print("Error in casm.query")
    print("  proj:", proj.path)
    print("  Attempted to execute: '" + args + "'")
    print("---- stdout: ---------------------")
    print(stdout)
    print("---- stderr: ---------------------")
    print(stderr)
    print("----------------------------------")
    raise


def _query_args(proj, columns, selection=None, verbatim=True, all=False, api=False):
  """
  Args:
       columns: iterable of strings corresponding to 'casm query -k' args
       selection: a Selection to query (default is "MASTER" selection)
       verbatim: if True, use 'casm query --verbatim' option (default is True)
       all: if True, use 'casm query --all' option (default is False)
       api: if True, args string as if for query_via_api, else as if for query_via_cli
  """
  if selection == None:
    selection = casm.project.Selection(proj)
  elif not isinstance(selection, casm.project.Selection):
    raise Exception("Error, argument 'selection' must be None or a Selection")

  args = "query -k "
  if api:
    args += "'"
  for k in columns:
    args += k + " "
  if api:
    args += "'"
  if selection.path != "MASTER":
    args += " -c " + selection.path
  if verbatim == True:
    args += " -v"
  if all and (selection.path not in ["CALCULATED", "ALL"]):
    args += " -a"
  args += " -o STDOUT"
  return args
