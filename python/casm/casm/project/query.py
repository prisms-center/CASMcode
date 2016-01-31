import casm
from casm.project import lib_ccasm
import pandas, StringIO, ctypes

def query(proj, columns, selection=None, verbatim=True):
    """Return a pandas DataFrame object containing the output of a 
       'casm query' command.
       
       Args:
         proj: Project to query (default is CASM project containing the current working directory)
         columns: iterable of strings corresponding to 'casm query -k' args
         selection: a Selection to query (default is "MASTER" selection)
         verbatim: if True, use 'casm query --verbatim' option (default is True)
         
       Returns:
         data: a pandas DataFrame containing the query results
    """
    return query_via_capi(proj, columns, selection, verbatim)
    
    
def query_via_cli(proj, columns, selection=None, verbatim=True):
  if selection == None:
    selection = casm.project.Selection()
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
  
  (stdout, stderr) = proj.command(args)
  
  try:
    return pandas.read_csv(StringIO.StringIO(stdout[1:]), sep=' *', engine='python')
  except:
    print "Error in casm.query"
    print "  proj:", proj.path
    print "  executable:", proj.casm_exe
    print "  Attempted to execute: '" + args + "'"
    print "---- stdout: ---------------------"
    print stdout
    print "---- stderr: ---------------------"
    print stderr
    print "----------------------------------"
    raise

def query_via_capi(proj, columns, selection=None, verbatim=True):
    
  if selection == None:
    selection = casm.project.Selection()
  elif not isinstance(selection, casm.project.Selection):
    raise Exception("Error, argument 'selection' must be None or a Selection")
  
  args = "query -k '"
  for k in columns:
    args += k + " "
  args = args[:-1] + "'"
  if selection.path != "MASTER":
    args += " -c " + selection.path
  if verbatim == True:
    args += " -v"
  args += " -o STDOUT"
  
  # construct stringstream objects
  ss = lib_ccasm.ostringstream_new()
  ss_err = lib_ccasm.ostringstream_new()
  
  lib_ccasm.primclex_check(proj.data())
  
  res = lib_ccasm.query(args, proj.data(), ss, ss_err)
  
  # copy string and delete stringstream
  qstr = ctypes.create_string_buffer(lib_ccasm.ostringstream_size(ss))
  lib_ccasm.ostringstream_strcpy(ss, qstr)
  lib_ccasm.ostringstream_delete(ss)
  
  # copy string and delete stringstream
  qstr_err = ctypes.create_string_buffer(lib_ccasm.ostringstream_size(ss_err))
  lib_ccasm.ostringstream_strcpy(ss_err, qstr_err)
  lib_ccasm.ostringstream_delete(ss_err)
  
  try:
    return pandas.read_csv(StringIO.StringIO(qstr.value[1:]), sep=' *', engine='python')
  except:
    print "Error in casm.query"
    print "  proj:", proj.path
    print "  Attempted to execute: '" + args + "'"
    print "---- stdout: ---------------------"
    print qstr.value
    print "---- stderr: ---------------------"
    print qstr_err.value
    print "----------------------------------"
    raise
    