import casm
import pandas, StringIO

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
  
  stdout, stderr, returncode = proj.command(args)
  
  try:
    return pandas.read_csv(StringIO.StringIO(stdout[1:]), sep=' *', engine='python')
  except:
    print "Error in casm.query"
    print "  proj:", proj.path
    print "  Attempted to execute: '" + args + "'"
    print "---- stdout: ---------------------"
    print stdout
    print "---- stderr: ---------------------"
    print stderr
    print "----------------------------------"
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
  args += " -t " + selection.type
  if verbatim == True:
    args += " -v"
  if all and (selection.path not in ["CALCULATED", "ALL"]):
    args += " -a"
  args += " -o STDOUT"
  return args
  

