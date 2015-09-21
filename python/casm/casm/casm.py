import os, json

def jobname(configdir):
    """Return a name for PBS jobs for configuration in 'configdir'
       Returns "SCEL_name.config_name".
    """
    tmp = os.path.split(os.path.abspath(configdir))
    return os.path.split(tmp[0])[1] + "." + tmp[1]

def casm_settings(dir=None):
    """
    Crawl up from dir to find '.casm'.  Read '.casm/project_settings.json' as json dict.
    If not found, return None.
    
    - Currently prepends "calctype." to the found calctype, so that it is of the form "calctype.X" as in 
      previous CASM versions. 
    
    
    If dir == None, set to os.getcwd()
    """
    if dir == None:
      dir = os.getcwd()
    curr = dir
    cont = True
    while cont == True:
        test_path = os.path.join(curr,".casm")
        if os.path.isdir(test_path):
            input = json.load( open(os.path.join(test_path, "project_settings.json")))
            input["curr_calctype"] = "calctype." + input["curr_calctype"]
            return input
        elif curr == os.path.dirname(curr):
            return None
        else:
            curr = os.path.dirname(curr)
    return None

def settings_path(name, calctype, configdir=None):
    """
    Crawl casm directory structure starting at configdir and moving upwards 
    towards the casm project root directory to find the most relevant settings file or directory.
    
    Looks for: ".../settings/calctype." + calctype + "/name'"

    If configdir == None, gets set to os.getcwd()
    
    Return None if name not found
    """
    if configdir == None:
        configdir = os.getcwd()
    curr = configdir
    cont = True
    while cont == True:
        check = os.path.join(curr,"settings", calctype, name)
        if os.path.exists(check):
            return check
        if os.path.exists(os.path.join(curr,".casm")):
            return None
        elif curr == os.path.dirname(curr):
            return None
        else:
            curr = os.path.dirname(curr)
    return None
