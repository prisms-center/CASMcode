import json

def write_eci(proj, eci, proj_settings=None):
    """
    Write eci.json
    
    Args:
      proj: a casm.Project object
      eci: an iterable of ECI values
      proj_settings: a casm.ProjectSettings object with settings used write eci.json file.
        Default=None uses proj.settings.
    
    """
    dir = proj.dir
    if proj_settings == None:
      proj_settings = proj.settings
    j = json.load(dir.basis(proj_settings))
    if len(j["cluster_functions"]) != len(eci):
      print "Error writing eci.json"
      print "Read:", dir.basis(proj_settings)
      print "Expected", len(j["cluster_functions"]), "ECI"
      print "Received", len(eci), "ECI"
      raise ValueError("ECI length mismatch")
    for index, value in enumerate(eci):
      j["cluster_functions"][index]["eci"] = value
    f = open(dir.eci(proj_settings), 'w')
    json.dump(j, f)
    f.close()