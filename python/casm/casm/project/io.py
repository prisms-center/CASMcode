import casm
import json

def write_eci(proj, eci, proj_settings=None):
    """
    Write eci.json
    
    Args:
      proj: a casm.Project object
      eci: an iterable of tuple of (index, values)
        index: linear index of basis function
        value: ECI value
      proj_settings: a casm.ProjectSettings object with settings used write eci.json file.
        Default=None uses proj.settings.
    
    """
    dir = proj.dir
    if proj_settings == None:
      proj_settings = proj.settings
    
    # read basis.json
    filename = dir.basis(proj_settings.bset())
    f = open(filename,'r')
    j = json.load(f)
    f.close()
    
    # edit to add eci
    for index, value in eci:
      j["cluster_functions"][index]["eci"] = value
    
    # pretty printing
    for entry in j["site_functions"]:
      if entry["basis"] != None:
        basis = entry["basis"]
        for key, val in basis.items():
          basis[key] = casm.NoIndent(val)
    for entry in j["cluster_functions"]:
      entry["orbit"] = casm.NoIndent(entry["orbit"])
      sites = entry["prototype"]["sites"]
      for i in range(len(sites)):
        sites[i] = casm.NoIndent(sites[i])
    
    # write eci.json
    filename = dir.eci(proj_settings.clex(), 
                       proj_settings.calctype(), 
                       proj_settings.ref(), 
                       proj_settings.bset(), 
                       proj_settings.eci())
    
    f = open(filename, 'w')
    f.write(json.dumps(j, indent=2, cls=casm.NoIndentEncoder))
    f.close()