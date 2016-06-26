import casm
import json

def write_eci(proj, eci, fit_details=None, proj_settings=None, verbose=False):
    """
    Write eci.json
    
    Args:
      proj: a casm.Project object
      eci: an iterable of tuple of (index, values)
        index: linear index of basis function
        value: ECI value
      fit_details: description of the fitting method used to generate the ECI,
        usually as output by casm.learn.to_json
      proj_settings: a casm.ProjectSettings object with settings used write eci.json file.
        Default=None uses proj.settings.
    
    """
    dir = proj.dir
    if proj_settings == None:
      proj_settings = proj.settings
    
    # read basis.json
    filename = dir.basis(proj_settings.bset())
    with open(filename,'r') as f:
      j = json.load(f)
    
    # edit to add fitting settings
    j["fit"] = fit_details
    
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
    
    if verbose:
      print "Writing:", filename
    with open(filename, 'w') as f:
      f.write(json.dumps(j, indent=2, cls=casm.NoIndentEncoder))

