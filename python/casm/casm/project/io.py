import casm
import json

def write_eci(proj, eci, fit_details=None, clex=None, verbose=False):
    """
    Write eci.json
    
    Arguments
    ---------
    
      proj: casm.project.Project instance
        The CASM project
      
      eci: an iterable of tuple of (index, values)
        index: linear index of basis function
        value: ECI value

      fit_details: Dict
        Description of the fitting method used to generate the ECI,
        usually as output by casm.learn.to_json
      
      clex: ClexDescription instance, optional, default=proj.settings.default_clex
        Specifies where to write the ECI
    
    """
    dir = proj.dir
    if clex is None:
      clex = proj.settings.default_clex
    
    # read basis.json
    filename = dir.basis(clex)
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
    filename = dir.eci(clex)
    
    if verbose:
      print "Writing:", filename, "\n"
    with open(filename, 'w') as f:
      f.write(json.dumps(j, indent=2, cls=casm.NoIndentEncoder))
    
    # refresh proj to reflect new eci
    proj.refresh(clear_clex=True)

