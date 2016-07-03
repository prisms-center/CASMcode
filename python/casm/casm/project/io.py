import casm
import json

def write_eci(proj, eci, clex=None, verbose=False):
    """
    Write eci.json
    
    Arguments
    ---------
    
      proj: casm.Project instance
        The CASM project
      
      eci: an iterable of tuple of (index, values)
        index: linear index of basis function
        value: ECI value
      
      clex: ClexDescription instance, optional, default=proj.settings.clex_default
        Specifies where to write the ECI
    
    """
    dir = proj.dir
    if clex is None:
      clex = proj.settings.default_clex
    
    # read basis.json
    filename = dir.basis(clex)
    with open(filename,'r') as f:
      j = json.load(f)
    
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
      print "Writing:", filename
    with open(filename, 'w') as f:
      f.write(json.dumps(j, indent=2, cls=casm.NoIndentEncoder))

