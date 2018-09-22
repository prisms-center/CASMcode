"""casm.project file io"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import six
from casm.misc import compat, noindent

def write_eci(proj, eci, fit_details=None, clex=None, verbose=False):
    """
    Write eci.json
    
    Arguments
    ---------
    
      proj: casm.project.Project instance
        The CASM project
      
      eci: List[(index, value)]
        index (int): linear index of basis function
        value (float): ECI value

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
    with open(filename, 'rb') as f:
        j = json.loads(f.read().decode('utf-8'))
    #print(json.dumps(j, indent=2))
    
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
          basis[key] = noindent.NoIndent(val)
    for entry in j["cluster_functions"]:
      entry["orbit"] = noindent.NoIndent(entry["orbit"])
      sites = entry["prototype"]["sites"]
      for i in range(len(sites)):
        sites[i] = noindent.NoIndent(sites[i])
    
    # write eci.json
    filename = dir.eci(clex)
    
    if verbose:
      print("Writing:", filename, "\n")
    with open(filename, 'wb') as f:
      f.write(six.u(json.dumps(j, indent=2, cls=noindent.NoIndentEncoder)).encode('utf-8'))
    
    # refresh proj to reflect new eci
    proj.refresh(clear_clex=True)
