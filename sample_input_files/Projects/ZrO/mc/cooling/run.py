import subprocess
import json
import os
import numpy as np
import copy
from casm.project import Project

d = {
  "mode" : "incremental",
  "dependent_runs" : True, 
  "motif" : {
    "configname" : "restricted_auto",
  },
  "initial_conditions" : {
    "param_chem_pot" : {
      "a" : -3.00
    },
    "temperature" : 1000.0,
    "tolerance" : 0.001
  },
  "final_conditions" : {
    "param_chem_pot" : {
      "a" : -3.00
    },
    "temperature" : 200.0,
    "tolerance" : 0.001
  },
  "incremental_conditions" : {
    "param_chem_pot" : {
      "a" : 0.0
    },
    "temperature" : -100.0,
    "tolerance" : 0.001
  }
}

with open('../example_grand_canonical/metropolis_grand_canonical.json','r') as f:
  input = json.load(f)

cwd = os.getcwd()
proj = Project()

index = 0
for xi in np.arange(-3.,0.,0.1):
  path = 'path.' + str(index)
  
  tinput = copy.deepcopy(input)
  tinput['driver'] = d
  tinput['driver']['initial_conditions']['param_chem_pot']['a'] = xi
  tinput['driver']['final_conditions']['param_chem_pot']['a'] = xi
  
  os.mkdir(path)
  os.chdir(path)
  
  with open('input.json','w') as f:
    json.dump(tinput, f)
  
  stdout, stderr, returncode = proj.command('monte -s input.json')
  
  print stdout
  
  os.chdir(cwd)
  index += 1




