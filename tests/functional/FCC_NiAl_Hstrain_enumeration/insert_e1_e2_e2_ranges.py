import json
import sys
import numpy as np


def main():
    f=open(sys.argv[1])
    strain_dof=json.load(f)

    strain_dof["min"]=[-0.8,0,0,0,0,0]
    strain_dof["max"]=[0.8,0.8,0.8,0,0,0]
    strain_dof["increment"]=[0.2,0.2,0.2,0,0,0]
    strain_dof["confignames"]=[strain_dof["initial_configuration"]["identifier"]]

    strain_dofs=json.dumps(strain_dof,indent=4)

    with open(sys.argv[1], "w") as outfile:
        outfile.write(strain_dofs)

if __name__=="__main__":
    main()
