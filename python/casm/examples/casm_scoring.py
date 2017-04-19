import pandas

def score(sel):
    for key in ['basis_deformation', 'lattice_deformation']:
        sel.data.loc[:,key] = pandas.to_numeric(sel.data[key], errors='coerce')
    return sel.data['basis_deformation']*0.5 + sel.data['lattice_deformation']*0.5
