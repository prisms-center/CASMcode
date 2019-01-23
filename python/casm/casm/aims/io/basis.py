import os
import re
import glob
import casm
import casm.project
import casm.aimswrapper

"""
Sample basis File:
BASIS_DIR_PATH = /absolute/path/to/species_defaults
SPECIES    ALIAS    init_mom
Mn3        Mn         3
Mn4        Mn         4
"""


class BasisError(Exception):
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return self.msg


class IndividualBasis:
    """
        The IndividualBasis class contains:
            self.name: the name as listed in the POS file
            self.alias: the species file lists the name, for convenience only
            self.tags: (dict) the tags that need to be modified in the geometry.in for this specie 
                       (i.e. adding initial_moment)
                       All values are stored as strings.
            self.basisdir_base: common directory for all basis inputs
            self.basis_location: location of basis directory relative to self.basisdir_base
            self.basisdir: directory containing particular basis (self.basisdir_base joined with self.basis_location)
            self.basis_file: the filename of the basis set for the specie
    """

    def __init__(self, values, tags, basisdir_base):
        """ Construct an IndividualBasis.
            
            Args:
                values: (str list) entries in basis file row
                tags: (str list) column names 4+ in basis file, basis tags that need to be modified
                basisdir_base: (str) common directory for all basis inputs
        """

        configdir = os.getcwd()
        _res = os.path.split(configdir)
        cfgname = os.path.split(_res[0])[1] + "/" + _res[1]
        casm_dirs = casm.project.DirectoryStructure(configdir)
        casm_sets = casm.project.ProjectSettings(configdir)
        clex = casm_sets.default_clex
        setfile = casm_dirs.settings_path_crawl("relax.json", cfgname, clex)
        settings = casm.aimswrapper.read_settings(setfile)

        if len(values) != (len(tags)):
            raise BasisError("Length of values != length of tags.\nvalues = " + str(values) + "\ntags = " + str(tags))
        self.name = values[0]
        self.alias = values[1]
        self.init_mom = values[2]

        try:
            self.write_basis = not (float(values[2]) == 0)
        except ValueError:
            raise BasisError("Could not read basis: " + str(values))
        self.basisdir_base = os.path.join(basisdir_base, settings["basis"])

        file_pattern = str(self.basisdir_base)+"/*_"+str(values[1]) + "_*"
        for basis_file in glob.glob(file_pattern):
            self.filename = basis_file

        self.basisdir = self.basisdir_base
        self.tags = dict()
        for i, key in enumerate(tags):
            self.tags[key] = values[i]


def basis_settings(filename):
    """ Returns a dict of IndividualBasis objects, with keys equal to their names. """
    try:
        file = open(filename)
    except IOError:
        raise BasisError("Could not open: '" + filename + "'")
    
    # Read BASIS_DIR_PATH from first line
    line = file.readline()
    m = re.match("BASIS_DIR_PATH\s*=\s*(.*)", line)
    if not m:
        err_str = 'Could not read BASIS_DIR_PATH.\n'
        err_str += 'Expected: BASIS_DIR_PATH = /path/to/location/of/species_defaults/\n'
        err_str += 'Found: "' + line + '"'
        raise BasisError(err_str)
    basis_dir_path = m.group(1)
    
    # Parsing the header
    header = file.readline().strip()
    column_names = header.split()
    if len(column_names) < 3:
        raise BasisError("Insufficient number of columns in basis file")
    tags = column_names[:3]
    basis_settings_loc = dict()
    for line in file:
        if line.strip():
            values = line.strip().split()
            basis_settings_loc[values[1]] = IndividualBasis(values, tags, basis_dir_path)
    file.close()

    return basis_settings_loc


def write_basis(filename, basis):
    file = open(filename, 'a')
    for s in basis:
        with open(basis[s].filename, 'r') as bf:
            species_data = bf.readlines()
        for line in species_data:
            file.write(line)
    file.close()
