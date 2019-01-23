import os
import gzip


class AimsRunError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class AimsRun:
    """ An object

        The AimsRun class contains:
            self.total_energy: final total energy
            self.forces: final forces on atoms (2d list of double)
            self.atom_type: list of atom types (list of str)
            self.atoms_per_type: list of number of atoms for each atom type (list of int)
            self.lattice:  final lattice (2d list of double)
            self.rec_lat: final reciprocal lattice (2d list of double)
            self.basis: final basis (2d list of double)
            self.coord_mode: coordinate mode for basis (always 'direct')
            self.is_complete: did the calculation run to completion? (bool)
            self.eigenvalues: eigenvalues and energies for doing band structure plots (list of 2d list of double)
            self.all_e_0: energy (e_0_energy) for each electronic step of each ionic step
            self.nelm: NELM (max. number of electronic steps)
    """
    def __init__(self, filename):
        """ Create a run object """
        self.filename = filename
        self.lattice = []
        self.total_energy = None
        self.forces = []
        self.atom_type = []
        self.atoms_per_type = []
        self.lattice = []
        self.rec_lat = []
        self.basis = []
        self.coord_mode = None
        self.is_complete = False
        self.efermi = None
        self.all_e_0 = []

        # Parse geometry first
        geofile = self.filename.replace("std.out", "geometry.in.next_step")

        if not os.path.isfile(geofile):
            err_str = 'Initial geometry is already relaxed, this is not what you want\n'
            err_str += 'Quitting, please perturb input geometry to force relaxation.'
            raise AimsRunError(err_str)

        with open(geofile, 'rb') as f:
            atom = []
            name = []
            for line in f:
                if len(line) > 1:
                    if line[0] == 'lattice_vector':
                        self.lattice.append([float(x) for x in line.split()[1:4]])
                    if line[0] == 'atom_frac':
                        self.coord_mode = 'direct'
                        atom.append([float(x) for x in line.split()[1:4]])
                        name.append(line.split()[4])
                    if line[0] == 'atom':
                        self.coord_mode = 'cartesian'
                        atom.append([float(x) for x in line.split()[1:4]])
                        name.append(line.split()[4])

        symbols = [name[0]]
        for i in range(len(name)):
            tmp = []
            for j in range(len(symbols)):
                tmp.append(symbols[j])
            if tmp.count(name[i]) == 0:
                symbols.append(name[i])

        nums = {}
        for i in range(len(symbols)):
            k = 0
            for j in range(len(atom)):
                if symbols[i] == name[j]:
                    k += 1
            nums[symbols[i]] = k

        for i in range(len(symbols)):
            self.atom_type.append(symbols[i])
            self.atoms_per_type.append(nums[symbols[i]])
            for j in range(len(name)):
                if name[j] == symbols[i]:
                    self.basis.append(atom[j])

        #  Parse output now
        if os.path.isfile(self.filename):
            if self.filename.split('.')[-1].lower() == 'gz':
                f = gzip.open(self.filename, 'rb')
            else:
                f = open(self.filename)
        elif os.path.isfile(self.filename + '.gz'):
            f = gzip.open(self.filename + '.gz', 'rb')
        else:
            raise AimsRunError('file not found: ' + self.filename)

        print('read file', self.filename)

        k = 0
        read_forces = False
        for line in f:
            if b'Total atomic forces' in line:
                read_forces = True
                self.forces = []
                k = 0
            if read_forces:
                line = line.split()
                if k == 0:  # skip line with text trigger
                    k += 1
                    continue
                self.forces.append([float(x) for x in line[2:5]])
                k += 1
                if k >= len(atom):
                    read_forces = False

            if b'| Total energy                  :' in line:
                self.all_e_0.append(float(line.split()[11]))

            if b'Total energy of the DFT / Hartree-Fock s.c.f. calculation      :' in line:
                self.total_energy = float(line.split()[11])

    def is_complete(self):
        """ Return True if FHI-aims ran to completion """
        return self.is_complete
