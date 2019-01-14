from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import numpy as np
import re,copy,math

##############################################NAMELIST BASE CLASS TO SHARE COMMON METHODS#########################

#These dictionaries represent the recognized options read from quantum espresso namelists
QUANTUM_ESPRESSO_INT_LIST ={
            'CONTROL' : ['nstep','iprint','nberrycyc','gdir','nppstr'],
            'SYSTEM' :['ibrav','nat','ntyp','nbnd','nr1','nr2','nr3','nr1s','nr2s','nr3s','nspin',\
                    'nqx1','nqx2','nqx3','lda_plus_u_kind','edir','report','esm_nfit','space_group','origin_choice'],
            'ELECTRONS' :['electron_maxstep','mixing_ndim','mixing_fixed_ns','ortho_para','diago_cg_maxiter','diago_david_ndim'],
            'IONS' : ['nraise','bfgs_ndim'],
            'CELL' : []
}
QUANTUM_ESPRESSO_FLOAT_LIST ={
            'CONTROL' : ['dt','max_seconds','etot_conv_thr','forc_conv_thr'],
            'SYSTEM' : ['celldm(1)','celldm(2)','celldm(3)','celldm(4)','celldm(5)','celldm(6)',\
                        'a','b','c','cosab','cosac','cosbc','tot_charge','ecutwfc','ecutrho',\
                        'degauss','ecfixed','qcutz','q2sigma','exx_fraction','screening_parameter',\
                        'ecutvcut','emaxpos','eopreg','eamp','fixed_magnetization(1)','fixed_magnetization(2)','fixed_magnetization(3)'\
                        'lambda','esm_w','esm_efield','fcp_mu','london_s6','london_rcut','xdm_a1','xdm_a2'],
            'ELECTRONS' :['conv_thr','conv_thr_init','conv_thr_multi','mixing_beta','diago_thr_init','efield','efield_cart(1)','efield_cart(2)','efield_cart(3)'],
            'IONS' : ['tempw','tolp','delta_t','upscale','trust_radius_max','trust_radius_max','trust_radius_ini','w_1','w_2'],
            'CELL' : ['press','wmass','cell_factor','press_conv_thr']
}
QUANTUM_ESPRESSO_BOOL_LIST ={
            'CONTROL' : ['wf_collect','tstress','tprnfor','lkpoint_dir','tefield','dipfield',\
                        'lelfield','lorbm','lberry','lfcpopt'],
            'SYSTEM' : ['nosym','nosym_evc','noinv','no_t_rev','force_symmorphic','use_all_frac',\
                        'one_atom_occupations','starting_spin_angle','noncolin','x_gamma_extrapolation',\
                        'lda_plus_u','lspinorb','london','xdm','uniqueb','rhombohedral'],
            'ELECTRONS' : ['scf_must_converge','adaptive_thr','diago_full_acc','tqr'],
            'IONS' : ['remove_rigid_rot','refold_pos'],
            'CELL' : []
}
QUANTUM_ESPRESSO_STR_LIST ={
            'CONTROL' : ['calculation','title','verbosity','restart_mode','outdir','wfcdir',\
                        'prefix','disk_io','pseudo_dir'],
            'SYSTEM' : ['occupations','smearing','input_dft','exxdiv_treatment','u_projection_type','constrained_magnetization'\
                        'assume_isolated','esm_bc','vdw_corr'],
            'ELECTRONS' :['mixing_mode','diagonalization','startingpot','startingwfc'],
            'IONS' : ['ion_dynamics','ion_positions','pot_extrapolation','wfc_extrapolation','ion_temperature'],
            'CELL' : ['cell_dynamics','cell_dofree']
}

QUANTUM_ESPRESSO_TOTAL_LIST ={
            'CONTROL' : QUANTUM_ESPRESSO_INT_LIST['CONTROL'] + QUANTUM_ESPRESSO_FLOAT_LIST['CONTROL'] + QUANTUM_ESPRESSO_BOOL_LIST['CONTROL'] + QUANTUM_ESPRESSO_STR_LIST['CONTROL'],
            'SYSTEM' : QUANTUM_ESPRESSO_INT_LIST['SYSTEM'] + QUANTUM_ESPRESSO_FLOAT_LIST['SYSTEM'] + QUANTUM_ESPRESSO_BOOL_LIST['SYSTEM'] + QUANTUM_ESPRESSO_STR_LIST['SYSTEM'],
            'ELECTRONS' :QUANTUM_ESPRESSO_INT_LIST['ELECTRONS'] + QUANTUM_ESPRESSO_FLOAT_LIST['ELECTRONS'] + QUANTUM_ESPRESSO_BOOL_LIST['ELECTRONS'] + QUANTUM_ESPRESSO_STR_LIST['ELECTRONS'],
            'IONS' : QUANTUM_ESPRESSO_INT_LIST['IONS'] + QUANTUM_ESPRESSO_FLOAT_LIST['IONS'] + QUANTUM_ESPRESSO_BOOL_LIST['IONS'] + QUANTUM_ESPRESSO_STR_LIST['IONS'],
            'CELL' : QUANTUM_ESPRESSO_INT_LIST['CELL'] + QUANTUM_ESPRESSO_FLOAT_LIST['CELL'] + QUANTUM_ESPRESSO_BOOL_LIST['CELL'] + QUANTUM_ESPRESSO_STR_LIST['CELL']
}




QUANTUM_ESPRESSO_CONTROL_LIST = QUANTUM_ESPRESSO_INT_LIST['CONTROL'] + QUANTUM_ESPRESSO_FLOAT_LIST['CONTROL'] + QUANTUM_ESPRESSO_BOOL_LIST['CONTROL'] + QUANTUM_ESPRESSO_STR_LIST['CONTROL']
QUANTUM_ESPRESSO_SYSTEM_LIST = QUANTUM_ESPRESSO_INT_LIST['SYSTEM'] + QUANTUM_ESPRESSO_FLOAT_LIST['SYSTEM'] + QUANTUM_ESPRESSO_BOOL_LIST['SYSTEM'] + QUANTUM_ESPRESSO_STR_LIST['SYSTEM']
QUANTUM_ESPRESSO_ELECTRONS_LIST = QUANTUM_ESPRESSO_INT_LIST['ELECTRONS'] + QUANTUM_ESPRESSO_FLOAT_LIST['ELECTRONS'] + QUANTUM_ESPRESSO_BOOL_LIST['ELECTRONS'] + QUANTUM_ESPRESSO_STR_LIST['ELECTRONS']
QUANTUM_ESPRESSO_IONS_LIST = QUANTUM_ESPRESSO_INT_LIST['IONS'] + QUANTUM_ESPRESSO_FLOAT_LIST['IONS'] + QUANTUM_ESPRESSO_BOOL_LIST['IONS'] + QUANTUM_ESPRESSO_STR_LIST['IONS']
QUANTUM_ESPRESSO_CELL_LIST = QUANTUM_ESPRESSO_INT_LIST['CELL'] + QUANTUM_ESPRESSO_FLOAT_LIST['CELL'] + QUANTUM_ESPRESSO_BOOL_LIST['CELL'] + QUANTUM_ESPRESSO_STR_LIST['CELL']

class NamelistError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class NamelistBase:
    """
    The Namelist class contains:
        title: header of the namelist
        tags: a dict of all Namelist tags

    All namelist tags and associated values are stored as key-value pairs in the dictionary called 'tags'.
    """
    def __init__(self, title, nameliststring):
        """ Construct a Namelist object from 'nameliststring'"""
        self.title = title
        self.read(nameliststring)

    def read(self, nameliststring):
        """ Read a  namelist """
        self.tags=dict()
        line_segments=re.split('\n',nameliststring)
        # parse nameliststringlist into self.tags
        for line in line_segments:
                commasplit = re.split(',',line)
                if len(commasplit) > 0:
                    line = re.split('=',commasplit[0])
                    if len(line) == 2:
                        if line[0].strip()[0]!= '!':
                            self.tags[line[0].strip()] =  line[1].strip()
        self._verify_tags()
        self._make_natural_type()

    def _make_natural_type(self):
        """ Convert self.tags values from strings into their 'natural type' (int, float, etc.) """
        for tag in self.tags:
            if self.tags[tag] == None or str(self.tags[tag]).strip() == "":
                self.tags[tag] = None
            else:
                if tag.lower() in QUANTUM_ESPRESSO_INT_LIST[self.title]:
                    try:
                        self.tags[tag] = int(self.tags[tag])
                    except ValueError:
                        raise NamelistError("Could not convert '" + tag + "' : '" + self.tags[tag] + "' to int")
                elif tag.lower() in QUANTUM_ESPRESSO_FLOAT_LIST[self.title]:
                    try:
                        self.tags[tag] = float(self.tags[tag].lower().replace('d','e'))
                    except ValueError:
                        raise NamelistError("Could not convert '" + tag + "' : '" + self.tags[tag] + "' to float")
                elif tag.lower() in QUANTUM_ESPRESSO_BOOL_LIST[self.title]:
                    if not self.tags[tag].lower() in ['.true.','.false.']:
                        raise NamelistError("Could not find '" + tag + "' : '" + self.tags[tag].lower() + "' in ['.true.','.false.']")
                    else:
                        self.tags[tag] = (self.tags[tag].lower() == '.true.')
                elif tag.lower() in QUANTUM_ESPRESSO_STR_LIST[self.title]:
                    self._check_string_tag(tag,self.tags[tag])

    def _verify_tags(self):
        """ Check that only allowed Namelist tags are in self.tags """
        TOT_LIST = QUANTUM_ESPRESSO_INT_LIST[self.title] + QUANTUM_ESPRESSO_FLOAT_LIST[self.title] + QUANTUM_ESPRESSO_BOOL_LIST[self.title] + QUANTUM_ESPRESSO_STR_LIST[self.title]
        for tag in self.tags:
            if tag in TOT_LIST:
                continue
            else:
                print("Warning: unknown "+ self.title + " tag '" + tag + "' with value '" + str(self.tags[tag]) + "'")

    def make_string(self):
        """Convert namelist to string for writing"""
        nameliststr=""
        for tag in self.tags:
            if self.tags[tag] == None or str(self.tags[tag]).strip() == "":
                pass
            elif type(self.tags[tag])==type(True):
                nameliststr= nameliststr + " " + tag + " = ." + str(self.tags[tag]).lower() + ".,\n"
            else:
                nameliststr= nameliststr + " " + tag + " = " + str(self.tags[tag]) + ",\n" 
        return nameliststr

    def _check_string_tag(self,tag,value):
        pass

##############################################CONTROL MODULE BEGINS HERE############################################
class ControlError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class Control(NamelistBase):
    """
    The Control class contains:
        tags: a dict of all Control tags

    All namelist tags and associated values are stored as key-value pairs in the dictionary called 'tags'.
   """
    def __init__(self, title, nameliststring):
        """ Construct a Control object from 'nameliststring'"""
        #self.read(nameliststring)
        NamelistBase.__init__(self,title,nameliststring)

    def _check_string_tag(self,tag,value):
        """ Check that string-valued tags are allowed values """
        if tag.lower() == 'calculation':
            if value.lower() not in ["'scf'" ,"'nscf'" ,"'bands'" ,"'relax'" ,"'md'" ,"'vc-relax'","'vc-md'"]:
                raise ControlError("Unknown 'calculation' value: '" + value)
        elif tag.lower() == 'verbosity':
            if value.lower() not in ["'low'","'medium'","'debug'","'high'","'default'","'minimal'"]:
                raise ControlError("Unknown 'verbosity' value: '" + value)
        elif tag.lower() == 'restart_mode':
            if value.lower() not in ["'from_scratch'","'restart'"]:
                raise ControlError("Unknown 'restart_mode' value: '" + value)
        elif tag.lower() == 'disk_io':
            if value.lower() not in ["'high'","'medium'","'low'","'none'"]:
                raise ControlError("Unknown 'disk_io' value: '" + value)

##############################################END CONTROL MODULE#############################################################



##############################################SYSTEM MODULE BEGINS HERE############################################
class SystemError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class System(NamelistBase):
    """
    The Sys class contains:
        tags: a dict of all System tags

    All namelist tags and associated values are stored as key-value pairs in the dictionary called 'tags'.
   """
    def __init__(self, title, nameliststring):
        """ Construct a Sys object from 'nameliststring'"""
        NamelistBase.__init__(self,title,nameliststring)

    def _check_string_tag(self,tag,value):
        """ Check that string-valued tags are allowed values """
        if tag.lower() == 'occupations':
            if value.lower() not in ["'smearing'","'tetrahedra'","'fixed'","'from_input'"]:
                raise SysError("Unknown 'occupations' value: '" + value)
        elif tag.lower() == 'smearing':
            if value.lower() not in ["'gaussian'","'gauss'","'methfessel-paxton'","'m-p'","'mp'","'marzari-vanderbilt'","'cold'","'m-v'","'mv'","'fermi-dirac'","'f-d'","'fd'"]:
                raise SysError("Unknown 'smearing' value: '" + value)
        elif tag.lower() == 'exxdiv_treatment':
            if value.lower() not in ["'gygi-baldereschi'","'vcut_spherical'","'vcut_ws'","'none'"]:
                raise SysError("Unknown 'exxdiv_treatment' value: '" + value)
        elif tag.lower() == 'u_projection_type':
            if value.lower() not in ["'atomic'","'ortho-atomic'","'norm-atomic'","'file'","'pseudo'"]:
                raise SysError("Unknown 'U_projection_type' value: '" + value)
        elif tag.lower() == 'constrained_magnetization':
            if value.lower() not in ["'none'","'total'","'atomic'","'total direction'","'atomic direction'"]:
                raise SysError("Unknown 'constrained_magnetization' value: '" + value)
        elif tag.lower() == 'assume_isolated':
            if value.lower() not in ["'none'","'makov-payne'","'m-p'","'mp'","'martyna-tuckerman'","'m-t'","'mt'","'esm'"]:
                raise SysError("Unknown 'assumed_isolated' value: '" + value)
        elif tag.lower() == 'esm_bc':
            if value.lower() not in ["'pbc'","'bc1'","'bc2'","'bc3'"]:
                raise SysError("Unknown 'esm_bc' value: '" + value)
        elif tag.lower() == 'vdw_corr':
            if value.lower() not in ["'grimme-d2'","'dft-d'","'ts'","'ts-vdw'","'xdm'"]:
                raise SysError("Unknown 'esm_bc' value: '" + value)

##############################################END SYSTEM MODULE#############################################################



##############################################ELECTRONS MODULE BEGINS HERE############################################
class ElectronsError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class Electrons(NamelistBase):

    """
    The Electrons class contains:
        tags: a dict of all Electrons tags

    All namelist tags and associated values are stored as key-value pairs in the dictionary called 'tags'.
   """
    def __init__(self,title, nameliststring):
        """ Construct an Electrons object from 'nameliststring'"""
        NamelistBase.__init__(self,title,nameliststring)

    def _check_string_tag(self,tag,value):
        """ Check that string-valued tags are allowed values """
        if tag.lower() == 'mixing_mode':
            if value.lower() not in ["'plain'","'tf'" ,"'local-tf'"]:
                raise ElectronsError("Unknown 'mixing_mode' value: '" + value)
        elif tag.lower() == 'diagonalization':
            if value.lower() not in ["'david'","'cg'","'cg-serial'"]:
                raise ElectronsError("Unknown 'diagonalization' value: '" + value)
        elif tag.lower() == 'startingpot':
            if value.lower() not in ["'atomic'","'file'"]:
                raise ElectronsError("Unknown 'startingpot' value: '" + value)
        elif tag.lower() == 'startingwfc':
            if value.lower() not in ["'atomic'","'atomic+random'","'random'","'file'"]:
                raise ElectronsError("Unknown 'startingwfc' value: '" + value)
                
##############################################END ELECTRONS MODULE#############################################################



##############################################IONS MODULE BEGINS HERE############################################
class IonsError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class Ions(NamelistBase):
    """
    The Ions class contains:
        tags: a dict of all Ions tags

    All namelist tags and associated values are stored as key-value pairs in the dictionary called 'tags'.
   """

    def __init__(self,title, nameliststring):
        """ Construct an Ions object from 'nameliststring'"""
        NamelistBase.__init__(self,title,nameliststring)
        
    def _check_string_tag(self,tag,value):
        """ Check that string-valued tags are allowed values """
        if tag.lower() == 'ion_dynamics':
            if value.lower() not in ["'bfgs'","'damp'" ,"'verlet'","'langevin'","'langevin-smc'","'beeman'"]:
                raise IonsError("Unknown 'ion_dynamics' value: '" + value)
        elif tag.lower() == 'ion_positions':
            if value.lower() not in ["'default'","'from_input'"]:
                raise IonsError("Unknown 'ion_positions' value: '" + value)
        elif tag.lower() == 'pot_extrapolation':
            if value.lower() not in ["'none'","'atomic'","'first_order'","'second_order'"]:
                raise IonsError("Unknown 'pot_extrapolation' value: '" + value)
        elif tag.lower() == 'wfc_extrapolation':
            if value.lower() not in ["'none'","'first_order'","'second_order'"]:
                raise IonsError("Unknown 'wfc_extrapolation' value: '" + value)
        elif tag.lower() == 'ion_temperature':
            if value.lower() not in ["'rescaling'","'rescale-v'","'rescale-t'","'reduce-t'","'berendsen'","'andersen'","'initial'","'not_controlled'"]:
                raise IonsError("Unknown 'ion_temperature' value: '" + value)

##############################################END IONS MODULE#############################################################

##############################################CELL MODULE BEGINS HERE############################################

class CellError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Cell(NamelistBase):
    """
    The Cell class contains:
        tags: a dict of all Cell tags

    All namelist tags and associated values are stored as key-value pairs in the dictionary called 'tags'.
   """
    def __init__(self, title, nameliststring):
        """ Construct a Cell object from 'nameliststring'"""
        NamelistBase.__init__(self,title,nameliststring)

    def _check_string_tag(self,tag,value):
        """ Check that string-valued tags are allowed values """
        if tag.lower() == 'cell_dynamics':
            if value.lower() not in ["'bfgs'","'damp-pr'" ,"'damp-w'","'none'","'sd'","'pr'","'w'"]:
                raise CellError("Unknown 'cell_dynamics' value: '" + value)
        elif tag.lower() == 'cell_dofree':
            if value.lower() not in ["'all'","'x'","'y'","'z'","'xy'","'xz'","'yz'","'xyz'","'shape'","'volume'","'2dxy'","'2dshape'"]:
                raise CellError("Unknown 'cell_dofree' value: '" + value)

##############################################END CELL MODULE#############################################################

##############################################ATOMIC_SPECIES MODULE BEGINS HERE############################################
class AtomicSpeciesError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class AtomicSpecies:
    """
    The AtomicSpecies class contains:
        species: A list of all atomic species
        masses: A list of all atomic mass of the corresponding species
        pseudos: A list of all pseudo potential files for the corresponding species
   """
    def __init__(self,cardstring):
        """ Construct an AtomicSpecies object from 'cardstring'"""
        self.read(cardstring)

    def read(self, cardstring):
        """ Read an AtomicSpecies card """
        self.species=[]
        self.masses=[]
        self.pseudos=[]
        line_segments=re.split('\n',cardstring)
        # parse card into self.species, self.masses, self.pseudos
        for line in line_segments:
                line = line.strip().split()
                if len(line) == 3:
                    self.species = self.species + [line[0].strip()]
                    try:
                        self.masses = self.masses + [float(line[1].strip().lower().replace('d','e'))]
                    except ValueError:
                        raise AtomicSpeciesError("Could not convert mass to float")
                    self.pseudos = self.pseudos + [line[2].strip()]
    
    def make_string(self):
        """Convert card to string for writing"""
        cardstr=""
        for specie,mass,pseudo in zip(self.species,self.masses,self.pseudos):
                cardstr= cardstr + " " + specie + " " + str(mass) + " " + pseudo + "\n" 
        return cardstr
##############################################END ATOMIC_SPECIES MODULE#############################################################


##############################################ATOMIC_POSITIONS MODULE BEGINS HERE############################################
class AtomicPositionsError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class AtomicPositions:
    """
    The AtomicPositions class contains:
        units: The type of units for positions
        coords: A list of all atomic and coordinates for each species along with fixing switches
   """
    def __init__(self,cardstring,units=""):
        """ Construct an AtomicSpecies object from 'cardstring'"""
        if units=="":
            units="alat" 
        self.units=units
        self.read(cardstring)

        self._check_units(self.units)

    def read(self, cardstring):
        """ Read an AtomicSpecies card """
        self.coords=[]
        line_segments=re.split('\n',cardstring)
        # parse card into self.coords
        # NOTE: This function has trouble because it assumes integer and float arithmetic are separate unlike the documentation of Quantum Espresso indicates.
        for line in line_segments:
                line = line.strip().split()
                if len(line) == 4:
                    try:
                        self.coords = self.coords + [(line[0].strip(),[float(eval(line[1].strip().replace('d','e'))),float(eval(line[2].strip().replace('d','e'))),float(eval(line[3].strip().replace('d','e')))])]
                    except ValueError:
                        raise AtomicPositionsError("Could not convert coordinates to floats")
                elif len(line) == 7:
                    try:
                        self.coords = self.coords + [(line[0].strip(),[float(eval(line[1].strip().replace('d','e'))),float(eval(line[2].strip().replace('d','e'))),float(eval(line[3].strip().replace('d','e')))],[int(line[4].strip()),int(line[5].strip()),int(line[6].strip())])]
                    except ValueError:
                        raise AtomicPositionsError("Could not convert coordinates to floats or if_pos to ints")

    def _check_units(self,units):
        """ Check that units is one of allowed values"""
        if units.lower() not in ['alat','bohr' ,'angstrom','crystal','crystal_sg']:
            raise AtomicPositionsError("Unknown coordinate units:'" + units)
        
    
    def make_string(self):
        """Convert card to string for writing"""
        cardstr=""
        for coord in self.coords:
            if len(coord)==3:
                cardstr= cardstr + " " + coord[0] + " " + str(coord[1][0]) + " " + str(coord[1][1]) + " " + str(coord[1][2]) + " " + str(coord[2][0]) + " " + str(coord[2][1]) + " " + str(coord[2][2]) + "\n" 
            else:
                cardstr= cardstr + " " + coord[0] + " " + str(coord[1][0]) + " " + str(coord[1][1]) + " " + str(coord[1][2]) + "\n" 
        return cardstr
##############################################END ATOMIC_POSITIONS MODULE#############################################################

##############################################CELL_PARAMETERS MODULE BEGINS HERE############################################
class CellParametersError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class CellParameters:
    """
    The CellParameters class contains:
        units: The type of units for vectors
        vectors: A list of all lattice vectors in the above units
   """
    def __init__(self,cardstring,units=""):
        """ Construct a CellParameters object from 'cardstring'"""
        if units=="":
            units="alat" 
        self.units=units
        self.read(cardstring)

        self._check_units(self.units)

    def read(self, cardstring):
        """ Read a CellParameters card """
        self.vectors=[]
        line_segments=re.split('\n',cardstring)
        # parse card into self.vectors
        for line in line_segments:
                line = line.strip().split()
                if len(line) == 3:
                    try:
                        self.vectors = self.vectors + [[float((line[0].strip().replace('d','e'))),float(line[1].strip().replace('d','e')),float(line[2].strip().replace('d','e'))]]
                    except ValueError:
                        raise CellParametersError("Could not convert vector coordinates to floats")
        self.vectors=np.array(self.vectors)
        if self.vectors.shape != (3,3):
            raise CellParametersError("Lattice Shape Error:" + np.array_str(self.vectors))

    def _check_units(self,units):
        """ Check that units is one of allowed values"""
        if units.lower() not in ['alat','bohr' ,'angstrom']:
            raise CellParametersError("Unknown coordinate units:'" + units)
        
    
    def make_string(self):
        """Convert card to string for writing"""
        cardstr=""
        for vector in self.vectors:
            cardstr= cardstr + " " + str(vector[0]) + " " + str(vector[1]) + " " + str(vector[2]) + "\n"
        return cardstr
##############################################END CELL_PARAMETERS MODULE#############################################################


##############################################K_POINTS MODULE BEGINS HERE############################################
class KPointsError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class KPoints:
    """
    The KPoints class contains:
        units: The type of K-point meshing
        coords: A list of all kpoints
        nks: number of kpoints (if the kpoints are given as a list)
   """
    def __init__(self,cardstring,units=""):
        """ Construct a KPoints object from 'cardstring'"""
        if units=="":
            units="tpiba"
        self.units=units
        self.nks=0
        self._check_units(self.units)
        self.read(cardstring)

    def read(self, cardstring):
        """ Read a KPoints card """
        self.coords=[]
        line_segments=re.split('\n',cardstring)
        # parse card into self.species, self.masses, self.pseudos
        if self.units == "gamma":
            pass
        elif self.units == "automatic":
            line=line_segments[0].strip().split()
            if len(line)==6:
                try: #automatic mesh settings are stored as [nk1,nk2,nk3,sk1,sk2,sk3]
                    self.coords=[int(line[0].strip()),int(line[1].strip()),int(line[2].strip()),int(line[3].strip()),int(line[4].strip()),int(line[5].strip())]
                except ValueError:
                    raise KPointsError("Could not convert meshsizing to ints")
            else:
                raise KPointsError("Incorrect format for automatic K_POINT Meshing")
        else:
            try:
                self.nks=int(line_segments[0].strip())
            except ValueError:
                raise KPointsError("Non integer number of kpoints")
            line_segments=line_segments[1:]
            for line in line_segments:

                line = line.strip().split()
                if len(line) == 4:
                    try:
                        self.coords = self.coords + [((float(line[0].strip()),float(line[1].strip()),float(line[2].strip())),float(line[3].strip()))]
                    except ValueError:
                        raise KPointsError("Could not convert coordinates or weights to floats")

    def _check_units(self,units):
        """ Check that units is one of allowed values"""
        if units.lower() not in ['tpiba','crystal' ,'tpiba_b','crystal_b','tpiba_c','crystal_c','automatic','gamma']:
            raise KPointsError("Unknown units:'" + units)
    
    def super_kpoints(self, prim, super):
        """ Assuming 'self' is the kpoints associated with a PRIM, it uses a scaling method to calculate
                  the kpoint-mesh for a supercell, such that it has a equal or greater kpoint
                  density than the prim. 
            
            Returns:
                super_kpoints: a Kpoints object for the supercell
            
            Args:
                prim: Poscar object for the prim OR None
                super: a Poscar object for the supercell 
        """
        backup_prim=copy.deepcopy(prim)
        backup_super=copy.deepcopy(super)
        super_kpoints = copy.deepcopy(self)
        
        
        if prim == None:
            raise KpointsError("No POSCAR was provided for the PRIM, so the PRIM KPOINTS could not be scaled!")

        try: 
            float(super.scaling)
            super.scaling="angstrom"
        except ValueError:
            pass
        if prim.scaling != super.scaling:
            print ("WARNING: The reference lattice vectors and newly scaled lattice vectors are not of same units!!!!")
            if super.scaling != "angstrom":
                if super.scaling == "alat":
                    super._lattice = super.celldm*super._lattice
                    super.scaling = "bohr"
                if super.scaling == "bohr":
                    super._lattice = 0.52918*super._lattice
                    super.scaling = "angstrom"
                super._reciprocal_lattice = 2.0*math.pi*np.linalg.inv(np.transpose(super._lattice))
            if prim.scaling != "angstrom":
                if prim.scaling == "alat":
                    prim._lattice = prim.celldm*prim._lattice
                    prim.scaling = "bohr"
                if prim.scaling == "bohr":
                    prim._lattice = 0.52918*prim._lattice
                    prim.scaling = "angstrom"
                prim._reciprocal_lattice = 2.0*math.pi*np.linalg.inv(np.transpose(prim._lattice))



        super_kpoints.coords = [1, 1, 1, self.coords[3], self.coords[4], self.coords[5]] #4th 5th and 6th numbers are shifts
        
        # calculate prim volumetric kpoint densities
        prim_density = self.density(prim)
        
        # calculate recip lattice vector lengths
        super_recip_vec_lengths = [np.linalg.norm(super.reciprocal_lattice(x)) for x in range(3)]
        
        # while supercell kpoint density is less than prim kpoint density
        while super_kpoints.density(super) < prim_density:
            
            # increase the number of subdivisions along the least dense super recip vector
            linear_density = [super_kpoints.coords[x]/super_recip_vec_lengths[x] for x in range(3)]
            min_index = linear_density.index(min(linear_density))
            super_kpoints.coords[min_index] += 1
            
            # set all subdivisions to be at similar linear density 
            scale = super_kpoints.coords[min_index] / super_recip_vec_lengths[min_index]
            for i in range(3):
                super_kpoints.coords[i] = int(math.ceil(scale * super_recip_vec_lengths[i]-0.1))
        # end while
        prim = backup_prim
        super = backup_super
        return super_kpoints
    
    def density(self, poscar):
        """ Return the kpoint density with respect to a Poscar.
            
            Args:
                poscar: a Poscar object
        """
        return (self.coords[0] * self.coords[1] * self.coords[2]) / poscar.reciprocal_volume() 
    
    def make_string(self):
        """Convert card to string for writing"""
        cardstr=""
        if self.units == "gamma":
            cardstr="\n"
        elif self.units == "automatic":
            cardstr=" " + str(self.coords[0]) + " " + str(self.coords[1]) + " " + str(self.coords[2]) + " " + str(self.coords[3]) + " " + str(self.coords[4]) + " " + str(self.coords[5]) + "\n"
        else:
            cardstr= str(self.nks) + "\n"    
            for coord in self.coords:
                if len(coord)==2:
                    cardstr= cardstr + " " + str(coord[0][0]) + " " + str(coord[0][1]) + " " + str(coord[0][2]) + " " + str(coord[1]) + "\n" 
        return cardstr
##############################################END K_POINTS MODULE#############################################################

########################################################INFILE MODULE BEGINS HERE############################################################


# List of namelists in Quantum espresso 
QUANTUM_ESPRESSO_NAMELIST_LIST = ['CONTROL','SYSTEM','ELECTRONS','IONS','CELL','INPUTPH','PATH',\
                                'INPUTCOND','PRESS_AI','WANNIER','INPUTPP','PLOT','DOS','BANDS',\
                                'PROJWFC']

#List of namelists which have associated Classes
QUANTUM_ESPRESSO_NAMELIST_OBJ_LIST = ['CONTROL','SYSTEM','ELECTRONS','IONS','CELL']

#List of Cards in Quantum espresso
QUANTUM_ESPRESSO_CARD_LIST = ['ATOMIC_SPECIES','ATOMIC_POSITIONS','K_POINTS','CELL_PARAMETERS',\
                                'CONSTRAINTS','OCCUPATIONS','ATOMIC_FORCES','qPointsSpecs',\
                                'CLIMBING_IMAGES','K_and_Energy_Points','ATOMIC_VELOCITIES','REF_CELL_PARMETERS',\
                                'PLOT_WANNIER','AUTOPILOT']

#List of Cards in Quantum espresso which have associated Classes 
QUANTUM_ESPRESSO_CARD_OBJ_LIST = ['ATOMIC_SPECIES','ATOMIC_POSITIONS','K_POINTS','CELL_PARAMETERS']

QUANTUM_ESPRESSO_BLOCK_LIST= QUANTUM_ESPRESSO_NAMELIST_LIST + QUANTUM_ESPRESSO_CARD_LIST

class InfileError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class Infile:
    """
    The Infile class contains:
        namelists: a dict of all Infile namelists (which contain settings)
        cards: a dict of all Infile cards

    All input namelists and associated tag lists are stored as key-value pairs in the dictionary called 'namelists'.
   """
    def __init__(self,filename, species=None, poscar=None, sort=True):
        """ Construct an Infile object from 'filename'"""
        self.read(filename, species, poscar, sort)

    def read(self, filename, species=None, poscar=None, sort=True):
        """ Read an .in file """
        self.namelists = dict()
        self.cards = dict()
        namelist_open=False
        curr_namelist=""
        nameliststring=""
        card_open=False
        curr_card=""
        cardstring=""
        extratag=""
        try:
            file = open(filename,'r')
        except:
            raise InfileError("Could not open file: '" + filename + "'")

        # parse Infile into self.namelist and self.cards
        for line in file:
            line=line.strip()
            if len(line)>0 and line[0] == '&':
                #The check a new namelist is a line starting with &
                if namelist_open:
                    raise InfileError("namelist embedding...Please check '" + filename + "' for proper formatting")
                namelist_open = True
                nameliststring=""
                line = re.split('&',line)
                curr_namelist = line[1].strip()
            elif len(line)>0 and line[0] == '/':
                #The end of a namelist is denoted by /
                if not namelist_open:
                    raise InfileError("namelist embedding...Please check '" + filename + "' for proper formatting")
                namelist_open = False

                namelistobj={
                  'CONTROL': lambda x : Control(curr_namelist,x),
                  'SYSTEM':lambda x : System(curr_namelist,x),
                  'ELECTRONS':lambda x : Electrons(curr_namelist,x),
                  'IONS':lambda x :Ions(curr_namelist,x),
                  'CELL':lambda x :Cell(curr_namelist,x),
                }.get(curr_namelist, lambda x: x)(nameliststring)

                self.namelists[curr_namelist]=namelistobj
            elif namelist_open:
              # white space in name list is removed later on
                    nameliststring=nameliststring + line + '\n'
            else:
                if not card_open:
                    line = line.split()
                    #if a card isn't open for reading and empty line is found skip to next line
                    if len(line) > 0: 
                        #this check for the start of a card is a Capitalized phrase
                        if line[0].isupper():
                            card_open = True
                            cardstring=""
                            curr_card = line[0].strip()
                            extratag=""
                            if curr_card in ["ATOMIC_POSITIONS","K_POINTS","CELL_PARAMETERS"] and len(line)>=2:
                                extratag=line[1].strip()
                elif line.strip()=='':
                    #the only way to check for end of card is an empty line
                    # new line type white space in card is detrimental because it ends the card prematurely
                    card_open = False
                    cardobj={
                      'ATOMIC_SPECIES': lambda x : AtomicSpecies(x),
                      'ATOMIC_POSITIONS':lambda x : AtomicPositions(x,extratag),
                      'K_POINTS':lambda x : KPoints(x,extratag),
                      'CELL_PARAMETERS': lambda x: CellParameters(x,extratag),
                    }.get(curr_card, lambda x: x)(cardstring)
                    self.cards[curr_card]=cardobj
                else:
                    cardstring=cardstring + line + '\n'
        self._verify_blocks()

        file.close()

    def _verify_blocks(self):
        """ Check that only allowed Infile namelists and cards are in self.namelists and self.cards """
        for namelist in self.namelists:
            if namelist in QUANTUM_ESPRESSO_NAMELIST_LIST:
                continue
            else:
                print("Warning: unknown Infile namelist: " + namelist)
        for card in self.cards:
            if card in QUANTUM_ESPRESSO_CARD_LIST:
                continue
            else:
                print("Warning: unknown Infile card: " + card)

    def rewrite_poscar_info(self,poscar,species=None):
        """Update Infile object to a new set of coordinates and lattice vectors"""
        ## This function is called on subsequent runs of quantum espresso i.e. moving from run.0 to run.1 etc.
        ## Knowing what is going on here is dependent on knowledge of quantumespresso/qeio/poscar.py
        ## Knowing what is going on here is dependent on knowledge of quantumespresso/qeio/species.py
        if "CELL_PARAMETERS" not in self.cards.keys():
            self.cards["CELL_PARAMETERS"]=CellParameters("0.0 0.0 0.0 \n  0.0 0.0 0.0 \n 0.0 0.0 0.0")
        if poscar.scaling in ['alat','bohr','angstrom']:
            self.cards["CELL_PARAMETERS"].units=poscar.scaling
        else:
            self.cards["CELL_PARAMETERS"].units="angstrom"
        self.cards["CELL_PARAMETERS"].vectors=poscar._lattice
        if "SYSTEM" in self.namelists.keys():
            self.namelists["SYSTEM"].tags["ibrav"]=0
            if "celldm(1)" in self.namelists["SYSTEM"].tags.keys():
                del self.namelists["SYSTEM"].tags["celldm(1)"]
        if "ATOMIC_POSITIONS" in self.cards.keys():
            if poscar.coord_mode in ['alat','bohr','angstrom','crystal']:
                self.cards["ATOMIC_POSITIONS"].units=poscar.coord_mode
            elif poscar.coord_mode[0].lower()=='d':
                self.cards["ATOMIC_POSITIONS"].units="crystal"
            elif poscar.coord_mode[0].lower()=='c':
                self.cards["ATOMIC_POSITIONS"].units="angstrom"
            coords=[]
            for specie in poscar.basis:
                if specie.SD_FLAG !="":
                    flags=specie.SD_FLAG.split()
                    int_flags=[]
                    for flag in flags:
                        int_flags+=[int(flag[0]=='T' or flag[0]=='1')]
                    coords+= [(specie.occ_alias,specie.position,int_flags)]
                elif species!=None:
                    if "if_pos" in species[specie.occupant].tags.keys() and species[specie.occupant].tags["if_pos"] !="":
                        flags=species[specie.occupant].tags["if_pos"].split(',')
                        int_flags=[]
                        for flag in flags:
                            int_flags+=[int(flag[0]=='T' or flag[0]=='1')]
                        coords+= [(specie.occ_alias,specie.position,int_flags)]
                    else:
                        coords+= [(specie.occ_alias,specie.position)]
                else:
                    coords+= [(specie.occ_alias,specie.position)]
            self.cards["ATOMIC_POSITIONS"].coords=coords
            if "SYSTEM" in self.namelists.keys():
                self.namelists["SYSTEM"].tags["nat"]=len(poscar.basis)
                self.namelists["SYSTEM"].tags["ntyp"]=len(poscar.basis_dict().keys())
            if "ATOMIC_SPECIES" in self.cards.keys():
                for specie in poscar.basis_dict().keys():
                    if specie not in self.cards["ATOMIC_SPECIES"].species:
                        print("WARNING: " + specie + " in ATOMIC_POSITIONS does not have corresponding entry in ATOMIC_SPECIES")
                        if species!=None:
                            print("Attempting to retrieve data for " + specie + " from SPECIES file....")
                            for entry in species.keys():
                                if species[entry].alias == specie:
                                    self.cards["ATOMIC_SPECIES"].species+= [specie]
                                    self.cards["ATOMIC_SPECIES"].masses+= [0.0] #look up masses here????
                                    if self.namelists["CONTROL"].tags["pseudo_dir"] == ("'" + species[entry].pseudo_base + "'"):
                                        self.cards["ATOMIC_SPECIES"].pseudos+= [species[entry].pseudo_location]
                                        print("Found! " + specie + " in SPECIES file with pseudopotential at " + species[entry].pseudodir)
                                    else: 
                                        print("WARNING: PSEUDO_DIR_PATH in SPECIES file does not match 'pseudo_dir' tag in &CONTROL namelist of Infile")
                                        print("WILL ASSUME 'pseudo_dir' tag is correct")
                                        self.cards["ATOMIC_SPECIES"].pseudos+= [species[entry].pseudo_location]
                                else:
                                    print("Search failed! No data for " + specie + " to be found") 
                for specie in self.cards["ATOMIC_SPECIES"].species:
                    if specie not in poscar.basis_dict().keys():
                        index= self.cards["ATOMIC_SPECIES"].species.index(specie)
                        self.cards["ATOMIC_SPECIES"].species.pop(index)
                        self.cards["ATOMIC_SPECIES"].masses.pop(index)
                        self.cards["ATOMIC_SPECIES"].pseudos.pop(index)


    def write(self, filename):
        try:
            infile_write = open(filename,'w')
        except IOError as e:
            raise e
        for namelist in sorted(self.namelists.keys(),key=lambda x: QUANTUM_ESPRESSO_NAMELIST_LIST.index(x)):
            if self.namelists[namelist] == None or str(self.namelists[namelist]).strip() == "":
                pass
            else:
                if namelist in QUANTUM_ESPRESSO_NAMELIST_OBJ_LIST:
                    infile_write.write('&{}\n'.format(namelist))
                    infile_write.write('{}'.format(self.namelists[namelist].make_string()))
                    infile_write.write('/\n')
                else:
                    #This branch is for namelists unrecognized by the qe wrapper
                    infile_write.write('&{}\n'.format(namelist))
                    infile_write.write('{}'.format(self.namelists[namelist]))
                    infile_write.write('/\n')
        for card in sorted(self.cards.keys(),key=lambda x: QUANTUM_ESPRESSO_CARD_LIST.index(x)):
            if self.cards[card] == None or str(self.cards[card]).strip() == "":
                pass
            else:
                if card in QUANTUM_ESPRESSO_CARD_OBJ_LIST:
                    if card in ["ATOMIC_POSITIONS","K_POINTS","CELL_PARAMETERS"]:
                        infile_write.write('{} {}\n'.format(card,self.cards[card].units))
                    else:
                        infile_write.write('{}\n'.format(card))
                    infile_write.write('{}'.format(self.cards[card].make_string()))
                    infile_write.write('\n')
                else:
                    #This branch is for cards unrecognized by the qe wrapper
                    infile_write.write('{}\n'.format(card))
                    infile_write.write('{}'.format(self.cards[card]))
                    infile_write.write('\n')
        infile_write.close()

