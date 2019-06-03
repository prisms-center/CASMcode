from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import re
from casm.wrapper.misc import remove_chars

# List of tags in VASP sorted by the data type associated with it
VASP_TAG_INT_LIST = ['ialgo','ibrion','icharg','images','ismear','ispin',\
                     'istart','isym','lorbit','nbands','ndav','ngx','ngxf',\
                     'ngy','ngyf','ngz','ngzf','npar','ncore','spind','nsw',\
                     'isif', 'kpar', 'voskown', 'nsim', 'nedos', 'lmaxfock',\
                     'lmaxmix', 'nkred','ivdw','nelmin', 'nelm', 'nelmdl',\
                     'ldautype','ldauprint', 'ldauprint']
VASP_TAG_FLOAT_LIST = ['ediff','ediffg','emax','emin','encut','potim','sigma',\
                     'enmax','symprec', 'time', 'hfscreen','amix','bmix',\
                     'amix_mag', 'bmix_mag']
VASP_TAG_BOOL_LIST = ['lcharg','lsorbit','lwave','lscalapack', 'lscalu',\
                     'lplane', 'lhfcalc', 'shiftred', 'evenonly', 'oddonly',\
                     'addgrid', 'ldau', 'lasph']
# Site-wise list of arrays of FLOAT
VASP_TAG_SITEF_LIST = ['magmom','rwigs']
# Species-wise list of arrays of FLOAT
VASP_TAG_SPECF_LIST = ['ldauu', 'ldauj']
# Site-wise list of arrays of INT
VASP_TAG_SPECI_LIST = ['ldaul']
VASP_TAG_STRING_LIST = ['algo','prec','system', 'precfock','lreal']

# The master list of VASP tags is a union of the above -> need to allow for 'miscellaneous' ?
VASP_TAG_LIST = VASP_TAG_INT_LIST + VASP_TAG_SITEF_LIST + VASP_TAG_SPECI_LIST + VASP_TAG_BOOL_LIST + VASP_TAG_FLOAT_LIST + VASP_TAG_STRING_LIST + VASP_TAG_SPECF_LIST


class IncarError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg

class Incar(object):
    """
    The INCAR class contains:
        tags: a dict of all INCAR settings

    All input tags and associated values are stored as key-value pairs in the dicionary called 'tags'.
   """
    def __init__(self,filename, species=None, poscar=None, sort=True):
        """ Construct an Incar object from 'filename'"""
        self.read(filename, species, poscar, sort)

    def read(self, filename, species=None, poscar=None, sort=True):
        """ Read an INCAR file """
        self.tags = dict()
        try:
            file = open(filename,'r')
        except:
            raise IncarError("Could not open file: '" + filename + "'")

        # parse INCAR into self.tags
        for line in file:
            line = re.split('=',re.split('#',line)[0])
            if len(line) == 2:
                self.tags[line[0].strip()] = line[1].strip()
        self._verify_tags()
        self._make_natural_type()

        if species != None:
            self.update(species, poscar, sort)

        file.close()


    def _make_natural_type(self):
        """ Convert self.tags values from strings into their 'natural type' (int, float, etc.) """
        for tag in self.tags:
            if self.tags[tag] is None or str(self.tags[tag]).strip() == "":
                self.tags[tag] = None
            else:
                if tag.lower() in VASP_TAG_INT_LIST:
                    try:
                        self.tags[tag] = int(self.tags[tag])
                    except ValueError:
                        raise IncarError("Could not convert '" + tag + "' : '" + self.tags[tag] + "' to int")
                elif tag.lower() in VASP_TAG_FLOAT_LIST:
                    try:
                        self.tags[tag] = float(self.tags[tag].lower().replace('d','e'))
                    except ValueError:
                        raise IncarError("Could not convert '" + tag + "' : '" + self.tags[tag] + "' to float")
                elif tag.lower() in VASP_TAG_BOOL_LIST:
                    if not self.tags[tag].lower() in ['.true.','.false.']:
                        raise IncarError("Could not find '" + tag + "' : '" + self.tags[tag].lower() + "' in ['.true.','.false.']")
                    else:
                        self.tags[tag] = (self.tags[tag].lower() == '.true.')
                elif tag.lower() in VASP_TAG_SITEF_LIST + VASP_TAG_SPECF_LIST:
                    temp = []
                    for value in self.tags[tag].split():
                        try:
                            item=value.split('*')
                            if len(item)==1:
                                temp.append(float(value))
                            else:
                                if item[0] != 0:
                                    temp.append(str(item[0])+'*'+str(float(item[1])))
                        except ValueError:
                            raise IncarError("Could not convert '" + tag + "' : '" + self.tags[tag] + "' to float list")
                    self.tags[tag] = temp
                elif tag.lower() in VASP_TAG_SPECI_LIST:
                    temp = []
                    for value in self.tags[tag].split():
                        try:
                            temp.append(int(value))
                        except ValueError:
                            raise IncarError("Could not convert '" + tag + "' : '" + self.tags[tag] + "' to int list")
                    self.tags[tag] = temp
                elif tag.lower() in VASP_TAG_STRING_LIST:
                    self._check_string_tag(tag,self.tags[tag])


    def _check_string_tag(self,tag,value):
        """ Check that string-valued tags are allowed values """
        if tag.lower() == 'prec':
            if value.lower() not in ['low' ,'medium' ,'high' ,'normal' ,'single' ,'accurate']:
                raise IncarError("Unknown 'prec' value: '" + value)
        elif tag.lower() == 'algo':
            if value.lower() not in ['normal','veryfast','fast','conjugate','all','damped','subrot','eigenval','none','nothing','chi','gw0','gw','scgw0','scgw']:
                raise IncarError("Unknown 'algo' value: '" + value)

    def _verify_tags(self):
        """ Check that only allowed INCAR tags are in self.tags """
        for tag in self.tags:
            if tag.lower() in VASP_TAG_LIST:
                continue
            else:
                print(("Warning: unknown INCAR tag '" + tag + "' with value '" + str(self.tags[tag]) + "'"))


    def update(self, species, poscar, sort=True):
        """ Update Incar object to reflect Species settings """

        if sort == False:
            # for each 'tag' in the IndividualSpecies, create a list in self.tags
            for key in list(species.values())[0].tags.keys():
                if key.lower() in (VASP_TAG_INT_LIST + VASP_TAG_FLOAT_LIST):
                    self.tags[key] = 0.
                    for site in poscar.basis:
                        self.tags[key] += float(species[site.occupant].tags[key])
                    if key.lower() in VASP_TAG_INT_LIST:
                        self.tags[key] = int(self.tags[key])
                else:
                    self.tags[key] = []
                    if key.lower() in (VASP_TAG_SPECF_LIST + VASP_TAG_SPECI_LIST):
                        # add the value of the 'tag' for each species into the self.tags list
                        for spec in poscar.type_atoms_alias:
                            self.tags[key].append(species[spec].tags[key])
                    else:
                        # add the value of the 'tag' for each atom into the self.tags list
                        for site in poscar.basis:
                            self.tags[key].append( species[site.occupant].tags[key] )
        else:
            pos = poscar.basis_dict()
            # for each 'tag' in the IndividualSpecies, create a list in self.tags
            for key in list(species.values())[0].tags.keys():
            # for key in species[species.keys()[0]].tags.keys():
                if key.lower() in (VASP_TAG_INT_LIST + VASP_TAG_FLOAT_LIST):
                    self.tags[key] = 0.
                    for site in poscar.basis:
                        self.tags[key] += float(species[site.occupant].tags[key])
                    if key.lower() in VASP_TAG_INT_LIST:
                        self.tags[key] = int(self.tags[key])
                else:
                    self.tags[key] = []
                    # add the value of the 'tag' for each atom into the self.tags list
                    for alias in sorted(pos.keys()):
                        if key.lower() in (VASP_TAG_SPECF_LIST + VASP_TAG_SPECI_LIST):
                          # for species-specific tags, use the value specified for the
                          # species whose pseudopotential is being used for this alias
                          for name in species.keys():
                            if species[name].alias == alias and species[name].write_potcar:
                              self.tags[key].append(species[name].tags[key])
                              break
                        else:
                            for name in species.keys():
                                count=0
                                for site in pos[alias]:
                                    if site.occupant == name:
                                        count += 1
                                if species[name].alias == alias:
                                    if count > 0:
                                        self.tags[key].append( str(count) + "*" + str(species[name].tags[key]) )

    def write(self, filename):
        try:
            incar_write = open(filename,'w')
        except IOError as e:
            raise e
        for tag in self.tags:
            if self.tags[tag] is None or str(self.tags[tag]).strip() == "":
                pass
            else:
                if tag.lower() in VASP_TAG_SITEF_LIST + VASP_TAG_SPECF_LIST:
                    incar_write.write('{} = {}\n'.format(tag.upper(),remove_chars(self.tags[tag], "[\[\],']")))
                elif tag.lower() in VASP_TAG_SPECI_LIST:
                    incar_write.write('{} = {}\n'.format(tag.upper(),remove_chars(self.tags[tag], "[\[\],']")))
                elif tag.lower() in VASP_TAG_BOOL_LIST:
                    if self.tags[tag] == True:
                        incar_write.write('{} = .TRUE.\n'.format(tag.upper()))
                    else:
                        incar_write.write('{} = .FALSE.\n'.format(tag.upper()))
                else:
                    incar_write.write('{} = {}\n'.format(tag.upper(),self.tags[tag]))
        incar_write.close()

