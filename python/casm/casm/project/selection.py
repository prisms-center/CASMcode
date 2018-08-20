from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

# conda's current version of pandas raises these warnings, but they are safe
# see: https://stackoverflow.com/questions/40845304
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import copy
from io import StringIO
import json
import os
import subprocess

import numpy as np
import pandas
import six

from casm.project.project import Project
from casm.project.query import query
from casm.misc import compat

class Selection(object):
    """
    A Selection object contains information about a CASM project

    Attributes
    ----------

      proj: casm.Project, optional, default=Project containing the current working directory
        the CASM project the selection belongs to

      path: string, optional, default="MASTER"
        path to selection file, or "MASTER" (Default="MASTER")

      all: bool, optional, default=True
        if True, self.data will include all configurations, whether selected or
        not. If False, only selected configurations will be included.

      data: pandas.DataFrame
        A pandas.DataFrame describing the selected configurations. Has at least
        'configname' and 'selected' (as bool) columns.

    """
    def __init__(self, proj=None, path="MASTER", all=True):
        """
        Construct a CASM Project representation.

        Arguments
        ---------

          proj: casm.Project, optional, default=Project containing the current working directory
            the CASM project the selection belongs to

          path: string, optional, default="MASTER"
            path to selection file, or "MASTER" (Default="MASTER")

          all: bool, optional, default=True
            if True, self.data will include all configurations, whether selected or
            not. If False, only selected configurations will be included.


        """
        if proj == None:
          proj = Project()
        elif not isinstance(proj, Project):
          raise Exception("Error constructing Selection: proj argument is not a CASM project")
        self.proj = proj

        self.path = path
        if os.path.isfile(path):
          self.path = os.path.abspath(path)

        self.all = all

        self._data = None

        # reserved for use by casm.plotting
        self.src = None


    @property
    def data(self):
        """
        Get Selection data as a pandas.DataFrame

        If the data is modified, 'save' must be called for CASM to use the modified selection.
        """
        if self._data is None:
          if self.path in ["MASTER", "ALL", "CALCULATED"]:
            self._data = query(self.proj, ['configname', 'selected'], self, all=self.all)
          elif self._is_json():
            self._data = pandas.read_json(self.path, 'r', orient='records')
          else:
            with open(self.path, compat.pandas_rmode()) as f:
                if compat.peek(f) == '#':
                    f.read(1)
                self._data = pandas.read_csv(f, sep=compat.str(' +'), engine='python')

          self._clean_data()

          if not self.all:
            self._data = self._data[self._data['selected']==True]

        return self._data


    def save(self, data=None, force=False):
        """
        Save the current selection. Also allows completely replacing the 'data'
        describing the selected configurations.

        Args:
          data: None (default), or pandas.DataFrame describing the Selection with
            'configname' and 'selected' columns. If path=="MASTER", Configurations
            not included in 'data' will be set to not selected.
          force: Boolean, force overwrite existing files
        """
        if self.path == "MASTER":

          if data is not None:
            self._data = data.copy()
            self._clean_data()

          if self._data is None:
            return

          clist = self.proj.dir.config_list()
          backup = clist + ".tmp"
          if os.path.exists(backup):
            raise Exception("File: " + backup + " already exists")

          # read
          with open(clist, 'rb') as f:
              j = json.loads(f.read().decode('utf-8'))

          for sk, sv in six.iteritems(j["supercells"]):
            for ck, cv in six.iteritems(sv):
              sv[ck]["selected"] = False

          # set selection
          for index, row in self._data.iterrows():
            scelname, configid = row["configname"].split('/')
            j["supercells"][scelname][configid]["selected"] = row["selected"]

          # write
          with open(backup, 'wb') as f:
              f.write(six.u(json.dumps(j, indent=2)).encode('utf-8'))
          os.rename(backup, clist)

          # refresh proj config list
          self.proj.refresh(read_configs=True)

        elif self.path in ["ALL", "CALCULATED"]:
          raise Exception("Cannot save the '" + self.path + "' Selection")

        else:

          if data is not None:
            self._data = data.copy()
            self._clean_data()

          if os.path.exists(self.path) and not force:
            raise Exception("File: " + self.path + " already exists")

          backup = self.path + ".tmp"
          if os.path.exists(backup):
            raise Exception("File: " + backup + " already exists")

          if self._is_json():
            self._data.to_json(backup, orient='records')
          else:
            self.data.loc[:,"selected"] = self.data.loc[:,"selected"].astype(np.int_)
            with open(backup, compat.pandas_wmode()) as f:
                f.write('# ')  # will make this optional in a future version
                self._data.to_csv(f, sep=compat.str(' '), index=False)
          os.rename(backup, self.path)


    def saveas(self, path, force=False):
        """
        Create a new Selection from this one, save and return it

        Args:
          path: path to selection file (Default="MASTER")
          force: Boolean, force overwrite existing files

        Returns:
          sel: the new Selection created from this one
        """
        sel = Selection(self.proj, path, all=self.all)
        sel._data = self.data.copy()
        sel.save(force=force)
        return sel


    def _is_json(self):
        return self.path[-5:].lower() == ".json"


    def _clean_data(self):
        self._data.loc[:,'selected'] = self._data.loc[:,'selected'].astype(bool)


    def query(self, columns, force=False, verbose=False):
        """
        Query requested columns and store them in 'data'. Will not overwrite
        columns that already exist, unless 'force'==True.

        Will query data for all configurations, whether selected or not, if
        self.all == True.
        """

        if force == False:
          _col = [x for x in columns if x not in self.data.columns]
        else:
          _col = columns

        if verbose:
          print("# Query requested:", columns)
          if force == False:
            print("# Use existing:", [x for x in columns if x in self.data.columns])
          else:
            print("# Overwrite existing:", [x for x in columns if x in self.data.columns])
          if len(_col) == 0:
            print("# No query necessary")
          else:
            print("# Querying:", _col)

        if len(_col) == 0:
          return

        df = query(self.proj, _col, self, all=self.all)

        if verbose:
          print("#   DONE\n")

        msg = "querying different numbers of records: {0}, {1}".format(
          self.data.shape, df.shape)
        assert self.data.shape[0] == df.shape[0], msg

        for c in df.columns:
          self.data.loc[:,c] = df.loc[:,c].values


    def write_pos(self, all=False):
        """
        Write POS file for configurations

        Arguments
        ---------

          all: bool, optional, default=False
            if True, will write POS file for all configurations in the selection
            whether selected or not. If False, only write POS file for selected
            configurations.
        """
        self.proj.capture("query -c " + self.path + " --write-pos")


    def add_data(self, name, data=None, force=False):
        """
        Equivalent to:
        if name not in sel.data.columns or force == True:
          if data is None:
            sel.query([name], force)
          else:
            sel.data.loc[:,name] = data
        """
        if name not in self.data.columns or force == True:
          if data is None:
            self.query([name], force)
          else:
            self.data.loc[:,name] = data
