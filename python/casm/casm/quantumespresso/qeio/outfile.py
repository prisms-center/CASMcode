from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import os, re, gzip

class OutfileError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Outfile(object):
    """Parse Outfiles.

       Currently, just contains:
           self.complete = True/False
           self.slowest_loop = float
           self.E = list of E0 values (in units of Rydberg)
    """
    def __init__(self,filename):
        self.filename = filename
        self.complete = False
        self.slowest_loop = None
        self.E=[]
        self.read()


    def read(self):
        """Parse Outfile file.  Currently just checks for completion, total cpu time, and E0."""
        if os.path.isfile(self.filename):
            if self.filename.split(".")[-1].lower() == "gz":
                f = gzip.open(self.filename)
            else:
                f = open(self.filename)
        elif os.path.isfile(self.filename+".gz"):
            f = gzip.open(self.filename+".gz")
        else:
            raise OutfileError("file not found: " + self.filename)
        prev_time=0.0
        for line in f:
            try:
                if re.search("End final coordinates",line):
                    self.complete = True
            except:
                raise OutfileError("Error reading 'End final coordinates' from line: '" + line + "'\nIn file: '" + self.filename + "'")

            try:
                if re.search("total cpu time spent up to now is",line):
                    t = float(line.split()[-2])
                    if self.slowest_loop == None:
                        self.slowest_loop = t
                    elif t - prev_time > self.slowest_loop:
                        self.slowest_loop = t - prev_time
                    prev_time=t
            except:
                pass

            try:
                if re.search("!\s*total energy",line):
                    self.E.append(float(line.split()[-2]))
            except:
                pass

        f.close()


