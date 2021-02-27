import os
import gzip


class StepsError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Steps(object):
    """Parse run file for ionic steps

       Currently, just contains:
           self.E = list of ionic step Total energy values
    """
    def __init__(self, filename):
        self.filename = filename
        self.E = []
        self.read()

    def read(self):
        """Parse file. Currently just collects E0 as list 'self.E'"""
        if os.path.isfile(self.filename):
            if self.filename.split(".")[-1].lower() == "gz":
                f = gzip.open(self.filename)
            else:
                f = gzip.open(self.filename + ".gz")
        else:
            raise StepsError("file not found: " + self.filename)

        self.E = []
        for line in f:
            if b"| Total energy                  :" in line:
                self.E.append(float(line.split()[6]))
        f.close()
