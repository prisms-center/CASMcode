import os, re, gzip

class OutcarError(Exception):
    def __init__(self,msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Outcar(object):
    """Parse OUTCAR files.

       Currently, just contains:
           self.complete = True/False
           self.slowest_loop = float
           self.kpoints = list of int, or none
    """
    def __init__(self,filename):
        self.filename = filename
        self.complete = False
        self.slowest_loop = None
        self.read()


    def read(self):
        """Parse OUTCAR file.  Currently just checks for completion, loop time, and kpoints."""
        self.kpts = None
        if os.path.isfile(self.filename):
            if self.filename.split(".")[-1].lower() == "gz":
                f = gzip.open(self.filename)
            else:
                f = open(self.filename)
        elif os.path.isfile(self.filename+".gz"):
            f = gzip.open(self.filename+".gz")
        else:
            raise OutcarError("file not found: " + self.filename)

        for line in f:
            try:
               if re.search("generate k-points for:",line):
                   self.kpts = map(int, line.split()[-3:])
            except:
               pass

            try:
                if re.search("Total CPU time used",line):
                    self.complete = True
            except:
                raise OutcarError("Error reading 'Total CPU time used' from line: '" + line + "'\nIn file: '" + self.filename + "'")

            try:
                if re.search("LOOP",line):
                    t = float(line.split()[-1])
                    if self.slowest_loop is None:
                        self.slowest_loop = t
                    elif t > self.slowest_loop:
                        self.slowest_loop = t
            except:
                pass

        f.close()


