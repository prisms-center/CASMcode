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
        self.ngx = None
        self.ngy = None
        self.ngz = None
        self.found_ngx = False
        self.forces = []

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

        force_index = False
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

            try:
                r = re.match(r"\s*dimension x,y,z\s*NGX\s*=\s*([0-9]*)\s*NGY\s*=\s*([0-9]*)\s*NGZ\s*=\s*([0-9]*)\s*", line)
                # if r:
                if r and not self.found_ngx:
                    self.ngx = int(r.group(1))
                    self.ngy = int(r.group(2))
                    self.ngz = int(r.group(3))
                    self.found_ngx = True
            except:
                pass

            if force_index:
                if '--' in line and len(self.forces) > 0:
                    force_index = False
                elif '--' not in line:
                    self.forces.append(map(float, line.split()[-3:]))

            try:
                if re.search("TOTAL-FORCE", line):
                    force_index = True
            except:
                pass

        f.close()


