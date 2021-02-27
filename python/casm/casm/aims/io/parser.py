import os
import gzip


class ParserError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Parser(object):
    """Parse FHI-aims output files.

       Currently, just contains:
           self.complete = True/False
           self.kpoints = list of int, or none
    """
    def __init__(self, filename):
        self.complete = False
        self.read(filename)

    def read(self, file):
        if os.path.isfile(file):
            if file.split(".")[-1].lower() == "gz":
                f = gzip.open(file)
            else:
                f = open(file)
        elif os.path.isfile(file + ".gz"):
            f = gzip.open(file + ".gz")
        else:
            raise ParserError("file not found: " + file)

        for line in f:
            try:
                if "Have a nice day." in line:
                    self.complete = True
                    return True
            except IOError:
                raise ParserError("Error reading 'Have a nice day.' from line: '" + line + "'\n"
                                  "NOT CONVERGED ERROR In file: '" + file + "'")

        f.close()
