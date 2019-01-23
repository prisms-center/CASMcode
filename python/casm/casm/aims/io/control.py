class ControlError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class Control:
    """
    The Control class contains:
        tags: a dict of all settings
   """

    def __init__(self, filename):
        """ Construct a Control object from 'filename'"""

        self.ctrl_content = []

        try:
            f = open(filename, 'r')
        except IOError:
            raise ControlError("Could not open file: '" + filename + "'")

        # read control.in as is
        for line in f:
            self.ctrl_content.append(line.strip())
        f.close()

    def write(self, filename):
        try:
            outfile = open(filename, 'w')
        except IOError as e:
            raise e
        for line in self.ctrl_content:
            line += '\n'
            outfile.write(line)

        outfile.close()
