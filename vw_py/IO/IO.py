import os

class IO(object):
    def __init__(self, infile=None, outfile=None):
        self.infile = infile
        self.outfile = outfile
        if infile is not None:
            self.setInFile(infile)
        if outfile is not None:
            self.setOutFile(outfile)
        return

    def setInFile(self, infile):
        # Setter method for input file
        self.infile = infile
        return

    def setOutFile(self, outfile):
        # Setter method for output file
        self.outfile = outfile
        return

    def ReadFile(self):
        Filename = self.infile
        if os.path.isfile(Filename):
            inFile = open(Filename, 'r')
        elif Filename is None:
            raise IOError("File not set")
        else:
            raise IOError("No %s file found!" % Filename)
        return inFile

    def WriteFile(self):
        Filename = self.outfile
        if Filename is not None:
            outFile = open(Filename, 'w')
        else:
            raise IOError("File not set")
        return outFile