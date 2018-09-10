import numpy as np
from vw_py.IO.IO import IO
from vw_py.Parsers.outcarhandler import OutcarHandler

class EigenParserV:
    """
    Class object to parse EIGENVAL file from VASP.

    """
    def __init__(self, filename='EIGENVAL'):
        self.io = IO(filename)
        self.outcar = OutcarHandler()
        self.numkp = None
        self.numstates = None
        self.states = None
        self.path = None
        self.wgt = None

        self.Parser()
        return

    def Parser(self):
        print("Reading VASP EIGENVAL file...")
        fileStr = self.io.ReadFile()

        print("Parsing VASP EIGENVAL file...")
        eig = fileStr.readlines()

        # Reading the header part, only number of kps and number of bands
        # Other header parts are removed
        self.numkp = int(eig[5].split()[1])
        self.numstates = int(eig[5].split()[2])
        num_lines = self.numstates + 2

        for i in range(7):
            del(eig[0])

        states_array = []
        path_array = []
        wgt_array = []

        # Parse the kp path and band data, then appends to the numpy array
        for i in range(self.numkp):
            path_array.append(eig[i * num_lines].split()[0:3])
            wgt_array.append(eig[i * num_lines].split()[3])
            for j in range(self.numstates):
                states_array.append(eig[i * num_lines + j + 1].split()[1])

        # Unifying the array data type to the float, and then reshape arrays to an appropriate shape
        states_array = np.array(states_array, dtype='d')
        path_array = np.array(path_array, dtype='d')
        wgt_array = np.array(wgt_array, dtype='d')

        states_array = np.reshape(states_array, (self.numkp, self.numstates))
        path_array = np.reshape(path_array, (self.numkp, 3))
        wgt_array = np.reshape(wgt_array, (self.numkp, 1))

        self.path = path_array
        self.states = states_array
        self.wgt = wgt_array

        print("Done!")
        return