import numpy as np

from Parsers.outcarhandler import OutcarHandler as outcar


class ProcarParserV:
    """
    Class object to parse PROCAR file from VASP.

    """
    def __init__(self, filename=None, spin=False):
        self.numkp = None
        self.numstates = None
        self.numions = None
        self.states = None
        self.path = None
        self.wgt = None
        self.proj = None
        self.proj_spin = None
        self.parsed = False
        self.spin = spin

        if filename is None:
            self.filename = 'PROCAR'
        else:
            self.filename = filename

        return

    def parser(self, spin):
        print("Reading VASP PROCAR file...")
        with open(self.filename, 'r') as filestr:
            print("Parsing VASP PROCAR file...")
            eig = filestr.readlines()

            # Reading the header part and remove
            self.numkp = int(eig[1].split()[3])
            self.numstates = int(eig[1].split()[7])
            self.numions = int(eig[1].split()[11])

            for i in range(3):
                del(eig[0])

            states_array = []
            path_array = []
            wgt_array = []
            proj_array = []

            # Parse the kp path and band data, then appends to the numpy array
            if spin is False:
                if outcar.param_from_outcar('LSORBIT') == 'F':
                    for i in range(self.numkp):
                        if len(eig[i * ((self.numions + 5) * self.numstates + 3)].split()) != 9:
                            point = []
                            for x in eig[i * ((self.numions + 5) * self.numstates + 3)].split(':')[1].split()[0:-3]:
                                point.append(x.replace("-", " "))
                            path_array.append(" ".join(point).split())
                            wgt_array.append(eig[i * ((self.numions + 5)
                                                      * self.numstates + 3)].split(':')[1].split()[-1])
                            # wgt_array.append(eig[i * (((self.numions + 1) * 4 + 4)
                            #                           * self.numstates + 3)].split(':')[1].split()[-1])

                        else:
                            path_array.append(eig[i * ((self.numions + 5) * self.numstates + 3)].split()[3:6])
                            wgt_array.append(eig[i * ((self.numions + 5) * self.numstates + 3)].split()[8])
                        for j in range(self.numstates):
                            states_array.append(eig[(i * ((self.numions + 5) * self.numstates + 3)
                                                     + 2 + j * (self.numions + 5))].split()[4])
                            for k in range(self.numions + 1):
                                proj_array.append(eig[(i * ((self.numions + 5) * self.numstates + 3)
                                                       + 2 + j * (self.numions + 5)) + 3 + k].split()[:])

                elif outcar.param_from_outcar('LSORBIT') == 'T':
                    for i in range(self.numkp):
                        if len(eig[i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)].split()) != 9:
                            point = []
                            for x in eig[i * (((self.numions + 1) * 4 + 4)
                                           * self.numstates + 3)].split(':')[1].split()[0:-3]:
                                point.append(x.replace("-", " "))
                            path_array.append(" ".join(point).split())
                            wgt_array.append(eig[i * (((self.numions + 1) * 4 + 4)
                                                      * self.numstates + 3)].split(':')[1].split()[-1])

                        else:
                            path_array.append(eig[i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)].split()[3:6])
                            wgt_array.append(eig[i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)].split()[8])

                        for j in range(self.numstates):
                            states_array.append(eig[(i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)
                                                     + 2 + j * ((self.numions + 1) * 4 + 4))].split()[4])
                            for k in range(self.numions + 1):
                                proj_array.append(eig[(i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)
                                                     + 2 + j * ((self.numions + 1) * 4 + 4)) + 3 + k].split()[:])

                # Unifying the array data type to the float, and then reshape arrays to an appropriate shape
                states_array = np.array(states_array, dtype='d')
                path_array = np.array(path_array, dtype='d')
                wgt_array = np.array(wgt_array, dtype='d')
                proj_array = np.array(proj_array)

                states_array = np.reshape(states_array, (self.numkp, self.numstates))
                path_array = np.reshape(path_array, (self.numkp, 3))
                wgt_array = np.reshape(wgt_array, (self.numkp, 1))
                proj_array = np.reshape(proj_array, (self.numkp, self.numstates, self.numions + 1, 5))

            if spin is True:
                proj_tot = []
                proj_mx = []
                proj_my = []
                proj_mz = []

                for i in range(self.numkp):
                    if len(eig[i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)].split()) != 9:
                        point = []
                        for x in eig[i * (((self.numions + 1) * 4 + 4)
                                              * self.numstates + 3)].split(':')[1].split()[0:-3]:
                            point.append(x.replace("-", " "))
                        path_array.append(" ".join(point).split())
                        wgt_array.append(eig[i * (((self.numions + 1) * 4 + 4)
                                               * self.numstates + 3)].split(':')[1].split()[-1])

                    else:
                        path_array.append(eig[i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)].split()[3:6])
                        wgt_array.append(eig[i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)].split()[8])
                    for j in range(self.numstates):
                        states_array.append(eig[(i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)
                                                 + 2 + j * ((self.numions + 1) * 4 + 4))].split()[4])
                        for k in range(self.numions + 1):
                            proj_tot.append(eig[(i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)
                                                   + 2 + j * ((self.numions + 1) * 4 + 4)) + 3 + k].split()[:])
                            proj_mx.append(eig[(i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)
                                                 + 2 + j * ((self.numions + 1) * 4 + 4)) + 3 + k
                                                 + (self.numions + 1)].split()[:])
                            proj_my.append(eig[(i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)
                                                 + 2 + j * ((self.numions + 1) * 4 + 4)) + 3 + k
                                                 + (self.numions + 1) * 2].split()[:])
                            proj_mz.append(eig[(i * (((self.numions + 1) * 4 + 4) * self.numstates + 3)
                                                 + 2 + j * ((self.numions + 1) * 4 + 4)) + 3 + k
                                                 + (self.numions + 1) * 3].split()[:])

                states_array = np.array(states_array, dtype='d')
                path_array = np.array(path_array, dtype='d')
                wgt_array = np.array(wgt_array, dtype='d')
                proj_array = np.array(proj_tot)
                proj_spin_array = np.array([proj_mx, proj_my, proj_mz])

                states_array = np.reshape(states_array, (self.numkp, self.numstates))
                path_array = np.reshape(path_array, (self.numkp, 3))
                wgt_array = np.reshape(wgt_array, (self.numkp, 1))
                proj_array = np.reshape(proj_array, (self.numkp, self.numstates, self.numions + 1, 5))
                proj_spin_array = np.reshape(proj_spin_array, (3, self.numkp, self.numstates, self.numions + 1, 5))

                self.proj_spin = proj_spin_array

            self.path = path_array
            self.states = states_array
            self.wgt = wgt_array
            self.proj = proj_array
            self.parsed = True
            print("Done!")
        return

    def as_dict(self):
        if self.parsed is False:
            self.parser(self.spin)
        else:
            pass

        dic = {'path': self.path,
               'states': self.states,
               'weight': self.wgt,
               'proj': self.proj,
               'proj_spin': self.proj_spin,
               'numbands': self.numstates,
               'numions': self.numions,
               'numkps': self.numkp
               }

        return dic
