import numpy as np
from collections import OrderedDict
from vw_py.generalutils.operator import StrucOperators
from vw_py.IO.IO import IO


class ContcarHandler(object):
    @staticmethod
    def unitvec():
        with open('CONTCAR', 'r') as cont:
            f = cont.readlines()
            unit_vec = np.array([])
            factor = float(f[1].split()[0])
            for i in range(2, 5):
                unit_vec = np.append(unit_vec, f[i].split())
            unit_vec = np.reshape(np.array(unit_vec, dtype='d'), (3, 3)) * float(factor)
            return unit_vec

    @staticmethod
    def atominfo():
        with open('CONTCAR', 'r') as cont:
            dic = OrderedDict()
            lines = cont.readlines()
            for i in range(len(lines[5].split())):
                dic[lines[5].split()[i]] = lines[6].split()[i]
            return dic


class StrucParser(object):
    """
    Parsing the opened / imported VASP structure files, into one list.
    Parsed data includes every information stored in VASP structure files.

    Members:
    - Parser

    # Output format:
    # 1. file index (str)
    # 2. unit vector (np array)
    # 3. atom info (ordered dic)
    # 4. Cartesian coord? (bool)
    # 5. selective dyn (bool)
    # 6. coord (np array)
    # 7. dyn (np array)
    # 8. vel (np array)
    # 9. atom list (np array)
    """

    def __init__(self, filename=None):
        self.fileindex = None
        self.unitvec = None
        self.atominfo = None
        self.cart = None
        self.seldyn = None
        self.coord = None
        self.dyn = None
        self.vel = None
        self.atomlist = None
        self.matrix = None
        self.angles = None
        self.length = None

        # self.parsed = None
        self.utils = StrucOperators()

        if filename is not None:
            pass
        else:
            filename = 'POSCAR'

        print("Reading " + str(filename) + "...")
        self.io = IO(filename)
        self.parser()
        self.cellmatrix()
        return

    def parser(self):
        # print("Reading VASP structure file...")
        f = self.io.ReadFile()

        # print("Parsing VASP structure file...")
        self.fileindex = f.readline().strip(' \n')
        scalfac = float(f.readline().strip(' \n'))

        unitvec = []
        for i in range(3):
            unitvec.append(f.readline().split())
        self.unitvec = np.array(unitvec, dtype='d') * scalfac

        atom_species = f.readline().split()
        atom_counts = f.readline().split()
        atominfo = OrderedDict()
        for i in range(len(atom_species)):
            atominfo[atom_species[i]] = int(atom_counts[i])
        self.atominfo = atominfo
        count = 0
        for x in atominfo:
            count += atominfo[x]

        determ_line = f.readline()
        if str(determ_line[0]) == "S":
            self.seldyn = True
            coordtype = f.readline().strip(' \n')
        else:
            self.seldyn = False
            coordtype = determ_line.strip(' \n')

        if coordtype[0] == 'D':
            self.cart = False
        elif coordtype[0] == 'd':
            self.cart = False
        elif coordtype[0] == 'C':
            self.cart = True
        elif coordtype[0] == 'c':
            self.cart = True

        tmp = f.read()

        coord = []
        dyn = []
        vel = []

        if not self.seldyn:
            tmp = np.array(tmp.split(), dtype='d')
            tmp = np.reshape(tmp, (int(len(tmp) / 3), 3))
            self.dyn = np.array([], dtype='str')
            for i in range(count):
                coord.append(tmp[i])
                if len(tmp) == count * 2:
                    vel.append(tmp[i + count])
            self.coord = np.array(coord)
            self.vel = np.array(vel)

        elif self.seldyn:
            tmp = np.array(tmp.split())
            tmpc = tmp[:(count * 6)]
            if len(tmp) == count * 9:
                vel = np.reshape(tmp[(count * 6):], (count, 3))
            self.vel = np.array(vel, dtype='d')
            for i in range(count):
                for j in range(3):
                    coord.append(tmpc[i * 6 + j])
                    dyn.append(tmpc[i * 6 + j + 3])
            self.coord = np.reshape(np.array(coord, dtype='d'), (count, 3))
            self.dyn = np.reshape(np.array(dyn, dtype='str'), (count, 3))

        atomlist = []
        for i in atominfo:
            for j in range(atominfo[i]):
                atomlist.append(i)
        self.atomlist = np.vstack(atomlist)

        return

    def cellmatrix(self):
        alpha = self.utils.angle(self.unitvec[1], self.unitvec[2])
        beta = self.utils.angle(self.unitvec[2], self.unitvec[0])
        gamma = self.utils.angle(self.unitvec[0], self.unitvec[1])

        ang = [alpha, beta, gamma]
        self.angles = np.array(ang)

        a = self.utils.vectorlength(self.unitvec[0])
        b = self.utils.vectorlength(self.unitvec[1])
        c = self.utils.vectorlength(self.unitvec[2])

        # ca = math.cos(math.radians(alpha))
        # cb = math.cos(math.radians(beta))
        # cg = math.cos(math.radians(gamma))
        # sg = math.sin(math.radians(gamma))
        #
        # vol = a * b * c * math.sqrt(1 - ((ca ** 2) + (cb ** 2) + (cg ** 2)) + (2 * ca * cb * cg))

        # mat = [[a, (b * cg), (c * cb)], [0, (b * sg), (c * (ca - (cb * cg)) / sg)], [0, 0, (vol / (a * b * sg))]]
        mat = self.unitvec
        self.matrix = np.array(mat)
        self.length = np.array([a, b, c])
        return
