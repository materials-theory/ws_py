import copy
import numpy as np
from vw_py.Parsers.outcarhandler import OutcarHandler
from vw_py.Parsers.structV import ContcarHandler
from vw_py.IO.IO import IO


class DosParserV(object):
    """
    Class object to parse DOSCAR file from VASP.

    """

    def __init__(self, filename='DOSCAR'):
        self.io = IO(filename)
        self.cont = ContcarHandler()
        self.outcar = OutcarHandler()

        self.numtdos = None
        self.numpdos = None
        self.numsumdos = None
        self.tdos = None
        self.pdos = None
        self.sumdos = None
        self.emax = None
        self.emin = None
        self.nedos = None
        self.atoms = None
        self.numpdos_add = 1

        self.ispin = int(self.outcar.param_from_outcar('ISPIN'))
        self.lorbit = int(self.outcar.param_from_outcar('LORBIT'))
        self.lsorbit = str(self.outcar.param_from_outcar('LSORBIT'))

        self.dosline()
        self.parser()
        return

    def dosline(self):

        orbit_mult = None
        # mag_mult = None

        tdos_avail = [3, 5]
        pdos_avail = [4, 7, 10, 13, 19, 37]

        # tDOS
        # ISPIN = 2, LSORBIT = F
        # 5 - E tu td inu ind
        # Others
        # 3 - E t in
        if self.ispin == 2 and self.lsorbit == "F":
            self.numtdos = 5
            self.numpdos_add = 2
        else:
            self.numtdos = 3

        # pDOS
        # ISPIN=1, LORBIT=10
        # 4  - E s p d
        # ISPIN=2, LORBIT=10
        # 7  - E su sd pu pd du dd
        # ISPIN=1, LORBIT=11,12
        # 10 - E s px py pz dxy dyz dz2 dxz dx2
        # ISPIN=1, LORBIT=10, LSORBIT=T
        # 13 - E st smx smy smz pt pmx pmy pmz dmt pmx pmy pmz

        # ISPIN=2, LORBIT=11,12
        # 19 - E su sd pxu pxd pyu pyd pzu pzd dxyu dxyd dyzu dyzd dz2u dz2d dxzu dxzd dx2u dx2d (????)
        # ISPIN=1,2, LORBIT=11,12, LSORBIT=T
        # 37 - E ???
        if self.lorbit == 10:
            orbit_mult = 3
        elif self.lorbit == 11 or self.lorbit == 12:
            orbit_mult = 9

        if self.lsorbit == "T":
            mag_mult = 4
        else:
            mag_mult = 1

        self.numpdos = (self.ispin * orbit_mult * mag_mult) + 1

        if self.numtdos not in tdos_avail or self.numpdos not in pdos_avail:
            raise NotImplementedError("Can't parse this DOSCAR file!")
        else:
            pass

        return

    def parser(self):
        print("Reading VASP DOSCAR file...")
        filestr = self.io.ReadFile()

        print("Parsing VASP DOSCAR file...")
        dos = filestr.readlines()

        # Reading the header part, only number of kps and number of bands
        # Other header parts are removed
        self.emax = float(dos[5].split()[0])
        self.emin = float(dos[5].split()[1])
        self.nedos = int(dos[5].split()[2])

        for i in range(5):
            del(dos[0])

        tdos_array = []
        pdos_array = []
        sumdos_array = []

        # Writing tDOS part to an array
        del dos[0]
        for i in range(self.nedos):
            tdos_array.append(dos.pop(0).split())
        tdos_array = np.reshape(np.array(tdos_array, dtype='d'), (self.nedos, self.numtdos))

        # Multiplying -1 to down-spin tdos
        if self.numpdos_add == 2:
            tdos_array[:, 2] = tdos_array[:, 2] * -1
            tdos_array[:, 4] = tdos_array[:, 4] * -1

        # Writing pDOS part to an array
        self.atoms = self.cont.atominfo()
        totalcount = 0
        for x in self.atoms.values():
            totalcount += int(x)

        for i in range(totalcount):
            del dos[0]
            tmp_array = []
            for j in range(self.nedos):
                tmp_array.append(dos.pop(0).split())
            tmp_array = np.array(tmp_array, dtype='d')

            if self.numpdos_add == 2:
                up = tmp_array[:, 1::2].sum(1)
                down = tmp_array[:, 2::2].sum(1)
                tmp_array = np.column_stack((tmp_array, up, down))
            else:
                sum = tmp_array[:, 1:].sum(1)
                tmp_array = np.column_stack((tmp_array, sum))

            pdos_array.append(tmp_array)

        # Multiplying -1 to down-spin pdos
        if self.numpdos_add == 2:
            for i in range(totalcount):
                for j in range(self.numpdos + self.numpdos_add):
                    if j == 0:
                        pass
                    else:
                        if np.mod(j, 2) == 0:
                            pdos_array[i][:, j] = pdos_array[i][:, j] * -1
                        else:
                            pass

        # Summing up pDOS to generate sumDOS array
        count = 0
        self.numsumdos = len(self.atoms)
        for x in self.atoms.keys():
            tmp_array = []
            for y in range(int(self.atoms[x])):
                if y == 0:
                    tmp_array = copy.deepcopy(pdos_array[y + count])
                else:
                    tmp_array += pdos_array[y + count]

            tmp_array[:, 0] = tmp_array[:, 0] / int(self.atoms[x])
            count += int(self.atoms[x])
            sumdos_array.append(tmp_array)

        sumdos_array = np.reshape(sumdos_array, (self.numsumdos, self.nedos, self.numpdos + self.numpdos_add))

        self.tdos = tdos_array
        self.pdos = pdos_array
        self.sumdos = sumdos_array

        print("Done!")
        return
