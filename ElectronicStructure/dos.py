import copy
from Parsers.outcarhandler import OutcarHandler
from Parsers.doscarV import DosParserV


class Densityofstates(object):

    def __init__(self):
        self.outcar = OutcarHandler()
        self.doscar = DosParserV()
        self.tdos = None
        self.pdos = None
        self.sumdos = None
        self.numtdos = None
        self.numpdos = None
        self.numsumdos = None
        self.atominfo = None
        self.nedos = None

        return

    def dos(self, fermi):
        tdos = copy.deepcopy(self.doscar.tdos)
        pdos = copy.deepcopy(self.doscar.pdos)
        sumdos = copy.deepcopy(self.doscar.sumdos)

        if fermi == 0.0:
            fermi_e = float(self.outcar.param_from_outcar('E-fermi'))
        else:
            fermi_e = fermi

        # Shifting for the fermi level
        for i in range(len(tdos)):
            tdos[i, 0] -= fermi_e

        for i in range(len(pdos)):
            for j in range(len(pdos[i])):
                pdos[i][j, 0] -= fermi_e

        for i in range(len(sumdos)):
            for j in range(len(sumdos[i])):
                sumdos[i][j, 0] -= fermi_e

        self.tdos = tdos
        self.pdos = pdos
        self.sumdos = sumdos
        self.numtdos = self.doscar.numtdos
        self.numpdos = self.doscar.numpdos
        self.numsumdos = self.doscar.numsumdos
        self.atominfo = self.doscar.atoms
        self.nedos = self.doscar.nedos

        return
