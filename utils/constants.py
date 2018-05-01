import scipy.constants as constants
from ase import units


class Unitconverter(object):
    def __init__(self):
        return

    @staticmethod
    def conversion_factor(unitin, unitout):
        conv = {'Ry': units.Ry, 'eV': 1.0, 'J': units._e}
        try:
            return conv[unitin] / conv[unitout]
        except KeyError:
            print('Error: Unknown unit.')


class ScientificConstants(object):
    def __init__(self):
        return

    # def getunit(self, ):
    #
    #     return unit