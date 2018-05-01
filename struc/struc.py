import numpy as np
from Parsers.structV import StrucParser
from utils.operator import StrucOperators


class Structurehandler(object):
    def __init__(self, structurefile):
        self.structure = StrucParser(structurefile)
        return

    def rotation(self, rotmatrix):
        a1 = np.dot(self.structure.matrix[0], rotmatrix[0][0]) + np.dot(self.structure.matrix[1], rotmatrix[1][0]) + np.dot(self.structure.matrix[2], rotmatrix[2][0])
        a2 = np.dot(self.structure.matrix[0], rotmatrix[0][1]) + np.dot(self.structure.matrix[1], rotmatrix[1][1]) + np.dot(self.structure.matrix[2], rotmatrix[2][1])
        a3 = np.dot(self.structure.matrix[0], rotmatrix[0][2]) + np.dot(self.structure.matrix[1], rotmatrix[1][2]) + np.dot(self.structure.matrix[2], rotmatrix[2][2])

        newbasis = np.array([a1, a2, a3])
        return newbasis

    def translation(self, transmatrix, cart=None):
        if self.structure.cart is False:
            if cart is False:
                newcoord = self.structure.coord + transmatrix
            elif cart is True:
                newcoord = self.structure.coord + np.dot(self.structure.matrix, transmatrix)

        elif self.structure.cart is True:
            if cart is True:
                newcoord = self.structure.coord + transmatrix
            elif cart is False:
                newcoord = self.structure.coord + np.dot(self.structure.matrix, np.linalg.inv(transmatrix))

        return newcoord
