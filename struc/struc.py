import numpy as np
from Parsers.structV import StrucParser
from generalutils.operator import StrucOperators


# TODO: Make module consider PBC and repeat atoms based on spg
class Structurehandler(object):
    def __init__(self, structurefile):
        self.structure = StrucParser(structurefile)
        return

    def rotation(self, rotmatrix):
        a1 = np.dot(self.structure.matrix[0], rotmatrix[0][0]) + np.dot(self.structure.matrix[1], rotmatrix[1][0]) \
             + np.dot(self.structure.matrix[2], rotmatrix[2][0])
        a2 = np.dot(self.structure.matrix[0], rotmatrix[0][1]) + np.dot(self.structure.matrix[1], rotmatrix[1][1]) \
             + np.dot(self.structure.matrix[2], rotmatrix[2][1])
        a3 = np.dot(self.structure.matrix[0], rotmatrix[0][2]) + np.dot(self.structure.matrix[1], rotmatrix[1][2]) \
             + np.dot(self.structure.matrix[2], rotmatrix[2][2])

        newbasis = np.array([a1, a2, a3])
        return newbasis

    def translation(self, transmatrix, cart=None):
        if self.structure.cart is False:
            if cart is False:
                self.structure.coord = self.structure.coord + transmatrix
            elif cart is True:
                self.structure.coord = self.structure.coord + np.dot(self.structure.matrix, transmatrix)

        elif self.structure.cart is True:
            if cart is True:
                self.structure.coord = self.structure.coord + transmatrix
            elif cart is False:
                self.structure.coord = self.structure.coord + np.dot(self.structure.matrix, np.linalg.inv(transmatrix))

        return

    def cartdirconvert(self):
        if self.structure.cart is False:
            self.structure.coord = StrucOperators.directtocartesian(self.structure.matrix, self.structure.coord)
            self.structure.cart = True

        else:
            self.structure.coord = StrucOperators.cartesiantodirect(self.structure.matrix, self.structure.coord)
            self.structure.cart = False

        return

    # TODO: Test working on non-orthogonal unit cell
    def find_current_vac(self, direction, tolerance):
        direction = direction - 1 # x = 1, y = 2, z = 3
        nondir = [0, 1, 2]
        # nondir = del(nondir[direction])

        if self.structure.cart is False:
            self.cartdirconvert()

        tmp_shift = np.zeros(3)

        if tolerance is None:
            tmp_shift[direction] = 3.0
        else:
            tmp_shift[direction] = tolerance

        coord = self.structure.coord + tmp_shift

        for x in coord:
            tmp = 0
            for y in nondir:
                tmp += x[y] * self.structure.matrix[y]

            tmp += x[direction]

            # if tmp >=

            if x[direction] >= 1.0:
                x[direction] -= 1.0
            else:
                pass

        top_atom = coord[0]
        bottom_atom = coord[0]

        # for x in coord:
        #     if x[direction] >= top_atom[direction]:
        #         top_atom = x
        #     elif x[direction] <= bottom_atom[direction]:
        #         bottom_atom = x
        #     else:
        #         pass

        for x in coord:
            if np.dot(x, self.structure.matrix)[direction] >= top_atom[direction]:
                top_atom = x
            elif x[direction] <= bottom_atom[direction]:
                bottom_atom = x
            else:
                pass
        # tmp_shift = np.zeros(3)
        #
        # direction = direction - 1 # x = 1, y = 2, z = 3
        #
        # if tolerance is None:
        #     if self.structure.cart is True:
        #         tmp_shift[direction] = 3.0
        #     else:
        #         tmp_shift[direction] = 0.1
        # else:
        #     tmp_shift[direction] = tolerance
        #
        # # tmp_shift = np.dot(tmp_shift, self.structure.unitvec)
        # # shift_norm = np.linalg.norm(tmp_shift)
        #
        # coord = self.structure.coord + tmp_shift
        #
        # for x in coord:
        #     if x[direction] >= 1.0:
        #         x[direction] -= 1.0
        #     else:
        #         pass
        #
        # top_atom = coord[0]
        # bottom_atom = coord[0]
        #
        # # for x in coord:
        # #     if x[direction] >= top_atom[direction]:
        # #         top_atom = x
        # #     elif x[direction] <= bottom_atom[direction]:
        # #         bottom_atom = x
        # #     else:
        # #         pass
        #
        # for x in coord:
        #     if np.dot(x, self.structure.matrix)[direction] >= top_atom[direction]:
        #         top_atom = x
        #     elif x[direction] <= bottom_atom[direction]:
        #         bottom_atom = x
        #     else:
        #         pass

        return

    def addvec(self, direction, vecheight, tolerance):
        return

    def as_dict(self):
        dic = {"index": self.structure.fileindex,
               "unitvec": self.structure.unitvec,
               "atoms": self.structure.atominfo,
               "cart": self.structure.cart,
               "seldyn": self.structure.seldyn,
               "coord": self.structure.coord,
               "dyn": self.structure.dyn,
               "vel": self.structure.vel,
               }

        return dic
