import numpy as np
import itertools
import math

from utils.operator import StrucOperators
from Parsers.structV import StrucParser


class PolyhedronV:

    def __init__(self, structurefile=None):

        self.poly_level = None
        self.bondmax = None
        self.able_poly = [4, 6, 8, 12, 20]
        self.neighbors = None
        self.ideal_angle = None
        self.vertices = None
        # self.center_atom = None
        # self.vertex_atom = None
        self.center_atom = []
        self.vertex_atom = []
        # self.poly_pos = None

        # Assigning other classes
        self.utils = StrucOperators()
        self.struc = StrucParser(structurefile)
        return

    def determpolylevel(self, polylevel):
        # Checks the level of polyhedron and determine whether it is a Platonic solids
        if polylevel not in self.able_poly:
            raise IOError("Only Platonic solids can be calculated!")
        else:
            self.poly_level = polylevel

        if self.poly_level == 4:
            self.vertices = 4
            self.ideal_angle = 109.47
            self.neighbors = 3
        elif self.poly_level == 6:
            self.vertices = 8
            self.ideal_angle = 70.53
            self.neighbors = 3
        elif self.poly_level == 8:
            self.vertices = 6
            self.ideal_angle = 90.0
            self.neighbors = 4
        elif self.poly_level == 12:
            self.vertices = 20
            self.ideal_angle = 53.55
            self.neighbors = 3
        elif self.poly_level == 20:
            self.vertices = 12
            self.ideal_angle = 64.27
            self.neighbors = 5
        else:
            raise IOError("Only Platonic solids can be calculated!")

        return

    def setpoly(self, polylevel, bondmax=3.0, centeratom=None, vertexatom=None):
        # Setter method to store polyhedron informations
        self.determpolylevel(polylevel)
        self.bondmax = bondmax

        if vertexatom is None or len(vertexatom) == 0:
            raise IOError("There must be a vertex atom to define the polyhedron")
        elif centeratom is None or len(centeratom) == 0:
            raise IOError("There must be a center atom to define the polyhedron")

        for x in vertexatom:
            if x not in self.struc.atominfo:
                raise IOError("No such atom found in a structure file")
        for x in centeratom:
            if x not in self.struc.atominfo:
                raise IOError("No such atom found in a structure file")

        self.center_atom.append(centeratom)
        self.vertex_atom.append(vertexatom)

        return

    # TODO: Handling polyhedron without a center atom
    # TODO: polypos format change
    def polyposition(self, naming=False):
        # Method for calculating the atomic position of each polyhedron, with given center atom and vertex atom
        # Output atomic positions are all in the Cartesian coordinate type
        center_pos = []
        vertex_pos = []
        poly_pos = []
        center_element = []
        vertex_element = []
        poly_element = []

        # Finding center and vertex atom information from concatenated dataset, and appends them into a separated list
        for x in enumerate(self.struc.atomlist):
            if x[1][0] in self.center_atom[0]:
                center_pos.append(self.struc.coord[x[0]])
                center_element.append(x[1][0])
            if x[1][0] in self.vertex_atom[0]:
                vertex_pos.append(self.struc.coord[x[0]])
                vertex_element.append(x[1][0])

        center_pos = np.reshape(center_pos, (len(center_pos), 3))
        vertex_pos = np.reshape(vertex_pos, (len(vertex_pos), 3))

        for i in range(len(center_pos)):
            # poly_pos.append(self.utils.get_points_in_sphere(
            #      self.struc.unitvec, self.struc.matrix, vertex_pos, center_pos[i], self.bondmax, self.struc.cart))
            poly_pos.append(self.utils.atomsinsphere(
               self.struc.unitvec, self.struc.matrix, vertex_pos, center_pos[i], self.bondmax, self.struc.cart))

        # Append the info of center atom at the last column of array
        poly_pos = np.append(poly_pos, center_pos[:, None, :], axis=1)

        poly_pos = np.array(poly_pos)

        # Rounding for prevent confusion
        poly_pos = np.round(poly_pos, 10)

        if not self.struc.cart:
            center_pos = self.utils.directtocartesian(self.struc.matrix, center_pos)
            vertex_pos = self.utils.directtocartesian(self.struc.matrix, vertex_pos)
            poly_pos = self.utils.directtocartesian(self.struc.matrix, poly_pos)

        if not naming:
            return poly_pos

        if naming:
            return poly_pos, center_element

    @staticmethod
    def polyvolume(poly):
        from scipy.spatial import ConvexHull
        vol = ConvexHull(np.array(poly)[0:-1, :], qhull_options='Pp').volume
        return vol

    def quadraticelongation(self, poly):
        # Method to calculate Quadratic elongation of a given polyhedron
        # K. Robinson et al., Science 172, 567-570 (1971)
        # <lambda> = 1/n * sum(1 to n, (l_i / l_0) ** 2)

        if self.poly_level != 8:
            raise NotImplementedError("Only octahedron calculation available")
        else:
            pass

        # l_zero calculated from the polyhedral volume
        volume = self.polyvolume(poly)
        l_zero = (math.pow((3.0 * volume) / 4.0, 1.0 / 3.0))

        l_sum = 0.0
        for i in range(self.vertices):
            l_sum += math.pow((self.utils.vectordistance(poly[-1], poly[i]) / l_zero), 2)

        qe = l_sum / self.vertices

        return qe

    def distortionindex(self, poly):
        # Method to calculate Distortion index of a given polyhedron
        # K. Robinson et al., Science 172, 567-570 (1971)
        # <D> = 1/n * sum(1 to n, (|l_i - l_avg| / l_avg) )
        l_list = []
        l_sum = 0.0
        di = 0.0

        for x in poly[:-1]:
            l_sum += self.utils.vectordistance(x, poly[-1])
            l_list.append(self.utils.vectordistance(x, poly[-1]))

        l_list = np.array(l_list)
        l_avg = l_sum / len(l_list)

        for x in l_list:
            di += (abs(x - l_avg) / l_avg)
        di /= len(l_list)

        return di

    def bondanglevariance(self, poly):
        # Method to calculate Bond angle variance of a given polyhedron
        # W. H. Baur, Acta Crystallogr. Sect. B-Struct. Sci. 30, 1195-1215 (1974)
        # (sigma) ** 2 = 1 / (m - 1) * sum(1 to m, (phi_i - phi_0) ** 2)

        # Setting the phi_0 and m value corresponding to each polyhedron
        phi_zero = self.ideal_angle
        points = poly[:-1]
        center = poly[-1]

        phi = []
        diags = []

        for x in points:
            diags.append(self.polyvertexdiag(poly, x))
        diags = np.array(diags)

        # TODO: How to solve double counting?
        for i in range(len(points)):
            p1 = points[i]
            for x in points:
                for y in diags[i]:
                    if not np.allclose(x, y):
                        if not np.allclose(x, p1):
                            vec1 = p1 - center
                            vec2 = x - center
                            phi.append(self.utils.angle(vec1, vec2))

        phi = np.array(phi)

        bav = np.sum((phi - phi_zero) ** 2) / (len(phi) - 2)

        return bav

    def polycentertovertex(self, poly):
        cvdist = np.array([])
        for i in range(len(poly)):
            cvdist = np.append(cvdist, self.utils.vectordistance(poly[-1], poly[i]))
        cvdist = np.vstack(cvdist[:-1])
        return cvdist

    def polyvertextovertex(self, poly):
        vvdist = []
        for x in itertools.combinations(poly[:-1], 2):
            vvdist.append(self.utils.vectordistance(x[0], x[1]))
        vvdist = np.vstack(vvdist)
        return vvdist

    @staticmethod
    def polyfacetindices(poly):
        from scipy.spatial import ConvexHull
        indices = ConvexHull(poly, qhull_options='Pp').simplices
        return indices

    @staticmethod
    def polyconvexpoints(poly):
        from scipy.spatial import ConvexHull
        points = ConvexHull(poly, qhull_options='Pp').points
        return points

    def polyvertexdiag(self, poly, point, index=False):
        simp = self.polyfacetindices(poly[:-1])
        points = self.polyconvexpoints(poly[:-1])
        ind = int()

        for i in range(len(points)):
            if np.allclose(points[i], point):
                # if i < len(poly) - 1:
                ind = i

        tmp = []
        for x in simp:
            if ind in x:
                tmp.append(x)
        tmp = np.array(tmp).flatten()

        diag_ind = []
        for i in range(len(points)):
            if i not in tmp:
                diag_ind.append(i)

        diag_pos = []
        for i in range(len(diag_ind)):
            diag_pos.append(points[diag_ind[i]])

        diag_pos = np.array(diag_pos)

        if index:
            return diag_ind

        if not index:
            return diag_pos

    def polyshift(self, poly, cellmatrix, dim=None, cart=True, sets=False):
        # Method which imports the positions of polyhedron and shifts it to a certain directions
        # When the option 'set' is turned on, then the set of 9 polyhedron shifted along x,y,z direction is provided
        shift_matrix = []

        if not sets:
            if type(dim) == list:
                shift_matrix = [int(dim[0]), int(dim[1]), int(dim[2])]
            elif type(dim) == str:
                if ',' in dim:
                    shift_matrix = [int(dim.split(',')[0]), int(dim.split(',')[1]), int(dim.split(',')[2])]
                else:
                    shift_matrix = [int(dim.split()[0]), int(dim.split()[1]), int(dim.split()[2])]
            elif dim is None:
                shift_matrix = [1, 1, 1]

            if cart:
                shift_matrix = cellmatrix.dot(shift_matrix)

            shifted = poly + shift_matrix

            return shifted

        elif sets:
            # TODO: Do we need to consider 27 shifts?
            shiftedset = []
            dimrange = np.arange(-1, 2)

            for x in itertools.product(dimrange, dimrange, dimrange):
                shiftedset.append(self.polyshift(poly, cellmatrix, list(x), cart))

            shiftedset = np.array(shiftedset)

            return shiftedset

    @staticmethod
    def sharingatom(poly1, poly2):
        # Currently only can handle a sharing atom within one periodic image
        sharing = []
        for x in poly1:
            for y in poly2:
                if np.allclose(x, y):
                    sharing.append(x)

        sharing = np.array(sharing)

        return sharing

    def polyrelation(self, poly1, poly2):
        sharing = self.sharingatom(poly1, poly2)
        rel = None

        if np.allclose(poly1, poly2):
            rel = 'I'

        elif len(sharing) == 0:
            rel = None

        elif len(sharing) == 1:
            rel = 'C'

        elif len(sharing) == 2:
            rel = 'E'

        elif len(sharing) >= 3:
            rel = 'F'

        return rel

    def octacvcangle(self, poly1, poly2, vertexonly):
        # Method to measure the centroid-sharing vertex-centroid angle between two ocathedra
        # vertexonly boolean toggles to calculate centroid without considering center atom
        p1 = None
        p2 = None

        sharing = self.sharingatom(poly1, poly2)
        cvc = []

        if vertexonly:
            p1 = poly1[:-1]
            p2 = poly2[:-1]
        elif not vertexonly:
            p1 = poly1
            p2 = poly2

        c1 = self.utils.centroid(p1)
        c2 = self.utils.centroid(p2)

        for x in sharing:
            vec1 = c1 - x
            vec2 = c2 - x
            cvc.append(self.utils.angle(vec1, vec2))

        return cvc

    def octamvmangle(self, poly1, poly2):
        # Method to measure the center atom-sharing vertex-center atom angle between two ocathedra
        c1 = poly1[-1]
        c2 = poly2[-1]
        sharing = self.sharingatom(poly1, poly2)
        mvm = []

        for x in sharing:
            vec1 = c1 - x
            vec2 = c2 - x
            mvm.append(self.utils.angle(vec1, vec2))

        return mvm

    def octavvvangle(self, poly1, poly2):
        # Method to measure the diagonal vertex-sharing vertex-diagonal vertex angle between two ocathedra
        sharing = self.sharingatom(poly1, poly2)
        vvv = []

        for x in sharing:
            diag1 = self.polyvertexdiag(poly1, x)
            diag2 = self.polyvertexdiag(poly2, x)
            vec1 = diag1 - x
            vec2 = diag2 - x
            vvv.append(self.utils.angle(vec1[0], vec2[0]))

        return vvv

    def octacvcoangle(self, poly1, poly2, centroid, vertexonly):
        # Method to measure the angle between the plane (defined by Centroid1, Centroid2, sharing vertex)
        # and the diagonal vertex from sharing vertex
        # vertexonly boolean toggles to calculate centroid without considering center atom
        # centroid boolean toggles to consider the position of center atom as a centroid
        c1 = None
        c2 = None
        p1 = None
        p2 = None

        cvco = []
        sharing = self.sharingatom(poly1, poly2)

        if vertexonly:
            p1 = poly1[:-1]
            p2 = poly2[:-1]
        elif not vertexonly:
            p1 = poly1
            p2 = poly2

        if centroid:
            c1 = self.utils.centroid(p1)
            c2 = self.utils.centroid(p2)
        elif not centroid:
            c1 = poly1[-1]
            c2 = poly2[-1]

        for x in sharing:
            vec1 = c1 - x
            vec2 = c2 - x
            diag1 = self.polyvertexdiag(poly1, x)
            diag2 = self.polyvertexdiag(poly2, x)
            diagvec1 = diag1[0] - x
            diagvec2 = diag2[0] - x
            plane = self.utils.planefromvecs(vec1, vec2)
            a1 = self.utils.angle(plane, diagvec1)
            a2 = self.utils.angle(plane, diagvec2)
            if a1 <= 90:
                cvco.append(a1)
            elif 180 >= a1 > 90:
                cvco.append(180 - a1)

            if a2 <= 90:
                cvco.append(a2)
            elif 180 >= a2 > 90:
                cvco.append(180 - a2)

        return cvco

    def bruteangles(self, poly1, poly2, centroid, vertexonly, edge):
        cvc_res = []
        mvm_res = []
        vvv_res = []
        cvco_res = []

        p1 = poly1
        p2_range = self.polyshift(poly2, self.struc.matrix, None, True, True)

        for x in p2_range:
            if not edge:
                if self.polyrelation(p1, x) == 'C':
                    cvc_res.append(self.octacvcangle(p1, x, vertexonly))
                    mvm_res.append(self.octamvmangle(p1, x))
                    vvv_res.append(self.octavvvangle(p1, x))
                    cvco_res.append(self.octacvcoangle(p1, x, centroid, vertexonly))

            if edge:
                if self.polyrelation(p1, x) == 'C' or self.polyrelation(p1, x) == 'E':
                    cvc_res.append(self.octacvcangle(p1, x, vertexonly))
                    mvm_res.append(self.octamvmangle(p1, x))
                    vvv_res.append(self.octavvvangle(p1, x))
                    cvco_res.append(self.octacvcoangle(p1, x, centroid, vertexonly))

        result = [cvc_res, mvm_res, vvv_res, cvco_res]

        return result

    @staticmethod
    def tiltanglevariance(angles, ideal, avg=False):
        # Method to calculate Tilt angle variance between given polyhedra
        # W. Jang et al., ~~~
        # <theta> = 1/m * sum(1 to m, (theta_i / theta_0) ** 2)

        tav = 0
        count = 0
        if avg:
            for x in angles:
                if len(x) != 0:
                    for y in x:
                        tav += math.pow((y[0] / ideal), 2)
                        count += 1
        elif not avg:
            for x in angles:
                if len(x) != 0:
                    for y in x:
                        tav += y[0]
                        count += 1

        if count != 0:
            tav = tav / count
        elif count == 0:
            tav = None

        return tav

    def polyanalyze(self):
        poly = self.polyposition()
        index = []
        volume = []
        qe = []
        di = []
        bav = []
        parameters = [self.center_atom, self.vertex_atom, self.bondmax]

        a = 0
        for x in poly:
            a += 1
            index.append(a)
            volume.append(self.polyvolume(x))
            qe.append(self.quadraticelongation(x))
            di.append(self.distortionindex(x))
            bav.append(self.bondanglevariance(x))

        return [index, volume, qe, di, bav, parameters]

    def polytiltanalyze(self, centroid, vertexonly, ideal, edge, tav):
        poly = self.polyposition()
        angles = []
        cvc = []
        mvm = []
        vvv = []
        cvco = []

        for x in itertools.combinations_with_replacement(poly, 2):
            angles.append(self.bruteangles(x[0], x[1], centroid, vertexonly, edge))

        for x in angles:
            cvc.append(x[0])
            mvm.append(x[1])
            vvv.append(x[2])
            cvco.append(x[3])

        # if edge:
        #     cvc = np.array(cvc).tolist()
        #     mvm = np.array(mvm).tolist()
        #     vvv = np.array(vvv).tolist()
        #     cvco = np.array(cvco).tolist()

        res_cvc = self.tiltanglevariance(cvc, ideal, tav)
        res_mvm = self.tiltanglevariance(mvm, ideal, tav)
        res_vvv = self.tiltanglevariance(vvv, ideal, tav)
        res_cvco = self.tiltanglevariance(cvco, 90.0, tav)

        return [res_cvc, res_mvm, res_vvv, res_cvco]
