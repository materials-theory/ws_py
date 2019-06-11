import math
import numpy as np
from scipy.special import erf
# import scipy.constants as constants
# from scipy import interpolate, ndimage


def vectorlength(vec):
    """
    Method to calculate unit vector length / the distance of a point from the origin

    Argument : vec (list or np array) (e.g. ['0', '0', '0'] or [0, 0, 0])

    Returns : float
    """
    len_vec = math.sqrt(math.pow(float(vec[0]), 2.0) + math.pow(float(vec[1]), 2.0) + math.pow(float(vec[2]), 2.0))

    return len_vec


def vectordistance(vec1, vec2):
    """
    Method to calculate the distance between two vector components

    Argument : vec (list or np array) (e.g. ['0', '0', '0'] or [0, 0, 0])

    Returns : float
    """
    v1 = None
    v2 = None

    if type(vec1) == np.ndarray:
        v1 = np.reshape(vec1, (3,))
    if type(vec2) == np.ndarray:
        v2 = np.reshape(vec2, (3,))

    dist_vec = math.sqrt(math.pow(float(v1[0]) - float(v2[0]), 2.0) +
                         math.pow(float(v1[1]) - float(v2[1]), 2.0) +
                         math.pow(float(v1[2]) - float(v2[2]), 2.0))
    return dist_vec


def cross(vec1, vec2):
    """
    Method to calculate the cross product of two vector components

    Argument : vec1, vec2 (list or np array) (e.g. ['0', '0', '0'] or [0, 0, 0])

    Returns : np array
    """

    cross_prod = np.cross(np.array(vec1, dtype='d'), np.array(vec2, dtype='d'))

    return cross_prod


def dot(vec1, vec2):
    """
    Method to calculate the dot product of two vector components

    Argument : vec1, vec2 (list or np array) (e.g. ['0', '0', '0'] or [0, 0, 0])

    Returns : np array
    """

    dot_prod = np.dot(np.array(vec1, dtype='d'), np.array(vec2, dtype='d'))

    return dot_prod


def volume(unitvec):
    """
    Method calculating the volume from three given vectors

    Argument : unitvec (list or np array) consists of three vectors

    Returns : float
    """

    vol = dot(unitvec[0], cross(unitvec[1], unitvec[2]))

    return vol


def recvec(unitvec):
    """
    Calculating the reciprocal vector with a given set of unit vectors

    Argument : unitvec (list or np array) consists of three vectors

    Returns : np array
    """

    astar = cross(unitvec[1], unitvec[2]) / volume(unitvec)
    bstar = cross(unitvec[2], unitvec[0]) / volume(unitvec)
    cstar = cross(unitvec[0], unitvec[1]) / volume(unitvec)

    return np.array([astar, bstar, cstar])


def directtocartesian(unitvec, posdata):
    """
    Converts the position data from Direct coordinate to Cartesian coordinate.

    Argument : unitvec and posdata (list or np array)

    Returns : np array
    """
    return np.dot(posdata, unitvec)


def cartesiantodirect(unitvec, posdata):
    """
    Converts the position data from Direct coordinate to Cartesian coordinate.

    Argument : unitvec and posdata (list or np array)

    Returns : np array
    """
    return np.dot(posdata, np.linalg.inv(unitvec))


def angle(vec1, vec2):
    """
    Calculates the angle between two vectors

    Argument : two vectors (list or np array)

    Returns : float
    """
    # noinspection PyTypeChecker
    ang = math.degrees(math.acos(round(dot(vec1, vec2) / (vectorlength(vec1) * vectorlength(vec2)), 10)))
    return ang


# def centroid(atom_pos, atom_weight=None):
#     pos = np.array(atom_pos, dtype='d')
#
#     if atom_weight is not None:
#         wgt = np.vstack(np.array(atom_pos, dtype='d'))
#         # wgt = np.concatenate((wgt,wgt,wgt),axis=1)
#     else:
#         wgt = np.vstack(np.ones([len(pos), 3]))
#
#     weighted_pos = pos * wgt
#     x_wgt_sum = wgt[:, 0].sum()
#     y_wgt_sum = wgt[:, 1].sum()
#     z_wgt_sum = wgt[:, 2].sum()
#
#     x_pos = (weighted_pos[:, 0].sum() / x_wgt_sum)
#     y_pos = (weighted_pos[:, 1].sum() / y_wgt_sum)
#     z_pos = (weighted_pos[:, 2].sum() / z_wgt_sum)
#
#     result = np.array([x_pos, y_pos, z_pos])
#
#     # return result.tolist()
#     return result
#
#
# # TODO: unitvec and matrix conflict solve
# def atomsinsphere(self, unitvec, matrix, pos, centerpoint, radius, cart=True, tolerance=0.01):
#     """
#     Method to pick out the points existing in the sphere with a center and radius value given.
#     First, it calculates the number of periodic images required to fit in the sphere with a radius.
#
#     """
#     center_cart = []
#     center_frac = []
#     vert_cart = []
#     vert_frac = []
#
#     recvec = np.array(self.recvec(unitvec))
#     recvec_len = np.array([self.vectorlength(recvec[0]),
#                            self.vectorlength(recvec[1]),
#                            self.vectorlength(recvec[2])])
#
#     unit_vec = np.array(unitvec)
#     unitvec_len = np.array([self.vectorlength(unit_vec[0]),
#                             self.vectorlength(unit_vec[1]),
#                             self.vectorlength(unit_vec[2])])
#
#     image_calls = radius * recvec_len + tolerance
#
#     if cart:
#         center_cart.append(centerpoint)
#         center_frac.append(self.cartesiantodirect(matrix, centerpoint))
#
#         if len(np.shape(pos)) != 1:
#             for x in pos:
#                 vert_cart.append(x)
#                 vert_frac.append(self.cartesiantodirect(matrix, x))
#
#         else:
#             vert_cart.append(pos)
#             vert_frac.append(self.cartesiantodirect(matrix, pos))
#
#     elif not cart:
#         center_frac.append(centerpoint)
#         center_cart.append(self.directtocartesian(matrix, centerpoint))
#
#         if len(np.shape(pos)) != 1:
#             for x in pos:
#                 vert_frac.append(x)
#                 vert_cart.append(self.directtocartesian(matrix, x))
#
#         else:
#             vert_frac.append(pos)
#             vert_cart.append(self.directtocartesian(matrix, pos))
#
#     center_cart = np.array(center_cart)
#     center_frac = np.array(center_frac)
#     vert_cart = np.array(vert_cart)
#     vert_frac = np.array(vert_frac)
#
#     min_unit = np.floor(center_frac - image_calls)
#     max_unit = np.ceil(center_frac + image_calls)
#
#     a_imag = np.arange(min_unit[0][0], max_unit[0][0])
#     b_imag = np.arange(min_unit[0][1], max_unit[0][1])
#     c_imag = np.arange(min_unit[0][2], max_unit[0][2])
#
#     imag_list = []
#     for x in a_imag:
#         for y in b_imag:
#             for z in c_imag:
#                 imag_list.append([x, y, z])
#
#     images = np.array(imag_list)
#
#     tmp_list = []
#     for x in images:
#         for y in vert_frac:
#             tmp_list.append(x + y)
#
#     atom_all_frac = np.array(tmp_list)
#     atom_all_cart = []
#     for x in atom_all_frac:
#         atom_all_cart.append(self.directtocartesian(matrix, x))
#
#     atom_all_cart = np.array(atom_all_cart)
#
#     dist_info = []
#     # tmp_list = list(atom_all_cart.tolist())
#     for x in atom_all_cart:
#         dist_info.append(self.vectordistance(x, center_cart))
#
#     dist_info = np.vstack(dist_info)
#
#     result_index = np.where(np.logical_and(np.greater(dist_info, 0.0), np.greater_equal(radius, dist_info)))
#
#     result = []
#     for i in range(len(result_index[0])):
#         if cart:
#             result.append(atom_all_cart[result_index[0][i]])
#         elif not cart:
#             result.append(atom_all_frac[result_index[0][i]])
#
#     return result
#
#
# def planefromvecs(self, vec1, vec2, coeff=False):
#     normalvec = self.cross(vec1, vec2)
#     if all(elem == 0 for elem in normalvec):
#         nonzero = [i for i, e in enumerate(vec1) if e != 0]
#         zero = [i for i, e in enumerate(vec1) if e == 0]
#         guess = (vec1[zero[0]] + vec1[zero[1]]) / vec1[nonzero[0]]
#         testvec = [1, 1, 1]
#         testvec[nonzero[0]] = guess
#         testvec = np.array(testvec)
#         if np.dot(vec2, testvec) == 0:
#             normalvec = testvec / self.vectorlength(testvec)
#         else:
#             raise ValueError("Plane normal vector is zero vector")
#
#     a = normalvec[0]
#     b = normalvec[1]
#     c = normalvec[2]
#     d = self.dot(normalvec, vec1)
#
#     planecoeff = [a, b, c, d]
#
#     if not coeff:
#         result = np.array(normalvec)
#     elif coeff:
#         result = np.array(planecoeff)
#
#     return result
#
#
# def planefitting(dataset, fitorder=1, grid=0.1, coeff=False):
#     import scipy.optimize
#     import functools
#
#     def plane(x, y, params):
#         avec = params[0]
#         bvec = params[1]
#         cvec = params[2]
#         z = avec * x + bvec * y + cvec
#         return z
#
#     def error(params, points):
#         result = 0
#         for (x, y, z) in points:
#             plane_z = plane(x, y, params)
#             diff = abs(plane_z - z)
#             result += diff ** 2
#         return result
#
#     # def cross(a, b):
#     #     return [a[1] * b[2] - a[2] * b[1],
#     #             a[2] * b[0] - a[0] * b[2],
#     #             a[0] * b[1] - a[1] * b[0]]
#
#     fun = functools.partial(error, points=dataset)
#     params0 = [0, 0, 0]
#     res = scipy.optimize.minimize(fun, params0)
#
#     a = res.x[0]
#     b = res.x[1]
#     c = res.x[2]
#
#     # point = np.array([0.0, 0.0, c])
#     # normal = np.array(cross([1, 0, a], [0, 1, b]))
#     # d = -point.dot(normal)
#     # xx, yy = np.meshgrid([-5, 10], [-5, 10])
#     # z = (-normal[0] * xx - normal[1] * yy - d) * 1. / normal[2]
#     #
#     # plt.show()
#
#     return [a, b, c]
#
#
# def get_points_in_sphere(self, unitvec, matrix, points, center, r, cart=True):
#     """
#     Function from Pymatgen project (needs citation?)
#     failure check in triclinic structures (?)
#
#     """
#     center_cart = []
#     center_frac = []
#     vert_cart = []
#     vert_frac = []
#
#     if cart:
#         center_cart.append(center)
#         center_frac.append(self.cartesiantodirect(matrix, center))
#
#         if len(np.shape(points)) != 1:
#             for x in points:
#                 vert_cart.append(x)
#                 vert_frac.append(self.cartesiantodirect(matrix, x))
#
#         else:
#             vert_cart.append(points)
#             vert_frac.append(self.cartesiantodirect(matrix, points))
#
#     elif not cart:
#         center_frac.append(center)
#         center_cart.append(self.directtocartesian(matrix, center))
#
#         if len(np.shape(points)) != 1:
#             for x in points:
#                 vert_frac.append(x)
#                 vert_cart.append(self.directtocartesian(matrix, x))
#
#         else:
#             vert_frac.append(points)
#             vert_cart.append(self.directtocartesian(matrix, points))
#
#     center_cart = np.array(center_cart)
#     center_frac = np.array(center_frac)
#     vert_cart = np.array(vert_cart)
#     vert_frac = np.array(vert_frac)
#
#     recpvec = np.array(self.recvec(unitvec))
#     recp_len = np.array([self.vectorlength(recpvec[0]),
#                            self.vectorlength(recpvec[1]),
#                            self.vectorlength(recpvec[2])])  # / (2 * np.pi)
#
#     # recp_len = np.array(self.reciprocal_lattice.abc) / (2 * np.pi())
#     nmax = float(r) * recp_len + 0.01
#
#     pcoords = center_frac
#     center = center_cart
#
#     n = len(vert_frac)
#     fcoords = np.array(vert_frac) % 1
#     indices = np.arange(n)
#
#     mins = np.floor(pcoords - nmax)[0]
#     maxes = np.ceil(pcoords + nmax)[0]
#     arange = np.arange(start=mins[0], stop=maxes[0])
#     brange = np.arange(start=mins[1], stop=maxes[1])
#     crange = np.arange(start=mins[2], stop=maxes[2])
#     arange = arange[:, None] * np.array([1, 0, 0])[None, :]
#     brange = brange[:, None] * np.array([0, 1, 0])[None, :]
#     crange = crange[:, None] * np.array([0, 0, 1])[None, :]
#     images = arange[:, None, None] + brange[None, :, None] +\
#         crange[None, None, :]
#
#     shifted_coords = fcoords[:, None, None, None, :] + \
#         images[None, :, :, :, :]
#
#     cart_coords = self.directtocartesian(matrix, fcoords)
#     cart_images = self.directtocartesian(matrix, images)
#     # cart_coords = []
#     # cart_images = [[[]]]
#     # for x in fcoords:
#     #     cart_coords.append(self.directtocartesian(matrix, x))
#
#     # cart_coords = np.array(cart_coords)
#     # cart_images = np.array(cart_images)
#
#     coords = cart_coords[:, None, None, None, :] + \
#         cart_images[None, :, :, :, :]
#     coords -= center[None, None, None, :]
#     coords **= 2
#     d_2 = np.sum(coords, axis=4)
#
#     within_r = np.where(d_2 <= r ** 2)
#     if cart:
#         result = self.directtocartesian(matrix, shifted_coords[within_r])
#     elif not cart:
#         result = shifted_coords[within_r]
#
#     return result
#         # np.sqrt(d_2[within_r]),
#         # indices[within_r[0]]


class StatisticFuncs(object):
    def __init__(self):
        return

    @staticmethod
    def gaussiandist(x, mu, sigma):
        func = 0.5 * (1 - erf((x - mu) / sigma))
        return func

    @staticmethod
    def fermidiracdist(x, mu, sigma):
        # boltzmann = constants.value('Boltzmann constant in eV/K')
        # sigma = boltzmann * temperature
        func = 1 / (np.exp((x - mu) / sigma) + 1)
        return func

    @staticmethod
    def boseeinsteindist(x, mu, sigma):
        # boltzmann = constants.value('Boltzmann constant in eV/K')
        # sigma = boltzmann * temperature
        func = 1 / (np.exp((x - mu) / sigma) - 1)
        return func

    @staticmethod
    def maxwellboltzmanndist(x, mu, sigma):
        # boltzmann = constants.value('Boltzmann constant in eV/K')
        # sigma = boltzmann * temperature
        func = 1 / (np.exp((x - mu) / sigma))
        return func

    @staticmethod
    def gaussianprob(x, mu, sigma):
        prob = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp((-1 * ((x - mu) ** 2)) / (2 * (sigma ** 2)))
        return prob
