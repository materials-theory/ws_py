import copy
import math
import numpy as np
import scipy.constants as constants

from scipy import interpolate
from vw_py.Parsers.eigenparserV import EigenParserV
from vw_py.Parsers.procarV import ProcarParserV as Procar
from vw_py.Parsers.outcarhandler import OutcarHandler as Outcar
from vw_py.generalutils.operator import StrucOperators as utils
from vw_py.generalutils.operator import StatisticFuncs as stats
from vw_py.Parsers.structV import ContcarHandler as cont


class Bandstructure(object):

    def __init__(self):
        self.eigen = None
        self.numband = None
        self.numkp = None
        self.band = None
        self.proj = None
        self.proj_mz = None
        self.kvec = None
        self.vbm = None
        self.cbm = None
        self.gap = None
        self.fermi = None
        self.emass = None
        self.path = None
        self.atominfo = None
        self.traj = None
        self.emass_residuals = None

        return

    def band_data(self, fermi, fakeweight, shift):

        self.eigen = EigenParserV()
        self.numband = self.eigen.numstates

        recvec = np.array(utils().recvec(cont.unitvec()))
        path = np.array([])
        band = np.array([])
        numkp = 0

        # Handling the fake-weight type bandstructure calculations
        if fakeweight is True:
            print("Fake-weight bandstructure calculation...")
            for i in range(len(self.eigen.path)):
                if self.eigen.wgt[i][0] != 0.0:
                    pass
                else:
                    numkp += 1
                    path = np.append(path, self.eigen.path[i])
                    band = np.append(band, self.eigen.states[i])

            path = np.reshape(path, (numkp, 3))
            band = np.reshape(band, (numkp, self.numband))

        else:
            print("Normal bandstructure calculation...")
            numkp = copy.deepcopy(self.eigen.numkp)
            path = copy.deepcopy(self.eigen.path)
            band = copy.deepcopy(self.eigen.states)

        if fermi == 0.0:
            fermi_e = float(Outcar.param_from_outcar('E-fermi'))
        else:
            fermi_e = fermi

        # Shifting for the fermi level value
        # Extra shifting of vbm only if it is turned on
        if shift is True:
            self.bandedge(fermi_e)
            band -= self.vbm[2]
        else:
            band -= fermi_e

        # Converting the path trajectory with 2pi scaling and reciprocal vector
        path_conv = np.array([])

        for i in range(len(path)):
            for j in range(3):
                path_conv = np.append(path_conv, (path[i] * recvec * 2 * math.pi)[:, j].sum())

        path_conv = np.array(path_conv, dtype='d')
        path_conv = np.reshape(path_conv, (len(path), 3))

        # Assigning k vector using the distance between k-points
        k_vec = np.array([])
        for i in range(len(path_conv)):
            if i == 0:
                k_vec = np.append(k_vec, 0.0)
            else:
                k_vec = np.append(k_vec, k_vec[i - 1] + (utils.vectordistance(path_conv[i], path_conv[i - 1])))

        k_vec = np.reshape(np.array(k_vec, dtype='d'), (numkp, 1))

        self.numkp = numkp
        self.band = band
        self.kvec = k_vec
        self.path = path
        self.traj = path_conv

        return

    @staticmethod
    def rawband():
        eigen = EigenParserV()
        numband = eigen.numstates
        recvec = np.array(utils().recvec(cont.unitvec()))

        numkp = copy.deepcopy(eigen.numkp)
        path = copy.deepcopy(eigen.path)
        band = copy.deepcopy(eigen.states)

        path_conv = []

        for i in range(len(path)):
            for j in range(3):
                path_conv.append((path[i] * recvec * 2 * math.pi)[:, j].sum())

        path_conv = np.array(path_conv, dtype='d')
        path_conv = np.reshape(path_conv, (len(path), 3))

        # Assigning k vector using the distance between k-points
        k_vec = []
        for i in range(len(path_conv)):
            if i == 0:
                k_vec.append(0.0)
            else:
                k_vec.append(k_vec[i - 1] + (utils.vectordistance(path_conv[i], path_conv[i - 1])))
        k_vec = np.reshape(np.array(k_vec, dtype='d'), (numkp, 1))

        return path, band, k_vec, numband

    def projectedband_data(self, fermi, fakeweight, shift, spin, atom, orbital):
        dic = Procar(None, spin).as_dict()

        numbands = dic['numbands']
        numkps = dic['numkps']
        # numions = dic['numions']
        path_orig = dic['path']
        states = dic['states']
        weight = dic['weight']
        proj = dic['proj']
        proj_spin = dic['proj_spin']

        recvec = np.array(utils().recvec(cont.unitvec()))
        path = []
        band = []
        numkp = 0

        atomnum = []
        if atom is not None:
            for x in atom:
                tmp = []
                tmp2 = []
                if '-' in x:
                    for y in range(int(x.split('-')[0]), int(x.split('-')[1])):
                        tmp.append(y)
                    atomnum.append(tmp)
                else:
                    tmp2.append(x)

                if len(tmp2) != 0:
                    atomnum.append(tmp2)
                else:
                    pass

        else:
            pass

        # Handling the fake-weight type bandstructure calculations
        if fakeweight is True:
            print("Fake-weight bandstructure calculation...")
            for i in range(len(path_orig)):
                if weight[i][0] != 0.0:
                    pass
                else:
                    numkp += 1
                    path.append(path_orig[i])
                    band.append(states[i])

        else:
            print("Normal bandstructure calculation...")
            numkp = copy.deepcopy(numkps)
            path = copy.deepcopy(path_orig)
            band = copy.deepcopy(states)

        numkp = np.array(numkp)
        path = np.array(path)
        band = np.array(band)

        if fermi == 0.0:
            fermi_e = float(Outcar.param_from_outcar('E-fermi'))
        else:
            fermi_e = fermi

        self.fermi = fermi_e

        # Shifting for the fermi level value
        # Extra shifting of vbm only if it is turned on
        if shift is True:
            self.bandedge(fermi_e)
            band -= self.vbm[2]
        else:
            band -= fermi_e

        # Converting the path trajectory with 2pi scaling and reciprocal vector
        path_conv = []

        for i in range(len(path)):
            for j in range(3):
                path_conv.append((path[i] * recvec * 2 * math.pi)[:, j].sum())

        path_conv = np.array(path_conv, dtype='d')
        path_conv = np.reshape(path_conv, (len(path), 3))

        # Assigning k vector using the distance between k-points
        k_vec = []
        for i in range(len(path_conv)):
            if i == 0:
                k_vec.append(0.0)
            else:
                k_vec.append(k_vec[i - 1] + (utils.vectordistance(path_conv[i], path_conv[i - 1])))

        k_vec = np.reshape(np.array(k_vec, dtype='d'), (numkp, 1))

        proj_atom = []
        proj_orb = []
        atom_data = []

        for x in cont.atominfo().items():
            atom_data.append(x)

        self.atominfo = atom_data
        numatoms = len(atom_data)

        if atom is None:
            for i in range(numbands):
                for j in range(numkp):
                    proj_atom.append(proj[j][i][-1])

            numatoms = 1

            # proj_atom = np.reshape(np.array(proj_atom), (1, self.eigen.numstates, self.eigen.numkp, 5))

        elif atom is not None:
            if len(atomnum) == 0:
                for i in range(len(atom_data)):
                    proj_atom.append([])

                for i in range(numbands):
                    for j in range(numkp):
                        count = 0
                        for k in range(len(atom_data)):
                            tmp = []
                            for l in range(int(atom_data[k][1])):
                                tmp.append(proj[j][i][l + count])
                            count += int(atom_data[k][1])
                            tmp = np.array(tmp, dtype='d')
                            proj_atom[k].append(np.sum(tmp, axis=0))

                # proj_atom = np.reshape(np.array(proj_atom), (numatoms, self.eigen.numstates, self.eigen.numkp, 5))

            else:
                for i in range(len(atomnum)):
                    proj_atom.append([])

                for i in range(numbands):
                    for j in range(numkp):
                        for k in range(len(atomnum)):
                            tmp = []
                            for l in atomnum[k]:
                                tmp.append(proj[j][i][int(l)])
                            tmp = np.array(tmp, dtype='d')
                            proj_atom[k].append(np.sum(tmp, axis=0))

                numatoms = len(atomnum)
                # proj_atom = np.reshape(np.array(proj_atom), (len(atomnum), self.eigen.numstates, self.eigen.numkp, 5))

        proj_atom = np.reshape(np.array(proj_atom), (numatoms, numbands, numkp, 5))

        # Ion type - # of band - # of kp
        if orbital is False:
            for i in range(len(proj_atom)):
                for j in range(numbands):
                    proj_orb.append(proj_atom[i, j, :, 4])

            proj_orb = np.reshape(np.array(proj_orb), (numatoms, 1, numbands, numkp))

        # Ion type - type of orbital - # of band - # of kp
        elif orbital is True:
            for i in range(len(proj_atom)):
                for k in range(4):
                    for j in range(numbands):
                        proj_orb.append(proj_atom[i, j, :, k + 1])

            proj_orb = np.reshape(np.array(proj_orb), (numatoms, 4, numbands, numkp))

        self.path = path_conv
        self.proj = proj_orb
        self.numkp = numkp
        self.band = band
        self.kvec = k_vec
        self.numband = numbands

        if spin is True:
            proj_atom = []
            proj_orb = []

            if atom is None:
                for i in range(numbands):
                    for j in range(numkp):
                        proj_atom.append(proj_spin[2][j][i][-1])
                numatoms = 1

            elif atom is not None:
                if len(atomnum) == 0:
                    for i in range(len(atom_data)):
                        proj_atom.append([])

                    for i in range(numbands):
                        for j in range(numkp):
                            count = 0
                            for k in range(len(atom_data)):
                                tmp = []
                                for l in range(int(atom_data[k][1])):
                                    tmp.append(proj_spin[2][j][i][l + count])
                                count += int(atom_data[k][1])
                                tmp = np.array(tmp, dtype='d')
                                proj_atom[k].append(np.sum(tmp, axis=0))

                else:
                    for i in range(len(atomnum)):
                        proj_atom.append([])

                    for i in range(numbands):
                        for j in range(numkp):
                            for k in range(len(atomnum)):
                                tmp = []
                                for l in atomnum[k]:
                                    tmp.append(proj_spin[2][j][i][int(l)])
                                tmp = np.array(tmp, dtype='d')
                                proj_atom[k].append(np.sum(tmp, axis=0))

                    numatoms = len(atomnum)

            proj_atom = np.reshape(np.array(proj_atom), (numatoms, numbands, numkp, 5))

            # Ion type - # of band - # of kp
            if orbital is False:
                for i in range(len(proj_atom)):
                    for j in range(numbands):
                        proj_orb.append(proj_atom[i, j, :, 4])

                proj_orb = np.reshape(np.array(proj_orb), (numatoms, 1, numbands, numkp))

            # Ion type - type of orbital - # of band - # of kp
            elif orbital is True:
                for i in range(len(proj_atom)):
                    for k in range(4):
                        for j in range(numbands):
                            proj_orb.append(proj_atom[i, j, :, k + 1])

                proj_orb = np.reshape(np.array(proj_orb), (numatoms, 4, numbands, numkp))

            self.proj_mz = proj_orb
        return

    def temperature_band(self, fermi, fakeweight, temperature,
                         emin, emax, egrid, kgrid, truncate_e,
                         smearmin, smearmax, smeargrid, method):

        self.projectedband_data(fermi, fakeweight, None, False, '', False)
        print("Adding temperature smearing...")
        dic = self.as_dict()
        kvec = dic['kvec']
        atom = dic['atom']
        # proj_spin = dic['pband_spin']
        canvas = np.zeros((kgrid, egrid))
        ewindow = np.linspace(emin, emax, egrid)  # - self.fermi
        dense_kvec = np.linspace(np.min(kvec), np.max(kvec), kgrid)
        projected = []

        for i in range(len(atom)):
            for j in range(dic['numband']):
                band = dic['band'].T[j]
                proj = dic['pband'][i][-1][j]  # #ion, #orb(-1 = total) #band

                band_interp = interpolate.interp1d(kvec[:, 0], band)
                proj_interp = interpolate.interp1d(kvec[:, 0], proj)

                dense_band = band_interp(dense_kvec)
                dense_proj = proj_interp(dense_kvec)

                for k in range(kgrid):
                    p = self.temperature_convolution(dense_band[k], dense_proj[k], temperature,
                                                     ewindow, truncate_e, smearmin, smearmax, smeargrid, method)
                    projected.append(p)
                    canvas[k] += p

        axes = [dense_kvec, ewindow]
        return canvas, axes

    @staticmethod
    def temperature_convolution(bandeigen, bandweight, temperature,
                                ewindow, truncate_e, smearmin, smearmax, smeargrid, method):
        smear_func = None
        kb = constants.value('Boltzmann constant in eV/K')
        smearwindow = np.linspace(smearmin, smearmax, smeargrid)
        smear_dist = stats.gaussianprob(smearwindow, 0.0, (temperature * kb))
        smear_dist = smear_dist / np.sum(smear_dist)

        if method[0] is 'G':
            # print("Using Gaussian smearing...\n")
            smear_func = stats.gaussiandist(smearwindow, 0.0, (temperature * kb))
        elif method[0] is 'F':
            # print("Using Fermi-Dirac smearing...\n")
            smear_func = stats.fermidiracdist(smearwindow, 0.0, (temperature * kb))
        elif method[0] is 'B':
            # print("Using Bose-Einstein smearing...\n")
            smear_func = stats.boseeinsteindist(smearwindow, 0.0, (temperature * kb))
        elif method[0] is 'M':
            # print("Using Maxwell-Boltzmann smearing...\n")
            smear_func = stats.maxwellboltzmanndist(smearwindow, 0.0, (temperature * kb))

        convoluted = smear_dist * smear_func
        convoluted = convoluted / np.sum(convoluted)

        bandrange = smearwindow + bandeigen
        distribution = convoluted * bandweight

        # Truncating both ends
        bandrange = np.insert(bandrange, 0, (ewindow[0] - truncate_e), axis=0)
        bandrange = np.insert(bandrange, bandrange.shape[0], (ewindow[-1] + truncate_e), axis=0)
        distribution = np.insert(distribution, 0, 0.0, axis=0)
        distribution = np.insert(distribution, distribution.shape[0], 0.0, axis=0)

        proj_interp = interpolate.interp1d(bandrange, distribution, kind='linear')
        projected = proj_interp(ewindow)

        return projected

    def bandedge(self, fermi):
        # Calculating the band edge position of semiconductor
        if self.eigen is None:
            self.eigen = EigenParserV()

        vbm = [0, 0, -10.0]
        cbm = [0, 0, 10.0]

        if fermi == 0.0:
            fermi_e = float(Outcar.param_from_outcar('E-fermi'))
        else:
            fermi_e = fermi

        for i, x in enumerate(self.eigen.states):
            for j, y in enumerate(x):
                if vbm[2] <= y <= fermi_e:
                    vbm = [i, j, y]
                elif cbm[2] >= y >= fermi_e:
                    cbm = [i, j, y]

        self.vbm = vbm
        self.cbm = cbm
        self.fermi = fermi_e
        self.gap = self.cbm[2] - self.vbm[2]

        return

    def effectivemass(self, band, kpoint, points, show_residuals):
        ev_joule_conv = 1.60217733 * 10 ** (-19)
        plank_const = 6.626070040 * 10 ** (-34)
        h_bar = plank_const / (2 * math.pi)
        m_zero = 9.1093821545 * 10 ** (-31)
        recvec = np.array(utils().recvec(cont.unitvec()))

        emass = []
        residuals = []

        if band is 0:
            kp = self.vbm[0]
            # kp2 = self.cbm[0]

            band = self.vbm[1]
            # band2 = self.cbm[1]

        else:
            kp = kpoint - 1
            band = band - 1

        # Left direction
        if kp - points < 0:
            pass
        else:
            eig = []
            path = []
            count = 0
            for i in range(points):
                if i == 0:
                    eig.append(self.eigen.states[kp - count][band])
                    path.append(self.eigen.path[kp - count])
                else:
                    if self.eigen.states[kp - count][band] == eig[i - 1]:
                        count += 1
                        eig.append(self.eigen.states[kp - count][band])
                        path.append(self.eigen.path[kp - count])
                    else:
                        eig.append(self.eigen.states[kp - count][band])
                        path.append(self.eigen.path[kp - count])
                count += 1

            tmp = []
            for x in path:
                for j in range(3):
                    tmp.append((x * recvec * 2 * math.pi)[:, j].sum())
            path = np.reshape(tmp, (len(path), 3))

            kvec = np.array([])
            for i in range(len(path)):
                if i == 0:
                    kvec = np.append(kvec, 0.0)
                else:
                    kvec = np.append(kvec, kvec[i - 1] + (utils.vectordistance(path[i], path[i - 1])))

            # print(kvec, eig)
            poly = np.polyfit(np.asfarray(kvec), np.asfarray(eig), 2, full=True)

            fit_grad = poly[0][0]
            curvature_quad = fit_grad * 2 * ev_joule_conv

            emass_raw = (h_bar ** 2) / (curvature_quad * 10 ** (-20))
            emass.append([emass_raw / m_zero, "Left"])
            residuals.append([poly[1][0], "Left"])

        # Right direction
        if kp + points > len(self.eigen.path):
            pass
        else:
            eig = []
            path = []
            count = 0
            for i in range(points):
                if i == 0:
                    eig.append(self.eigen.states[kp + count][band])
                    path.append(self.eigen.path[kp + count])
                else:
                    if self.eigen.states[kp + count][band] == eig[i - 1]:
                        count += 1
                        eig.append(self.eigen.states[kp + count][band])
                        path.append(self.eigen.path[kp + count])
                    else:
                        eig.append(self.eigen.states[kp + count][band])
                        path.append(self.eigen.path[kp + count])
                count += 1

            tmp = []
            for x in path:
                for j in range(3):
                    tmp.append((x * recvec * 2 * math.pi)[:, j].sum())
            path = np.reshape(tmp, (len(path), 3))

            kvec = np.array([])
            for i in range(len(path)):
                if i == 0:
                    kvec = np.append(kvec, 0.0)
                else:
                    kvec = np.append(kvec, kvec[i - 1] + (utils.vectordistance(path[i], path[i - 1])))

            # print(kvec, eig)
            poly = np.polyfit(np.asfarray(kvec), np.asfarray(eig), 2, full=True)

            fit_grad = poly[0][0]
            curvature_quad = fit_grad * 2 * ev_joule_conv

            emass_raw = (h_bar ** 2) / (curvature_quad * 10 ** (-20))
            emass.append([emass_raw / m_zero, "Right"])
            residuals.append([poly[1][0], "Right"])

        self.emass = np.array(emass)
        self.emass_residuals = np.array(residuals)

        return

    def as_dict(self):
        dic = {'numband': self.numband,
               'numkp': self.numkp,
               'band': self.band,
               'pband': self.proj,
               'pband_mz': self.proj_mz,
               'path': self.path,
               'atom': self.atominfo,
               'kvec': self.kvec
               }

        return dic
