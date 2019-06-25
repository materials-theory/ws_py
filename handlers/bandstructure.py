import re
import numpy as np

from collections import defaultdict


class Bandstructure(object):

    def __init__(self, eigenvalues: dict = None, projection: dict = None, ktrace: list = None,
                 structure: object = None, efermi: float = None, shift="vbm",
                 kweight: list = None, atom: list = None, orbital: list = None, ymin: float = None, ymax: float = None):
        self.eigenvalues = eigenvalues
        self.projection = projection
        self.ktrace = ktrace
        self.structure = structure
        self.efermi = efermi
        self.shift = shift
        self.kweight = kweight
        self.atom = atom
        self.orbital = orbital
        self.ymin = ymin
        self.ymax = ymax
        self.metal = False
        self.shifted = False
        self.projected = {}
        self.proj_info = {}

        self.shift_band()

        if ymin or ymax:
            self.apply_window(ymin, ymax)

        if projection is not None:
            self.projected = self.projection_orbital(self.projection_atom(self.projection, self.atom), self.orbital)

    def nested_dict(self):
        return defaultdict(self.nested_dict)

    def projection_atom(self, target, atom_group=None, proj_others=True):
        to_sum = []
        projected = self.nested_dict()

        if atom_group is not None:
            for x in atom_group:
                if type(x) is list:
                    to_sum.append([int(y) - 1 for y in x])
                else:
                    if "-" in x:
                        to_sum.append([a for a in range(int(x.split("-")[0]) - 1, int(x.split("-")[1]))])
                    else:
                        # to_sum.append(re.findall(r"[\w']+", x))
                        to_sum.append([int(x) - 1])

            if proj_others is True:
                tmp = []
                for i in range(1, len(self.structure.elements)):
                    if i not in [a for b in to_sum for a in b]:
                        tmp.append(i)
                if tmp:
                    to_sum.append(tmp)
                    atom_group.append("others")

        elif atom_group is None:
            to_sum = [[idx for idx in range(len(self.structure.elements))]]
            atom_group = [""]

        # TODO: simplifying nested for loops

        for spin, x in target["ion"].items():
            for kp, y in x.items():
                for band, z in y.items():
                    for atom_idx, occ in enumerate(z):
                        for group_idx, group in enumerate(to_sum):
                            if atom_idx in group:
                                if len(to_sum) == 1:
                                    if len(projected[""][spin][kp][band]) == 0:
                                        projected[""][spin][kp][band] = np.array(occ)
                                    else:
                                        projected[""][spin][kp][band] += np.array(occ)
                                else:
                                    if len(projected["atom_group_" + str(group_idx)][spin][kp][band]) == 0:
                                        projected["atom_group_" + str(group_idx)][spin][kp][band] = np.array(occ)
                                    else:
                                        projected["atom_group_" + str(group_idx)][spin][kp][band] += np.array(occ)

        for group, x in projected.items():
            for spin, y in x.items():
                for kp, z in y.items():
                    for band, items in z.items():
                        projected[group][spin][kp][band] = items.tolist()

        for i in range(len(atom_group)):
            if type(atom_group[i]) is str:
                atom_group[i] = atom_group[i].replace("-", "to")
            elif type(atom_group[i]) is list:
                tmp = ""
                for x in atom_group[i]:
                    tmp = tmp + x + "_"
                atom_group[i] = tmp.rstrip("_")
        self.proj_info["atom_group"] = atom_group
        return projected

    # TODO: need to simplify...

    def projection_orbital(self, target, orbital=None, proj_others=True, total=True):
        to_sum = []
        projected = self.nested_dict()

        if orbital is not None:
            for x in orbital:
                if type(x) is list:
                    to_sum.append(x)
                else:
                    to_sum.append(re.findall(r"[\w']+", x))

            if proj_others is True:
                tmp = []
                for x in self.projection["field"]:
                    if x[0] not in [a for b in to_sum for a in b]:
                        tmp.append(x[0])
                if tmp:
                    to_sum.append(tmp)
                    orbital.append(tmp)

        if total is True:
            to_sum.append([x for orb in self.projection["field"] for x in orb[0]])
            orbital.append([["total"]])

        for i in range(len(to_sum)):
            for j in range(len(to_sum[i])):
                for idx, orb in enumerate(self.projection["field"]):
                    if to_sum[i][j] == orb[0]:
                        to_sum[i][j] = [orb[0], idx]

        for group, x in target.items():
            for spin, y in x.items():
                for kp, z in y.items():
                    for band, occ in z.items():
                        for orb_group_idx, orb_group in enumerate(to_sum):
                            projected[group]["orb_group_" + str(orb_group_idx)][spin][kp][band] = 0
                            for orb in orb_group:
                                projected[group]["orb_group_" + str(orb_group_idx)][spin][kp][band] += occ[orb[1]]

        for i in range(len(orbital)):
            tmp = ""
            for x in orbital[i]:
                tmp += x[0]
            orbital[i] = tmp
        self.proj_info["orbital_group"] = orbital
        return projected

    def ktrace_convert(self, ktrace):
        kline = []
        for i in range(len(ktrace)):
            if i == 0:
                kline.append(0.0)
            else:
                kline.append(kline[i - 1] + np.linalg.norm((np.array(ktrace[i]) - np.array(ktrace[i - 1])) *
                                                           np.array(self.structure.get_reciprocal_vec())) * 2 * np.pi)
        return kline

    def metal(self, tol_efermi=1e-4, tol_spillover=0.2):
        for x in self.eigenvalues["band"]:
            for i, y in enumerate(self.eigenvalues["band"][x]):
                for j, z in enumerate(self.eigenvalues["band"][x][y]):
                    if z[0] > self.efermi + tol_efermi and z[1] >= tol_spillover:
                        return True
        return False

    def vbm(self):
        vbm = (-10.0, None, None, None)
        for spin, x in self.eigenvalues.items():
            for kp, y in x.items():
                for band, z in y.items():
                    if vbm[0] <= z < self.efermi:
                        vbm = (z, int(re.findall(r"[\d']+", band)[0]), int(re.findall(r"[\d']+", kp)[0]), spin)
        return vbm

    def cbm(self):
        cbm = (10.0, None, None, None)
        for spin, x in self.eigenvalues.items():
            for kp, y in x.items():
                for band, z in y.items():
                    if cbm[0] >= z > self.efermi:
                        cbm = (z, int(re.findall(r"[\d']+", band)[0]), int(re.findall(r"[\d']+", kp)[0]), spin)
        return cbm

    def gap(self):
        if self.metal is True:
            return (0.0, None)
        vbm = self.vbm()
        cbm = self.cbm()
        gap = (self.cbm()[0] - self.vbm()[0], [vbm[2], cbm[2]])
        return gap

    def shift_band(self):
        to_shift = 0.0
        if self.shift is None:
            pass
        elif self.shift in "vbmVBM":
            to_shift = self.vbm()[0]
        elif self.shift in "fermiF":
            to_shift = self.efermi

        for spin, x in self.eigenvalues.items():
            for kp, y in x.items():
                for band, z in y.items():
                    self.eigenvalues[spin][kp][band] -= to_shift
        self.shifted = True

    def apply_window(self, ymin, ymax):
        if ymin is None:
            ymin = -30.0
        if ymax is None:
            ymax = 30.0

        for spin, x in self.eigenvalues.items():
            band_to_remove = []
            for i in range(len(x["kpoint_1"])):
                tmp = []
                for kp, y in x.items():
                    tmp.append(y["band_" + str(i + 1)])
                if min(tmp) <= ymin:
                    band_to_remove.append("band_" + str(i + 1))
                elif ymax <= max(tmp):
                    band_to_remove.append("band_" + str(i + 1))
            for r in band_to_remove:
                for kp, y in x.items():
                    del(y[r])

            for proj_spin, p in self.projection["ion"].items():
                for r in band_to_remove:
                    for kp, y in p.items():
                        del (y[r])

    def get_plot_dict(self, eig, proj):
        to_plot = {"band": eig,
                   "projection": proj,
                   "proj_info": self.proj_info,
                   "kline": self.ktrace_convert(self.ktrace)}

        return to_plot
