import re
import numpy as np

from collections import defaultdict


class Bandstructure(object):

    def __init__(self, eigenvalues: dict = None, projection: dict = None, ktrace: list = None,
                 structure: object = None, efermi: float = None, shift="vbm",
                 kweight: list = None, atom: list = None, orbital: list = None):
        self.eigenvalues = eigenvalues
        self.projection = projection
        self.ktrace = ktrace
        self.structure = structure
        self.efermi = efermi
        self.shift = shift
        self.kweight = kweight
        self.atom = atom
        self.orbital = orbital
        self.metal = False
        self.shifted = False
        self.projected = {}
        self.proj_info = {}

        self.shift_band()

        if projection is not None:
            self.projected = self.projection_orbital(self.projection_atom(self.projection, self.atom), self.orbital)
            # if atom is not None and orbital is not None:
            #     self.projected = self.projection_orbital(self.projection_atom(self.projection, self.atom), self.orbital)
            # elif atom is not None and orbital is None:
            #     self.projected = self.projection_atom(self.projection, self.atom)
            # elif atom is None and orbital is not None:
            #     self.projected = self.projection_orbital(self.projection, self.orbital)

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

        elif proj_others is True and atom_group is not None:
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

        self.proj_info["atom"] = to_sum
        self.proj_info["atom_input"] = atom_group
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

        if proj_others is True and orbital is not None:
            tmp = []
            for x in self.projection["field"]:
                if x[0] not in [a for b in to_sum for a in b]:
                    tmp.append(x[0])
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
                            for orb in orb_group:
                                if len(to_sum) == 1:
                                    if type(projected[group][""][spin][kp][band]) is not np.ndarray:
                                        projected[group][""][spin][kp][band] = np.array(occ[orb[1]])
                                    else:
                                        projected[group][""][spin][kp][band] += np.array(occ[orb[1]])
                                else:
                                    if type(projected[group]["orb_group_" + str(orb_group_idx)][spin][kp][band]) is not np.ndarray:
                                        projected[group]["orb_group_" + str(orb_group_idx)][spin][kp][band] = np.array(occ[orb[1]])
                                    else:
                                        projected[group]["orb_group_" + str(orb_group_idx)][spin][kp][band] += np.array(occ[orb[1]])

        for atom_group, w in projected.items():
            for orb_group, x in w.items():
                for spin, y in x.items():
                    for kp, z in y.items():
                        for band, items in z.items():
                            projected[atom_group][orb_group][spin][kp][band] = items.tolist()

        self.proj_info["orbital"] = to_sum
        self.proj_info["orbital_input"] = orbital
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
        vbm = (-10.0, None, None)
        for x in self.eigenvalues["band"]:
            for i, y in enumerate(self.eigenvalues["band"][x]):
                for j, z in enumerate(self.eigenvalues["band"][x][y]):
                    if vbm[0] <= z[0] < self.efermi:
                        vbm = (z[0], j, i)
        return vbm

    def cbm(self):
        cbm = (10.0, None, None)
        for x in self.eigenvalues["band"]:
            for i, y in enumerate(self.eigenvalues["band"][x]):
                for j, z in enumerate(self.eigenvalues["band"][x][y]):
                    if cbm[0] >= z[0] > self.efermi:
                        cbm = (z[0], j, i)
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

    def get_plot_dict(self, eig, proj):
        to_plot = {"band": eig,
                   "projection": proj,
                   "proj_info": self.proj_info,
                   "kline": self.ktrace_convert(self.ktrace)}

        return to_plot
