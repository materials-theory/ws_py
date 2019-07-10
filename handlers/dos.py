import re
import numpy as np

from collections import defaultdict


class Dos(object):
    def __init__(self, tdos: dict = None, pdos: dict = None, structure: object = None, efermi: float = None,
                 shift: bool = True, others: bool = False, atom: list = None, orbital: list = None):
        self.tdos = tdos
        self.pdos = pdos
        self.structure = structure
        self.efermi = efermi
        self.shifted = False
        self.projected = {}
        self.proj_info = {}

        if shift is True:
            self.shift_dos()

        self.projected = self.projection_orbital(self.projection_atom(self.pdos, atom, others), orbital, others)

    def nested_dict(self):
        return defaultdict(self.nested_dict)

    def shift_dos(self):
        to_shift = self.efermi
        for spin, x in self.tdos.items():
            for i in range(len(x["energy"])):
                x["energy"][i] -= to_shift
        self.shifted = True

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

        for spin, x in target["ion"].items():
            for egrid, y in x.items():
                for atom_idx, occ in enumerate(y.values()):
                    for group_idx, group in enumerate(to_sum):
                        if atom_idx in group:
                            if len(to_sum) == 1:
                                if len(projected[""][spin][egrid]) == 0:
                                    projected[""][spin][egrid] = np.array(occ)
                                else:
                                    projected[""][spin][egrid] += np.array(occ)
                            else:
                                if len(projected["atom_group_" + str(group_idx)][spin][egrid]) == 0:
                                    projected["atom_group_" + str(group_idx)][spin][egrid] = np.array(occ)
                                else:
                                    projected["atom_group_" + str(group_idx)][spin][egrid] += np.array(occ)

        for group, x in projected.items():
            for spin, y in x.items():
                for egrid, z in y.items():
                    projected[group][spin][egrid] = z.tolist()

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

    def projection_orbital(self, target, orbital=None, proj_others=True):
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
                for x in self.pdos["field"]:
                    if x[0] not in [a for b in to_sum for a in b]:
                        tmp.append(x[0])
                if tmp:
                    to_sum.append(tmp)
                    orbital.append([["others"]])

        if orbital is None:
            orbital = []
            for x in self.pdos["field"]:
                to_sum.append([x[0]])
                orbital.append([[x[0]]])

        for i in range(len(to_sum)):
            for j in range(len(to_sum[i])):
                for idx, orb in enumerate(self.pdos["field"]):
                    if to_sum[i][j] == orb[0]:
                        to_sum[i][j] = [orb[0], idx]

        for group, x in target.items():
            for spin, y in x.items():
                for egrid, occ in y.items():
                    for orb_group_idx, orb_group in enumerate(to_sum):
                        projected[group]["orb_group_" + str(orb_group_idx)][spin][egrid] = 0
                        for orb in orb_group:
                            projected[group]["orb_group_" + str(orb_group_idx)][spin][egrid] += occ[orb[1]]

        for i in range(len(orbital)):
            tmp = ""
            for x in orbital[i]:
                tmp += x[0]
            orbital[i] = tmp

        self.proj_info["orbital_group"] = orbital
        return projected

    def get_plot_dict(self, tdos, pdos):
        energy = self.nested_dict()

        for spin, x in self.tdos.items():
            energy[spin] = x["energy"]

        to_plot = {"tdos": tdos,
                   "pdos": pdos,
                   "proj_info": self.proj_info,
                   "energy": energy}

        return to_plot
