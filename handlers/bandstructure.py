import numpy as np


class Bandstructure(object):

    def __init__(self, eigenvalues: dict = None, projected: dict = None, ktrace: list = None,
                 structure: object = None, efermi: float = None, shift="vbm",
                 kweight: list = None, atom: list = None, orbital: list = None):
        self.eigenvalues = eigenvalues
        self.projected = projected
        self.ktrace = ktrace
        self.structure = structure
        self.efermi = efermi
        self.shift = shift
        self.kweight = kweight
        self.atom = atom
        self.orbital = orbital
        self.metal = False

        a = self.get_plot_dict("band")
        return

    def eigenval_reformat(self):
        reformatted = []
        for x in self.eigenvalues["band"].keys():
            for y in self.eigenvalues["band"][x].keys():
                tmp = []
                for z in self.eigenvalues["band"][x][y]:
                    tmp.append(z[0])
                reformatted.append(tmp)
        return list(zip(*reformatted))

    def projection_atom(self, atom, sum_others=True):
        to_sum = []
        for x in atom:
            if type(x) is list:
                to_sum.append(x)
            else:
                to_sum.append([a for a in range(int(x.split("-")[0]), int(x.split("-")[1]) + 1)])

        if sum_others is True:
            tmp = []
            for i in range(1, len(self.structure.elements) + 1):
                if i not in [a for b in to_sum for a in b]:
                    tmp.append(i)
            to_sum.append(tmp)
        return

    def projection_orbital(self, orbital):
        return

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

    def get_plot_dict(self, plot):
        count = 0
        to_plot = {}
        to_shift = 0.0

        if self.shift is None:
            pass
        elif "v" in self.shift:
            to_shift = self.vbm()[0]
        elif "f" in self.shift:
            to_shift = self.efermi

        if plot == "band":
            to_plot["total"] = {"spin_1": {}}
            for x in self.eigenval_reformat():
                count += 1
                to_plot["total"]["spin_1"]["band_" + str(count)] = tuple(np.subtract(x, to_shift))

        if plot == "pband":
            return

        to_plot["kline"] = self.ktrace_convert(self.ktrace)
        return to_plot
