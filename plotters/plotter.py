import re
import numpy as np
from plotters.programs import IgorPlot


class BandPlot(object):
    def __init__(self, data: dict = None, outfile: str = None, prefix: str = "", program: str = "igor"):
        self.data = data
        self.outfile = outfile
        self.prefix = prefix
        self.program = program

        if self.outfile is None:
            if self.program == "igor":
                self.outfile = "band.itx"

    def add_guide(self, points=None):
        guideline = []
        kline = self.data["kline"]
        if points is None:
            for i in range(len(kline)):
                if i == 0:
                    guideline.append(kline[0])
                elif kline[i] == kline[i - 1]:
                    guideline.append(kline[i])
            guideline.append(kline[-1])

        else:
            slices = len(kline) / points
            for i in range(slices):
                guideline.append(kline[slices * points])

    def guide_symbols(self):
        return

    def plot(self):
        if self.program == "igor":
            plot = IgorPlot(self.outfile, presettype="band")
            proj_z = None
            bs = [z for spin, x in self.data["band"].items() for kp, y in x.items() for band, z in y.items()]
            band_label = [(self.prefix + "_" + band + "_" + spin).lstrip("_").rstrip("_")
                          for spin, x in self.data["band"].items() for band in x[list(x)[0]]]

            plot.add_wave([bs, band_label])

            if self.data["projection"] != {}:
                p = []
                for k in range(len(self.data["kline"])):
                    for atom_grp, v in self.data["projection"].items():
                        for orb_grp, w in v.items():
                            for spin, x in w.items():
                                for band, occ in x["kpoint_" + str(k + 1)].items():
                                    p.append(occ)

                proj_label = [(self.prefix + "_proj_" + self.data["proj_info"]["atom_group"][atom_idx] + "_" +
                               self.data["proj_info"]["orbital_group"][orb_idx] + "_" + band + "_" + spin)
                                  .lstrip("_").rstrip("_").replace("__", "_").replace("x2-y2", "x2my2")
                              for atom_idx, v in enumerate(self.data["projection"].items())
                              for orb_idx, w in enumerate(v[1].items())
                              for spin, x in w[1].items() for band in x[list(x)[0]]]

                plot.add_wave([p, proj_label])
                proj_z = np.reshape(proj_label, (int(len(proj_label) / len(band_label)), len(band_label))).tolist()

            plot.add_wave([self.data["kline"], [(self.prefix + "_band_kv").lstrip("_")]])
            plot.to_itx()
            plot.make_graph([band_label, (self.prefix + "_band_kv").lstrip("_")])

            if proj_z:
                for x in proj_z:
                    plot.make_graph([band_label, (self.prefix + "_band_kv").lstrip("_")],
                        title=x[0].replace("proj_", "").replace(re.findall(r"band_[\d']+", x[0])[0], "").rstrip("_"))
                    plot.plot_z(band_label, x)
        return


class DosPlot(object):
    def __init__(self, data: dict = None, outfile: str = None, prefix: str = "", program: str = "igor",
                 split_file: bool = False):
        self.data = data
        self.outfile = outfile
        self.prefix = prefix
        self.program = program

        if self.outfile is None:
            if self.program == "igor":
                self.outfile = "dos.itx"

        return

    def plot(self):
        if self.program == "igor":
            plot = IgorPlot(self.outfile, presettype="dos")

            tdos = [values for spin, x in self.data["tdos"].items() for key, y in x.items() for values in y]
            tdos = np.reshape(tdos, (3, int(len(tdos) / 3))).T.flatten().tolist()
            tdos_label = [(self.prefix + "_" + key + "_" + spin).lstrip("_").rstrip("_")
                          for spin, x in self.data["tdos"].items() for key in x.keys()]

            plot.add_wave([tdos, tdos_label])

            pdos = [values for atom_grp, v in self.data["pdos"].items() for orb_grp, w in v.items()
                    for spin, x in w.items() for values in x.values()]
            pdos_label = [(self.prefix + self.data["proj_info"]["atom_group"][atom_idx] + "_" +
                           self.data["proj_info"]["orbital_group"][orb_idx] + "_" + spin)
                              .lstrip("_").rstrip("_").replace("__", "_").replace("x2-y2", "x2my2")
                          for atom_idx, v in enumerate(self.data["pdos"].items())
                          for orb_idx, w in enumerate(v[1].items()) for spin in w[1].keys()]

            plot.add_wave([pdos, pdos_label])

            plot.to_itx()

            to_plot = pdos_label
            to_plot.insert(0, tdos_label[1])

            plot.make_graph([to_plot, tdos_label[0]], title=self.prefix + "_dos")
