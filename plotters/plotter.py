import re
import numpy as np
from plotters.programs import IgorPlot


class BandPlot(object):
    def __init__(self, data: dict = None, outfile: str = "band.itx", prefix: str = "", program: str = "igor"):
        self.data = data
        self.outfile = outfile
        self.prefix = prefix
        self.program = program

        self.itx_plot()
        return

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

    def itx_plot(self):
        plot = IgorPlot(self.outfile)
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
                              .lstrip("_").rstrip("_").replace("__", "_")
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
                                title=x[0].replace("proj_", "").replace(re.findall(r"band_[\d']+_", x[0])[0], ""))
                plot.plot_z(band_label, x)
        return
