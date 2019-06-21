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
        bs = [z for spin, x in self.data["band"].items() for kp, y in x.items() for band, z in y.items()]
        band_label = [(self.prefix + "_" + band + "_" + spin).lstrip("_").rstrip("_")
                      for spin, x in self.data["band"].items() for band in x[list(x)[0]]]

        plot.add_wave([bs, band_label])

        if self.data["projection"] is not {}:
            p = [occ for atom_grp, v in self.data["projection"].items() for orb_grp, w in v.items() for spin, x in
                  w.items() for kp, y in x.items() for band, occ in y.items()]

            proj_label = [(self.prefix + "_proj_" + atom_grp + "_" + orb_grp + "_" + band + "_" + spin)
                              .lstrip("_").rstrip("_").replace("__", "_")
                          for atom_grp, v in self.data["projection"].items()
                          for orb_grp, w in v.items() for spin, x in w.items() for band in x[list(x)[0]]]

            plot.add_wave([p, proj_label])

        plot.add_wave([self.data["kline"], [(self.prefix + "_band_kv").lstrip("_")]])
        plot.to_itx()
        return
