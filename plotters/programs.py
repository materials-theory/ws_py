import json
import os

class IgorPlot(object):
    def __init__(self, outfile, options=None):
        self.outfile = outfile
        self.options = options
        self.wavelist = []
        self.wavecount = 0

        with open(os.path.join(os.path.dirname(__file__), "igor_preset.json"), "r") as preset:
            self.options = json.load(preset)

    # [[value lists], [wavename]]
    def add_wave(self, wave):
        self.wavelist.append(wave)

    def plot_z(self):
        return

    def write_wave(self, outfile):
        for x in self.wavelist:
            outfile.write("WAVES/D ")
            for y in x[1]:
                outfile.write(y + " ")
            outfile.write("\n")

            outfile.write("BEGIN\n")
            self.data_write(outfile, x)
            outfile.write("END\n")

        return

    def data_write(self, outfile, data):
        line = int(len(data[0]) / len(data[1]))
        count = 0
        for i in range(line):
            for j in range(len(data[1])):
                outfile.write("{: >14.10f}".format(data[0][len(data[1]) * count + j]) + " ")
            outfile.write("\n")
            count += 1

    def to_itx(self):
        with open(self.outfile, "w") as out:
            out.write("IGOR\n")
            self.write_wave(out)
            return
