import json
import os


class IgorPlot(object):
    def __init__(self, outfile, options=None):
        self.outfile = outfile
        self.options = options
        self.wavelist = []

        with open(os.path.join(os.path.dirname(__file__), "igor_preset.json"), "r") as preset:
            self.preset = json.load(preset)

    # [[value lists], [wavename]]
    def add_wave(self, wave):
        self.wavelist.append(wave)

    def plot_z(self, targetwave, zwave, zmin=0.00, zmax=1.00, color="Red", reverse=True):
        if len(targetwave) != len(zwave):
            raise IndexError("Length of target wave and z wave not matching!")

        if reverse:
            reverse = 1
        else:
            reverse = 0

        with open(self.outfile, "a") as out:
            for i in range(len(targetwave)):
                out.write("X ModifyGraph lsize=2, mode(" + targetwave[i] + ")=0, zColor(" + targetwave[i] + ")={"
                          + zwave[i] + "," + str(zmin) + "," + str(zmax) + "," + color + "," + str(reverse) + "}\n")

    def write_wave(self, outfile):
        for x in self.wavelist:
            outfile.write("WAVES/D ")
            for y in x[1]:
                outfile.write(y + " ")
            outfile.write("\n")

            outfile.write("BEGIN\n")
            self.data_write(outfile, x)
            outfile.write("END\n")

    @staticmethod
    def data_write(outfile, data):
        line = int(len(data[0]) / len(data[1]))
        count = 0
        for i in range(line):
            for j in range(len(data[1])):
                outfile.write("{: >14.10f}".format(data[0][len(data[1]) * count + j]) + " ")
            outfile.write("\n")
            count += 1

    def make_graph(self, xypair, preset=True, title=None):
        with open(self.outfile, "a") as out:
            out.write("\nX Display ")

            for i in range(len(xypair[0])):
                if i != len(xypair[0]) - 1:
                    out.write(xypair[0][i] + ",")
                else:
                    if title is not None:
                        out.write(xypair[0][i] + " vs " + xypair[1] + " as \"" + title + "\"\n")
                    else:
                        out.write(xypair[0][i] + " vs " + xypair[1] + "\n")

            if preset:
                for key, value in self.preset.items():
                    if hasattr(value, "keys"):
                        for subkey, subvalues in value.items():
                            if type(subvalues) is str:
                                out.write("X " + key + " " + subkey + " " + str(self.preset[key][subkey]) + "\n")
                            else:
                                out.write("X " + key + " " + subkey + "=" + str(self.preset[key][subkey]) + "\n")
                    else:
                        out.write("X " + key + " " + value + "\n")

    def to_itx(self, xdata=None, ydata=None):
        with open(self.outfile, "w") as out:
            out.write("IGOR\n")
            self.write_wave(out)
        if xdata is not None and ydata is not None:
            for i in range(len(xdata)):
                self.make_graph([xdata[i], ydata[i]])
            return
