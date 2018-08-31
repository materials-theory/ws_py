import numpy as np
import math
from IO.IO import IO

# TODO: Moving out the import statement

class Plotter(object):
    def __init__(self, outfile, inputname=None):
        self.outfile = outfile
        self.banddata = None
        self.dosdata = None
        self.io = None
        if inputname is None:
            self.input_name = input("Please input the name of system  :  ")
        else:
            self.input_name = inputname
        return

    def boltzplot(self, boltzdic, shift):
        if self.outfile is None:
            outfile = 'boltzplot.itx'
        else:
            outfile = self.outfile

        print("-----------------------------------------------------------------------------------------------")
        print("Plotting BoltzTraP data...")

        params = ['N', 'DOS', 'S', 'st', 'R', 'k', 'c', 'chi']
        temps = sorted(list(boltzdic.keys()))

        with open(outfile, 'w') as out:
            out.write("IGOR\n")

            for t in temps:
                data = []
                prefix = str(t) + "K"
                wavehead = wavename(params, self.input_name, prefix)
                wavehead.append(str(self.input_name) + "_" + prefix + "_energy")
                waveheader(out, wavehead)

                energies = sorted(list(boltzdic[t].keys()))
                for y in params:
                    data.append(self.boltzdata(boltzdic[t], y))

                for i in range(len(data[0])):
                    for lines in data:
                        out.write(str(lines[i]))
                        out.write(" ")
                    out.write(str(energies[i] - shift))
                    out.write("\n")
                out.write("END\n")

            fermi = (energies[0] + energies[-1]) / 2.0
            for x in params:
                out.write("X Display" + " as " + '"' + str(self.input_name) + "_" + str(x) + '"' + "\n")
                for t in temps:
                    prefix = str(t) + "K"
                    out.write("X AppendToGraph " + str(self.input_name) + "_" + prefix + "_" + str(x)
                              + " vs " + str(self.input_name) + "_" + prefix + "_energy" + "\n")

                graphpreset(out, None, None, (fermi - shift - 1), (fermi - shift + 1))

    @staticmethod
    def boltzdata(dic, index):
        data = []

        energies = sorted(list(dic.keys()))
        for energy in energies:
            data.append(dic[energy][index])

        return data

    def plotband(self, fermi, fakeweight, shift, guide, ktraj):
        from ElectronicStructure.bandstructure import Bandstructure

        if self.outfile is None:
            self.outfile = 'band.itx'
        else:
            self.outfile = self.outfile + '.itx'

        self.io = IO(None, self.outfile)
        self.banddata = Bandstructure()
        self.banddata.band_data(fermi, fakeweight, shift)

        out = self.io.WriteFile()

        print("-----------------------------------------------------------------------------------------------")
        print("Plotting bandstructure...")
        input_name = self.input_name

        # Writing the header part of .itx
        out.write("IGOR\n")
        out.write("WAVES/D ")
        out.write(input_name + "_Band_kv")

        tmp_wavename = []
        for i in range(self.banddata.numband):
            tmp_wavename.append(" " + input_name + "_Band_" + str(i + 1))
            out.write(" " + input_name + "_Band_" + str(i + 1))
        out.write("\nBEGIN\n")

        # Writing the wave part of .itx
        for i in range(self.banddata.numkp):
            out.write(str(self.banddata.kvec[i][0]))
            for j in range(self.banddata.numband):
                out.write(" " + str(self.banddata.band[:, j][i]))
            out.write("\n")
        out.write("END\n")

        # Line for the high-symmetry point guideline
        if guide is True:
            guideline = []
            for i in range(self.banddata.numkp):
                if i == 0:
                    guideline.append(self.banddata.kvec[i][0])
                elif self.banddata.kvec[i][0] == self.banddata.kvec[i - 1][0]:
                    guideline.append(self.banddata.kvec[i][0])
                elif i == (self.banddata.numkp - 1):
                    guideline.append(self.banddata.kvec[i][0])
            guideline = np.array(guideline)
            out.write("WAVES/D ")
            out.write(" " + input_name + "_guide " + input_name + "_guide_y1 " + input_name + "_guide_y2\n")
            out.write("BEGIN\n")
            for i in range(len(guideline)):
                out.write(str(guideline[i]) + " -30.00 30.00\n")
            out.write("END\n")

        if ktraj is True:
            out.write("WAVES/D ")
            out.write(" " + input_name + "_ktraj_x " + input_name + "_ktraj_y " + input_name + "_ktraj_z\n")
            out.write("BEGIN\n")
            for i in range(len(self.banddata.traj)):
                out.write(" " + str(self.banddata.traj[:, 0][i]) + " " + str(self.banddata.traj[:, 1][i]) + " " + str(self.banddata.traj[:, 2][i]) + "\n")
            out.write("END\n")

        # Writing the graph plotting part of .itx
        tmp = "X Display" + "".join(tmp_wavename)
        limchar = 15

        if len(tmp_wavename) <= limchar:
            out.write(str(tmp))
            out.write(" vs " + input_name + "_Band_kv")
            out.write(" as " + ' \"' + "band_" + input_name + '\"' + "\n")

        else:
            numline = math.ceil(len(tmp) / limchar)
            numpart = math.ceil(len(tmp_wavename) / numline)

            for i in range(int(numline)):
                if i == 0:
                    out.write("X Display" + "".join(tmp_wavename[:numpart]))
                    out.write(" vs " + input_name + "_Band_kv")
                    out.write(" as " + '\"' + "band_" + input_name + '\"' + "\n")
                elif len("".join(tmp_wavename[(i * numpart):(i * numpart + numpart)])) == 0:
                    pass
                else:
                    out.write("X AppendToGraph" + "".join(tmp_wavename[(i * numpart):(i * numpart + numpart)]))
                    out.write(" vs " + input_name + "_Band_kv\n")

        preset = ("X DefaultFont/U \"Times New Roman\"\n"
                  "X ModifyGraph width=255.118,height=340.157\n"
                  "X ModifyGraph marker=19\n"
                  "X ModifyGraph lSize=1.5\n"
                  "X ModifyGraph tick(left)=2,tick(bottom)=3,noLabel(bottom)=2\n"
                  "X ModifyGraph mirror=1\n"
                  "X ModifyGraph zero(left)=8\n"
                  "X ModifyGraph fSize=28\n"
                  "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                  "X ModifyGraph standoff=0\n"
                  "X ModifyGraph axThick=1.5\n"
                  "X ModifyGraph axisOnTop=1\n"
                  "X Label left \"\Z28 Energy (eV)\"\n"
                  "X ModifyGraph zero(bottom)=0;DelayUpdate\n"
                  "X SetAxis left -3,3\n"
                  "X ModifyGraph zeroThick(left)=2.5")

        preset_guide = ("\nX AppendToGraph " + input_name + "_guide_y1 " + input_name + "_guide_y2 vs " +
                        input_name + "_guide\n"
                        "X ModifyGraph mode(" + input_name + "_guide_y1)=1,rgb(" + input_name + "_guide_y1)=(0,0,0)\n"
                        "X ModifyGraph mode(" + input_name + "_guide_y2)=1,rgb(" + input_name + "_guide_y2)=(0,0,0)\n"
                        "X SetAxis left -3,3"
                        )

        # Writing the graph layout part of .itx
        out.write(preset)
        if guide is True:
            out.write(preset_guide)

        return

    def plotprojectedband(self, fermi, fakeweight, shift, guide, spin, atom, orbital):
        # TODO: Orbital-specific .itx generation
        # TODO: Energy window
        # TODO: Changing the method names for easy reading
        # TODO: Wave name entry change for easy use of wave browser in Igor

        from ElectronicStructure.bandstructure import Bandstructure

        if self.outfile is None:
            self.outfile = 'pband.itx'
        else:
            self.outfile = self.outfile + '.itx'

        self.io = IO(None, self.outfile)
        self.banddata = Bandstructure()
        self.banddata.projectedband_data(fermi, fakeweight, shift, spin, atom, orbital)

        out = self.io.WriteFile()

        print("-----------------------------------------------------------------------------------------------")
        print("Plotting projected bandstructure...")
        input_name = self.input_name

        def bandname():
            tmp_wavename = []
            for i in range(self.banddata.numband):
                tmp_wavename.append(" " + input_name + "_Band_" + str(i + 1))

            return tmp_wavename

        def pbandname():
            wavename_prefix = []
            if atom is not None:
                if atom == []:
                    atominfo = self.banddata.atominfo
                    tmp = []
                    for x in atominfo:
                        tmp.append(x[0])
                    wavename_prefix.append(tmp)

                else:
                    tmp = []
                    for x in atom:
                        tmp2 = []
                        if '-' in x:
                            tmp.append(x.replace('-', 'to'))
                        else:
                            tmp2.append(x)
                            if len(tmp2) != 0:
                                tmp.append("_".join(tmp2))

                    wavename_prefix.append(tmp)

            elif atom is None:
                wavename_prefix.append([''])

            if orbital is True:
                wavename_prefix.append(['s', 'p', 'd', 't'])

            elif orbital is False:
                wavename_prefix.append([''])

            return wavename_prefix

        def pwavename(pbandname):
            proj_wavename = []
            for i in range(len(pbandname[0])):
                for j in range(len(pbandname[1])):
                    for k in range(self.banddata.numband):
                        tmp = ("Proj_" + input_name + "_Band_" + str(k + 1) + "_" + str(pbandname[0][i])
                               + "_" + str(pbandname[1][j])).replace("__", "_")
                        if tmp[-1] == "_":
                            tmp = tmp[:-1]
                        proj_wavename.append(tmp)
            return proj_wavename

        def figurename(atom, orbital):
            figname = input_name + "_" + str(atom) + "_" + str(orbital)
            figname = figname.replace("__", "_")
            if figname[-1] == "_":
                figname = figname[:-1]
            return figname

        def itxheader(bandname):
            # Writing the header part of .itx
            out.write("IGOR\n")
            out.write("WAVES/D ")
            out.write(input_name + "_Band_kv")
            for x in bandname:
                out.write(x)
            out.write("\nBEGIN\n")

        def bandwave():
            # Writing the wave part of .itx
            for i in range(self.banddata.numkp):
                out.write(str(self.banddata.kvec[i][0]))
                for j in range(self.banddata.numband):
                    out.write(" " + str(self.banddata.band[:, j][i]))
                out.write("\n")
            out.write("END\n")

        def guidewave():
            # Line for the high-symmetry point guideline
            if guide is True:
                guideline = []
                for i in range(self.banddata.numkp):
                    if i == 0:
                        guideline.append(self.banddata.kvec[i][0])
                    elif self.banddata.kvec[i][0] == self.banddata.kvec[i - 1][0]:
                        guideline.append(self.banddata.kvec[i][0])
                    elif i == (self.banddata.numkp - 1):
                        guideline.append(self.banddata.kvec[i][0])
                guideline = np.array(guideline)
                out.write("WAVES/D ")
                out.write(" " + input_name + "_guide " + input_name + "_guide_y1 " + input_name + "_guide_y2\n")
                out.write("BEGIN\n")
                for i in range(len(guideline)):
                    out.write(str(guideline[i]) + " -30.00 30.00\n")
                out.write("END\n")

        def pbandwave(pbandname, pwavename):
            # Wave for the projected bands
            out.write("\nWAVES/D ")

            for x in pwavename:
                out.write(x + " ")
            out.write("\nBEGIN\n")

            for i in range(self.banddata.numkp):
                for j in range(len(pbandname[0])):
                    for k in range(len(pbandname[1])):
                        for l in range(self.banddata.numband):
                            out.write(" " + str(self.banddata.proj[j][k][l][i]))
                out.write("\n")
            out.write("END\n")

        def spinmzwave(pbandname):
            mz_wavename = []
            for i in range(len(pbandname[0])):
                for j in range(len(pbandname[1])):
                    for k in range(self.banddata.numband):
                        tmp = (("Proj_" + input_name + "_Band_"
                                + str(k + 1) + "_" + str(pbandname[0][i]) + "_" + str(pbandname[1][j] + "_mz")
                                ).replace("__", "_"))
                        if tmp[-1] == "_":
                            tmp = tmp[:-1]
                        mz_wavename.append(tmp)

            # Wave for the projected bands
            out.write("\nWAVES/D ")

            for x in mz_wavename:
                out.write(x + " ")
            out.write("\nBEGIN\n")

            for i in range(self.banddata.numkp):
                for j in range(len(pbandname[0])):
                    for k in range(len(pbandname[1])):
                        for l in range(self.banddata.numband):
                            out.write(" " + str(self.banddata.proj_mz[j][k][l][i]))
                out.write("\n")
            out.write("END\n")

        def displayband(bandname, figurename=None):
            # Writing the graph plotting part of .itx
            tmp = "\nX Display" + "".join(bandname)
            limchar = 15
            if figurename is None:
                if len(bandname) <= limchar:
                    out.write(str(tmp))
                    out.write(" vs " + input_name + "_Band_kv")
                    out.write(" as " + ' \"' + "band_" + input_name + '\"' + "\n")

                else:
                    numline = math.ceil(len(tmp) / limchar)
                    numpart = math.ceil(len(bandname) / numline)

                    for i in range(int(numline)):
                        if i == 0:
                            out.write("\nX Display" + "".join(bandname[:numpart]))
                            out.write(" vs " + input_name + "_Band_kv")
                            out.write(" as " + '\"' + "band_" + input_name + '\"' + "\n")
                        elif len("".join(bandname[(i * numpart):(i * numpart + numpart)])) == 0:
                            pass
                        else:
                            out.write("X AppendToGraph" + "".join(bandname[(i * numpart):(i * numpart + numpart)]))
                            out.write(" vs " + input_name + "_Band_kv\n")

            else:
                if len(bandname) <= limchar:
                    out.write(str(tmp))
                    out.write(" vs " + input_name + "_Band_kv")
                    out.write(" as " + ' \"' + "band_" + figurename + '\"' + "\n")

                else:
                    numline = math.ceil(len(tmp) / limchar)
                    numpart = math.ceil(len(bandname) / numline)

                    for i in range(int(numline)):
                        if i == 0:
                            out.write("\nX Display" + "".join(bandname[:numpart]))
                            out.write(" vs " + input_name + "_Band_kv")
                            out.write(" as " + '\"' + "band_" + figurename + '\"' + "\n")
                        elif len("".join(bandname[(i * numpart):(i * numpart + numpart)])) == 0:
                            pass
                        else:
                            out.write("X AppendToGraph" + "".join(bandname[(i * numpart):(i * numpart + numpart)]))
                            out.write(" vs " + input_name + "_Band_kv\n")

        def displayprojection(bandname, pbandname, minE=0, maxE=1, color='Red'):
            preset_proj = ("\nX ModifyGraph mode(%s)=0,lsize=2,zColor(%s)={%s,%f,%f,%s,1}"
                               % (bandname, bandname, pbandname, minE, maxE, color))
            # if downspin is True:
            #     preset_proj = ("\nX ModifyGraph mode(%s)=0,lsize=2,zColor(%s)={%s,-1,0,%s,0}"
            #                    % (bandname, bandname, pbandname, color))

            out.write(preset_proj)

        def displaypreset():
            preset = ("X DefaultFont/U \"Times New Roman\"\n"
                      "X ModifyGraph width=255.118,height=340.157\n"
                      "X ModifyGraph marker=19\n"
                      "X ModifyGraph lSize=1.5\n"
                      "X ModifyGraph tick(left)=2,tick(bottom)=3,noLabel(bottom)=2\n"
                      "X ModifyGraph mirror=1\n"
                      "X ModifyGraph zero(left)=8\n"
                      "X ModifyGraph fSize=28\n"
                      "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                      "X ModifyGraph standoff=0\n"
                      "X ModifyGraph axThick=1.5\n"
                      "X ModifyGraph axisOnTop=1\n"
                      "X Label left \"\Z28 Energy (eV)\"\n"
                      "X ModifyGraph zero(bottom)=0;DelayUpdate\n"
                      "X SetAxis left -3,3\n"
                      "X ModifyGraph zeroThick(left)=2.5\n")

            preset_guide = ("\nX AppendToGraph " + input_name + "_guide_y1 " + input_name + "_guide_y2 vs " +
                            input_name + "_guide\n"
                            "X ModifyGraph mode(" + input_name + "_guide_y1)=1,rgb(" + input_name + "_guide_y1)=(0,0,0)\n"
                            "X ModifyGraph mode(" + input_name + "_guide_y2)=1,rgb(" + input_name + "_guide_y2)=(0,0,0)\n"
                            "X SetAxis left -3,3"
                            )

            out.write(preset)
            if guide is True:
                out.write(preset_guide)

        # Writing the itx file
        wavename = bandname()
        pbandkind = pbandname()
        pwavetitle = pwavename(pbandkind)
        itxheader(wavename)
        bandwave()
        guidewave()
        # displayband(wavename)
        # displaypreset()

        if spin is True:
            spinmzwave(pbandkind)
            for i in range(len(pbandkind[0])):
                for j in range(len(pbandkind[1])):
                    figname = figurename(pbandkind[0][i], pbandkind[1][j])
                    displayband(wavename, str(figname + "_mz"))
                    displaypreset()
                    for k in range(len(wavename)):
                        mzwave = "Proj_" + input_name + "_Band_" + str(wavename[k].split('_')[-1]) + "_" \
                                 + str(pbandkind[0][i]) + "_" + str(pbandkind[1][j] + "_mz")
                        mzwave = mzwave.replace("__", "_")
                        displayprojection(wavename[k].strip(), mzwave, -0.5, 0.5, 'RedWhiteBlue')

        elif spin is False:
            pbandwave(pbandkind, pwavetitle)
            for i in range(len(pbandkind[0])):
                for j in range(len(pbandkind[1])):
                    figname = figurename(pbandkind[0][i], pbandkind[1][j])
                    displayband(wavename, figname)
                    displaypreset()
                    for k in range(len(wavename)):
                        pwave = "Proj_" + input_name + "_Band_" + str(wavename[k].split('_')[-1]) + "_" \
                                + str(pbandkind[0][i]) + "_" + str(pbandkind[1][j])
                        pwave = pwave.replace("__", "_")
                        if pwave[-1] == "_":
                            pwave = pwave[:-1]

                        displayprojection(wavename[k].strip(), pwave)

        return

    def temp_pband(self, canvas, axes):

        if self.outfile is None:
            outfile = "tband.itx"
        else:
            outfile = self.outfile + ".itx"

        canvas = canvas.T
        kgrid = np.insert(axes[0], np.shape(axes[0]), (2 * axes[0][-1] - axes[0][-2]))
        egrid = np.insert(axes[1], np.shape(axes[1]), (2 * axes[1][-1] - axes[1][-2]))

        name = []
        for i in range(len(canvas)):
            name.append(self.input_name + "_grid_" + str(i + 1))

        nameheader = (self.input_name + "_grid_")
        figuretitle = ("tBand_" + self.input_name)

        with open(outfile, 'w') as out:
            waveheader(out, name, True)
            for i in range(len(canvas[0])):
                for j in range(len(canvas)):
                    # out.write(':.6e'.format(canvas[j][i]))
                    out.write("%.6e" % (canvas[j][i]))
                    out.write(" ")
                out.write("\n")
            out.write("END\n")

            out.write("WAVES/D ")
            out.write(self.input_name + "_kvec\n")
            out.write("BEGIN\n")
            for i in range(len(kgrid)):
                out.write(str(kgrid[i]) + "\n")
            out.write("END\n")

            out.write("WAVES/D ")
            out.write(self.input_name + "_energy\n")
            out.write("BEGIN\n")
            for i in range(len(egrid)):
                out.write(str(egrid[i]) + "\n")
            out.write("END\n")

            out.write("X SetDataFolder root:\n")
            out.write("X Concatenate WaveList(\"%s*\", \";\", \"\"), %s\n" % (nameheader, figuretitle))
            # out.write("X NewImage/K=0  root:%s" % figuretitle)
            out.write("X Display as \"%s\"\n" % figuretitle)
            out.write("X AppendImage %s vs {%s, %s}\n" %
                      (figuretitle, (self.input_name + "_kvec"), (self.input_name + "_energy")))
            out.write("X ModifyImage %s ctab={0,0.05,Grays,1}\n" % figuretitle)

            graphpreset(out, None, None, None, None, 255.118, 340.157)


    def plotdielectric(self, direction, toplot, drude, plasmasq, tau):
        from Parsers import xmlParserV as xml

        np.seterr(divide='ignore', invalid='ignore')

        if self.outfile is None:
            self.outfile = 'optics.itx'
        else:
            self.outfile = self.outfile + '.itx'

        self.io = IO(None, self.outfile)
        out = self.io.WriteFile()

        if toplot is None:
            toplot = 'eI'
        else:
            pass

        def wavename(waves, suffix=''):
            wavename = []
            for x in waves:
                wavename.append(str(x + "_" + suffix + "_" + self.input_name + " ").replace("__", "_"))
            return wavename

        def itxheader(wavename):
            # Writing the header part of .itx
            out.write("WAVES/D ")
            for x in wavename:
                out.write(str(x))
            out.write("\nBEGIN\n")

        def itxwave(data):
            for x in data.T:
                for y in x:
                    out.write(str(y) + " ")
                out.write("\n")
            out.write("END\n")

        def itxdisplay(parameter, wavename, xAxis=0):
            if parameter == "eI":
                param = 0
            elif parameter == "eR":
                param = 1
            elif parameter == "n":
                param = 2
            elif parameter == "k":
                param = 3
            elif parameter == "alpha":
                param = 4
            elif parameter == "ELS":
                param = 5
            elif parameter == "R":
                param = 6
            elif parameter == "T":
                param = 7
            elif parameter == "A":
                param = 8
            else:
                param = None

            if param is not None:
                out.write("X Display " + str(wavename[param]).rstrip(" ") + " vs " + str(wavename[xAxis + 9].rstrip(" ")
                          + " as " + '"' + str(wavename[param]).rstrip(" ") + '"' + "\n"))

                if xAxis == 0:
                    label_bot = ("X Label bottom " + '"' + "\Z28E (eV)" + '"' + "\n")
                elif xAxis == 1:
                    label_bot = ("X Label bottom " + '"' + "\Z28\F'Symboll\F'Times New Roman (nm)" + '"' + "\n")

                label_left = ("X Label left " + '"' + "\Z28" + str(wavename[param]).rstrip(" ") + '"' + "\n")

                preset = ("X DefaultFont/U \"Times New Roman\"\n"
                          "X ModifyGraph width=340.157,height=226.772\n"
                          "X ModifyGraph marker=19\n"
                          "X ModifyGraph lSize=1.5\n"
                          "X ModifyGraph tick=2\n"
                          "X ModifyGraph mirror=1\n"
                          "X ModifyGraph zero(bottom)=8\n"
                          "X ModifyGraph fSize=28\n"
                          "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                          "X ModifyGraph standoff=0\n"
                          "X ModifyGraph axThick=1.5\n"
                          "X ModifyGraph axisOnTop=1\n"
                          "X Setaxis bottom 0,5\n"
                          )
                out.write(preset)
                out.write(label_left)
                out.write(label_bot)

            else:
                pass

        def coeff_calc(imag, real, freq):
            n_index = np.sqrt((np.sqrt(imag ** 2 + real ** 2) + real) / 2.0)
            k_coeff = np.sqrt((np.sqrt(imag ** 2 + real ** 2) - real) / 2.0)
            alpha = (4 * np.pi * k_coeff) / freq
            ELS = imag / (real ** 2 + imag ** 2)
            R = ((n_index - 1) ** 2 + k_coeff ** 2) / ((n_index + 1) ** 2 + k_coeff ** 2)
            T = 1 - R
            A = np.log10(1 / T)

            return [n_index, k_coeff, alpha, ELS, R, T, A]

        data = xml.xmlParserV().dielectric

        print("-----------------------------------------------------------------------------------------------")
        print("Plotting dielectric response...")
        # input_name = input("Please input the name of system  :  ")

        di_imag = np.array(data[0])
        di_real = np.array(data[1])

        # Dielectric function is arranged like:
        # E / XX / YY / ZZ / XY / YZ / ZX

        freq = di_imag[:, 0]
        freq_nm = 1240.0 / freq

        if drude:
            print("Calculating Drude contributions in dielectric responses...")

            eI_drude = (tau * plasmasq) / (freq * (freq ** 2 + tau ** 2))
            eR_drude = 1 - plasmasq / (freq * (freq ** 2 + tau ** 2))

            for i in range(3):
                di_imag[:, i+1] = di_imag[:, i+1] + eI_drude
                di_real[:, i+1] = di_real[:, i+1] + eR_drude

        if direction is False:
            eI = np.average(di_imag[:, [1, 2, 3]], 1)
            eR = np.average(di_real[:, [1, 2, 3]], 1)

            coeffs = coeff_calc(eI, eR, freq)
            waves = ['eI', 'eR', 'n', 'k', 'alpha', 'ELS', 'R', 'T', 'A', 'freq', 'freq_nm']

            toitx = [eI, eR]

            for x in coeffs:
                toitx.append(x)

            toitx.append(freq)
            toitx.append(freq_nm)

            wave = wavename(waves)
            out.write("IGOR\n")
            itxheader(wave)
            itxwave(np.array(toitx))

            if len(toplot) > 1:
                for x in toplot:
                    itxdisplay(x, wave)
            else:
                for x in toplot:
                    itxdisplay(x, wave)

        else:
            eI_trans = np.average(di_imag[:, [1, 2]], 1)
            eI_longi = np.average(di_imag[:, [3]], 1)
            eR_trans = np.average(di_real[:, [1, 2]], 1)
            eR_longi = np.average(di_real[:, [3]], 1)

            coeffs_trans = coeff_calc(eI_trans, eR_trans, freq)
            coeffs_longi = coeff_calc(eI_longi, eR_longi, freq)

            waves = ['eI', 'eR', 'n', 'k', 'alpha', 'ELS', 'R', 'T', 'A', 'freq', 'freq_nm']

            toitx = [eI_trans, eR_trans]

            for x in coeffs_trans:
                toitx.append(x)

            toitx.append(freq)
            toitx.append(freq_nm)

            wave = wavename(waves, 'trans')
            out.write("IGOR\n")
            itxheader(wave)
            itxwave(np.array(toitx))

            for x in toplot.split():
                itxdisplay(x, wave)

            toitx = [eI_longi, eR_longi]

            for x in coeffs_trans:
                toitx.append(x)

            toitx.append(freq)
            toitx.append(freq_nm)

            wave = wavename(waves, 'longi')
            itxheader(wave)
            itxwave(np.array(toitx))

            if len(toplot) > 1:
                for x in toplot:
                    itxdisplay(x, wave)
            else:
                for x in toplot:
                    itxdisplay(x, wave)

        print("Done!")
        return

    def plotdos(self, fermi, combine):

        from ElectronicStructure.dos import Densityofstates
        self.input_name = "_" + self.input_name
        self.dosdata = Densityofstates()
        self.dosdata.dos(fermi)
        atom = self.dosdata.atominfo

        print("-----------------------------------------------------------------------------------------------")
        print("Plotting density-of-states...")

        # Defining the wave name
        tdosname = None
        tdosdisplay = None
        tdosdisplayforcomb = None
        pdosname = []
        pdosdisplay = []
        sumdosname = []
        sumdosdisplay = []
        sumdosappend = []

        if self.dosdata.numtdos == 3:
            tdosname = ("E_tDOS" + self.input_name
                        + " tDOS" + self.input_name
                        + " intDOS" + self.input_name)
            tdosdisplay = ("X Display"
                           + " tDOS" + self.input_name + " vs"
                           + " E_tdos" + self.input_name + " as"
                           + ' \"' + "tDOS" + self.input_name + '\"' + "\n")
            tdosdisplayforcomb = ("X Display"
                                  + " tDOS" + self.input_name + " vs"
                                  + " E_tdos" + self.input_name + " as"
                                  + ' \"' + "DOS" + self.input_name + '\"' + "\n")

        elif self.dosdata.numtdos == 5:
            tdosname = ("E_tDOS" + self.input_name
                        + " tDOS_up" + self.input_name
                        + " tDOS_dw" + self.input_name
                        + " intDOS_up" + self.input_name
                        + " intDOS_dw" + self.input_name)
            tdosdisplay = ("X Display"
                           + " tDOS_up" + self.input_name
                           + " tDOS_dw" + self.input_name + " vs"
                           + " E_tdos" + self.input_name + " as"
                           + ' \"' + "tDOS" + self.input_name + '\"' + "\n")
            tdosdisplayforcomb = ("X Display"
                                  + " tDOS_up" + self.input_name
                                  + " tDOS_dw" + self.input_name + " vs"
                                  + " E_tdos" + self.input_name + " as"
                                  + ' \"' + "DOS" + self.input_name + '\"' + "\n")

        if self.dosdata.numpdos == 4:
            for x in atom.keys():
                sumname = ("E_sumDOS_" + x + "" + self.input_name + " "
                           + "sumDOS_" + x + "_s" + self.input_name + " "
                           + "sumDOS_" + x + "_p" + self.input_name + " "
                           + "sumDOS_" + x + "_d" + self.input_name + " "
                           + "sumDOS_" + x + "" + self.input_name + " "
                           )
                sumdis = ("X Display "
                          + "sumDOS_" + x + "_s" + self.input_name + " "
                          + "sumDOS_" + x + "_p" + self.input_name + " "
                          + "sumDOS_" + x + "_d" + self.input_name + " "
                          + "sumDOS_" + x + "" + self.input_name + " "
                          + "vs "
                          + "E_sumDOS_" + x + "" + self.input_name + " "
                          + "as "
                          + ' \"' + "sumDOS_" + x + "" + self.input_name + '\"' + "\n"
                          )
                sumapp = ("X AppendToGraph "
                          + "sumDOS_" + x + "_s" + self.input_name + " "
                          + "sumDOS_" + x + "_p" + self.input_name + " "
                          + "sumDOS_" + x + "_d" + self.input_name + " "
                          + "sumDOS_" + x + "" + self.input_name + " "
                          + "vs "
                          + "E_sumDOS_" + x + "" + self.input_name + "\n"
                          )

                for y in range(int(atom[x])):
                    tmpname = ("E_pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_s" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_p" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_d" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                               )
                    tmpdisplay = ("X Display "
                                  + "pDOS_" + x + str(y + 1) + "_s" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_p" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_d" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                                  + "vs "
                                  + "E_pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                                  + "as "
                                  + ' \"' + "pDOS_" + x + str(y + 1) + "" + self.input_name + '\"' + "\n"
                                  )
                    pdosname.append(tmpname)
                    pdosdisplay.append(tmpdisplay)
                sumdosname.append(sumname)
                sumdosdisplay.append(sumdis)
                sumdosappend.append(sumapp)

        elif self.dosdata.numpdos == 7:
            for x in atom.keys():
                sumname = ("E_sumDOS_" + x + "" + self.input_name + " "
                           + "sumDOS_" + x + "_s_up" + self.input_name + " "
                           + "sumDOS_" + x + "_s_dw" + self.input_name + " "
                           + "sumDOS_" + x + "_p_up" + self.input_name + " "
                           + "sumDOS_" + x + "_p_dw" + self.input_name + " "
                           + "sumDOS_" + x + "_d_up" + self.input_name + " "
                           + "sumDOS_" + x + "_d_dw" + self.input_name + " "
                           + "sumDOS_" + x + "" + self.input_name + " "
                           )
                sumdis = ("X Display "
                          + "sumDOS_" + x + "_s_up" + self.input_name + " "
                          + "sumDOS_" + x + "_s_dw" + self.input_name + " "
                          + "sumDOS_" + x + "_p_up" + self.input_name + " "
                          + "sumDOS_" + x + "_p_dw" + self.input_name + " "
                          + "sumDOS_" + x + "_d_up" + self.input_name + " "
                          + "sumDOS_" + x + "_d_dw" + self.input_name + " "
                          + "sumDOS_" + x + "" + self.input_name + " "
                          + "vs "
                          + "E_sumDOS_" + x + "" + self.input_name + " "
                          + "as "
                          + ' \"' + "sumDOS_" + x + "" + self.input_name + '\"' + "\n"
                          )
                sumapp = ("X AppendToGraph "
                          + "sumDOS_" + x + "_s_up" + self.input_name + " "
                          + "sumDOS_" + x + "_s_dw" + self.input_name + " "
                          + "sumDOS_" + x + "_p_up" + self.input_name + " "
                          + "sumDOS_" + x + "_p_dw" + self.input_name + " "
                          + "sumDOS_" + x + "_d_up" + self.input_name + " "
                          + "sumDOS_" + x + "_d_dw" + self.input_name + " "
                          + "sumDOS_" + x + "" + self.input_name + " "
                          + "vs "
                          + "E_sumDOS_" + x + "" + self.input_name + "\n"
                          )

                for y in range(int(atom[x])):
                    tmpname = ("E_pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_s_up" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_s_dw" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_p_up" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_p_dw" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_d_up" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_d_dw" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                               )
                    tmpdisplay = ("X Display "
                                  + "pDOS_" + x + str(y + 1) + "_s_up" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_s_dw" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_p_up" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_p_dw" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_d_up" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_d_dw" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                                  + "vs "
                                  + "E_pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                                  + "as "
                                  + ' \"' + "pDOS_" + x + str(y + 1) + "" + self.input_name + '\"' + "\n"
                                  )
                    pdosname.append(tmpname)
                    pdosdisplay.append(tmpdisplay)
                sumdosname.append(sumname)
                sumdosdisplay.append(sumdis)
                sumdosappend.append(sumapp)

        elif self.dosdata.numpdos == 10:
            for x in atom.keys():
                sumname = ("E_sumDOS_" + x + "" + self.input_name + " "
                           + "sumDOS_" + x + "_s" + self.input_name + " "
                           + "sumDOS_" + x + "_px" + self.input_name + " "
                           + "sumDOS_" + x + "_py" + self.input_name + " "
                           + "sumDOS_" + x + "_pz" + self.input_name + " "
                           + "sumDOS_" + x + "_dxy" + self.input_name + " "
                           + "sumDOS_" + x + "_dyz" + self.input_name + " "
                           + "sumDOS_" + x + "_dz2" + self.input_name + " "
                           + "sumDOS_" + x + "_dxz" + self.input_name + " "
                           + "sumDOS_" + x + "_dx2" + self.input_name + " "
                           + "sumDOS_" + x + "" + self.input_name + " "
                           )
                sumdis = ("X Display "
                          + "sumDOS_" + x + "_s" + self.input_name + " "
                          + "sumDOS_" + x + "_px" + self.input_name + " "
                          + "sumDOS_" + x + "_py" + self.input_name + " "
                          + "sumDOS_" + x + "_pz" + self.input_name + " "
                          + "sumDOS_" + x + "_dxy" + self.input_name + " "
                          + "sumDOS_" + x + "_dyz" + self.input_name + " "
                          + "sumDOS_" + x + "_dz2" + self.input_name + " "
                          + "sumDOS_" + x + "_dxz" + self.input_name + " "
                          + "sumDOS_" + x + "_dx2" + self.input_name + " "
                          + "sumDOS_" + x + "" + self.input_name + " "
                          + "vs "
                          + "E_sumDOS_" + x + "" + self.input_name + " "
                          + "as "
                          + ' \"' + "sumDOS_" + x + "" + self.input_name + '\"' + "\n"
                          )
                sumapp = ("X AppendToGraph "
                          + "sumDOS_" + x + "_s" + self.input_name + " "
                          + "sumDOS_" + x + "_px" + self.input_name + " "
                          + "sumDOS_" + x + "_py" + self.input_name + " "
                          + "sumDOS_" + x + "_pz" + self.input_name + " "
                          + "sumDOS_" + x + "_dxy" + self.input_name + " "
                          + "sumDOS_" + x + "_dyz" + self.input_name + " "
                          + "sumDOS_" + x + "_dz2" + self.input_name + " "
                          + "sumDOS_" + x + "_dxz" + self.input_name + " "
                          + "sumDOS_" + x + "_dx2" + self.input_name + " "
                          + "sumDOS_" + x + "" + self.input_name + " "
                          + "vs "
                          + "E_sumDOS_" + x + "" + self.input_name + "\n"
                          )
                for y in range(int(atom[x])):
                    tmpname = ("E_pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_s" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_px" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_py" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_pz" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_dxy" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_dyz" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_dz2" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_dxz" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_dx2" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                               )
                    tmpdisplay = ("X Display "
                                  + "pDOS_" + x + str(y + 1) + "_s" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_px" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_py" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_pz" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_dxy" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_dyz" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_dz2" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_dxz" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_dx2" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                                  + "vs "
                                  + "E_pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                                  + "as "
                                  + ' \"' + "pDOS_" + x + str(y + 1) + "" + self.input_name + '\"' + "\n"
                                  )
                    pdosname.append(tmpname)
                    pdosdisplay.append(tmpdisplay)
                sumdosname.append(sumname)
                sumdosdisplay.append(sumdis)
                sumdosappend.append(sumapp)

        elif self.dosdata.numpdos == 13:
            for x in atom.keys():
                sumname = ("E_sumDOS_" + x + "" + self.input_name + " "
                           + "sumDOS_" + x + "_s_t" + self.input_name + " "
                           + "sumDOS_" + x + "_s_mx" + self.input_name + " "
                           + "sumDOS_" + x + "_s_my" + self.input_name + " "
                           + "sumDOS_" + x + "_s_mz" + self.input_name + " "
                           + "sumDOS_" + x + "_p_t" + self.input_name + " "
                           + "sumDOS_" + x + "_p_mx" + self.input_name + " "
                           + "sumDOS_" + x + "_p_my" + self.input_name + " "
                           + "sumDOS_" + x + "_p_mz" + self.input_name + " "
                           + "sumDOS_" + x + "_d_t" + self.input_name + " "
                           + "sumDOS_" + x + "_d_mx" + self.input_name + " "
                           + "sumDOS_" + x + "_d_my" + self.input_name + " "
                           + "sumDOS_" + x + "_d_mz" + self.input_name + " "
                           + "sumDOS_" + x + "" + self.input_name + " "
                           )
                sumdis = ("X Display "
                          + "sumDOS_" + x + "_s_t" + self.input_name + " "
                          + "sumDOS_" + x + "_s_mx" + self.input_name + " "
                          + "sumDOS_" + x + "_s_my" + self.input_name + " "
                          + "sumDOS_" + x + "_s_mz" + self.input_name + " "
                          + "sumDOS_" + x + "_p_t" + self.input_name + " "
                          + "sumDOS_" + x + "_p_mx" + self.input_name + " "
                          + "sumDOS_" + x + "_p_my" + self.input_name + " "
                          + "sumDOS_" + x + "_p_mz" + self.input_name + " "
                          + "sumDOS_" + x + "_d_t" + self.input_name + " "
                          + "sumDOS_" + x + "_d_mx" + self.input_name + " "
                          + "sumDOS_" + x + "_d_my" + self.input_name + " "
                          + "sumDOS_" + x + "_d_mz" + self.input_name + " "
                          + "sumDOS_" + x + "" + self.input_name + " "
                          + "vs "
                          + "E_sumDOS_" + x + "" + self.input_name + " "
                          + "as "
                          + ' \"' + "sumDOS_" + x + "" + self.input_name + '\"' + "\n"
                          )
                sumapp = ("X Display "
                          + "sumDOS_" + x + "_s_t" + self.input_name + " "
                          + "sumDOS_" + x + "_s_mx" + self.input_name + " "
                          + "sumDOS_" + x + "_s_my" + self.input_name + " "
                          + "sumDOS_" + x + "_s_mz" + self.input_name + " "
                          + "sumDOS_" + x + "_p_t" + self.input_name + " "
                          + "sumDOS_" + x + "_p_mx" + self.input_name + " "
                          + "sumDOS_" + x + "_p_my" + self.input_name + " "
                          + "sumDOS_" + x + "_p_mz" + self.input_name + " "
                          + "sumDOS_" + x + "_d_t" + self.input_name + " "
                          + "sumDOS_" + x + "_d_mx" + self.input_name + " "
                          + "sumDOS_" + x + "_d_my" + self.input_name + " "
                          + "sumDOS_" + x + "_d_mz" + self.input_name + " "
                          + "sumDOS_" + x + "" + self.input_name + " "
                          + "vs "
                          + "E_sumDOS_" + x + "" + self.input_name + "\n"
                          )

                for y in range(int(atom[x])):
                    tmpname = ("E_pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_s_t" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_s_mx" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_s_my" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_s_mz" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_p_t" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_p_mx" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_p_my" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_p_mz" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_d_t" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_d_mx" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_d_my" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "_d_mz" + self.input_name + " "
                               + "pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                               )
                    tmpdisplay = ("X Display "
                                  + "pDOS_" + x + str(y + 1) + "_s_t" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_s_mx" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_s_my" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_s_mz" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_p_t" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_p_mx" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_p_my" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_p_mz" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_d_t" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_d_mx" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_d_my" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "_d_mz" + self.input_name + " "
                                  + "pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                                  + "vs "
                                  + "E_pDOS_" + x + str(y + 1) + "" + self.input_name + " "
                                  + "as "
                                  + ' \"' + "pDOS_" + x + str(y + 1) + "" + self.input_name + '\"' + "\n"
                                  )
                    pdosname.append(tmpname)
                    pdosdisplay.append(tmpdisplay)
                sumdosname.append(sumname)
                sumdosdisplay.append(sumdis)
                sumdosappend.append(sumapp)

        # Defining the graph preset
        preset = ("X ModifyGraph width=340.157,height=170.0785\n"
                  "X ModifyGraph marker=19\n"
                  "X ModifyGraph lSize=1.5\n"
                  "X ModifyGraph tick=2\n"
                  "X ModifyGraph mirror=1\n"
                  "X ModifyGraph zero(bottom)=8\n"
                  "X ModifyGraph fSize=28\n"
                  "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                  "X ModifyGraph standoff=0\n"
                  "X ModifyGraph axThick=1.5\n"
                  "X ModifyGraph axisOnTop=1\n"
                  "X Label left \"\\Z28 Density-of-states (arb. unit)\"\n"
                  "X Label bottom \"\\Z28 Energy (eV)\"\n"
                  "X SetAxis bottom -4,4\n"
                  "X DefaultFont/U \"Times New Roman\"\n"
                  "X ModifyGraph noLabel(left)=1,lblMargin(left)=30\n"
                  "X ModifyGraph axThick=2.5\n")

        # When writing in a separate itx file
        if not combine:

            # tDOS writing
            # Header part
            if self.outfile is None:
                outfile = 'tDOS.itx'
            else:
                outfile = 'tDOS_' + self.outfile + '.itx'
            self.io = IO(None, outfile)
            tdos = self.io.WriteFile()

            tdos.write("IGOR\n")
            tdos.write("WAVES/D ")
            tdos.write(tdosname)
            tdos.write("\nBEGIN\n")

            # Wave part
            for i in range(self.dosdata.nedos):
                tdos.write(str(self.dosdata.tdos[i, 0]))
                for j in range(self.dosdata.numtdos):
                    if j == 0:
                        pass
                    else:
                        tdos.write(" " + str(self.dosdata.tdos[i, j]))
                tdos.write("\n")
            tdos.write("END\n")

            tdos.write(tdosdisplay)
            tdos.write(preset)

            # pDOS writing
            count = 0
            for x in atom.keys():
                for y in range(int(atom[x])):
                    if self.outfile is None:
                        outfile = 'pDOS_' + x + str(y + 1) + ".itx"
                    else:
                        outfile = 'pDOS_' + x + str(y + 1) + "_" + self.outfile + '.itx'
                    self.io = IO(None, outfile)
                    pdos = self.io.WriteFile()

                    # Header part
                    pdos.write("IGOR\n")
                    pdos.write("WAVES/D ")
                    pdos.write(pdosname[count])
                    pdos.write("\nBEGIN\n")

                    # Wave part
                    for i in range(self.dosdata.nedos):
                        pdos.write(str(self.dosdata.pdos[count][i, 0]))
                        for j in range(self.dosdata.numpdos + 1):
                            if j == 0:
                                pass
                            else:
                                pdos.write(" " + str(self.dosdata.pdos[count][i, j]))
                        pdos.write("\n")
                    pdos.write("END\n")

                    pdos.write(pdosdisplay[count])
                    pdos.write(preset)

                    count += 1

            # sumDOS writing
            count = 0
            for x in atom.keys():
                if self.outfile is None:
                    outfile = 'sumDOS_' + x + ".itx"
                else:
                    outfile = 'sumDOS_' + x + "_" + self.outfile + '.itx'
                self.io = IO(None, outfile)
                sumdos = self.io.WriteFile()

                # Header part
                sumdos.write("IGOR\n")
                sumdos.write("WAVES/D ")
                sumdos.write(sumdosname[count])
                sumdos.write("\nBEGIN\n")

                # Wave part
                for i in range(self.dosdata.nedos):
                    sumdos.write(str(self.dosdata.sumdos[count][i, 0]))
                    for j in range(self.dosdata.numpdos + 1):
                        if j == 0:
                            pass
                        else:
                            sumdos.write(" " + str(self.dosdata.sumdos[count][i, j]))
                    sumdos.write("\n")
                sumdos.write("END\n")

                sumdos.write(sumdosdisplay[count])
                sumdos.write(preset)
                count += 1

        # When writing in a combined itx file
        else:
            if self.outfile is None:
                outfile = 'dos.itx'
            else:
                outfile = 'dos_' + self.outfile + '.itx'

            self.io = IO(None, outfile)
            dos = self.io.WriteFile()

            # # Header part
            # dos.write("IGOR\n")
            # dos.write("WAVES/D ")
            # dos.write(tdosname + " " + "".join(pdosname) + " " + " ".join(sumdosname))
            # dos.write("\nBEGIN\n")
            #
            # # Wave part
            # for i in range(self.dosdata.nedos):
            #     dos.write(str(self.dosdata.tdos[i, 0]))
            #     for j in range(self.dosdata.numtdos):
            #         if j == 0:
            #             pass
            #         else:
            #             dos.write(" " + str(self.dosdata.tdos[i, j]))
            #
            #     count = 0
            #     for x in atom.keys():
            #         for y in range(int(atom[x])):
            #             for j in range(self.dosdata.numpdos + 1):
            #                 dos.write(" " + str(self.dosdata.pdos[count][i, j]))
            #         count += 1
            #
            #     count = 0
            #     for k in range(len(atom.keys())):
            #         for j in range(self.dosdata.numpdos + 1):
            #             dos.write(" " + str(self.dosdata.sumdos[count][i, j]))
            #         count += 1
            #
            #     dos.write("\n")
            # dos.write("END\n")

            dos.write("IGOR\n")
            dos.write("WAVES/D ")
            dos.write(tdosname)
            dos.write("\nBEGIN\n")

            # Wave part
            for i in range(self.dosdata.nedos):
                dos.write(str(self.dosdata.tdos[i, 0]))
                for j in range(self.dosdata.numtdos):
                    if j == 0:
                        pass
                    else:
                        dos.write(" " + str(self.dosdata.tdos[i, j]))
                dos.write("\n")
            dos.write("END\n")

            # pDOS writing
            count = 0
            for x in atom.keys():
                for y in range(int(atom[x])):

                    # Header part
                    dos.write("WAVES/D ")
                    dos.write(pdosname[count])
                    dos.write("\nBEGIN\n")

                    # Wave part
                    for i in range(self.dosdata.nedos):
                        dos.write(str(self.dosdata.pdos[count][i, 0]))
                        for j in range(self.dosdata.numpdos + 1):
                            if j == 0:
                                pass
                            else:
                                dos.write(" " + str(self.dosdata.pdos[count][i, j]))
                        dos.write("\n")
                    dos.write("END\n")
                    count += 1

            # sumDOS writing
            for k in range(int(len(atom.keys()))):
                # Header part
                dos.write("WAVES/D ")
                dos.write(sumdosname[k])
                dos.write("\nBEGIN\n")

                # Wave part
                for i in range(self.dosdata.nedos):
                    dos.write(str(self.dosdata.sumdos[k][i, 0]))
                    for j in range(self.dosdata.numpdos + 1):
                        if j == 0:
                            pass
                        else:
                            dos.write(" " + str(self.dosdata.sumdos[k][i, j]))
                    dos.write("\n")
                dos.write("END\n")

            dos.write(tdosdisplayforcomb)
            for i in range(len(sumdosappend)):
                dos.write(sumdosappend[i])
            dos.write(preset)

        print("Done!")
        return


def waveheader(outfile, wavetitle, filestart=False):
    # Writing the header part of .itx
    if filestart is True:
        outfile.write("IGOR\n")
    outfile.write("WAVES/D ")
    for x in wavetitle:
        outfile.write(str(x))
        outfile.write(" ")
    outfile.write("\nBEGIN\n")


def wavename(waves, name, prefix='', suffix=''):
    wavetitle = []
    for x in waves:
        title = None
        if prefix is not '':
            if suffix is not '':
                title = str(name + "_" + prefix + "_" + x + "_" + suffix)
            elif suffix is '':
                title = str(name + "_" + prefix + "_" + x)
        elif prefix is '':
            if suffix is not '':
                title = str(name + "_" + x + "_" + suffix)
            elif suffix is '':
                title = str(name + "_" + x)
        wavetitle.append(title)
    return wavetitle


def graphpreset(outfile, leftmin, leftmax, botmin, botmax, width=340.157, height=226.772):
    # Writing the header part of .itx
    preset = ("X DefaultFont/U \"Times New Roman\"\n"
              "X ModifyGraph marker=19\n"
              "X ModifyGraph lSize=1.5\n"
              "X ModifyGraph tick=2\n"
              "X ModifyGraph mirror=1\n"
              "X ModifyGraph fSize=28\n"
              "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
              "X ModifyGraph standoff=0\n"
              "X ModifyGraph axThick=1.5\n"
              "X ModifyGraph axisOnTop=1\n")

    outfile.write(preset)

    outfile.write("X ModifyGraph width=%f,height=%f\n" % (width, height))

    if leftmin is not None:
        leftaxis = ("X Setaxis left %4f, %4f\n" % (leftmin, leftmax))
        outfile.write(leftaxis)
    if botmin is not None:
        botaxis = ("X Setaxis bottom %4f, %4f\n" % (botmin, botmax))
        outfile.write(botaxis)
