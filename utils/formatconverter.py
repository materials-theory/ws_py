from polyhedron.polyhandler import PolyhedronV as poly
from utils import operator
from ase import io as aio


class ConvertV(object):

    def __init__(self, infile=None, outfile=None):
        self.infile = infile
        self.outfile = outfile
        self.utils = operator.StrucOperators
        # self.IO = IO.IO.IO(infile, outfile)

    def cif2vasp(self):
        atoms = aio.read(self.infile)
        aio.write(self.outfile + ".vasp", atoms, format="vasp")
        with open(self.outfile + ".vasp", "r") as out:
            tmp = out.readlines()
        with open(self.outfile + ".vasp", "w") as out:
            tmp.insert(5, tmp[0])
            out.write("".join(tmp))
        return

    def vasp2qe(self):
        atoms = aio.read(self.infile, format="vasp")
        aio.write(self.outfile + ".in", atoms, format="espresso-in")
        return

    def vasp2cif(self):
        atoms = aio.read(self.infile, format="vasp")
        aio.write(self.outfile + "cif", atoms, format="cif")
        return

    def qe2vasp(self):
        atoms = aio.read(self.infile)
        # atoms = aio.read(self.infile, format="espresso-in")
        # atoms = aio.read(self.infile, format="espresso-out")
        aio.write(self.outfile + ".vasp", atoms, format="vasp")
        with open(self.outfile + ".vasp", "r") as out:
            tmp = out.readlines()
        with open(self.outfile + ".vasp", "w") as out:
            tmp.insert(5, tmp[0])
            out.write("".join(tmp))
        return

    # Converting each octahedron in VASP structure file to xyz file for SHAPE program
    # Only for test use!
    def polytoxyz(self, face, length, center, vertex, shape_feeder=False, centroid=False):
        p = poly(self.infile)
        p.setpoly(face, length, center, vertex)
        polydata, centername = p.polyposition(naming=True)

        if self.outfile is None:
            filename = (str(self.infile) + "_converted")
        else:
            filename = (str(self.outfile))

        if centroid:
            filename = filename + "_ideal"
        else:
            filename = filename + "_real"

        if not shape_feeder:
            filename = filename + ".xyz"
        else:
            filename = filename + ".dat"

        with open(filename, "w") as out:
            if not shape_feeder:
                print("Writing .xyz data...")
                out.write(str(len(polydata) * 7) + "\n")
                out.write(str(filename) + "\n")
                for i in range(len(polydata)):
                    for j in range(7):
                        if j != 6:
                            out.write(str(vertex[0]) + " ")
                            for x in polydata[i][j]:
                                out.write(str(x) + " ")
                        else:
                            out.write(str(centername[i]) + " ")
                            if not centroid:
                                for x in polydata[i][j]:
                                    out.write(str(x) + " ")
                            else:
                                centroid_pos = self.utils.centroid(polydata[i][:-1])
                                for x in centroid_pos:
                                    out.write(str(x) + " ")
                        out.write("\n")

            if shape_feeder:
                print("Writing .dat data...")
                out.write("$ " + str(filename) + "\n")
                out.write("%fullout" + "\n")
                out.write("6 7\n")
                out.write("3\n")

                for i in range(len(polydata)):
                    out.write("Oc-6\n")
                    for j in range(7):
                        if j != 6:
                            out.write(str(vertex[0]) + " ")
                            for x in polydata[i][j]:
                                out.write(str(x) + " ")
                        else:
                            out.write(str(centername[i]) + " ")
                            if not centroid:
                                for x in polydata[i][j]:
                                    out.write(str(x) + " ")
                            else:
                                centroid_pos = self.utils.centroid(polydata[i][:-1])
                                for x in centroid_pos:
                                    out.write(str(x) + " ")
                        out.write("\n")

        return
