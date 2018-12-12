#!/usr/local/bin/python

import argparse
import os
import sys
import numpy as np

installpath = '/Users/Woosun/Dropbox/Dev/ws_py'
sys.path.extend([installpath])

from vw_py.generalutils.formatconverter import ConvertV
from vw_py.boltztrap import boltztrap
from vw_py.Plotter import plotter
from vw_py.ElectronicStructure.bandstructure import Bandstructure
from vw_py.polyhedron.polyhandler import PolyhedronV
from vw_py.struc.struc import Structurehandler
from decimal import Decimal


def executepoly(args):
    p = PolyhedronV(args.input)
    p.setpoly(args.face, args.length, args.center, args.vertex)
    result = p.polyanalyze()

    if not args.silent:
        print("-----------------------------------------------------------------------------------------")
        print("Analyzing the polyhedra in %s...                                                     "
              % args.input)
        print("-----------------------------------------------------------------------------------------")
        print("  Center : %2s,  Vertex : %2s,  Maximum bond length : %2f                             "
              % (", ".join(result[5][0][0]), ", ".join(result[5][1][0]), result[5][2]))
        print("-----------------------------------------------------------------------------------------")
        print("   Index   |  Volume  |    QE    |    DI    |    BAV     |                               ")
        print("-----------------------------------------------------------------------------------------")

        for i in range(len(result[0])):
            print(" Poly %3i  | %8f | %8f | %8f | %8f |                               "
                  % (result[0][i], result[1][i], result[2][i], result[3][i], result[4][i]))

        print("-------------------------------------------------------------------------------------------")
        print("Please be aware that the definition of \"bond angle\" can be different!")

        if args.tilt:
            print("-------------------------------------------------------------------------------------------")
            print("Processing the octahedral tilt angle analysis...")
            print("-------------------------------------------------------------------------------------------")
            tiltresult = p.polytiltanalyze(args.centroid, args.vertexonly, args.ideal, args.edge, args.tav)
            if args.centroid:
                print("CENTROID of each octahedra is used to calculate CVCO angle")
            elif not args.centroid:
                print("CENTER METAL ATOM of each octahedra is used to calculate CVCO angle")
            if args.vertexonly:
                print("Center metal atom is EXCLUDED in the calculation of the centroid")
            elif not args.vertexonly:
                print("Center metal atom is INCLUDED in the calculation of the centroid")
            if args.edge:
                print("Edge-sharing octahedra are INCLUDED in the calculation of tilting angles")
            elif not args.edge:
                print("Edge-sharing octahedra are EXCLUDED in the calculation of tilting angles")
            if args.tav:
                print("Tilt angle analysis result is provided in the format of TILT ANGLE VARIANCE")
            elif not args.tav:
                print("Tilt angle analysis result is provided in the format of AVERAGE")
            print("-----------------------------------------------------------------------------------------")
            print("     C-V-C   |    M-V-M    |    V-V-V    |    CVC-V   |                                  ")
            print("-----------------------------------------------------------------------------------------")
            ang = []
            for x in tiltresult:
                if x is not None:
                    ang.append("%8f" % x)
                elif x is None:
                    ang.append("   None   ")

            print("  %8s |  %8s |  %8s |  %8s  |                                  "
                  % (ang[0], ang[1], ang[2], ang[3]))

    elif args.silent or args.output:
        import pandas as pd
        file = []
        index = []
        center = []
        vertex = []
        vol = []
        bav = []
        di = []
        qe = []
        face = []

        if args.output is None:
            outfile = 'polyresult.csv'
        else:
            outfile = args.output + '.csv'
        # out = IO(None, outfile).WriteFile()

        for i in range(len(result[0])):
            index.append(result[0][i])
            vol.append(result[1][i])
            qe.append(result[2][i])
            di.append(result[3][i])
            bav.append(result[4][i])
            file.append(args.input)
            center.append(", ".join(result[5][0][0]))
            vertex.append(", ".join(result[5][1][0]))
            face.append(args.face)

        index.append('average')
        file.append(args.input)
        center.append(", ".join(result[5][0][0]))
        vertex.append(", ".join(result[5][1][0]))
        face.append(args.face)
        qe.append(np.average(qe))
        di.append(np.average(di))
        bav.append(np.average(bav))
        vol.append(np.average(vol))

        if args.tilt:
            cvc = []
            mvm = []
            vvv = []
            cvcv = []
            tiltresult = p.polytiltanalyze(args.centroid, args.vertexonly, args.ideal, args.edge, args.tav)
            ang = []
            for x in tiltresult:
                if x is not None:
                    ang.append("%8f" % x)
                elif x is None:
                    ang.append("   None   ")

            for i in range(len(result[0])):
                cvc.append('NaN')
                mvm.append('NaN')
                vvv.append('NaN')
                cvcv.append('NaN')

            cvc.append(ang[0])
            mvm.append(ang[1])
            vvv.append(ang[2])
            cvcv.append(ang[3])

        if args.tilt:
            data = {'file': file,
                    'index': index,
                    'center': center,
                    'vertex': vertex,
                    'volume': vol,
                    'bav': bav,
                    'di': di,
                    'qe': qe,
                    'face': face,
                    'cvc': cvc,
                    'mvm': mvm,
                    'vvv': vvv,
                    'cvcv': cvcv,
                    }
            df = pd.DataFrame(data, columns=['file', 'index', 'center', 'vertex', 'face',
                                             'volume', 'bav', 'di', 'qe',
                                             'cvc', 'mvm', 'vvv', 'cvcv'])

        else:
            data = {'file': file,
                    'index': index,
                    'center': center,
                    'vertex': vertex,
                    'volume': vol,
                    'bav': bav,
                    'di': di,
                    'qe': qe,
                    'face': face,
                    }
            df = pd.DataFrame(data, columns=['file', 'index', 'center', 'vertex', 'face',
                                             'volume', 'bav', 'di', 'qe'])
        if not os.path.isfile(outfile):
            df.to_csv(outfile, mode='a', index=False)
            print(outfile + " Generated!")
        elif len(df.columns) != len(pd.read_csv(outfile, nrows=1).columns):
            raise Exception(
                "Columns do not match!! Dataframe has " + str(len(df.columns)) + " columns. CSV file has " + str(
                    len(pd.read_csv(outfile, nrows=1).columns)) + " columns.")
        elif not (df.columns == pd.read_csv(outfile, nrows=1).columns).all():
            raise Exception("Columns and column order of dataframe and csv file do not match!!")
        else:
            df.to_csv(outfile, mode='a', index=False, header=False)
            print("Found " + outfile + ", Appended to " + outfile + "!")

    return


def executeelec(args):
    e = Bandstructure()
    e.bandedge(args.fermi)
    print("E-Fermi  : %4f eV" % e.fermi)
    print("Band gap : %4f eV" % e.gap)
    if args.edge is True:
        print("VBM : %4f eV (k-point # %3i, band # %3i)" % (e.vbm[2], e.vbm[0] + 1, e.vbm[1] + 1))
        print("CBM : %4f eV (k-point # %3i, band # %3i)" % (e.cbm[2], e.cbm[0] + 1, e.cbm[1] + 1))

    if args.mass is True:
        e.effectivemass(args.band, args.kpoint, args.points, args.residuals)
        for x in e.emass:
            print("Effective mass of carrier at k-point # %3i, band # %3i to the %5s : %5f" %
                  (args.kpoint, args.band, x[1], float(x[0])))
        if args.residuals is True:
            for x in e.emass_residuals:
                print("Residuals of the least-square fit to the %5s : %4.4E" %
                      (x[1], Decimal(x[0])))
        print("-----------------------------------------------------------------------------------------")

    return


def executeplotband(args):
    p = plotter.Plotter(args.output, args.prefix)
    p.plotband(args.fermi, args.fake, args.shift, args.guide, args.ktraj)
    return


def executeplotpband(args):
    p = plotter.Plotter(args.output, args.prefix)
    p.plotprojectedband(args.fermi, args.fake, args.shift, args.guide, args.spin, args.atom, args.orbital)
    return


def executeplotdos(args):
    p = plotter.Plotter(args.output, args.prefix)
    p.plotdos(args.fermi, args.combine)
    return


def executeplotdielectric(args):
    p = plotter.Plotter(args.output, args.prefix)
    p.plotdielectric(args.direction, args.plot, args.drude, args.plasmasq, args.tau)
    return


def executecif2vasp(args):
    c = ConvertV(args.input, args.output)
    c.cif2vasp()
    return


def executecif2dmol(args):
    c = ConvertV(args.input, args.output)
    c.cif2dmol()
    return

def executecif2qe(args):
    c = ConvertV(args.input, args.output)
    c.cif2qe()
    return

def executevasp2qe(args):
    c = ConvertV(args.input, args.output)
    c.vasp2qe()
    return


def executevasp2cif(args):
    c = ConvertV(args.input, args.output)
    c.vasp2cif()
    return


def executeqe2vasp(args):
    c = ConvertV(args.input, args.output)
    c.qe2vasp()
    return


def executepolyxyz(args):
    c = ConvertV(args.input, args.output)
    c.polytoxyz(args.face, args.length, args.center, args.vertex, args.shapefeeder, args.centroid)
    return


def executeboltzprep(args):
    b = boltztrap.BoltztrapPreparation()
    b.vasp2boltz(args.output, args.input, args.efermi, args.ewindow, args.egrid,
                 args.tmax, args.tstep, args.murange, args.lpfac, True)
    return


def executeboltzpost(args):
    b = boltztrap.BoltztrapPost(args.input)
    # b.name(args.input)
    shift = 0.0
    if args.shift:
        if args.fermi is None:
            shift = b.intransparser()['fermi']
        else:
            shift = args.fermi
    p = plotter.Plotter(args.output)
    p.boltzplot(b.traceparser(), shift)
    return


def executetempband(args):
    b = Bandstructure().temperature_band(args.fermi, args.fake, args.temp, args.emin, args.emax, args.egrid,
                                         args.kgrid, args.trunc, args.smin, args.smax, args.sgrid, args.method)
    p = plotter.Plotter(args.output).temp_pband(b[0], b[1])
    return


def executestruccartdir(args):
    s = Structurehandler(args.input)

    if s.structure.cart is True:
        print("Converting Cartesian to Direct coordinate...")

    else:
        print("Converting Direct to Cartesian coordinate...")

    s.cartdirconvert()
    dic = s.as_dict()

    if args.dyn is False:
        dic['seldyn'] = False

    with open(args.output, 'w') as out:
        out.write(dic['index'] + "\n")
        out.write("1.0\n")

        for x in dic['unitvec']:
            for y in x:
                out.write("{:15.10f}".format(y) + "   ")
            out.write("\n")

        for x in dic['atoms']:
            out.write("{:>4}".format(x) + "   ")
        out.write("\n")

        for x in dic['atoms']:
            out.write("{:4d}".format(dic['atoms'][x]) + "   ")
        out.write("\n")

        if dic['seldyn'] is True:
            out.write("Selective Dynamics\n")

        if dic['cart'] is True:
            out.write("Cartesian\n")
        elif dic['cart'] is False:
            out.write("Direct\n")

        for i in range(len(dic['coord'])):
            for x in dic['coord'][i]:
                out.write("{:15.10f}".format(x) + "   ")
                # out.write(str(x) + " ")

            if dic['seldyn'] is True:
                for x in dic['dyn'][i]:
                    out.write(str(x) + "   ")
            out.write("\n")

        if args.vel is True:
            for x in dic['vel']:
                for y in x:
                    out.write("{:15.10f}".format(y) + "   ")
                out.write("\n")

    print("Done!")
    return


def main():
    description = """ws_py : a python library to handle the VASP input / output files.
Functions : poly, elec, plot, conv, boltz
Type vw.py [function] --help for more detailed informations.
"""

    descpoly = """
-------------------------------------------------------------------------------------------------------
"poly" : Tool for handling the polyhedron-related structures.
Usage example : .py poly -c W -v O [additional options]

Required arguments
  -c [center] : Setting the center atom of the polyhedron. 
  -v [vertex] : Setting the vertex atom of the polyhedron.
Options
  -i [input]  : Name of the input file. (Default : POSCAR)
  -o [output] : Name of the output file. Without this, it will not provide additional output file.
  -s [silent] : Silent option. Automatically turns on the -o option.
  -f [face]   : Setting the number of faces in polyhedron. (Default : 8)
  -l [length] : Setting the bond length limit in building a polyhedron. (Default : 3.0)
  -t          : Toggle on the octahedral tilt angle analysis.
Sub-options for -t
  -ideal [angle] : Setting the ideal angle in calculating tilt angle variance. (Default : 180.0)
  -centroid      : Use centroid in calculating CVC-O angle. Else, the center atom is used.
  -vo            : Exclude the center atom when calculating the centroid of octahedron.
  -edge          : Include edge-sharing octahedra in the tilt angle analysis.
  -tav           : Change the format of angle analysis in the TAV, rather than the average.
-------------------------------------------------------------------------------------------------------"""

    descelec = """
-------------------------------------------------------------------------------------------------------
"elec" : Tool for handling the electronic structure result.
Usage example : .py elec [additional options]

Options
  -fermi [energy] : Manually setting the fermi level.
  -m              : Toggle on the effective mass calculation tool.
  -e              : Toggle on the band edge position tool.
Sub-options for "-m" 
  -k              : Setting the k-point number to calculate carrier effective mass.
  -b              : Setting the band number to calculate carrier effective mass.
  -p              : Setting the number of eigenvalue points to interpolate.
  -r              : Displays the residuals of the least-square parabolic fit.
-------------------------------------------------------------------------------------------------------"""

    descplot = """
-------------------------------------------------------------------------------------------------------
"plot" : Tool for plotting.
Usage example : .py plot [band|dos|pband|temp|diel] [additional options]

Required arguments (But mutually exclusive)
  band  : Plot band structure.
  pband : Plot projected band structure.
  dos   : Plot density-of-states.
  diel  : Plot dielectric responses.

Sub-options for "band"
  -o [output]       : Setting the name of the output file. (Default : band.itx)
  -fermi [float]    : Manually setting the fermi level. (Default: E-fermi from OUTCAR)
  -s                : Toggle on the automatic shifting of VBM to 0
  -g                : Draw the guiding line for high-symmetry points.
  -fake             : Toggle fake-weight method.

Sub-options for "dos"
  -o [output]       : Setting the name of the output file. (Default : *dos.itx)
  -fermi [float]    : Manually setting the fermi level. (Default: E-fermi from OUTCAR)
  -c                : Merges every DOS data into a one single .itx file.

Sub-options for "pband"
  -o [output]       : Setting the name of the output file. (Default : pband.itx)
  -fermi [float]    : Manually setting the fermi level. (Default: E-fermi from OUTCAR)
  -s                : Toggle on the automatic shifting of VBM to 0
  -g                : Draw the guiding line for high-symmetry points.
  -fake             : Toggle fake-weight method.
  -atom [ints]      : Atom-projected bandstructure. 
  -orb              : Orbital-projected bandstructure.
  -spin             : Toggle spin(mz)-projected bandstructure mode in spin-polarized calculation.

Sub-options for "temp"
  -o [output]       : Setting the name of the output file. (Default : tband.itx)
  -fermi [float]    : Manually setting the fermi level. (Default: E-fermi from OUTCAR)
  -t [float]        : Temperature in Kelvin scale. (Default: 300K)
  -emin [float]     : Minimum value of eigenvalues. (Default: -20eV)
  -emax [float]     : Maximum value of eigenvalues. (Default: 15eV)
  -egrid [int]      : Number of eigenvalue grids (y-axis). (Default: 2000)
  -etrunc [float]   : Upper and lower truncation energy. (Default: 1eV)
  -kgrid [int]      : Number of k-grids (x-axis). (Default: 800)
  -smin [float]     : Lower limit of smearing window. (Default: -2eV)
  -smax [float]     : Upper limit of smearing window. (Default: 2eV)
  -sgrid [int]      : Number of smearing grids (local x-axis). (Default: 2000)
  -m [string]       : Type of smearing (G, FD, BE, MB available). (Default: Gaussian)
  -fake             : Toggle fake-weight method.
   
Sub-options for "diel"
  -o [output]       : Setting the name of the output file. (Default : dielec.itx)
  -d                : Toggle on the transversal and longitudinal direction separation.
  -p                : Selecting the variables to plot in Igor itx file. (eI, eR, n, k, alpha, ELS, R, T, A)
  -D                : Toggle on the manual calculation / inclusion of the Drude peaks.
  -P                : Plasma frequency squared value required in calculating the Drude peaks.
  -t                : Tau (relaxation time) term required in calculating the Drude peaks. (Default : 0.1)
-------------------------------------------------------------------------------------------------------
"""

    descconv = """
-------------------------------------------------------------------------------------------------------
"conv" : Tool for converting file formats.
Usage example : .py conv [cif2vasp|vasp2qe|polyxyz] -i [input] -o [output]

Required arguments (But mutually exclusive)
  cif2vasp  : Converts cif file to VASP structure format.
  vasp2qe   : Converts VASP structure file to QE input file format.
  vasp2cif  : Converts VASP structure file to cif file format.
  qe2vasp   : Extract the structure file from QE input file and writes to VASP structure file format. 
  polyxyz   : Extract the octahedron data in VASP and creates xyz files. 

Common options
  -i [input]  : Name of the input file.
  -o [output] : Name of the output file.
  
Sub-options for "polyxyz"
  -i [input]  : Name of the input file.
  -o [output] : Name of the output file.
  -c [center] : Setting the center atom of the polyhedron. 
  -v [vertex] : Setting the vertex atom of the polyhedron.
  -f [face]   : Setting the number of faces in polyhedron. (Default : 8)
  -l [length] : Setting the bond length limit in building a polyhedron. (Default : 3.0)
  -C          : Replaces the center atom position to the centroid of octahedron
  -S          : Toggles the output format in SHAPE program compatible
-------------------------------------------------------------------------------------------------------
"""

    descboltz = """
-------------------------------------------------------------------------------------------------------
"boltz" : Tool for pre / postprocessing files for BoltzTraP.
Usage example : .py boltz [prep|post] [additional options]

Required arguments (But mutually exclusive)
  prep  : Prepare files required to use BoltzTraP.
  post  : Postprocess BoltzTraP outputs.
  
Sub-options for "prep"
  -i [string]       : Name of file to extract structure data (Default : POSCAR)
  -o [string]       : Name of the output file (Default : boltztrap)
  -f [float]        : Fermi energy from scf run (Default : None, automatically reads from OUTCAR)
  -w [float]        : Energy window along Fermi E (Default : 5 eV)
  -g [float]        : Energy grid distance within given window (Default : 0.007 eV)
  -t [integer]      : Maximum temperature to consider (Default : 800 K)
  -step [integer]   : Temperature gradient (Default : 50 K)
  -mu [float]       : Chemical potential window (Default : 2 eV)
  -lpfac [integer]  : Number of lattice points per k-point (Default : 5)

Sub-options for "post"
  -o [string]       : Name of the output file (Default : boltzplot.itx)
  -s                : Toggles on the shifting Fermi level to zero

-------------------------------------------------------------------------------------------------------
    """

    descstruc = """
-------------------------------------------------------------------------------------------------------
"struc" : Tool to handle structure files.
Usage example : .py struc [additional options]

Required arguments (But mutually exclusive)
  trans    :  Translates atomic positions with given translation matrix (1x3)
  cartdir  :  Converts Cartesian coordinate system to Direct coordinate system, and vice versa.

Sub-options for "trans"
  -i [string]       : Name of input structure file (Default : POSCAR)
  -o [string]       : Name of output structure file (Default : POSCAR_conv)
  -t [list]         : Translation matrix (e.g. 0 0 1)
  -c                : Turn on if the given translation matrix is in the Cartesian coordinates

Sub-options for "cartdir"
  -i [string]       : Name of input structure file (Default : POSCAR)
  -o [string]       : Name of output structure file (Default : POSCAR_conv)
  -s                : Removes "Selective Dynamics" part from the file
  -v                : Includes the velocity part at the bottom of the file
-------------------------------------------------------------------------------------------------------
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Subparsers
    subparsers = parser.add_subparsers(title="Functions")

    # Polyhedron parsers
    parser_poly = subparsers.add_parser("poly", formatter_class=argparse.RawTextHelpFormatter, description=descpoly)
    parser_poly.add_argument("-i", dest="input", type=str, default="POSCAR")
    parser_poly.add_argument("-o", dest="output", type=str, default=None)
    parser_poly.add_argument("-s", dest="silent", action='store_true')
    parser_poly.add_argument("-c", dest="center", type=str, required=True, nargs='*')
    parser_poly.add_argument("-v", dest="vertex", type=str, required=True, nargs='*')
    parser_poly.add_argument("-f", dest="face", type=int, default="8")
    parser_poly.add_argument("-l", dest="length", type=float, default="3.0")
    parser_poly.add_argument("-t", dest="tilt", action='store_true')
    parser_poly.add_argument("-centroid", dest="centroid", action='store_true')
    parser_poly.add_argument("-vo", dest="vertexonly", action='store_true')
    parser_poly.add_argument("-edge", dest="edge", action='store_true')
    parser_poly.add_argument("-tav", dest="tav", action='store_true')
    parser_poly.add_argument("-ideal", dest="ideal", type=float, default="180.0")

    parser_poly.set_defaults(func=executepoly)

    # Electronic structure analysis parsers
    parser_elec = subparsers.add_parser("elec", formatter_class=argparse.RawTextHelpFormatter, description=descelec)

    # parser_elec.add_argument("-i", dest="input", type=str, default="None", help="Name of the input file.")
    parser_elec.add_argument("-fermi", dest="fermi", type=float, default="0.0")
    parser_elec.add_argument("-m", dest="mass", action='store_true')
    parser_elec.add_argument("-e", dest="edge", action='store_true')
    parser_elec.add_argument("-b", dest="band", type=int, default=0)
    parser_elec.add_argument("-k", dest="kpoint", type=int, default=0)
    parser_elec.add_argument("-p", dest="points", type=int, default=5)
    parser_elec.add_argument("-r", dest="residuals", action='store_true')
    parser_elec.set_defaults(func=executeelec)

    # Plotting tool parsers
    parser_plot = subparsers.add_parser("plot", formatter_class=argparse.RawTextHelpFormatter, description=descplot)

    plotsubparsers = parser_plot.add_subparsers()
    parser_plot_band = plotsubparsers.add_parser("band")
    parser_plot_band.add_argument("-o", dest="output", type=str, default=None)
    parser_plot_band.add_argument("-fermi", dest="fermi", type=float, default="0.0")
    parser_plot_band.add_argument("-s", dest="shift", action='store_true')
    parser_plot_band.add_argument("-fake", dest="fake", action='store_true')
    parser_plot_band.add_argument("-g", dest="guide", action='store_true')
    parser_plot_band.add_argument("-k", dest="ktraj", action='store_true')
    parser_plot_band.add_argument("-i", dest="prefix", type=str, default=None)
    parser_plot_band.set_defaults(func=executeplotband)

    parser_plot_pband = plotsubparsers.add_parser("pband")
    parser_plot_pband.add_argument("-o", dest="output", type=str, default=None)
    parser_plot_pband.add_argument("-fermi", dest="fermi", type=float, default="0.0")
    parser_plot_pband.add_argument("-s", dest="shift", action='store_true')
    parser_plot_pband.add_argument("-fake", dest="fake", action='store_true')
    parser_plot_pband.add_argument("-g", dest="guide", action='store_true')
    # parser_plot_pband.add_argument("-c", dest="combine", action='store_true')
    parser_plot_pband.add_argument("-spin", dest="spin", action='store_true')
    # parser_plot_pband.add_argument("-atom", dest="atom", action='store_true')
    parser_plot_pband.add_argument("-atom", dest="atom", type=str, default=None, nargs='*')
    parser_plot_pband.add_argument("-orb", dest="orbital", action='store_true')
    parser_plot_pband.add_argument("-i", dest="prefix", type=str, default=None)
    parser_plot_pband.set_defaults(func=executeplotpband)

    parser_plot_dos = plotsubparsers.add_parser("dos")
    parser_plot_dos.add_argument("-o", dest="output", type=str, default=None)
    parser_plot_dos.add_argument("-fermi", dest="fermi", type=float, default="0.0")
    parser_plot_dos.add_argument("-c", dest="combine", action='store_true')
    parser_plot_dos.add_argument("-i", dest="prefix", type=str, default=None)
    parser_plot_dos.set_defaults(func=executeplotdos)

    parser_plot_diel = plotsubparsers.add_parser("diel")
    parser_plot_diel.add_argument("-o", dest="output", type=str, default=None)
    parser_plot_diel.add_argument("-d", dest="direction", action='store_true')
    parser_plot_diel.add_argument("-p", dest="plot", type=str, default=None, nargs='*')
    parser_plot_diel.add_argument("-D", dest="drude", action='store_true')
    parser_plot_diel.add_argument("-P", dest="plasmasq", type=float, default=None)
    parser_plot_diel.add_argument("-t", dest="tau", type=float, default=0.1)
    parser_plot_diel.add_argument("-i", dest="prefix", type=str, default=None)
    parser_plot_diel.set_defaults(func=executeplotdielectric)

    parser_plot_temp = plotsubparsers.add_parser("temp")
    parser_plot_temp.add_argument("-o", dest="output", type=str, default=None)
    parser_plot_temp.add_argument("-fermi", dest="fermi", type=float, default="0.0")
    parser_plot_temp.add_argument("-t", dest="temp", type=float, default="300.0")
    parser_plot_temp.add_argument("-emin", dest="emin", type=float, default="-20.0")
    parser_plot_temp.add_argument("-emax", dest="emax", type=float, default="15.0")
    parser_plot_temp.add_argument("-egrid", dest="egrid", type=int, default="2000")
    parser_plot_temp.add_argument("-kgrid", dest="kgrid", type=int, default="800")
    parser_plot_temp.add_argument("-etrunc", dest="trunc", type=float, default="1")
    parser_plot_temp.add_argument("-smin", dest="smin", type=float, default="-2.0")
    parser_plot_temp.add_argument("-smax", dest="smax", type=float, default="2.0")
    parser_plot_temp.add_argument("-sgrid", dest="sgrid", type=int, default="2000")
    parser_plot_temp.add_argument("-m", dest="method", type=str, default='FD')
    parser_plot_temp.add_argument("-fake", dest="fake", action='store_true')
    parser_plot_temp.add_argument("-g", dest="guide", action='store_true')
    parser_plot_temp.set_defaults(func=executetempband)

    # Converting tool parsers
    parser_conv = subparsers.add_parser("conv", formatter_class=argparse.RawTextHelpFormatter, description=descconv)
    convsubparsers = parser_conv.add_subparsers()

    parser_conv_cif2vasp = convsubparsers.add_parser("cif2vasp")
    parser_conv_cif2vasp.add_argument("-i", dest="input", type=str, required=True)
    parser_conv_cif2vasp.add_argument("-o", dest="output", type=str, required=True)
    parser_conv_cif2vasp.set_defaults(func=executecif2vasp)

    parser_conv_cif2dmol = convsubparsers.add_parser("cif2dmol")
    parser_conv_cif2dmol.add_argument("-i", dest="input", type=str, required=True)
    parser_conv_cif2dmol.add_argument("-o", dest="output", type=str, required=True)
    parser_conv_cif2dmol.set_defaults(func=executecif2dmol)

    parser_conv_vasp2qe = convsubparsers.add_parser("vasp2qe")
    parser_conv_vasp2qe.add_argument("-i", dest="input", type=str, required=True)
    parser_conv_vasp2qe.add_argument("-o", dest="output", type=str, required=True)
    parser_conv_vasp2qe.set_defaults(func=executevasp2qe)

    parser_conv_qe2vasp = convsubparsers.add_parser("qe2vasp")
    parser_conv_qe2vasp.add_argument("-i", dest="input", type=str, required=True)
    parser_conv_qe2vasp.add_argument("-o", dest="output", type=str, required=True)
    parser_conv_qe2vasp.set_defaults(func=executeqe2vasp)

    parser_conv_vasp2cif = convsubparsers.add_parser("vasp2cif")
    parser_conv_vasp2cif.add_argument("-i", dest="input", type=str, required=True)
    parser_conv_vasp2cif.add_argument("-o", dest="output", type=str, required=True)
    parser_conv_vasp2cif.set_defaults(func=executevasp2cif)

    parser_conv_cif2qe = convsubparsers.add_parser("cif2qe")
    parser_conv_cif2qe.add_argument("-i", dest="input", type=str, required=True)
    parser_conv_cif2qe.add_argument("-o", dest="output", type=str, required=True)
    parser_conv_cif2qe.set_defaults(func=executecif2qe)

    parser_conv_polyxyz = convsubparsers.add_parser("polyxyz")
    parser_conv_polyxyz.add_argument("-i", dest="input", type=str, required=True)
    parser_conv_polyxyz.add_argument("-o", dest="output", type=str)
    parser_conv_polyxyz.add_argument("-c", dest="center", type=str, required=True, nargs='*')
    parser_conv_polyxyz.add_argument("-v", dest="vertex", type=str, required=True, nargs='*')
    parser_conv_polyxyz.add_argument("-f", dest="face", type=int, default="8")
    parser_conv_polyxyz.add_argument("-l", dest="length", type=float, default="3.0")
    parser_conv_polyxyz.add_argument("-C", dest="centroid", action='store_true')
    parser_conv_polyxyz.add_argument("-S", dest="shapefeeder", action='store_true')
    parser_conv_polyxyz.set_defaults(func=executepolyxyz)

    # Boltztrap related functions
    parser_boltz = subparsers.add_parser("boltz", formatter_class=argparse.RawTextHelpFormatter, description=descboltz)
    boltzsubparsers = parser_boltz.add_subparsers()

    parser_boltz_prep = boltzsubparsers.add_parser("prep")
    parser_boltz_prep.add_argument("-i", dest="input", type=str, default='POSCAR')
    parser_boltz_prep.add_argument("-o", dest="output", type=str, default='boltztrap')
    parser_boltz_prep.add_argument("-f", dest="efermi", type=float, default=None)
    parser_boltz_prep.add_argument("-w", dest="ewindow", type=float, default=5)
    parser_boltz_prep.add_argument("-g", dest="egrid", type=float, default=0.007)
    parser_boltz_prep.add_argument("-t", dest="tmax", type=int, default=800)
    parser_boltz_prep.add_argument("-step", dest="tstep", type=int, default=50)
    parser_boltz_prep.add_argument("-mu", dest="murange", type=float, default=2)
    parser_boltz_prep.add_argument("-lpfac", dest="lpfac", type=int, default=5)
    # parser_boltz_prep.add_argument("-g", dest="generic", action='store_false')
    parser_boltz_prep.set_defaults(func=executeboltzprep)

    parser_boltz_post = boltzsubparsers.add_parser("post")
    parser_boltz_post.add_argument("-i", dest="input", type=str, default='boltztrap')
    parser_boltz_post.add_argument("-o", dest="output", type=str, default=None)
    parser_boltz_post.add_argument("-s", dest="shift", action='store_true')
    parser_boltz_post.add_argument("-f", dest="fermi", type=float, default=0.0)
    parser_boltz_post.set_defaults(func=executeboltzpost)

    parser_struc = subparsers.add_parser("struc", formatter_class=argparse.RawTextHelpFormatter, description=descstruc)
    strucsubparsers = parser_struc.add_subparsers()

    parser_struc_cartdir = strucsubparsers.add_parser("cartdir")
    parser_struc_cartdir.add_argument("-i", dest="input", type=str, default='POSCAR')
    parser_struc_cartdir.add_argument("-o", dest="output", type=str, default='POSCAR_conv')
    parser_struc_cartdir.add_argument("-s", dest="dyn", action='store_false')
    parser_struc_cartdir.add_argument("-v", dest="vel", action='store_true')
    parser_struc_cartdir.set_defaults(func=executestruccartdir)

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()
