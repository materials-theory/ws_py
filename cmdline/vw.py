#!/usr/local/bin/python

import argparse
import sys

from parsers.vasprun import Vasprun
from handlers.dos import Dos
from handlers.bandstructure import Bandstructure
from plotters.plotter import *


def executeplotdos(args):
    if args.xml:
        dic = Vasprun(args.input).to_dic(band=False, pband=False, dos=True)
    else:
        return

    dos = Dos(tdos=dic["dos"]["tdos"], pdos=dic["dos"]["pdos"], structure=dic["structure"],
              efermi=dic["energy"]["e_fermi"], others=args.others, atom=args.atom, orbital=args.orb)

    plotter = DosPlot(dos.get_plot_dict(dos.tdos, dos.projected), args.output, args.prefix, args.program, args.split)
    return plotter.plot()


def executeplotband(args):
    if args.xml:
        dic = Vasprun(args.input).to_dic(band=True, pband=False, dos=False)
    else:
        return

    bs = Bandstructure(eigenvalues=dic["band"], ktrace=dic["kpoints"]["list"],
                      structure=dic["structure"], efermi=dic["energy"]["e_fermi"], shift=args.shift,
                      ymin=args.ymin, ymax=args.ymax)

    plotter = BandPlot(bs.get_plot_dict(bs.eigenvalues, bs.projected), args.output, args.prefix, args.program)
    return plotter.plot()


def executeplotpband(args):
    if args.xml:
        dic = Vasprun(args.input).to_dic(band=True, pband=True, dos=False)
    else:
        return

    bs = Bandstructure(eigenvalues=dic["band"], projection=dic["pband"], ktrace=dic["kpoints"]["list"],
                      structure=dic["structure"], efermi=dic["energy"]["e_fermi"], shift=args.shift,
                      atom=args.atom, orbital=args.orb, ymin=args.ymin, ymax=args.ymax)

    plotter = BandPlot(bs.get_plot_dict(bs.eigenvalues, bs.projected), args.output, args.prefix, args.program)
    return plotter.plot()


def main():
    description = ""
    descplot = ""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(title="Functions")

    parser_plot = subparsers.add_parser("plot", formatter_class=argparse.RawTextHelpFormatter, description=descplot)
    plotsubparsers = parser_plot.add_subparsers()

    parser_plot_dos = plotsubparsers.add_parser("dos")
    parser_plot_dos.add_argument("-x", dest="xml", action="store_true")
    parser_plot_dos.add_argument("-o", dest="output", type=str, default="dos.itx")
    parser_plot_dos.add_argument("-f", dest="fermi", type=float, default="0.0")
    parser_plot_dos.add_argument("-s", dest="shift", action="store_true")
    parser_plot_dos.add_argument("-O", dest="others", action="store_true")
    parser_plot_dos.add_argument("-atom", dest="atom", type=str, default=None, nargs='*')
    parser_plot_dos.add_argument("-orb", dest="orb", type=str, default=None, nargs='*')
    parser_plot_dos.add_argument("-spin", dest="spin", action="store_true")
    parser_plot_dos.add_argument("-p", dest="prefix", type=str, default=None)
    parser_plot_dos.add_argument("-prog", dest="program", type=str, default=None)
    parser_plot_dos.add_argument("-S", dest="split", action="store_true")
    parser_plot_dos.set_defaults(func=executeplotdos)

    parser_plot_band = plotsubparsers.add_parser("band")
    parser_plot_band.add_argument("-x", dest="xml", action="store_true")
    parser_plot_band.add_argument("-o", dest="output", type=str, default="band.itx")
    parser_plot_band.add_argument("-f", dest="fermi", type=float, default="0.0")
    parser_plot_band.add_argument("-ymin", dest="ymin", type=float, default=None)
    parser_plot_band.add_argument("-ymax", dest="ymax", type=float, default=None)
    parser_plot_band.add_argument("-s", dest="shift", type=str, default="fermi")
    parser_plot_band.add_argument("-fake", dest="fake", action="store_true")
    parser_plot_band.add_argument("-g", dest="guide", action="store_true")
    parser_plot_band.add_argument("-p", dest="prefix", type=str, default=None)
    parser_plot_band.add_argument("-prog", dest="program", type=str, default=None)
    parser_plot_band.set_defaults(func=executeplotband)

    parser_plot_pband = plotsubparsers.add_parser("pband")
    parser_plot_pband.add_argument("-x", dest="xml", action="store_true")
    parser_plot_pband.add_argument("-o", dest="output", type=str, default="pband.itx")
    parser_plot_pband.add_argument("-f", dest="fermi", type=float, default="0.0")
    parser_plot_pband.add_argument("-ymin", dest="ymin", type=float, default=None)
    parser_plot_pband.add_argument("-ymax", dest="ymax", type=float, default=None)
    parser_plot_pband.add_argument("-s", dest="shift", type=str, default="fermi")
    parser_plot_pband.add_argument("-fake", dest="fake", action="store_true")
    parser_plot_pband.add_argument("-g", dest="guide", action="store_true")
    parser_plot_pband.add_argument("-atom", dest="atom", type=str, default=None, nargs='*')
    parser_plot_pband.add_argument("-orb", dest="orb", type=str, default=None, nargs='*')
    parser_plot_pband.add_argument("-spin", dest="spin", action="store_true")
    parser_plot_pband.add_argument("-p", dest="prefix", type=str, default=None)
    parser_plot_pband.add_argument("-prog", dest="program", type=str, default=None)
    parser_plot_pband.set_defaults(func=executeplotpband)

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()
