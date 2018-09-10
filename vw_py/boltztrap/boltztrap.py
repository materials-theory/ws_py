import math
import spglib
import numpy as np

from collections import defaultdict
from ase import io

from vw_py.Parsers.outcarhandler import OutcarHandler
from vw_py.ElectronicStructure import bandstructure as bs
from vw_py.generalutils.constants import Unitconverter as Units


class BoltztrapPreparation(object):
    def __init__(self):
        return

    def vasp2boltz(self, filename='boltztrap',
                         structure='POSCAR',
                         efermi=None,
                         ewindow=5,
                         egrid=0.007,
                         tmax=800,
                         tstep=50,
                         murange=2,
                         lpfac=5,
                         boltz_generic=True):

        if efermi is None:
            efermi = float(OutcarHandler.param_from_outcar('E-fermi'))
        else:
            pass

        ev_to_ry = Units.conversion_factor('eV', 'Ry')
        fermi = efermi * ev_to_ry
        window = ewindow * ev_to_ry
        grid = egrid * ev_to_ry
        mu = murange * ev_to_ry

        self.write_intrans_boltztrap(filename, fermi, window, grid, tmax, tstep, mu, lpfac, boltz_generic)
        self.write_struct_boltztrap(structure, filename)
        self.write_energy_boltztrap(None, filename)
        return

    @staticmethod
    def write_intrans_boltztrap(filename, efermi, ewindow, egrid, tmax, tstep, murange, lpfac, boltz_generic):

        outfile = str(filename + ".intrans")
        n_elec = float(OutcarHandler.param_from_outcar('NELECT'))

        with open(outfile, 'w') as out:
            if boltz_generic:
                out.write("GENE          # use generic interface\n")
            else:
                out.write("WIEN          # use wien interface\n")
            out.write("0 0 0 0.0         # iskip (not presently used) idebug setgap shiftgap \n")
            out.write("%7.5f %7.5f %7.5f %i # E-fermi (Ry), energy grid, energy window, NELECT\n"
                      % (efermi, egrid, ewindow, int(n_elec)))
            out.write("CALC                    # CALC (calculate expansion coeff), NOCALC read from file\n")
            out.write("%i                      # lpfac, number of latt-points per k-point\n" % lpfac)
            out.write("BOLTZ                   # run mode (only BOLTZ is supported)\n")
            out.write("%7.5f                   # (efcut) energy range of chemical potential\n" % murange)
            out.write("%i %i                   # Tmax, temperature grid\n" % (tmax, tstep))
            out.write("-1.                     # energy range of bands in sig_xxx and dos_xxx (xxx is band number)\n")
            out.write("HISTO\n")

        return

    @staticmethod
    def write_struct_boltztrap(structure, filename):
        ao = io.read(structure)
        spacegroup = spglib.get_spacegroup(ao, symprec=1e-5)
        ao.info = {'spacegroup': spacegroup}

        filename = str(filename + ".struct")

        with open(filename, "w") as out:
            out.write('HTE output' + '\n')  # title
            latt = ao.get_cell() / 0.5291772083

            for i in range(3):
                line = ''
                for j in range(3):
                    line = line + "%12.5f" % latt[i][j]
                out.write(line + '\n')

            krot = get_kspace_operations(ao)
            out.write(str(len(krot)) + '\n')

            for iop in range(len(krot)):
                for i in range(3):
                    for j in range(3):
                        out.write(str(krot[iop][i][j]) + ' ')
                    out.write('\n')

        return

    @staticmethod
    def write_energy_boltztrap(bandstructure, filename):
        ev_to_ry = Units.conversion_factor('eV', 'Ry')
        filename = str(filename + ".energy")

        if bandstructure is None:
            banddata = bs.Bandstructure().rawband()
        else:
            banddata = bandstructure

        path = np.round(banddata[0], 8).tolist()
        eigen = np.round(banddata[1], 8).tolist()
        numband = int(banddata[3])

        for i in range(len(path)):
            for j in range(3):
                path[i][j] = abs(path[i][j])


        with open('test.txt', 'w') as p:
            for x in path:
                for y in x:
                    p.write(str(y) + " ")
                p.write("\n")

        iden_points = []
        numpoint = []

        for num, p in enumerate(path):
            if p not in iden_points:
                iden_points.append(p)
                numpoint.append(num)

        with open(filename, 'w') as out:
            out.write('HTE output\n')
            out.write('%i\n' % int(len(iden_points)))
            for num, line in enumerate(iden_points):
                out.write('%12.8f %12.8f %12.8f %d\n' % (line[0], line[1], line[2], numband))
                for band in eigen[numpoint[num]]:
                    out.write('%12.8f\n' % (band * ev_to_ry))

        return


#TODO: Complete tensor parsers
class BoltztrapPost(object):
    def __init__(self, filename='boltztrap'):
        self.filename = filename
        return

    # def name(self, filename='boltztrap'):
    #     self.filename = filename
    #     return

    def intransparser(self):
        intransname = str(self.filename + ".intrans")
        intrans = {}
        ry_to_ev = Units.conversion_factor('Ry', 'eV')

        with open(intransname, 'r') as inp:
            lines = inp.readlines()
            intrans['fermi'] = float(lines[2].split()[0]) * ry_to_ev
            intrans['grid'] = float(lines[2].split()[1]) * ry_to_ev
            intrans['window'] = float(lines[2].split()[2]) * ry_to_ev
            intrans['nelect'] = int(lines[2].split()[3])
            intrans['lpfac'] = int(lines[4].split()[0])
            intrans['murange'] = float(lines[6].split()[0]) * ry_to_ev
            intrans['tmax'] = int(lines[7].split()[0])
            intrans['tgrid'] = int(lines[7].split()[1])

        return intrans

    def structparser(self):
        structname = str(self.filename + ".struct")
        struct = {}

        with open(structname, 'r') as inp:
            inp.readline()
            struct['x'] = np.array(inp.readline().split())
            struct['y'] = np.array(inp.readline().split())
            struct['z'] = np.array(inp.readline().split())

        return struct

    def outputtransparser(self, filename='boltztrap'):
        return

    def condtensparser(self, filename='boltztrap'):
        return

    def halltensparser(self, filename='boltztrap'):
        return

    def transdosparser(self, filename='boltztrap'):
        return

    def traceparser(self):
        tracename = str(self.filename + ".trace")
        parsed = defaultdict(dict)
        ry_to_ev = Units.conversion_factor('Ry', 'eV')

        with open(tracename, 'r') as trace:
            trace.readline() # deleting header line
            lines = trace.readlines()
            for i in range(len(lines)):
                insidedict = {'N': float(lines[i].split()[2]),
                              'DOS': float(lines[i].split()[3]),
                              'S': float(lines[i].split()[4]),
                              'st': float(lines[i].split()[5]),
                              'R': float(lines[i].split()[6]),
                              'k': float(lines[i].split()[7]),
                              'c': float(lines[i].split()[8]),
                              'chi': float(lines[i].split()[9])}

                parsed[int(float(lines[i].split()[1]))][float(lines[i].split()[0]) * ry_to_ev] = insidedict

        return parsed

def get_kspace_operations(ao, symprec_spglib=1e-5):
    # returns k-space operations for the atoms object ao
    kops = None
    rot = spglib.get_symmetry_dataset(ao, symprec=symprec_spglib)['rotations']
    kops = []
    for iop in range(len(rot)):
        mat = rot[iop]
        newop = True
        for i in range(len(kops)):
            if cmp_mat(kops[i], mat):
                newop = False
                break
        if newop:
            kops.append(mat)
    return kops


def cmp_mat(mat1, mat2):
    absdiff = abs(mat1 - mat2)
    value = sum(sum(absdiff))
    if value > 1.0e-10:
        return False
    else:
        return True