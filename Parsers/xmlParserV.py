
# import os
# import itertools
# import math
# import numpy as np
from xml.etree import ElementTree as ET
# from core import IO


#TODO: Change the way calling parser -- do not call many times
class xmlParserV:
    def __init__(self):
        # self.IO = IO(filename)
        # self.generator = None
        # self.incar = None
        # self.kpoints = None
        # self.parameters = None
        # self.atominfo = None
        # self.initstruc = None
        # self.finstruc = None
        self.dielectric = None
        self.other_dielectric = None
        self.lfe = None
        # self.Parser(filename)

        self.parser()
        return

    def parser(self, filename='vasprun.xml'):
        print("Parsing vasprun.xml...")
        for event, elem in ET.iterparse(filename):
            tag = elem.tag
            # if not parsed_header:
            #     if tag == "generator":
            #         self.generator = self._parse_params(elem)
            #     elif tag == "incar":
            #         self.incar = self._parse_params(elem)
            #     elif tag == "kpoints":
            #         self.kpoints, self.actual_kpoints, \
            #         self.actual_kpoints_weights = self._parse_kpoints(
            #             elem)
            #     elif tag == "parameters":
            #         self.parameters = self._parse_params(elem)
            #     elif tag == "structure" and elem.attrib.get("name") == \
            #             "initialpos":
            #         self.initial_structure = self._parse_structure(elem)
            #     elif tag == "atominfo":
            #         self.atomic_symbols, self.potcar_symbols = \
            #             self._parse_atominfo(elem)
            #         self.potcar_spec = [{"titel": p,
            #                              "hash": None} for
            #                             p in self.potcar_symbols]

            # if tag == "calculation":
            #     parsed_header = True
            #     if not self.parameters.get("LCHIMAG", False):
            #         ionic_steps.append(self._parse_calculation(elem))
            #     else:
            #         ionic_steps.extend(self._parse_chemical_shift_calculation(elem))
            # elif parse_dos and tag == "dos":
            #     try:
            #         self.tdos, self.idos, self.pdos = self._parse_dos(elem)
            #         self.efermi = self.tdos.efermi
            #         self.dos_has_errors = False
            #     except Exception as ex:
            #         self.dos_has_errors = True
            # elif parse_eigen and tag == "eigenvalues":
            #     self.eigenvalues = self._parse_eigen(elem)
            # elif parse_projected_eigen and tag == "projected":
            #     self.projected_eigenvalues = self._parse_projected_eigen(
            #         elem)
            if tag == "dielectricfunction":
                if "comment" not in elem.attrib:
                    self.dielectric = self._parse_diel(elem)
                    self.lfe = False

                elif elem.attrib["comment"] == \
                        "HEAD OF MICROSCOPIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE)":
                    self.other_dielectric = self._parse_diel(elem)

                elif elem.attrib["comment"] == \
                        "1 + v P,  with REDUCIBLE POLARIZABILTY P=P_0 (1 -(v+f) P_0)^-1":
                    pass

                elif elem.attrib["comment"] == \
                        "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))":
                    self.dielectric = self._parse_diel(elem)
                    self.lfe = True

                elif elem.attrib["comment"] == \
                         "screened Coulomb potential":
                    pass


            # if ("comment" not in elem.attrib) or \
            #                 elem.attrib["comment"] == \
            #                 "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))":
            #     self.dielectric = self._parse_diel(elem)
            # else:
            #     self.other_dielectric[elem.attrib["comment"]] = self._parse_diel(elem)


            # elif tag == "structure" and elem.attrib.get("name") == \
            #         "finalpos":
            #     self.final_structure = self._parse_structure(elem)
            # elif tag == "dynmat":
            #     hessian, eigenvalues, eigenvectors = self._parse_dynmat(elem)
            #     natoms = len(self.atomic_symbols)
            #     hessian = np.array(hessian)
            #     self.force_constants = np.zeros((natoms, natoms, 3, 3), dtype='double')
            #     for i in range(natoms):
            #         for j in range(natoms):
            #             self.force_constants[i, j] = hessian[i * 3:(i + 1) * 3, j * 3:(j + 1) * 3]
            #     phonon_eigenvectors = []
            #     for ev in eigenvectors:
            #         phonon_eigenvectors.append(np.array(ev).reshape(natoms, 3))
            #     self.normalmode_eigenvals = np.array(eigenvalues)
            #     self.normalmode_eigenvecs = np.array(phonon_eigenvectors)

        print("Done!")

        return

    @staticmethod
    def _parse_diel(elem):
        imag = []
        real = []

        for x in elem.find("imag").find("array").find("set").findall("r"):
            tmp = []
            for y in x.text.split():
                tmp.append(float(y))
            imag.append(tmp)

        for x in elem.find("real").find("array").find("set").findall("r"):
            tmp = []
            for y in x.text.split():
                tmp.append(float(y))
            real.append(tmp)

        elem.clear()

        # return np.array(imag, )
        return [imag, real]

    # def SplitLine(self, line):
    #     elements = []
    #
    #     for element in line.strip().split():
    #         strippedElement = element.strip()
    #
    #         if strippedElement != "":
    #             elements.append(strippedElement)
    #
    #     return elements

# xmlParserV()
