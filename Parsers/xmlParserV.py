from xml.etree import ElementTree as ET


#TODO: Change the way calling parser -- do not call many times
#TODO: VASP 5.4.4 density-density and current-current dielectric function parsing problem
#Currently adopted from "Vasprun" class from Pymatgen (http://pymatgen.org)
#Need to be updated

class xmlParserV:
    def __init__(self):
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
