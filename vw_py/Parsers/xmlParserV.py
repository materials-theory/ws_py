from xml.etree import ElementTree as etree

# TODO: Change the way calling parsers -- do not call many times

# Modified version of "Vasprun" class in Pymatgen (http://pymatgen.org)
# Need to be updated


class xmlParserV:
    def __init__(self):
        self.version = None
        self.dielectric = {}
        self.other_dielectric = {}
        self.lfe = {}
        # self.Parser(filename)

        self.parser()
        return

    def parser(self, filename='vasprun.xml'):
        print("Parsing vasprun.xml...")
        for event, elem in etree.iterparse(filename):
            tag = elem.tag
            if tag == "generator":
                for x in elem.findall("i"):
                    if x.attrib["name"] == "version":
                        self.version = x.text.split()[0][0:5]

            if tag == "dielectricfunction":
                if "comment" not in elem.attrib:
                    if self.version >= "5.4.4":
                        if "density" not in self.dielectric:
                            self.dielectric["density"] = self._parse_diel(elem)
                        else:
                            self.dielectric["current"] = self._parse_diel(elem)
                    else:
                        self.dielectric["normal"] = self._parse_diel(elem)
                    self.lfe = False

                elif elem.attrib["comment"] == \
                        "HEAD OF MICROSCOPIC DIELECTRIC TENSOR (INDEPENDENT PARTICLE)":
                    if self.version >= "5.4.4":
                        if "density" not in self.dielectric:
                            self.other_dielectric["density"] = self._parse_diel(elem)
                        else:
                            self.other_dielectric["current"] = self._parse_diel(elem)
                    else:
                        self.other_dielectric["normal"] = self._parse_diel(elem)

                elif elem.attrib["comment"] == \
                        "1 + v P,  with REDUCIBLE POLARIZABILTY P=P_0 (1 -(v+f) P_0)^-1":
                    pass

                elif elem.attrib["comment"] == \
                        "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))":
                    if self.version >= "5.4.4":
                        if "density" not in self.dielectric:
                            self.dielectric["density"] = self._parse_diel(elem)
                        else:
                            self.dielectric["current"] = self._parse_diel(elem)
                    else:
                        self.dielectric["normal"] = self._parse_diel(elem)
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
