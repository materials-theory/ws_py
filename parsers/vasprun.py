import numpy as np
import xml.etree.cElementTree as et

from handlers.structure import AtomicStructure

# TODO: make all dicts to list for memory management -> maybe transfer back to dicts using to_dict


class Vasprun(object):

    def __init__(self, infile, parse_eig=True, parse_dos=True, parse_pband=True):
        self.infile = infile
        self.generator = None
        self.incar = None
        self.kpoints = None
        self.parameters = None
        self.atominfo = None
        self.init_structure = None
        self.fin_structure = None
        self.energy = None
        self.time = {}
        self.calculation = {}

        self._parse_xml(self.infile, parse_eig, parse_dos, parse_pband)

    def _parse_xml(self, infile, parse_eig, parse_dos, parse_pband):
        for event, elem in et.iterparse(infile):
            if elem.tag == "generator":
                self.generator = self._parse_elem(elem)
            elif elem.tag == "incar":
                self.incar = self._parse_elem(elem)
            elif elem.tag == "kpoints":
                self.kpoints = self._parse_elem(elem)
            elif elem.tag == "parameters":
                self.parameters = self._parse_elem(elem)
            elif elem.tag == "atominfo":
                self.atominfo = self._parse_elem(elem)
                simple = []
                for x in self.atominfo["atomtypes"]["type"]:
                    simple.append([x[1][0], x[0][0]])
                self.atominfo["simple"] = simple
            elif elem.tag == "structure":
                if elem.get("name") == "initialpos":
                    self.init_structure = self._parse_structure(elem)
                elif elem.get("name") == "finalpos":
                    self.fin_structure = self._parse_structure(elem)

            elif elem.tag == "calculation":
                for elem2 in elem:
                    if elem2.tag == "scstep":
                        pass
                    elif elem2.tag == "structure":
                        pass
                    elif elem2.tag == "energy":
                        self.energy = self._parse_elem(elem2)
                    elif elem2.tag == "time":
                        timetype = elem2.get("name")
                        self.time[timetype] = list(map(float, elem2.text.split()))
                    elif elem2.tag == "eigenvalues" and parse_eig is True:
                        self.calculation["eigenvalues"] = self._parse_elem(elem2)["array"]
                    elif elem2.tag == "separator":
                        name = elem2.get("name")
                        self.calculation[name] = self._parse_elem(elem2)
                    elif elem2.tag == "dos":
                        self.energy["e_fermi"] = float(elem2.find("i").text)
                        if parse_dos is True:
                            tmp = self._parse_elem(elem2)
                            self.calculation["total_dos"] = tmp["total"]["array"]
                            self.calculation["partial_dos"] = tmp["partial"]["array"]
                            # self._to_np_array(self.calculation["partial_dos"]["gridpoints"])
                            tmp = None
                    elif elem2.tag == "projected" and parse_pband is True:
                        self.calculation["pband"] = self._parse_elem(elem2)["array"]
                        # self._to_np_array(self.calculation["pband"]["ion"])

        return

    @staticmethod
    def _parse_keys(tag, keytype, key):
        if tag in ["v", "c"]:
            if keytype == "string":
                return list(map(str, key.split()))

            elif keytype == "logical":
                return list(map(bool, key.split()))

            elif keytype == "int":
                return list(map(int, key.split()))

            else:
                if tag in ["c"]:
                    return list(map(str, key.split()))
                else:
                    return list(map(float, key.split()))

        elif tag in ["r"]:
            return list(map(float, key.split()))

        elif tag in ["i"]:
            if keytype == "string":
                if key is None:
                    return None
                else:
                    return key.strip(" ")

            elif keytype == "logical":
                if "T" in key or ".t" in key:
                    return True
                elif "F" in key or ".f" in key:
                    return False
                else:
                    raise TypeError("Unknown logical key value...")

            elif keytype == "int":
                return int(key)

            else:
                return float(key)

        else:
            raise ValueError("Unknown input tag...")

    def _parse_elem(self, elem, recursive_root=None):
        parsed = {}
        if recursive_root is None:
            for x in elem:
                if x.tag not in ["i", "v", "c", "rc", "r"]:
                    if x.get("param") is not None:
                        root = x.get("param")
                    elif x.get("comment") is not None:
                        root = x.get("comment")
                    elif x.get("name") is not None:
                        root = x.get("name")
                    elif len(list(x)) == 0:
                        root = None
                    else:
                        root = x.tag

                    if root is not None:
                        parsed_recursive = self._parse_elem(x, root)
                        parsed[root] = parsed_recursive

                else:
                    elemtype = x.get("type")
                    if x.get("name") is not None:
                        name = x.get("name")
                    else:
                        name = x.tag

                    parsed[name] = self._parse_keys(x.tag, elemtype, x.text)

        else:
            if elem.tag == "varray":
                v = []
                for x in elem:
                    v.append(self._parse_keys(x.tag, x.get("type"), x.text))
                parsed = v

            elif elem.tag == "array":
                dim = []
                field = []
                for x in elem:
                    if x.tag == "dimension":
                        dim.append(x.text.strip())
                    elif x.tag == "field":
                        field.append([x.text.strip(), x.get("type")])
                    elif x.tag == "set":
                        parsed[dim[0]] = self._parse_set(x)
                        parsed["dimension"] = dim
                        parsed["field"] = field

            else:
                parsed = self._parse_elem(elem)
        elem.clear()
        return parsed

    def _parse_set(self, set_to_parse):
        parsed = {}
        s = []
        for x in set_to_parse:
            if x.tag == "set":
                root = x.get("comment").replace(" ", "_")
                parsed[root] = self._parse_set(x)
            elif x.tag in ["r", "c"]:
                s.append(self._parse_keys(x.tag, x.get("type"), x.text))
            elif x.tag in ["rc"]:
                rc = []
                for y in x:
                    rc.append(self._parse_keys(y.tag, y.get("type"), y.text))
                s.append(rc)

        if len(s) != 0:
            parsed = s

        return parsed

    def _parse_structure(self, elem):
        structure = self._parse_elem(elem)
        if "selective" not in structure.keys():
            selective = None
        else:
            selective = structure["selective"]

        parsed = AtomicStructure(None, structure["crystal"]["basis"], self.atominfo["atoms"]["ion"],
                                 self.atominfo["simple"], False, structure["positions"], selective)
        return parsed

    def _to_np_array(self, dict, recursive_root=None):
        for root, x in dict.items():
            if type(x) is list:
                dict[root] = np.array(dict[root])
            else:
                self._to_np_array(x, root)


    def to_dic(self, title="vaspcalc", initstruc=False, parameters=False, band=False, pband=False, dos=False):
        dic = {"title": title,
               "version": self.generator["version"],
               "incar": self.incar,
               "kpoints": {"list": self.kpoints["kpointlist"],
                           "weight": self.kpoints["weights"]},
               "atominfo": {"simple": self.atominfo["simple"],
                            "full": self.atominfo["atoms"]},
               "energy": self.energy,
               "time": self.time["totalsc"][-1],
               "structure": self.fin_structure,
               }

        potcar = {}
        for x in self.atominfo["atomtypes"]["type"]:
            potcar[x[1][0]] = [x[-1], x[2][0], x[3][0]]
        dic["potcar"] = potcar

        if initstruc is True:
            dic["structure_init"] = self.init_structure
        if parameters is True:
            dic["parameters"] = self.parameters
        if band is True:
            dic["band"] = self.calculation["eigenvalues"]
        if pband is True:
            dic["pband"] = self.calculation["pband"]
        if dos is True:
            dic["dos"] = {"tdos": self.calculation["total_dos"],
                          "pdos": self.calculation["partial_dos"]
                          }
        return dic
