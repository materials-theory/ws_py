import xml.etree.ElementTree as ET
import numpy as np


class Vasprun(object):

    def __init__(self, infile, skip_info=False):
        self.infile = infile
        self.generator = None
        self.incar = None
        self.kpoints = None
        self.parameters = None
        self.atominfo = None
        self.init_structure = None
        self.fin_structure = None
        self.prim_structure = None
        self.calculation = None

        # tree = ET.parse(infile)
        # root = tree.getroot()

        if skip_info is False:
            self._parse_info(self.infile)

    def _parse_info(self, infile):
        for event, elem in ET.iterparse(infile):
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
                # pass
            elif elem.tag == "structure":
                if elem.attrib["name"] == "initialpos":
                    self.init_structure = self._parse_elem(elem)
                elif elem.attrib["name"] == "finalpos":
                    self.fin_structure = self._parse_elem(elem)
                elif elem.attrib["name"] == "primitive_cell":
                    self.prim_structure = self._parse_elem(elem)
        return

    def _parse_calc(self):
        return

    def _parse_keys(self, tag, type, key):
        if tag in ["v"]:
            if type == "string":
                return list(map(str, key.split()))

            elif type == "logical":
                return list(map(bool, key.split()))

            elif type == "int":
                return list(map(int, key.split()))

            else:
                return list(map(str, key.split()))

        elif tag in ["r"]:
            return list(map(float, key.split()))

        elif tag in ["i"]:
            if type == "string":
                if key is None:
                    return None
                else:
                    return key.strip(" ")

            elif type == "logical":
                if "T" in key or ".t" in key:
                    return True
                elif "F" in key or ".f" in key:
                    return False
                else:
                    raise TypeError("Unknown logical key value...")

            elif type == "int":
                return int(key)

            else:
                return float(key)

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
                        # parsed[x.tag] = self._parse_values(x.type, x.text)
                        root = None

                        # if x.tag in ["dimension", "field"]:
                        #     pass
                        # elif x.tag in ["set"]:
                        #     parsed[x.tag] = self._parse_set(x)
                        # else:
                        #     parsed[x.tag] = self._parse_values(None, x.text)
                        # root = None

                    if root is not None:
                        parsed_recursive = self._parse_elem(x, root)
                        # parsed[root] = parsed_recursive[root]

                else:
                    type = x.get("type")
                    if x.get("name") is not None:
                        name = x.get("name")
                    else:
                        name = x.tag

                    if x.tag == "i":
                        parsed[name] = self._parse_values(type, x.text)
                    elif x.tag == "v":
                        parsed[name] = self._parse_array(type, x.text)

        else:
            if elem.tag == "varray":
                v = []
                for x in elem:
                    v.append(self._parse_array(x.get("type"), x.text))
                parsed[recursive_root] = v

            elif elem.tag == "array":
                a = []
                dim = []
                field = []
                for x in elem:
                    if x.tag == "dimension":
                        dim.append(x.text)
                    elif x.tag == "field":
                        field.append([x.text, x.get("type")])
                    elif x.tag == "set":
                        if x.get("comment") is None:
                            self._parse_set(x, None)
                        else:
                            root = x.get("comment")
                            self._parse_set(x, root)
                        # a.append(root)


            else:
                parsed[recursive_root] = self._parse_elem(elem)

        return parsed

    def _parse_elem_tag_exception(self, elem, tag):
        return

    def _parse_set(self, set, recursive_root):
        s = []
        if recursive_root is None:
            for x in set:
                if x.tag in ["r"]:
                    s.append(self._parse_array())
                elif x.tag in ["c"]:
                    s.append(self._parse_)
                elif x.tag in ["rc"]:
                    pass
                elif x.tag in ["set"]:
                    pass



        return


