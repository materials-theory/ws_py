import numpy as np


class Bandstructure(object):

    def __init__(self, eigenvalues: dict = None, projected: dict = None, ktrace: list = None,
                 efermi: float = 0.0, shift: bool = True, atom: list = None, orbital: list = None):
        self.eigenvalues = eigenvalues
        self.projected = projected
        self.ktrace = ktrace
        self.efermi = efermi
        self.shift = shift
        self.atom = atom
        self.orbital = orbital

    def projection_atom(self, atom):

        return

    def projection_orbital(self, orbital):
        return

