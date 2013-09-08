"""
Definition of monomers for HP Model
"""

from lib.Monomer import Monomer

class TypeH(Monomer):
    name = "H"
    interaction_table = {'H': -1, 'P': 0}

class TypeP(Monomer):
    name = "P"
    interaction_table = {'H': 0, 'P': 0}
