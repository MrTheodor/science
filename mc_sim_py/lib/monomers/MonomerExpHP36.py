"""
Definition of monomers for HP Model for temp = 36
"""

from lib.Monomer import Monomer

epsilon = 10**-23 

class TypeH(Monomer):
    name = "H"
    interaction_table = {'H': -epsilon, 'P': 0}

class TypeP(Monomer):
    name = "P"
    interaction_table = {'H': 0, 'P': 0}
