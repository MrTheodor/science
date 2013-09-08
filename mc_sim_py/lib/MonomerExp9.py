from Monomer import Monomer

class TypeA(Monomer):
    name = "A"
    interaction_table = {'A': 0, 'B': 0.5, '_s': 0}

class TypeB(Monomer):
    name = 'B'
    interaction_table = {'A': 0.5, 'B': 0, '_s': 0.5}
