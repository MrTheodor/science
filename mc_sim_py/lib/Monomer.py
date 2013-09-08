from Lattice import Lattice

from monomers import *

class Monomer(object):
    """
    Single monomer object
    """
    _position = (0, 0, 0)
    _open = (0, 0, 0)
    next_direction = None ## from this monomer to the next one
    neighbours_position = []
    neighbours_open_position = []
    energy = 0
    name = None
    interaction_table = {}

    chain = None
    prev_monomer = None
    next_monomer = None

    def __init__(self):
        if self.name is None:
            raise NotImplementedError()
        if self.name == '_s':
            raise Exception("_s is reserved name for solvent!")
    
    def get_position(self):
        return self._position

    def set_position(self, value):
        refresh_nb = False
        self.next_direction = None

        if self._position != value:
            refresh_nb = True
        self._position = tuple(value)
        if refresh_nb:
            self.neighbours_position = Lattice.get_neighbours(tuple(value))
    
    position = property(get_position, set_position)

    def get_open_position(self):
        return self._open
    
    def set_open_position(self, value):
        refresh_nb = False
        if self._open != value:
            refresh_nb = True

        self._open = tuple(value)
        if refresh_nb:
            self.neighbours_open_position = Lattice.get_neighbours_open(tuple(value))
    
    open_position = property(get_open_position, set_open_position)

    def __str__(self):
        return "Monomer type: %s, position: (%s)" % (self.name, ",".join(map(str, self._position))) 


## Monomer types for Experiment4 and Experiment8

class TypeC(Monomer):
    name = "C"
    interaction_table = {'C': -1, 'O': -0.1}

class TypeO(Monomer):
    name = 'O'
    interaction_table = {'O': -1, 'C': -0.1}

class TypeA(Monomer):
    name = "A"
#    interaction_table = {'A': 0, 'B': 1, '_s': 0}
    interaction_table = {'A': 0, 'B': 1, '_s': 0}

class TypeB(Monomer):
    name = "B"
    #interaction_table = {'A': 1, 'B': 0, '_s': 0} ## old one
    interaction_table = {'A': 1, 'B': 0, '_s': 1}

class TypeC2(Monomer):
    name = "C"
    interaction_table = {'C': -1, 'O': 0}

class TypeO2(Monomer):
    name = 'O'
    interaction_table = {'O': -1, 'C': 0}

##

class TypeA2(Monomer):
    name = "A2"
    interaction_table = {'A2': 0, 'B2': 1, '_s': 0 }

class TypeB2(Monomer):
    name = "B2"
    interaction_table = {'A2': 1, 'B2': 0, '_s': 1 }



