from Box import Lattice

from monomers import *

cdef class Monomer(object):
    """
    Single monomer object
    """
    _position = (0, 0, 0)
    _open = (0, 0, 0)
    cdef public int next_direction ## from this monomer to the next one
    neighbours_position = []
    cdef float energy
    cdef public name
    interaction_table = {}

    cdef public chain
    cdef public Monomer prev_monomer
    cdef public Monomer next_monomer

    def __init__(self):
        self.chain = None
        self.name = ""
        if self.name is None:
            raise NotImplementedError()
        if self.name == '_s':
            raise Exception("_s is reserved name for solvent!")

        self.next_direction = -1
    
    def get_position(self):
        return self._position
    """
    @property
    def next_direction(self):
        #if self._next_direction is None:
        #    if self.next_monomer is not None:
        #        self._next_direction = self.neighbours_position.index(self.next_monomer._position)
        #
        return self._next_direction

    @next_direction.setter
    def next_direction(self, value):
        self._next_direction = value
    """

    def set_position(self, value):
        self._next_direction = None
        self._position = tuple(value)
        self.neighbours_position = Lattice.get_neighbours(tuple(value))[0]
        self.neighbours_open_position = Lattice.get_neighbours(tuple(value))[1]
    
    position = property(get_position, set_position)

    def get_open_position(self):
        return self._open
    
    def set_open_position(self, value):
        self._open = tuple(value)
    
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

##

class TypeA2(Monomer):
    name = "A2"
    interaction_table = {'A2': 0, 'B2': 1, '_s': 0 }

class TypeB2(Monomer):
    name = "B2"
    interaction_table = {'A2': 1, 'B2': 0, '_s': 1 }


class TypeCys(Monomer):
    name = "Cys"
    interaction_table = {'Cys': -5.44, 'Ile': -5.5, 'Ser': -2.86, 'Val': -4.96, 'Gly': -3.16, 'Gln': -2.85, 'Pro': -3.07, 'Lys': -1.95, 'Thr': -3.11, 'Phe': -5.8, 'Ala': -3.57, 'Met': -4.99, 'Asp': -2.41, 'Leu': -5.83, 'His': -3.6, 'Arg': -2.57, 'Trp': -4.95, 'Glu': -2.27, 'Asn': -2.59, 'Tyr': -4.16}
    
class TypeIle(Monomer):
    name = "Ile"
    interaction_table = {'Cys': 0.49, 'Ile': -6.54, 'Ser': -3.52, 'Val': -6.05, 'Gly': -3.78, 'Gln': -3.67, 'Pro': -3.76, 'Lys': -3.01, 'Thr': -4.03, 'Phe': 0.06, 'Ala': -4.58, 'Met': -0.01, 'Asp': -3.17, 'Leu': -7.04, 'His': -4.14, 'Arg': -3.63, 'Trp': -5.78, 'Glu': -3.27, 'Asn': -3.24, 'Tyr': -5.25}
    
class TypeSer(Monomer):
    name = "Ser"
    interaction_table = {'Cys': 0.69, 'Ile': 0.59, 'Ser': -1.67, 'Val': 0.55, 'Gly': 0.14, 'Gln': -1.49, 'Pro': -1.57, 'Lys': -1.05, 'Thr': -0.06, 'Phe': 0.44, 'Ala': 0.18, 'Met': 0.53, 'Asp': -1.63, 'Leu': 0.6, 'His': -2.11, 'Arg': -1.62, 'Trp': 0.38, 'Glu': -1.48, 'Asn': -1.58, 'Tyr': 0.14}
    
class TypeVal(Monomer):
    name = "Val"
    interaction_table = {'Cys': 0.52, 'Ile': -0.01, 'Ser': -3.05, 'Val': -5.52, 'Gly': -3.38, 'Gln': -3.07, 'Pro': -3.32, 'Lys': -2.49, 'Thr': -3.46, 'Phe': 0.1, 'Ala': -4.04, 'Met': 0.18, 'Asp': -2.48, 'Leu': -0.04, 'His': -3.58, 'Arg': -3.07, 'Trp': -5.18, 'Glu': -2.67, 'Asn': -2.83, 'Tyr': -4.62}
    
class TypeGly(Monomer):
    name = "Gly"
    interaction_table = {'Cys': 0.68, 'Ile': 0.62, 'Ser': -1.82, 'Val': 0.51, 'Gly': -2.24, 'Gln': -1.66, 'Pro': -1.87, 'Lys': -1.15, 'Thr': -2.08, 'Phe': 0.62, 'Ala': 0.18, 'Met': 0.46, 'Asp': -1.59, 'Leu': 0.65, 'His': -2.15, 'Arg': -1.72, 'Trp': 0.24, 'Glu': -1.22, 'Asn': -1.74, 'Tyr': 0.2}
    
class TypeGln(Monomer):
    name = "Gln"
    interaction_table = {'Cys': 0.64, 'Ile': 0.37, 'Ser': 0.11, 'Val': 0.46, 'Gly': 0.24, 'Gln': -1.54, 'Pro': -1.73, 'Lys': -1.29, 'Thr': -0.08, 'Phe': 0.3, 'Ala': 0.24, 'Met': 0.2, 'Asp': -1.46, 'Leu': 0.42, 'His': -1.98, 'Arg': -1.8, 'Trp': 0.19, 'Glu': -1.42, 'Asn': -0.1, 'Tyr': -0.12}
    
class TypePro(Monomer):
    name = "Pro"
    interaction_table = {'Cys': 0.53, 'Ile': 0.39, 'Ser': 0.14, 'Val': 0.31, 'Gly': 0.13, 'Gln': -0.08, 'Pro': -1.75, 'Lys': -0.04, 'Thr': 0.04, 'Phe': 0.25, 'Ala': 0.2, 'Met': 0.16, 'Asp': 0.14, 'Leu': 0.35, 'His': 0.15, 'Arg': -0.05, 'Trp': -0.33, 'Glu': 0.07, 'Asn': 0.18, 'Tyr': -0.23}
    
class TypeLys(Monomer):
    name = "Lys"
    interaction_table = {'Cys': 0.83, 'Ile': 0.32, 'Ser': -0.15, 'Val': 0.33, 'Gly': 0.03, 'Gln': -0.46, 'Pro': -0.97, 'Lys': -0.12, 'Thr': -0.19, 'Phe': 0.33, 'Ala': 0.11, 'Met': 0.31, 'Asp': -1.01, 'Leu': 0.37, 'His': 0.23, 'Arg': 0.24, 'Trp': -0.1, 'Glu': -1.28, 'Asn': -0.3, 'Tyr': -0.46}
    
class TypeThr(Monomer):
    name = "Thr"
    interaction_table = {'Cys': 0.67, 'Ile': 0.3, 'Ser': -1.96, 'Val': 0.36, 'Gly': 0.1, 'Gln': -1.9, 'Pro': -1.9, 'Lys': -1.31, 'Thr': -2.12, 'Phe': 0.41, 'Ala': 0.1, 'Met': 0.28, 'Asp': -1.8, 'Leu': 0.4, 'His': -2.42, 'Arg': -1.9, 'Trp': 0.37, 'Glu': -1.74, 'Asn': -1.88, 'Tyr': 0.13}
    
class TypePhe(Monomer):
    name = "Phe"
    interaction_table = {'Cys': 0.54, 'Ile': -6.84, 'Ser': -4.02, 'Val': -6.29, 'Gly': -4.13, 'Gln': -4.1, 'Pro': -4.25, 'Lys': -3.36, 'Thr': -4.28, 'Phe': -7.26, 'Ala': -4.81, 'Met': -0.2, 'Asp': -3.48, 'Leu': -7.28, 'His': -4.77, 'Arg': -3.98, 'Trp': -6.16, 'Glu': -3.56, 'Asn': -3.75, 'Tyr': -5.66}
    
class TypeAla(Monomer):
    name = "Ala"
    interaction_table = {'Cys': 0.51, 'Ile': 0.05, 'Ser': -2.01, 'Val': 0.08, 'Gly': -2.31, 'Gln': -1.89, 'Pro': -2.03, 'Lys': -1.31, 'Thr': -2.32, 'Phe': 0.17, 'Ala': -2.72, 'Met': 0.15, 'Asp': -1.7, 'Leu': 0.13, 'His': -2.41, 'Arg': -1.83, 'Trp': 0.07, 'Glu': -1.51, 'Asn': -1.84, 'Tyr': 0.09}
    
class TypeMet(Monomer):
    name = "Met"
    interaction_table = {'Cys': 0.46, 'Ile': -6.02, 'Ser': -3.03, 'Val': -5.32, 'Gly': -3.39, 'Gln': -3.3, 'Pro': -3.45, 'Lys': -2.48, 'Thr': -3.51, 'Phe': -6.56, 'Ala': -3.94, 'Met': -5.46, 'Asp': -2.57, 'Leu': -6.41, 'His': -3.98, 'Arg': -3.12, 'Trp': -5.55, 'Glu': -2.89, 'Asn': -2.95, 'Tyr': -4.91}
    
class TypeAsp(Monomer):
    name = "Asp"
    interaction_table = {'Cys': 0.91, 'Ile': 0.71, 'Ser': -0.19, 'Val': 0.89, 'Gly': 0.13, 'Gln': -0.09, 'Pro': -1.33, 'Lys': -1.68, 'Thr': -0.14, 'Phe': 0.75, 'Ala': 0.26, 'Met': 0.77, 'Asp': -1.21, 'Leu': 0.89, 'His': -2.32, 'Arg': -2.29, 'Trp': 0.3, 'Glu': -1.02, 'Asn': -0.24, 'Tyr': -0.07}
    
class TypeLeu(Monomer):
    name = "Leu"
    interaction_table = {'Cys': 0.57, 'Ile': -0.08, 'Ser': -3.92, 'Val': -6.48, 'Gly': -4.16, 'Gln': -4.04, 'Pro': -4.2, 'Lys': -3.37, 'Thr': -4.34, 'Phe': 0.03, 'Ala': -4.91, 'Met': 0.01, 'Asp': -3.4, 'Leu': -7.37, 'His': -4.54, 'Arg': -4.03, 'Trp': -6.14, 'Glu': -3.59, 'Asn': -3.74, 'Tyr': -5.67}
    
class TypeHis(Monomer):
    name = "His"
    interaction_table = {'Cys': 0.65, 'Ile': 0.66, 'Ser': 0.26, 'Val': 0.7, 'Gly': 0.5, 'Gln': 0.31, 'Pro': -2.25, 'Lys': -1.35, 'Thr': 0.16, 'Phe': 0.39, 'Ala': 0.47, 'Met': 0.28, 'Asp': -0.19, 'Leu': 0.67, 'His': -3.05, 'Arg': -2.16, 'Trp': 0.08, 'Glu': -0.16, 'Asn': 0.29, 'Tyr': 0.09}
    
class TypeArg(Monomer):
    name = "Arg"
    interaction_table = {'Cys': 0.93, 'Ile': 0.41, 'Ser': -0.01, 'Val': 0.47, 'Gly': 0.18, 'Gln': -0.26, 'Pro': -1.7, 'Lys': -0.59, 'Thr': -0.07, 'Phe': 0.42, 'Ala': 0.3, 'Met': 0.38, 'Asp': -0.91, 'Leu': 0.43, 'His': 0.14, 'Arg': -1.55, 'Trp': -0.11, 'Glu': -1.04, 'Asn': -0.02, 'Tyr': -0.3}
    
class TypeTrp(Monomer):
    name = "Trp"
    interaction_table = {'Cys': 0.3, 'Ile': 0.02, 'Ser': -2.99, 'Val': 0.11, 'Gly': -3.42, 'Gln': -3.11, 'Pro': -3.73, 'Lys': -2.69, 'Thr': -3.22, 'Phe': 0.0, 'Ala': -3.82, 'Met': -0.29, 'Asp': -2.84, 'Leu': 0.08, 'His': -3.98, 'Arg': -3.41, 'Trp': -5.06, 'Glu': -2.99, 'Asn': -3.07, 'Tyr': -4.66}
    
class TypeGlu(Monomer):
    name = "Glu"
    interaction_table = {'Cys': 0.91, 'Ile': 0.46, 'Ser': -0.19, 'Val': 0.55, 'Gly': 0.36, 'Gln': -0.19, 'Pro': -1.26, 'Lys': -1.8, 'Thr': -0.22, 'Phe': 0.52, 'Ala': 0.3, 'Met': 0.3, 'Asp': 0.05, 'Leu': 0.55, 'His': -2.15, 'Arg': -2.27, 'Trp': 0.0, 'Glu': -0.91, 'Asn': -0.21, 'Tyr': -0.25}
    
class TypeAsn(Monomer):
    name = "Asn"
    interaction_table = {'Cys': 0.97, 'Ile': 0.87, 'Ser': 0.1, 'Val': 0.77, 'Gly': 0.22, 'Gln': -1.71, 'Pro': -1.53, 'Lys': -1.21, 'Thr': 0.02, 'Phe': 0.72, 'Ala': 0.36, 'Met': 0.62, 'Asp': -1.68, 'Leu': 0.79, 'His': -2.08, 'Arg': -1.64, 'Trp': 0.3, 'Glu': -1.51, 'Asn': -1.68, 'Tyr': 0.17}
    
class TypeTyr(Monomer):
    name = "Tyr"
    interaction_table = {'Cys': 0.64, 'Ile': 0.11, 'Ser': -2.78, 'Val': 0.23, 'Gly': -3.01, 'Gln': -2.97, 'Pro': -3.19, 'Lys': -2.6, 'Thr': -3.01, 'Phe': 0.05, 'Ala': -3.36, 'Met': -0.1, 'Asp': -2.76, 'Leu': 0.1, 'His': -3.52, 'Arg': -3.16, 'Trp': -0.04, 'Glu': -2.79, 'Asn': -2.76, 'Tyr': -4.17}
