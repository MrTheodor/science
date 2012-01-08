# 
# Monte carlo simulation for Polymer on lattice
# Author: Jakub Krajniak <jkrajniak@gmail.com>
#
# 

import random
import math

#SEED = 7777 
#random.seed(SEED)

class Monomer(object):
    _position = (0,0,0)
    id = ''

    _interaction_table = {}

    neighbours_table = []
    
    @property
    def position(self)
        return self._position

    @position.setter
    def position(self, value):
        pass


    def get_interaction(self, monomer):
        return self._interaction_table[monomer.id]

class A(Monomer):
    id = 'A'


class Lattice(object):
    """
    3-D lattice, cubic

    """
    box_size = 0
    p_size = 0
    energy = 0
   
    _polymer_chain = []
    
    def __init__(self, box_size, polymer_chain_size):
        self.box_size = box_size
        self.size = polymer_chain_size
        
        self.randomize()
    
    def randomize(self):
        """
        Randomize polymer

        """
        ##  pick random starting position
        pos = [ random.randrange(0, self.size) for tmp in range(3) ]
        
        number_units = 1
        while number_units < self.size:
            
            ## pickup random position
            nb_pos = self.neighbours(pos, random.randrange(6))
            print nb_pos

            if nb_pos in self._polymer_chain:
                continue

            self._polymer_chain.append(nb_pos)
            pos = nb_pos[:]

            number_units += 1



    def neighbours(self, position, number):

        positions = (self._get_left_x, 
            self._get_right_x, 
            self._get_top_y, 
            self._get_down_y, 
            self._get_top_z, 
            self._get_down_z)

        return positions[number](position)
    
    def snapshot(self):
        """
        Return snapshot of configuration to RasMol format
        """

        file = open('snapshot', 'w+')
        file.writelines("%d\n\n" % len(self._polymer_chain))
        [ file.writelines("C %s\n" % (" ".join(map(str, self._polymer_chain[idx])))) for idx in range(len(self._polymer_chain)) ]
        file.close()
        

    
    ## get positions
    def _get_left_x(self, cp):
        return ((cp[0] - 1) % self.box_size, cp[1], cp[2])

    def _get_right_x(self, cp):
        return ((cp[0] + 1) % self.box_size, cp[1], cp[2])

    def _get_top_y(self, cp):
        return (cp[0], (cp[1] + 1) % self.box_size, cp[2])

    def _get_down_y(self, cp):
        return (cp[0], (cp[1] - 1) % self.box_size, cp[2])
    
    def _get_top_z(self, cp):
        return (cp[0], cp[1], (cp[2] + 1) % self.box_size)
    
    def _get_down_z(self, cp):
        return (cp[0], cp[1], (cp[2] - 1) % self.box_size)

    
class MC(object):
    
    ## MC config
    sweeps = 0
    
    ## physical constants
    _T = None

    ## improvments, precompute exp vector
    _exp = [0]*9

    ## stats
    acceptance = 0

    def __init__(self, size, polymer_size, T, sweeps):
        """
        @parms
            size - size of box
            polymer_size - size of polymer chain
            T - temperature
            sweeps - number of sweeps
        """
        self.lattice_size = size
        self.polymer_size = polymer_size
        self.lattice_size_sq = size**3
        self.T = T

        self.lattice = Lattice(size, polymer_size)
        
        self.sweeps = sweeps

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, value):
        self._T = value

    def _step(self):
        """
        One step
        Return if accept or reject new configuration
        """

        pass


    def run(self):
        """
        Run simulation
        """
        
        pass

mc = MC(50, 30, 1, 100)
mc.lattice.snapshot()
    

## for doc test
#if __name__ == "__main__":
#    import doctest
#    doctest.testmod(extraglobals={'lattice': Lattice(10)})
