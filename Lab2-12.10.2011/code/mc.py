# 
# Monte carlo simulation for Ising model
# Author: Jakub Krajniak <jkrajniak@gmail.com>
#
# 

import random
import math
import itertools

#SEED = 7777 
#random.seed(SEED)

class Lattice(object):
    """
    2-D lattice

    """
    size = 0
    energy = 0
    _lattice = []
    _dictionary = []
    _all_coordinate = []
    
    def __init__(self, lattice_size):
        self.size = lattice_size
        self._dictionary = (-1, 1)
        
        self.randomize()
        
        self.total_energy()

        ## cross product to get all coordinate
        self._all_coordinate = list(itertools.product(range(self.size), range(self.size)))
    
    def randomize(self):
        """
        Randomize lattice

        >>> len(t._lattice)
        10
        >>> len(t._lattice[0])
        10

        """

        # initialize, clear lattice

        self._lattice = []
        for n in range(self.size):
            self._lattice.append([ random.choice(self._dictionary) for x in range(self.size) ])
    
    def total_energy(self):
        """
        Return total energy of lattice
        """

        e = 0.0

        for pos in self._all_coordinate:
            e += sum(self[pos] * self.get_nb_values(pos))

        self.energy = e

        return e


    def energy_difference(self, position):
        """
        Return energy change between two configuration that differ only by 
        one spin (on position 'position')

        """
        delta = 2 * self[position] * sum(self.get_nb_values(position))
        
        return delta

    def total_magnetization(self, per_spin = False):
        """
        Return total magnetization
        """
        M = sum([ self[pos] for pos in self._all_coordinate ])

        if per_spin:
            M = M / float(self.size**2)

        return M


    def flip_spin(self, position, random=False):
        """
        Flip spin on 'position'
        Return new spin value

        """

        if random:
            position = self.get_random_spin_position()

        self[position] *= -1

        return self[position]
    
    def get_nb_position(self, position):
        """
        Return neigbours spin value in the format of tuple:
         (top, right, down, left)
        >>> n = lattice.get_nb_position((0,0))
        >>> n == (1, 1, self.size - 1, self.size - 1)

        """

        order = (self._get_top, self._get_right, self._get_down, self._get_left)

        return tuple([ fun(position) for fun in order ])

    def get_nb_values(self, position):
        """
        Return value of neigbours spins around position
        
        """

        position = self.get_nb_position(position)
        values = tuple([ self[pos] for pos in position ])

        return values
    
    def get_random_spin_position(self):
        """
        Return random spin position
        """

        try:
            random_position = random.choice(list(self._all_coordinate))
        except:
            from pdb import set_trace;set_trace() ############################## Breakpoint ##############################

        return random_position
    
    ## get positions
    def _get_left(self, cp):
        return ((cp[0] - 1) % self.size, cp[1])

    def _get_right(self, cp):
        return ((cp[0] + 1) % self.size, cp[1])

    def _get_top(self, cp):
        return (cp[0], (cp[1] + 1) % self.size)

    def _get_down(self, cp):
        return (cp[0], (cp[1] - 1) % self.size)

    
    ## for containers
    def __getitem__(self, key):
        try:
            return self._lattice[key[0]][key[1]]
        except:
            from pdb import set_trace;set_trace() ############################## Breakpoint ##############################


    def __setitem__(self, key, value):
        self._lattice[key[0]][key[1]] = value

    def __delitem(self, key):
        self._lattice[key[0]][key[1]] = None

    def __contains__(self, key):
        x = key[0]
        y = key[1]
        size = self.size - 1
        if x > 0 and x < size and y > 0 and y < size:
            return True
        else:
            return False

    def copy(self):
        import copy
        
        data = copy.copy(self._lattice)
        dictionary = copy.copy(self._dictionary)
        self._lattice = None
        self._dictionary = []
        try:
            c = copy.copy(self)
            c._lattice = data
            c._dictionary = dictionary
        finally:
            self._lattice = data
            self._dictionary = dictionary

        return c



class MC(object):
    
    ## MC config
    sweeps = 0
    eq_time = 0
    
    # lattice
    lattice_size = 0
    lattice_size_sq = 0
    lattice = []
    
    ## physical constants
    _T = None
    J = 1

    ## improvments, precompute exp vector
    _exp = [0]*9

    ## stats
    acceptance = 0

    def __init__(self, size, T, sweeps, eq_time = 0):
        """
        @parms
            size - size of lattice
            T - temperature
            sweeps - number of sweeps
            eq_time - equilibrium time
        """
        self.lattice_size = size
        self.lattice_size_sq = size**2
        self.T = T

        self.lattice = Lattice(size)
        self.lattice.energy *= self.J

        self.sweeps = sweeps
        self.eq_time = eq_time

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, value):
        self._T = value

        ## precomputate exponential values for given temperatur
        self._exp[0] = 1
        self._exp[4] = math.exp( (-1.0/self._T)*4 )
        self._exp[8] = math.exp( (-1.0/self._T)*8 )


    def _step(self):
        """
        One step
        Flip random spin and check energy difference
        if is lower than 0, accept, otherwise compute exp(-dE/kbT) and check

        Return if accept or reject new configuration
        """

        spin_pos = self.lattice.get_random_spin_position()
        deltaE = self.J * self.lattice.energy_difference(spin_pos)
        
        accept = False

        if deltaE < 0:
            accept = True
        else:
            r = random.random()
            if r <= self._exp[deltaE]:
                accept = True

        if accept: ## flip spin
            self.lattice.energy += deltaE
            self.lattice.flip_spin(spin_pos)

        return accept

    def run(self):
        """
        Run simulation
        """

        magnetization_vector = []
        energy_vector = []

        for sweep in xrange(self.sweeps):
            
            for i in xrange(self.lattice_size_sq):
                step = self._step()
                self.acceptance += int(step)

            ##
            if sweep >= self.eq_time:
                magnetization_vector.append(self.lattice.total_magnetization())
                energy_vector.append(self.lattice.energy)
        
        return {
            'magnetization': magnetization_vector,
            'energy': energy_vector
        }


## for doc test
#if __name__ == "__main__":
#    import doctest
#    doctest.testmod(extraglobals={'lattice': Lattice(10)})
