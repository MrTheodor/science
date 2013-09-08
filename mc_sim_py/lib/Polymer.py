'''
Created on 2 Nov 2011

@author: Jakub Krajniak <jkrajniak@gmail.com>
'''

import copy
import itertools
import logging
import numpy as np
import random

#from random import randrange
from numpy.random import randint as randrange

from lib.Lattice import Lattice
from lib.Monomer import Monomer

import settings

import time

class Chain(object):
    """
    Polymer sequence of PolymerBlocks
    """
    ## private
    _total_length = 0
    
    positions_list = []
    _tmp = {}
    
    ## public
    chain = []
    _monomer_types = set()
    _chain_structure = {}
    box = None
    total_energy = 0
    name = ""
    idx = 0
    
    #
    last_values = None
    
    TYPE_OF_MOVE = ['begin_rotate', 'end_rotate','segment_rotate', 'end_snake', 'begin_snake', 'pull_move']
    MOVE_STATS = [0,0,0,0,0,0]
    B_ROTATE = 0
    E_ROTATE = 1
    S_ROTATE = 2
    E_SNAKE = 3
    B_SNAKE = 4
    PULL_MOVE = 5

    MOVE_INSIDE = [ B_ROTATE, E_ROTATE, S_ROTATE, ]
    MOVE_ONE_STEP = [ E_SNAKE, B_SNAKE, PULL_MOVE ]
    
    temperature_history = []
    move_name = None
    type_of_move = None

    ## simple stats
    min_energy = 0.0
    min_rg = 0.0
    
    def __init__(self, chain_sequence, name=""):
        self.chain = chain_sequence
        self._total_length = len(chain_sequence)
        
        self.name = name
        if self.name == "":
            self.name = [ x.name for x in chain_sequence ][0:10]
        
        self.positions_list = []
        for m in self.chain:
            self.positions_list.append(m._position)
            self._monomer_types.add(m.name)

        self.name = name
        if self.name == "":
            self.name = "".join([ x for x in self._monomer_types ][0:10])
 
        if settings.DEBUG or settings.PROFILER:
            random.seed(1)
            np.random.seed(1)

    def move(self, type_of_move):
        """
        @int type_of_move: type of move, if -1 then it is random, otherwise you can put values from TYPE_OF_MOVE
        """
        #type_of_move = self.PULL_MOVE

        self.move_name = self.TYPE_OF_MOVE[type_of_move]
        self.type_of_move = type_of_move

        ## performe move
        
        #fnc = getattr(self, self.move_name)
        #value = fnc()
        if type_of_move == self.B_ROTATE:
            value = self.begin_rotate()
        elif type_of_move == self.E_ROTATE:
            value = self.end_rotate()
        elif type_of_move == self.S_ROTATE:
            value = self.segment_rotate()
        elif type_of_move == self.E_SNAKE:
            value = self.end_snake()
        elif type_of_move == self.B_SNAKE:
            value = self.begin_snake()
        elif type_of_move == self.PULL_MOVE:
            value = self.pull_move()
            
        if value: self.MOVE_STATS[type_of_move] += 1

        return value
    
    def valid_chain(self, assertion=True):
        """
        Check if chain is valid
        """
        
        if self._total_length != len(self.chain):
            raise Exception("Length is not correct")
        
        ## check next_direction
        """
        for x in xrange(0, self._total_length - 1):
            curr = self.chain[x]
            next_m = self.chain[x+1]
            if curr.neighbours_position[curr.next_direction] != next_m._position:
                raise Exception("Problem with next_direction for %d" % x )
        """

        build_position_list = [ x._position for x in self.chain ]
        if build_position_list != self.positions_list:
            raise Exception('Chain positions_list problem')
        
        first_distance = np.array(self.chain[0]._open) - np.array(self.chain[1]._open)
        first_distance = first_distance.dot(first_distance)
        if first_distance != Lattice.neighbour_distance:
            raise Exception("Distance not valid, should be: %f, but is %f" % (Lattice.neighbour_distance, first_distance))
        for idx in range(1, self._total_length - 1):
            prev_monomer = self.chain[idx-1]
            curr_monomer = self.chain[idx]
            
            distance = np.array(prev_monomer._open) - np.array(curr_monomer._open)
            distance = distance.dot(distance)
            if first_distance != distance:
                raise Exception("Distance between (%d, %d, %d) and (%d, %d, %d) is not correct: %f for idx: %d" % (prev_monomer._open + curr_monomer._open + (distance,idx)))

            if prev_monomer.neighbours_position[prev_monomer.next_direction] != curr_monomer._position:
                raise Exception("Problem with next_direction for %d" % (idx-1))
            
        ## check if it is SAW
        ## build position_list of all occupied place
        ## check if every position appear only once on the list
        positions_list = []
        for c in self.box.chain_list:
            positions_list += c.positions_list

        max_count = max([ positions_list.count(e) for e in positions_list ])
        if max_count > 1:
            raise Exception("SAW problem!")

        return True
    
    def get_direction(self, monomer_b, monomer_a):
        """
        Return direction between two monomers, use to avoid some moves
        End point: b
        Start point: a
        """
        
        if not isinstance(monomer_a, Monomer):
            monomer_a = self[monomer_a]
        if not isinstance(monomer_b, Monomer):
            monomer_b = self[monomer_b]
        
        vector = np.array(monomer_b._open) - np.array(monomer_a._open)
        
        try:
            direction = Lattice.get_direction(tuple(vector))
        except Exception, e:
            #logging.debug("get_direction for (%d, %d, %d) - (%d, %d, %d): (%d, %d, %d)" % (monomer_b._open, monomer_a._open, vector))
            raise e
        
        return direction
    
    def end_rotate(self, monomer_idx = None):
        """
        rotate with begin
        """
        random_direction = randrange(0, Lattice.z)
        new_position = self.chain[-2].neighbours_position[random_direction]
        if self.box.is_occupied(new_position) or new_position in self.positions_list:
            #logging.info("Reject rotate; position: (%d, %d, %d) is occupied" % (new_position))
            return False
        
        monomer = self.chain[-1]
        # copy
        old_position = monomer._position
        old_open_position = monomer._open
        old_next_direction = self.chain[-2].next_direction
        old_energy = self.total_energy
        
        self.chain[-1].position = new_position
        self.chain[-1]._open = Lattice.get_open_coordinate(self.chain[-2]._open, random_direction)
        self.chain[-2].next_direction = random_direction
        self.positions_list[-1] = new_position

        
        perform_move = True
        new_energy = None
        dE = None

        if self.box.athermal_state == False:
            new_energy = self.box.calculate_chain_energy(self)
        
            dE = new_energy - old_energy
        
            if self.box.accept(dE, None):
                perform_move = True
            else:
                perform_move = False

        if perform_move:
            if not new_energy:
                new_energy = self.box.calculate_chain_energy(self)
            self.total_energy = new_energy
        else:
            monomer.position = old_position
            monomer._open = old_open_position
            self.chain[-2].next_direction = old_next_direction
            self.chain[-1] = monomer
            self.positions_list[-1] = old_position
            self.total_energy = old_energy

        return perform_move
    
    def begin_rotate(self, monomer_idx = None):
        """
        rotate with begin
        """
        random_direction = randrange(0, Lattice.z)
        #new_position = Lattice.get_coordinate(self.chain[1]._position, random_direction)
        new_position = self.chain[1].neighbours_position[random_direction]
        if self.box.is_occupied(new_position) or new_position in self.positions_list:
            #logging.info("Reject rotate; position: (%d, %d, %d) is occupied" % (new_position))
            return False
        
        new_monomer = self.chain[0]

        old_position = self.chain[0]._position
        old_open_position = self.chain[0]._open
        old_next_direction = self.chain[0].next_direction
        old_energy = self.total_energy

        self.chain[0].position = new_position
        self.chain[0]._open = Lattice.get_open_coordinate(self.chain[1]._open, random_direction)
        self.chain[0].next_direction = Lattice.TURN_DIRECTION[random_direction]
        self.positions_list[0] = new_position
        
        perform_move = True
        new_energy = None
        
        if self.box.athermal_state == False:
        
            new_energy = self.box.calculate_chain_energy(self)
        
            dE = new_energy - old_energy
        
            if self.box.accept(dE, None):
                perform_move = True
            else:
                perform_move = False

        if perform_move:
            if not new_energy:
                new_energy = self.box.calculate_chain_energy(self)
            self.total_energy = new_energy
        else:
            new_monomer.position = old_position
            new_monomer._open = old_open_position
            new_monomer.next_direction = old_next_direction
            self.chain[0] = new_monomer
            self.positions_list[0] = old_position
            # restore old energy
            self.total_energy = old_energy

        return perform_move
    
    def begin_snake(self, monomer_idx = None):
        return self._snake(0)
    
    def end_snake(self, monomer_idx = None):
        return self._snake(1)
    
    def _snake(self, move_type):
        """
        Snake move
        """
        
        #direction_to_choice = range(0, Lattice.z)
        if settings.DEBUG: self.valid_chain()    
        if move_type == 0:
            monomer = self.chain[0]
        elif move_type == 1:
            monomer = self.chain[-1]
            
        #random new position
        random_direction = randrange(0, Lattice.z)
        new_position = monomer.neighbours_position[random_direction]
        open_position = Lattice.get_open_coordinate(monomer._open, random_direction)

        if self.box.is_occupied(new_position) or new_position in self.positions_list:
            return False
        
        # copy
        old_position_list = self.positions_list[:]
        old_open_position_list = [ x._open for x in self.chain ]
        old_next_direction_list = [ x.next_direction for x in self.chain ]
        old_energy = self.total_energy
        
        self.positions_list = []
        
        """
        Update position of monomers in chain (copy)
        """
        if move_type == 0: # left
            self.positions_list = [new_position]
            
            for idx in range(1, self._total_length):
                self.chain[idx].position = old_position_list[idx - 1]
                self.chain[idx]._open = old_open_position_list[idx - 1]
                self.chain[idx].next_direction = old_next_direction_list[idx - 1]
            
            self.positions_list += old_position_list[0: self._total_length - 1]
            self.chain[0].position = new_position
            self.chain[0]._open = open_position
            self.chain[0].next_direction = self.box.lattice.TURN_DIRECTION[random_direction]
            
        elif move_type == 1: # right
            for idx in range(0, self._total_length - 1):
                self.chain[idx].position = old_position_list[idx + 1]
                self.chain[idx]._open = old_open_position_list[idx + 1]
                self.chain[idx].next_direction = old_next_direction_list[idx + 1]
                
            self.chain[-1].position = new_position
            self.chain[-1]._open = open_position
            self.chain[-2].next_direction = random_direction
            
            self.positions_list = old_position_list[1:] + [new_position]
       

        new_energy = self.box.calculate_chain_energy(self)
        
        dE = new_energy - old_energy
         
        if dE <= 0.0 or self.box.athermal_state:
            if settings.DEBUG: self.valid_chain()
            self.total_energy = new_energy

            return True
        elif self.box.accept(dE, None):
            if settings.DEBUG: self.valid_chain()
            self.total_energy = new_energy
            return True
        else:
            for idx in xrange(self._total_length):
                self.chain[idx].position = old_position_list[idx]
                self.chain[idx]._open = old_open_position_list[idx]
                self.chain[idx].next_direction = old_next_direction_list[idx]
            
            self.positions_list = old_position_list[:]
            self.total_energy = old_energy
            
            if settings.DEBUG: self.valid_chain()
            
        return False
    
    def segment_rotate(self, monomer_idx = None):
        """
        @Monomer|int monomer: monomer or monomer_idx on chain_list
        
        Rotate :monomer
        """
        if monomer_idx is None:
            monomer_idx = randrange(0, self._total_length - 1)
        
        if monomer_idx == 0 or monomer_idx == self._total_length - 1:
            #logging.debug("Tried to rotate monomer of id: %d which is in the end of polymer chain. It doesn't make sense" % monomer_idx)
            return False

        if monomer_idx + 1 >= self._total_length:
            return False

        prev_monomer = self.chain[monomer_idx - 1]
        rotate_monomer = self.chain[monomer_idx]
        next_monomer = self.chain[monomer_idx + 1]
        
        ## we can chouse only position from this set
        common_neighbours_positions = set(prev_monomer.neighbours_position).intersection(next_monomer.neighbours_position)
        number_of_positions = len(common_neighbours_positions)
        
        if number_of_positions == 0:
            return False
        
        new_position = list(common_neighbours_positions)[randrange(0, number_of_positions)]
        if not new_position or len(new_position) != 3:
            if __debug__: logging.debug("New_position %s:" % str(new_position))
            return False
        
        try:
            direction = prev_monomer.neighbours_position.index(new_position)
        except Exception, e:
            if __debug__: logging.error("Reject segment, error in direction, (%d, %d, %d)" % (new_position))
            return False
            
        open_position = Lattice.get_open_coordinate(prev_monomer._open, direction)
        
        if self.box.is_occupied(new_position):
            #logging.debug("Reject segment; Position: (%d, %d, %d) is occupied" % (new_position))
            return False
        
        # make copy of old monomer
        old_position = rotate_monomer._position
        old_open_position = rotate_monomer._open
        old_prev_next_direction = self.chain[monomer_idx - 1].next_direction
        old_next_direction = rotate_monomer.next_direction

        old_energy = self.total_energy

        self.chain[monomer_idx].position = new_position
        self.chain[monomer_idx]._open = open_position
        self.chain[monomer_idx].next_direction = self.chain[monomer_idx].neighbours_position.index(next_monomer._position)
        self.chain[monomer_idx - 1].next_direction = prev_monomer.neighbours_position.index(new_position)
        self.positions_list[monomer_idx] = new_position

        
        new_energy = self.box.calculate_chain_energy(self)
        
        #print self.total_energy - current_energy + new_energy, self.total_energy
        
        dE = new_energy - old_energy
        accept = False
        if dE <= 0.0 or self.box.athermal_state:
            accept = True
        elif self.box.accept(dE, new_position):
            accept = True
        
        if accept:
            #self.valid_chain()
            self.total_energy = new_energy
            return True
        else: # revert
            self.chain[monomer_idx].position = old_position
            self.positions_list[monomer_idx] = old_position
            self.chain[monomer_idx]._open = old_open_position
            self.chain[monomer_idx - 1].next_direction = old_prev_next_direction
            self.chain[monomer_idx].next_direction = old_next_direction
            self.total_energy = old_energy

        
        return False
    
    def pull_move(self, monomer_idx = None):
        """
        Pull move, N. Lesh, M. Mitzemacher, S. Whitesides
        """
        ## flags
        valid_configuration = False
        if settings.DEBUG: self.valid_chain()
        #if monomer_idx is None:
        monomer_idx = randrange(0, self._total_length - 1)
        ## two special cases at the begin/end
        case = 3
        old_energy = self.total_energy
        
        is_occupied = self.box.is_occupied
        is_neighbour = self.box.lattice.is_neighbour

        ##  trials with optimization
        last_affected_idx = 0

        # case
        # 0 - left
        # 1 - right
        # 3 - inside
        #  0 - right
        #  1 - left
        
        #mod_pos = 0
            
        if monomer_idx == 0: ## left
            dir_1 = randrange(0, Lattice.z)
            dir_2 = randrange(0, Lattice.z)
            pos_1 = self.chain[0].neighbours_position[dir_1]
            pos_2 = self.box.lattice.get_coordinate(pos_1, dir_2)
            
            if is_occupied(pos_1) or is_occupied(pos_2):
                return False

            #if not is_neighbour(pos_1, pos_2):
            #    return False

            m1 = self.chain[0]
            m2 = self.chain[1]
           
            open_pos_1 = self.box.lattice.get_open_coordinate(m1._open, dir_1)
            open_pos_2 = self.box.lattice.get_open_coordinate(open_pos_1, dir_2)

            old_position_list = self.positions_list[:]
            #old_monomer_position = [ m._position for m in self.chain ]
            old_monomer_position = old_position_list
            old_monomer_open_position = [ m._open for m in self.chain ]
            old_monomer_next_direction = [ m.next_direction for m in self.chain ]

            m1.position = pos_2
            m1._open = open_pos_2

            m2.position = pos_1
            m2._open = open_pos_1
            
            ## old energy
        
            #self[0] = m1
            m1.chain = self
            self.chain[0] = m1
            self.positions_list[0] = m1._position

            #self[1] = m2
            m2.chain = self
            self.chain[1] = m2
            self.positions_list[1] = m2._position

            #mod_pos = 2
    
            for idx in xrange(2, self._total_length):
                i_monomer = self.chain[idx]
                p_monomer = self.chain[idx - 1]

                diff = np.array(i_monomer._open) - np.array(p_monomer._open)
                diff = diff.dot(diff)

                if diff == 2:
                    p_monomer.next_direction = p_monomer.neighbours_position.index(i_monomer._position)
                    break
                
                #mod_pos += 1

                i_monomer.position = old_monomer_position[idx - 2]
                i_monomer._open = old_monomer_open_position[idx - 2]
                i_monomer.next_direction = old_monomer_next_direction[idx - 2]
                #self[idx] = i_monomer
                i_monomer.chain = self
                self.chain[idx] = i_monomer
                self.positions_list[idx] = i_monomer._position

                last_affected_idx = idx

            
            self.chain[0].next_direction = self.chain[0].neighbours_position.index(self.chain[1]._position)
            self.chain[1].next_direction = self.chain[1].neighbours_position.index(self.chain[2]._position)
            case = 0
            if settings.DEBUG: self.valid_chain()
        elif monomer_idx == self._total_length - 1: ## right
            dir_1 = randrange(0, Lattice.z)
            dir_2 = randrange(0, Lattice.z)
            pos_1 = self.chain[-1].neighbours_position[dir_1]
            pos_2 = self.box.lattice.get_coordinate(pos_1, dir_2)

            if is_occupied(pos_1) or is_occupied(pos_2):
                return False

            open_pos_1 = self.box.lattice.get_open_coordinate(self.chain[-1]._open, dir_1)
            open_pos_2 = self.box.lattice.get_open_coordinate(open_pos_1, dir_2)

            old_position_list = self.positions_list[:]
            #old_monomer_position = [ m._position for m in self.chain ]
            old_monomer_position = old_position_list
            old_monomer_open_position = [ m._open for m in self.chain ]
            old_monomer_next_direction = [ m.next_direction for m in self.chain ]
            
            m1 = self.chain[-1]
            m2 = self.chain[-2]

            m1.position = pos_2
            m1._open = open_pos_2

            m2.position = pos_1
            m2._open = open_pos_1

            #mod_pos = 2

            for idx in xrange(self._total_length - 2, -1, -1):
                i_monomer = self.chain[idx]
                n_monomer = self.chain[idx+1]

                diff = np.array(i_monomer._open) - np.array(n_monomer._open)
                diff = diff.dot(diff)

                if diff == 2:
                    i_monomer.next_direction = i_monomer.neighbours_position.index(n_monomer._position)
                    break
                
                #mod_pos += 1

                i_monomer.position = old_monomer_position[idx + 2]
                i_monomer._open = old_monomer_open_position[idx + 2]
                i_monomer.next_direction = old_monomer_next_direction[idx + 2]

                #self[idx] = i_monomer
                i_monomer.chain = self
                self.chain[idx] = i_monomer
                self.positions_list[idx] = i_monomer._position

                last_affected_idx = idx
            
            case = 1

            self.chain[-2].next_direction = self.chain[-2].neighbours_position.index(self.chain[-1]._position) 
            self.chain[-1].next_direction = 0

            if settings.DEBUG: self.valid_chain()

        if case == 3: 
            prev_monomer = self.chain[monomer_idx - 1]
            monomer = self.chain[monomer_idx]
            next_monomer = self.chain[monomer_idx + 1]

            random_neighbour_prev = randrange(0, Lattice.z)
            #if self.box.is_occupied(prev_monomer.neighbours_position[random_neighbour_prev]):
            #    return False

            #L can be adjenct to i-1 or i+1, randomly choose
            sub_case = randrange(0, 2) ## 0 - right, 1 - left
            
            if sub_case == 0:
                L_position = next_monomer.neighbours_position[random_neighbour_prev]
                L_open_position = Lattice.get_open_coordinate(next_monomer._open, random_neighbour_prev)
                Lnext_direction = self.box.lattice.get_neighbours(L_position).index(next_monomer._position)
            elif sub_case == 1:
                L_position = prev_monomer.neighbours_position[random_neighbour_prev]
                L_open_position = Lattice.get_open_coordinate(prev_monomer._open, random_neighbour_prev)
                Lprev_direction = self.box.lattice.get_neighbours(L_position).index(prev_monomer._position)

            C_position = monomer.neighbours_position[random_neighbour_prev]
            C_open_position = Lattice.get_open_coordinate(monomer._open, random_neighbour_prev)

            diff = np.array(L_position) - np.array(C_position)
            diff = diff.dot(diff)

            if diff != 2:
                return False

            if is_occupied(L_position):
                return False

            if is_occupied(C_position):
                if (sub_case == 0 and  C_position != prev_monomer._position) or (sub_case == 1 and C_position != next_monomer._position):
                    return False


            # old_energy
            #old_energy = self.box.calculate_chain_energy(self)
            
            # direction L->C
            C_nb = self.box.lattice.get_neighbours(C_position)
            CL_direction = C_nb.index(L_position)
            LC_direction = self.box.lattice.TURN_DIRECTION[CL_direction]

            # save old position list
            old_position_list = self.positions_list[:]
            #old_monomer_position = [ m._position for m in self.chain ]
            old_monomer_position = old_position_list
            old_monomer_open_position = [ m._open for m in self.chain ]
            old_monomer_next_direction = [ m.next_direction for m in self.chain ]

            if (sub_case == 0 and prev_monomer._position == C_position) or \
                (sub_case == 1 and next_monomer._position == C_position):
                valid_configuration = True
            

            monomer.position = L_position
            monomer._open = L_open_position
            if sub_case == 0:
                monomer.next_direction = Lnext_direction
            elif sub_case == 1:
                monomer.next_direction = LC_direction

            #self[monomer_idx] = monomer
            monomer.chain = self
            self.chain[monomer_idx] = monomer
            self.positions_list[monomer_idx] = monomer._position

            if sub_case == 0:
                prev_monomer.position = C_position
                prev_monomer._open = C_open_position
                prev_monomer.next_direction = CL_direction
                #self[monomer_idx - 1] = prev_monomer
                prev_monomer.chain = self
                self.chain[monomer_idx - 1] = prev_monomer
                self.positions_list[monomer_idx - 1] = prev_monomer._position

            elif sub_case == 1:
                if not valid_configuration:
                   next_monomer.position = C_position
                   next_monomer._open = C_open_position
                   next_monomer.next_direction = self.box.lattice.TURN_DIRECTION[random_neighbour_prev]
                   #self[monomer_idx + 1] = next_monomer
                   next_monomer.chain = self
                   self.chain[monomer_idx + 1] = next_monomer
                   self.positions_list[monomer_idx + 1] = next_monomer._position

                prev_monomer.next_direction = random_neighbour_prev
                #self[monomer_idx - 1] = prev_monomer
                prev_monomer.chain = self
                self.chain[monomer_idx - 1] = prev_monomer
                self.positions_list[monomer_idx - 1] = prev_monomer._position
            
            #mod_pos = 2

            if not valid_configuration:
                if sub_case == 0: ## right
                    for j in xrange(monomer_idx - 2, -1, -1):
                        j_monomer = self.chain[j]
                        i_monomer = self.chain[j+1]

                        diff = np.array(j_monomer._open) - np.array(i_monomer._open)
                        diff = diff.dot(diff)

                        if diff == 2:
                            j_monomer.next_direction = j_monomer.neighbours_position.index(i_monomer._position)
                            break

                        #mod_pos += 1

                        j_monomer.position = old_monomer_position[j+2]
                        j_monomer._open = old_monomer_open_position[j+2]
                        j_monomer.next_direction = j_monomer.neighbours_position.index(i_monomer._position)
                        #self[j] = j_monomer
                        j_monomer.chain = self
                        self.chain[j] = j_monomer
                        self.positions_list[j] = j_monomer._position

                elif sub_case == 1:
                    for j in xrange(monomer_idx + 2, self._total_length):
                        j_monomer = self.chain[j]
                        i_monomer = self.chain[j-1]

                        diff = np.array(j_monomer._open) - np.array(i_monomer._open)
                        diff = diff.dot(diff)

                        if diff == 2:
                            i_monomer.next_direction = i_monomer.neighbours_position.index(j_monomer._position)
                            break
                        
                        #mod_pos += 1

                        j_monomer.position = old_monomer_position[j-2]
                        j_monomer._open = old_monomer_open_position[j-2]
                        i_monomer.next_direction = i_monomer.neighbours_position.index(j_monomer._position)
                        #self[j] = j_monomer
                        j_monomer.chain = self
                        self.chain[j] = j_monomer
                        self.positions_list[j] = j_monomer._position
                
            if settings.DEBUG: self.valid_chain()
        
        ## end move procedure
        # old_energy, emergency case
        new_energy = self.box.calculate_chain_energy(self)
        dE = new_energy - old_energy

        accept_move = False

        if dE <= 0.0 or self.box.athermal_state:
            accept_move = True
        else:
            accept_move = self.box.accept(dE, None)
        
        if accept_move is True:
            self.total_energy = new_energy
            #self.positions_list = [ m._position for m in self.chain ]

            self.box.global_energy_update([self.idx])
            
            #f = open("mod_pos_%d" % self._total_length, "w+")
            #f.write("%d\n" % mod_pos)
            #f.close()

            #self.valid_chain()
        else: ## revert
            if case == 3:
                if sub_case == 0:
                    monomer_list = xrange(0, monomer_idx+2)
                elif sub_case == 1:
                    monomer_list = xrange(monomer_idx-1, self._total_length)
            elif case == 0:
                monomer_list = xrange(0, self._total_length)
            elif case == 1: ##
                monomer_list = xrange(0, self._total_length)
            #monomer_list = changed_monomer_idx
            for idx in monomer_list:
                m = self.chain[idx]
                old_pos = old_monomer_position[idx]
                old_dir = old_monomer_next_direction[idx]

                if m._position != old_pos or m.next_direction != old_dir:
                    m.position = old_monomer_position[idx]
                    m._open = old_monomer_open_position[idx]
                    m.next_direction = old_monomer_next_direction[idx]
                    #self[idx] = m
                    m.chain = self
                    self.chain[idx] = m
                    self.positions_list[idx] = m._position

            #self.positions_list = old_position_list[:]
            self.total_energy = old_energy
            #self.valid_chain()

        return accept_move

    def get_length(self):
        return self._total_length
    
    def set_length(self, value):
        """
        Prevent change of chain length
        """
        raise Exception("Chain length cannot be modified")
    
    length = property(get_length, set_length)

    def get_chain_structure(self):
        """
        Return chain structure
        """
        if not self._chain_structure:
            group_idx = 0
            last_type = self.chain[0].name
            self._chain_structure = { "%s_%d" % (last_type, group_idx): [0,0] }
            for idx in range(1, self._total_length):
                key = "%s_%d" % (self.chain[idx-1].name, group_idx)
                new_key = "%s_%d" % (self.chain[idx].name, group_idx + 1)
                if self.chain[idx].name != self.chain[idx-1].name:
                    self._chain_structure[key][1] = idx
                    self._chain_structure[new_key] = [idx,0]
                    group_idx += 1
            self._chain_structure[key][1] = self._total_length - 1
        return self._chain_structure
    
    chain_structure = property(get_chain_structure)

    def get_monomer_types(self):
        if self._monomer_types != set():
            self._monomer_types = set([ m.name for m in self.chain ])
        return self._monomer_types

    monomer_types = property(get_monomer_types)

    ########### Calculations ####################
    def get_values_file_header(self):
        if self.monomer_types == set():
            self.monomer_types.update([ m.name for m in self.chain ])
        
        
        return [ "rg_global" ] + map(lambda x: "rg_%s" % x, list(self.monomer_types)) + ['distance', 'total_energy']


    def rg(self, type=None):
        """
        @self self: self of polymer
        Calculate radius of gyration (actually square) for single self
        """
        if type:
            monomer_positions = [ np.array(x._open) for x in self.chain if x.name == type ]
        else:
            monomer_positions = [ np.array(x._open) for x in self.chain ]

        center_of_mass = list(np.average(monomer_positions, axis=0))  #reduce(lambda p1, p2: p1 + p2, monomer_positions) / float(len(monomer_positions))
        ret = sum([ np.dot(pos-center_of_mass, pos-center_of_mass) for pos in monomer_positions ]) / float(len(monomer_positions))
        
        """
        ret = sum([ math.pow(sum(map(lambda x: x**2, map(lambda z,p: z-p, pos, center_of_mass))),2) for pos in self._tmp['monomer_positions'] ])
        for pos in self._tmp['monomer_positions']:
            distance = tools.array_diff(pos, center_of_mass)
            ret += math.pow(sum(map(lambda x: x**2, distance)), 2)
        """
        return float(ret)
    
    
    def center_of_mass(self, monomer_type = None, positions_list = None):
        """
        @string|None monomer_type: Type of monomer, to calculate center of mass for its groups
        Return center of mass
        """
        
        if positions_list:
            monomer_positions = np.array(positions_list)
        else:
            if monomer_type:
                monomer_positions = [ np.array(x._open) for x in self.chain if x.name == monomer_type ]
            else:
                monomer_positions = [ np.array(x._open) for x in self.chain if x.name == monomer_type ]
        
        center_of_mass = list(np.average(monomer_positions, axis=0))  #reduce(lambda p1, p2: p1 + p2, monomer_positions) / float(len(monomer_positions))
        #ret = map(lambda x: float(x) / self._total_length, reduce(lambda p1, p2: map(lambda x,y: x+y, p1, p2), monomer_positions))

        return center_of_mass

    def calculate_center_of_mass(self):
        """
        Return center of mass for global and of each of groups in multiblock copolymer
        """
        monomer_groups = self.chain_structure
        ret = [self.center_of_mass()]
        for group, v in monomer_groups.iteritems():
            positions = self.positions_list[v[0]:v[1]]
            cmm = self.center_of_mass(positions_list = positions)
            ret += [cmm]

        return ret

    
    def end_to_end_distance(self):
        """
        Return normalized distance between two vectors
        """
        first_monomer = np.array(self.chain[0]._open)
        last_monomer = np.array(self.chain[-1]._open)
        
        distance = last_monomer - first_monomer
        distance = distance.dot(distance)

        return distance

    def end_to_end_vector(self):
        """
        Return end-to-end vector
        """
        first_monomer = np.array(self.chain[0]._open)
        last_monomer = np.array(self.chain[-1]._open)

        return last_monomer - first_monomer
    
    def calculate_values(self):
        """
        Save results to file, increment counter
        """

        self.valid_chain()
        
        if self.monomer_types == set():
            self._monomer_types = set([ m.name for m in self.chain ])

        rg = [self.rg()] # + [ self.rg(t) for t in self.monomer_types ]

        #distance = self.end_to_end_distance()
        #self.total_energy = self.box.calculate_chain_energy(self)
        
        self.last_values = rg + [self.total_energy]

        if self.min_energy > self.total_energy:
            self.min_energy = self.total_energy

        if self.min_rg > rg[0]:
            self.min_rg = rg[0]

        return True
    
    def calculate_rdf(self, monomer_type=None):
        """
        Calculate Radial Distribution Function in the range of 0->box_size + 1 from the center of mass
        Histogram
        """
        maxbin = 1000
        dr = float(min(self.box.box_size)/float(maxbin)) ## bin width
        hist = [0] * (maxbin+1)
        
        monomer_positions = [ np.array(m._open) for m in self.chain if monomer_type == None or m.name == monomer_type ]
        cmm = reduce(lambda p1, p2: p1 + p2, monomer_positions) / float(self._total_length) ## center of mass
        
        diff = [ np.math.sqrt((m - cmm).dot(m - cmm)) for m in monomer_positions ]
        
        bins = filter(lambda x: x <= maxbin, [ int(np.ceil(x/dr)) for x in diff ])
        hist = [ bins.count(x) for x in range(maxbin) ]

        return hist 

        ## calculate for each of monomer_position
        """
        for particle in self.chain:
            center_of_mass = np.array(particle._open)
            monomer_positions = [ x - center_of_mass for x in chain]
            
            square = [ sum(map(lambda x: x**2, x)) for x in monomer_positions ]

            for i in len(shells):
                count = len(filter(lambda x: x >= i - 0.25 and x < x + 0.25, square))
                shells[i] += count
        ##
        # http://www.physics.emory.edu/~weeks/idl/gofr2.html
        shells /= len(chain)
        shells = [ shells[i] / (4*np.pi* (i**2)*0.5) for i in range(len(shells)) ]
        """
    def calculate_rdf2(self, monomer_type=None):
        """
        Calculate Radial Distribution Function in the range of 0->box_size + 1 from the center of mass
        """
        monomer_positions = [ np.array(m._open) for m in self.chain if monomer_type == None or m.name == monomer_type ]
        cmm = reduce(lambda p1, p2: p1 + p2, monomer_positions) / float(self._total_length) ## center of mass
        
        diff = [ np.math.sqrt((m - cmm).dot(m - cmm)) for m in monomer_positions ]
        
        return diff

    def calculate_rdf_mean(self, monomer_type=None):
        """
        Calculate average 
        """
        monomer_positions = [ np.array(m._open) for m in self.chain if monomer_type == None or m.name == monomer_type ]
        
        ### so we calculate diff between each pair of monomers, not only difference between cmm and monomers
        monomer_positions_product = list(itertools.product(monomer_positions, monomer_positions))
        diff = map(lambda x: x[0].dot(x[1]), monomer_positions_product)

        return diff
    

    def __getitem__(self, key):
        """
        Return monomer
        """
        if isinstance(key, tuple):
            return self.chain[self.positions_list.index(key)]

        return self.chain[key]

    def __setitem__(self, key, value):
        """
        Set monomer in chain
        """
        value.chain = self
        self.chain[key] = value
        self.positions_list[key] = value._position
