'''
Created on 1 Nov 2011
@author: teodor
'''

import itertools
import logging
import math
import random
from mpi4py import MPI
import numpy as np
#import bigfloat
from math import exp

import settings
import tools
import numpy

import time

from lib.Polymer import Chain
from lib.Lattice import Lattice

class Box(object):

    box_size = (0,0,0)
    chain_list = []

    fixed_objects = []

    number_of_chains = 0
    lattice = None
    total_energy = 0.0

    athermal_state = False

    temperature_history = []

    acceptation_count = 0

    ## statistic information
    _localT = -1.0

    def __init__(self, box_size, localT = 0.0):
        """
        @tuple box_size
        @float localT

        Main class for Box

        """
        self.box_size = box_size
        self.chain_list = []
        self.disabled_chains = None
        self.localT = localT

        # global set for Lattice argh!
        Lattice.box_size = box_size
        self.lattice = Lattice(box_size)

        if self.box_size == (0,0,0) or not self.lattice:
            raise Exception("Problem with BOX initializing")

    def get_localT(self):
        return self._localT

    def set_localT(self, value):
        self._localT = value
        ## do some other fancy stuff
        Box._locaT = value
        self.temperature_history.append(value) ## move to chain
        for chain in self.chain_list:
            chain.temperature_history.append(value)

        ## load new value of exp
        if value in settings.EXP_TABLES:
            settings.EXP_TABLE = settings.EXP_TABLES[value]
            logging.info('change exp_table %f ' % value)
        else:
            settings.EXP_TABLE = tools.load_file(settings.ROOT_DIR +  \
                settings.FILES_TEMPLATE['exp'] % value, {})
            settings.EXP_TABLES[value] = settings.EXP_TABLE
            logging.info('load new exp_table %f' % value)

        logging.info("Set temperature to: %f" % value)


    localT = property(get_localT, set_localT)

    #########################################################################

    def random_strategy(self, chain):
        """
        Put chain in random mode
        """
        random_point = tuple(map(lambda x: x/2, self.box_size))
        while self.is_occupied(random_point) or not Lattice.is_valid_coordinate(random_point):
            random_point = tuple([ random.randrange(0, self.box_size[x]) for x in range(3) ])

        if __debug__: logging.info("Chain start point: %d, %d, %d" % random_point)

        f_monomer = chain.chain[0]
        f_monomer.position = random_point
        f_monomer._open = random_point
        chain[0] = f_monomer

        for idx in range(1, chain._total_length):
            prev_monomer = chain.chain[idx - 1]
            monomer = chain.chain[idx]

            not_found = True
            direction_list = range(0, Lattice.z)
            direction_list_idx = Lattice.z
            while not_found:
                try:
                    direction_idx = random.randrange(0, direction_list_idx)
                except:
                    return False
                direction = direction_list[direction_idx]

                pos = Lattice.get_coordinate(prev_monomer._position, direction)
                open_pos = Lattice.get_open_coordinate(prev_monomer._open, direction)

                if pos in chain.positions_list or  \
                    (self.number_of_chains > 1 and self.is_occupied(pos)):

                    not_found = True
                    direction_list.remove(direction)
                    direction_list_idx -= 1
                else:
                    not_found = False

            monomer.position = pos
            monomer.open_position = open_pos
            #monomer.next_direction = None
            prev_monomer.next_monomer = monomer
            prev_monomer.next_direction = direction

            chain[idx] = monomer
            chain.chain[idx-1] = prev_monomer

        ## rebuild directions
        for idx in range(1, chain._total_length):
            p = chain.chain[idx - 1]
            m = chain.chain[idx]
            d = p.neighbours_position.index(m._position)
            p.next_direction = d
            chain.chain[idx - 1].next_direction = d

        chain.chain[-1].next_direction = 0

        if __debug__: logging.info('Chain inside box, length=%d' % chain.length)
        return chain

    def greedy_strategy(self, chain):
        """
        Put chain in the box based on a greedy strategy
        """
        random_point = tuple(map(lambda x: x/2, self.box_size))
        while self.is_occupied(random_point) or not Lattice.is_valid_coordinate(random_point):
            random_point = tuple([ random.randrange(0, self.box_size[x]) for x in range(3) ])

        f_monomer0 = chain.chain[0]
        f_monomer0.position = random_point
        f_monomer0._open = random_point
        chain[0] = f_monomer0

        next_direction = random.randrange(0, Lattice.z)

        f_monomer1 = chain.chain[1]
        f_monomer1.position = Lattice.get_coordinate(random_point, next_direction)
        f_monomer1._open = Lattice.get_open_coordinate(random_point, next_direction)
        chain[1] = f_monomer1
        chain[0].next_direction = next_direction

        pos_list = [f_monomer0._position, f_monomer1._position]

        direction_list = range(0, Lattice.z)
        for idx in range(2, chain._total_length):
            prev_monomer = chain.chain[idx - 1]
            monomer = chain.chain[idx]
            min_energy = self.calculate_chain_energy(chain)
            min_pos = None
            min_open_pos = None
            min_dir = 0

            np.random.shuffle(direction_list)

            last_valid_pos = None
            last_valid_open_pos = None

            for dir in direction_list:
                pos = Lattice.get_coordinate(prev_monomer._position, dir)
                open_pos = Lattice.get_open_coordinate(prev_monomer._open, dir)
                if pos in chain.positions_list or \
                    (self.number_of_chains > 1 and self.is_occupied(pos)) or \
                    pos in pos_list:
                    pass
                else:
                    last_valid_pos = pos
                    last_valid_open_pos = open_pos
                    monomer.position = pos
                    monomer._open = open_pos
                    chain[idx] = monomer
                    m_e = self.calculate_chain_energy(chain, slice=[0, idx+1])
                    if min_energy > m_e:
                        min_energy = m_e
                        min_pos = pos
                        min_open_pos = open_pos
                        min_dir = dir

            monomer.position = min_pos if min_pos else last_valid_pos
            monomer._open = min_open_pos if min_open_pos else last_valid_open_pos
            chain[idx] = monomer

            prev_monomer.next_direction = min_dir
            chain[idx - 1] = prev_monomer

        ## rebuild directions
        for idx in range(1, chain._total_length):
            p = chain.chain[idx - 1]
            m = chain.chain[idx]
            d = p.neighbours_position.index(m._position)
            p.next_direction = d
            chain.chain[idx - 1].next_direction = d

        if __debug__: logging.info("Greedy strategy, chain with min energy: %f" % \
                self.calculate_chain_energy(chain))
        return chain


    def add_chain(self, chain, new = True):
        """
        @Polymer chain: chain object
        @boolean new: is this a new chain?
        Add chain to the box
        """
        #random_point = self.lattice.get_random_point()
        chain.box = self

        chain.idx = self.number_of_chains

        self.chain_list.append(chain)
        self.number_of_chains += 1

        if new == True:
            """
            New virgin chain, need to update position of every monomer
            """
            chain_inside = False
            while not chain_inside:
                return_chain = self.random_strategy(chain)
                if return_chain is False:
                    chain_inside = False
                    print "Put again"
                else:
                    chain_inside = True
                    chain = return_chain
            #chain = self.greedy_strategy(chain)

        """
        Only add chain to the chain_list, without update monomer positions
        """
        chain.valid_chain()
        ## update energy of every chain in the box, also every monomer has own energy
        self.global_energy_update()

        print "Chain in box with energy", chain.total_energy

        return len(self.chain_list) - 1

    def rebuild_direction(self, chain):
        ## rebuild directions
        for idx in range(1, chain._total_length):
            p = chain.chain[idx - 1]
            m = chain.chain[idx]
            d = p.neighbours_position.index(m._position)
            p.next_direction = d
            chain.chain[idx - 1].next_direction = d

        chain.chain[-1].next_direction = 0
        return chain


    def refresh_chains(self):
        """
        Refresh chain positions_list and direction list
        """
        for chain in self.chain_list:
            monomer0 = chain.chain[0]
            for idx in range(1, chain._total_length):
                pass


    def remove_chain(self, chain):
        chain_idx = tools.index(self.chain_list, chain)

        del self.chain_list[chain_idx]
        self.number_of_chains -= 1

    def accept(self, energy_change, new_position = None):
        """
        @float energy_change: change of energy
        @tuple new_position: new position of monomer

        Check if change of energy is acceptable (according to current temperature)
        Accept if energy_change is below
        """

        localT = float(self._localT)
        return_value = False

        ## first check if momoner doesn't overlap
        if new_position != None:
            if self.is_occupied(new_position):
                return False

        if localT <= 0.0 or self.athermal_state is True:
            return_value = True
        elif energy_change <= 0.0:
            return_value = True
        else:
            ## cache of exp function
            #logging.debug("Energy change: %f" % (energy_change))
            #if settings.EXP_TABLE.has_key(energy_change):
            #    exp_value = settings.EXP_TABLE[energy_change]
            #else:

            if settings.EXP_TABLE.has_key(localT):
                if settings.EXP_TABLE[localT].has_key(energy_change):
                    exp_value = settings.EXP_TABLE[localT][energy_change]
                else:
                    exp_value = exp((-1*energy_change) / localT)
                    settings.EXP_TABLE[localT][energy_change] = exp_value
            else:
                exp_value = exp((-1*energy_change) / localT)
                settings.EXP_TABLE[localT] = { energy_change: exp_value}

            #    settings.EXP_TABLE[energy_change] = exp_value
            ## monte carlo acceptance
            r = numpy.random.random()
            if r < exp_value:
                return_value = True

        if return_value == True:
            self.acceptation_count += 1

        return return_value

    def is_occupied(self, position, exclude_chain = None):
        """
        @tuple position: position to check
        @Polymer exclude_chain: chain which should not be checked

        Check if given :position: is occupied more than once by some monomers
        """

        if position is None:
            return False

        if isinstance(position, set):
            for chain in self.chain_list:
                if chain is exclude_chain:
                    continue

                if position.issubset(set(chain.positions_list)):
                    return True
        else:
            for chain in self.chain_list:
                if chain is exclude_chain:
                    continue
                if position in chain.positions_list:
                    return True

        return False

    def global_energy_update(self, exclude_chains = []):
        """
        Iterate over every chain in box and update its energy
        """

        self.total_energy = 0.0

        for chain_idx in range(self.number_of_chains):
            if chain_idx in exclude_chains:
                continue
            chain_energy = self.calculate_chain_energy(self.chain_list[chain_idx])
            self.chain_list[chain_idx].total_energy = chain_energy
            self.total_energy += chain_energy

        ## we also update energy of box, why not

    def update_chain_energy(self, chain):
        """
        @Polymer chain: chain of monomers

        Iterate over chain and update energy of every monomer in chain
        """

        chain.total_energy = self.calculate_chain_energy(chain) # / 2.0 ## cause we calculate energies twice

        return chain

    def calculate_chain_energy(self, chain, slice=[]):
        """
        @Chain chain: chain of monomers
        @list slice: [start, stop] slice of polymer chain

        """

        total_energy = 0.0

        if slice != []:
            start = slice[0]
            stop = slice[1]
        else:
            start = 0
            stop = chain._total_length

        for monomer_idx in xrange(start,stop):
            energy = self.calculate_monomer_energy(chain.chain[monomer_idx])
            chain.chain[monomer_idx].energy = energy
            total_energy += energy

        chain.total_energy = total_energy

        return total_energy

    def calculate_monomer_energy(self, monomer, exclude_positions = None):
        """
        @Monomer monomer: monomer object
        @tuple|list exclude_positions: list of exclude positions

        Calculate energy around :monomer on position :monomer
        So we take into account :monomers around this one and use interaction_table to
        get total energy.

        Exclude positions, if we want to exclude from list of interacted monomres, some positions
        that we now they should not be consider as an occupied
        Eg. in movement, we should not take into account old position of :monomer if we
        try to calculate new energy
        """

        # find monomer on this position
        # if monomer is an Monomer object, we don't search for it
        if isinstance(exclude_positions, tuple):
            exclude_positions = [exclude_positions]

        if exclude_positions:
            neighbours = [ nb for nb in monomer.neighbours_position if nb not in exclude_positions ]
        else:
            neighbours = monomer.neighbours_position ## save positions of neighbours, but we don't know if on each of them there is a monomer

        neighbours = set(neighbours)
        solvent_positions = Lattice.z
        neighbours_interaction_energy = []

        for chain in self.chain_list:
            #monomers_idx = filter(None, map(lambda x: chain.positions_list.index(x) if x in neighbours else None, chain.positions_list))
            #neighbours_interaction_energy += map(lambda x: chain.chain[x].interaction_table[chain.chain[x].name], monomers_idx)
            check_positions = set(chain.positions_list).intersection(neighbours)
            neighbours_interaction_energy += [ monomer.interaction_table[x.name] for x in
                                            chain.chain if x._position in check_positions
                                            ]
            solvent_positions -= len(neighbours_interaction_energy)

        ## calculate energy
        E = 0.0
        E += sum(neighbours_interaction_energy)

        if solvent_positions > 0 and monomer.interaction_table.has_key('_s'):
            if monomer.interaction_table['_s'] != 0.0:
                E += solvent_positions*monomer.interaction_table['_s']

        return E

    def get_box(self):
        """
        Return list of chains in box, valid for send via MPI
        """
        ret = []
        for c in self.chain_list:
            pos = [ m._position for m in c.chain ]
            open_pos = [ m._open for m in c.chain ]
            ret.append( (pos, open_pos))

        return ret

    def update_box(self, conformations):
        """
        """
        for idx in range(len(conformations)):
            pos, open_pos = conformations[idx]
            chain = self.chain_list[idx]
            for i in range(chain._total_length):
                chain[i].position = pos[i]
                chain[i]._open = open_pos[i]

            chain.positions_list = pos

            chain = self.rebuild_direction(chain)
            chain.valid_chain()

            self.chain_list[idx] = chain

    def swap(self):
        """
        Swap box according to PT algorithm
        """
        if not settings.PT:
            return False

        pt = settings.PT_MODULE

        pt_rank = pt.rank
        pt.switch_count += 1

        pt.comm.Barrier()

        if pt_rank == 0: ## only rank = 0, master node
            time_start = time.time()
            ## get temperature and energy from every child nodes
            swapT = np.zeros(pt.size)
            swapE = np.zeros(pt.size)

            self.global_energy_update()

            swapT[0] = self._localT
            swapE[0] = self.total_energy

            # swap_table
            index = 1

            # message format:
            # rank, step, localT, energy
            while index < pt.size:
                status = MPI.Status()
                val = pt.recv(rank = MPI.ANY_SOURCE, tag = settings.TAG, status = status)
                source = status.Get_source()
                if len(val) == 3 and source != 0:
                    swapT[source] = val[1]
                    swapE[source] = val[2]
                    index += 1

            tmp = [0]*pt.size
            for i in range(pt.size):
                for j in range(pt.size):
                    if pt.temp_list[i] == swapT[j] or \
                        ((pt.temp_list[i] >= swapT[j] - 0.00001) and \
                         (pt.temp_list[i] <= swapT[j] + 0.00001)):

                        tmp[i] = j
                        break

            ## swap random configuration
            rank_vector = np.random.permutation(pt.size)
            for rank_vector_idx in range(1, pt.size):
                swap = False
                selected_rank = tmp[rank_vector[rank_vector_idx]]

                if swapE[selected_rank] <= swapE[selected_rank-1]:
                    swap = True

                    if __debug__: logging.debug('energy[%d]:%f <= energy[%d]:%d' % (selected_rank,
                                    swapE[selected_rank],
                                    selected_rank-1,
                                    swapE[selected_rank-1]))

                else:
                    r = np.random.random()

                    value = ( (1.0/swapT[selected_rank - 1]) - (1.0/swapT[selected_rank]))
                    value *= (swapE[selected_rank] - swapE[selected_rank - 1])

                    exp_value = math.exp(-1*value)

                    if r < exp_value:
                        swap = True

                if swap:
                     tmpT = swapT[selected_rank]

                     if __debug__: logging.info("%d: Swap temp, to %f from %d on step %d" % (settings.RANK,
                                    swapT[selected_rank - 1],
                                    selected_rank, settings.CURRENT_STEP))

                     swapT[selected_rank] = swapT[selected_rank-1]
                     swapT[selected_rank-1] = tmpT

                        ## FOPT
                     pt.switch_histogram[rank_vector_idx - 1] += 1


                for ii in range(pt.size):
                    for jj in range(pt.size):
                        if pt.temp_list[ii] == swapT[jj] or \
                            ((pt.temp_list[ii] >= swapT[jj] - 0.00001) and \
                            (pt.temp_list[ii] <= swapT[jj] + 0.00001)):

                            tmp[ii] = jj
                            break

            ## FOPT, swaped
            tmpUpDown = 0
            for i in range(pt.size):
                tmpUpDown = pt.up_down_tab[i]
                if (swapT[i] == pt.temp_list[pt.size - 1]) or \
                    ((swapT[i] <= pt.temp_list[pt.size - 1] + 0.00001) and \
                    (swapT[i] >= pt.temp_list[pt.size - 1] - 0.00001)):
                    pt.up_down_tab[i] = 2
                if (swapT[i] == pt.temp_list[0]) or \
                    ((swapT[i] <= pt.temp_list[0] + 0.00001) and \
                    (swapT[i] >= pt.temp_list[0] - 0.00001)):
                    pt.up_down_tab[i] = 1

                if pt.up_down_tab[i] !=0 and tmpUpDown != pt.up_down_tab[i]:
                    pt.round_trip[i] += 1

            for i in range(pt.size):
                if pt.up_down_tab[i] == 1:
                    pt.up_histogram[i] += 1
                if pt.up_down_tab[i] == 2:
                    pt.down_histogram[i] += 1

            # predict if any of replicas is at the end of simulation
            # if CURRENT_STEP + FREQ_SWAP >= settings.MC_MODULE.stop_step
            #tmp_step = swapStep[:]
            #tmp_step += settings.FREQ_SWAP
            #check = filter(lambda x: x < settings.MC_MODULE.stop_step, tmp_step)
            #if len(check) != pt.size:
            #    max_current_step = max(swapStep)
            #    current_step = 1
            #    logging.info("Some of nodes are too fast...")
            #    settings.MC_MODULE.set_stop_simulation()
            #else:
            #    current_step = -1

            # synchronization current_step between replicas, send max value

            ## send new values
            for i in range(1, pt.size):
                data = (pt.rank, swapT[i], -1)
                pt.send(data, rank = i, tag = settings.TAG)

            if swapT[0] != 0.0:
                self.localT = swapT[0]
                settings.LOCAL_T = swapT[0]
            else:
                pass
            time_delta = time.time() - time_start
            if __debug__: logging.info("Swap time: %f" % time_delta)
        else:
            self.global_energy_update()
            ## send data to 0-node
            data = (pt.rank, self._localT, self.total_energy)
            pt.send(data, rank = 0, tag = settings.TAG) ## send

            ## get new value of temperature from 0-node
            status = MPI.Status()
            send_by, temperature, energy = pt.recv(rank = 0, tag = settings.TAG, status = status)

            if temperature != 0.0:
                self.localT = temperature
                settings.LOCAL_T = temperature

                if __debug__:logging.debug('%d: Set temperature to: %f from %d' % (settings.RANK,
                                temperature,
                                send_by))

                #self.global_energy_update()
            else:
                #logging.debug('temperature == 0.0!!')
                pass

        ## temperature histogram
        for i in range(pt.size):
            if (self._localT == pt.temp_list[i]) or \
                ((self._localT <= pt.temp_list[i] + 0.00001) and \
                 (self._localT >= pt.temp_list[i] - 0.00001)):
                pt.temperature_histogram[i] += 1
        ###

    def jmol_snapshot(self, file_name):
        """
        Return current state of box in xyz file format
        """
        f = open(file_name, "w+")

        coordinates = []

        for chain in self.chain_list:
            for monomer in chain.chain:
                coordinates.append([monomer.name] + list(np.array(monomer._open)*settings.JMOL_SCALE))

        f.writelines(["%d\n\n" % len(coordinates)])
        f.writelines([ "%s %s\n" % (d[0], " ".join(map(str, d[1:]))) for d in coordinates ])
        f.close()

    def jmol_snapshot_hdf(self, hdfTbl, hdfSet, set_name):
        coordinates = []
        for chain in self.chain_list:
            idx = 0
            for m in chain.chain:
                coordinates.append([idx] + list(m._open))
                idx+=1

        scipy_array = coordinates
        hdfTbl.createArray(hdfSet, set_name, scipy_array, set_name)

    def jmol_snapshot_string(self):
        """
        Return current state of box in xyz file format
        """

        coordinates = []

        for chain in self.chain_list:
            for monomer in chain.chain:
                coordinates.append([monomer.name ] + list(monomer._open))

        return coordinates


    def load_jmol_conf(self, file_name, number_of_chains=1):
        """
        Load chain configuration from JMOL file

        @string file_name: Name of jmol file
        @dict map_dict: dict with Monomer class
        @int number_of_chains: number of chains in jmol file

        """

        f = open(file_name, "r+")
        data = f.readlines()

        header = int(data[0])
        coordinate_data = map(lambda x: x.split(" "), data[2:])

        chains = []
        for chain_idx in range(number_of_chains):
            chain_seq = []
            slice_coordinate = coordinate_data[chain_idx * header:]
            self.chain_list[chain_idx].positions_list = []
            for c_idx in range(len(slice_coordinate)):
                c = slice_coordinate[c_idx]
                crr = tuple(map(lambda x: float(x)/settings.JMOL_SCALE, tuple(c[1:])))
                print "Load", crr
                self.chain_list[chain_idx].chain[c_idx]._open = crr
                self.chain_list[chain_idx].chain[c_idx].position = Lattice.convert_to_periodic_coordinate(crr)
                self.chain_list[chain_idx].positions_list.append(self.chain_list[chain_idx].chain[c_idx]._position)


            self.chain_list[chain_idx] = self.rebuild_direction(self.chain_list[chain_idx])
