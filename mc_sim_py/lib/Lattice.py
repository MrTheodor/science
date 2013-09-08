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
# import bigfloat
from math import exp

import settings
import tools
import numpy

import time

class Lattice(object):
    """
    FCC lattice
    """

    DIRECTION = [
              (-1, 0, -1),
              (1, 0, 1),
              (-1, 0, 1),
              (1, 0, -1),
              (0, 1, -1),
              (0, -1, 1),
              (0, 1, 1),
              (0, -1, -1),
              (-1, 1, 0),
              (1, -1, 0),
              (-1, -1, 0),
              (1, 1, 0)
    ]
    XZ_DIRECTION = [ 0, 1, 2, 3 ]
    YZ_DIRECTION = [ 4, 5, 6, 7 ]
    XY_DIRECTION = [ 8, 9, 10, 11 ]

    # # get index of turn direction
    # # 0 -> 1
    # # 1 -> 0
    TURN_DIRECTION = [1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10]

    # # cache
    NEIGHBOURS_POSITION = {}
    NEIGHBOURS_OPEN_POSITION = {}
    ROTATION_WAYS = {}

    box_size = (0, 0, 0)
    coordinate = None
    number_of_points = 0

    z = 12  # number of neighbours
    neighbour_distance = 2

    def __init__(self, size):
        self.box_size = size
        self.number_of_points = (size[0] * size[1] * size[2]) / 2

    @classmethod
    def is_valid_coordinate(cls, c):
        return ((sum(c) % 2) == 0)

    @classmethod
    def get_neighbours(cls, position):
        """
        @position: tuple
        @type: int - 1 -- bc positions, 2 - open positions

        Return list of position of neighbours
        """

        nb_open_positions = cls.DIRECTION + np.array(position)
        nb_positions = nb_open_positions % cls.box_size
        nb_positions = [ tuple(x) for x in nb_positions ]

        return nb_positions

    @classmethod
    def get_neighbours_open(cls, position):
        nb_open_positions = cls.DIRECTION + np.array(position)
        nb_open_positions = [ tuple(x) for x in nb_open_positions ]


        return nb_open_positions

    @classmethod
    def is_neighbour(cls, a, b):
        """
        Is position a neighbour of position b
        """
        # nb_list = cls.get_neighbours(b)
        # return a in nb_list
        diff = np.array(a) - np.array(b)
        # diff = list(diff)
        # if diff.count(1) == 2 and diff.count(0) == 1:
        #    return True
        diff = diff.dot(diff)
        if diff == cls.neighbour_distance:
            return True
        return False


    @classmethod
    def get_coordinate(cls, current_position, direction):
        """
        Return new coordinate, from current_position in some direction
        Applied periodic boundary conditions

        @tuple current_position: tuple - current position
        @int direction: - direction to move

        DIRECTION [x][y][z]
         0: -1  0 -1
         1:  1  0  1
         2: -1  0  1
         3:  1  0 -1
         4:  0  1 -1
         5:  0 -1  1
         6:  0  1  1
         7:  0 -1 -1
         8: -1  1  0
         9:  1 -1  0
        10: -1 -1  0
        11:  1  1  0
        """

        bs = cls.box_size  # # size of box
        coordinate = tuple((np.array(current_position) + np.array(cls.DIRECTION[direction])) % bs)
        return coordinate

    @classmethod
    def get_open_coordinate(cls, current_position, direction):
        coordinate = tuple(np.array(current_position) + np.array(cls.DIRECTION[direction]))

        return coordinate

    @classmethod
    def convert_to_open_coordinate(cls, prev_position, converted_position):
        """
        @tuple prev_position:
        @tuple converted_position:
        Convert periodic bounduary condition to open one
        """
        difference = tools.array_diff(converted_position, prev_position)
        # distance = tuple(map(lambda x: cmp(x, 0), difference))
        # distance_length =  tools.dot(distance,distance)
        # if distance_length > 2:
        #    import ipdb;ipdb.set_trace() #============ BREAKPOINT ============#
        direction = cls.get_direction(tuple(map(lambda x: cmp(x, 0), difference)))
        return cls.get_open_coordinate(prev_position, direction)

    @classmethod
    def convert_to_periodic_coordinate(cls, position):
        """
        Convert coordinate to periodic coordinate
        """
        print "box_size", cls.box_size
        return tuple(map(lambda x, box_size: x % box_size, position, cls.box_size))

    @classmethod
    def get_direction(cls, vector):
        """
        Return direction of given vector or exception
        """
        return cls.DIRECTION.index(tuple(map(lambda x: cmp(x, 0), vector)))


    @classmethod
    def get_square_distance(cls, monomer_a, monomer_b):
        """
        Return distance^2 between monomers.
        """
        distance = np.array(monomer_b) - np.array(monomer_a)
        distance = distance.dot(distance)

        return distance
