'''
@author: teodor
'''

import tools
import random
import time
from mpi4py import MPI
import numpy as np
import settings
import logging
import numpy as np

import os

class PT(object):
    comm = None
    rank = 0
    size = 0
    name = None

    min_temp = 0.0
    max_temp = 0.0

    number_of_nodes = 1
    temp_list = []
    temp_delta = 0.0
    _stepT = 0.0

    current_temperature_range = None
    down_histogram = []
    up_histogram = []
    up_down_tab = []
    round_trip = []
    switch_histogram = []
    switch_count = 0

    ##
    world_list = {}

    def __init__(self, temp_list = None):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        self.down_histogram = np.zeros(self.size)
        self.up_histogram = np.zeros(self.size)
        self.up_down_tab = np.zeros(self.size)
        self.temperature_histogram = np.zeros(self.size)
        self.round_trip = np.zeros(self.size)
        self.switch_histogram = np.zeros(self.size)
        self.switch_count = 0

        ## initialize random generator
        ## each of replics need own seed
        if not settings.DEBUG and not settings.PROFILER:
            new_seed = time.time() * (self.rank if self.rank > 1 else 1) + self.rank
            np.random.seed(int(new_seed))
            random.seed(int(new_seed))

        ## build linear temp list
        if self.size > 1:
            self.temp_delta = (self.max_temp - self.min_temp) / float(self.size)
            if self.temp_delta == 0.0:
                self.temp_delta = 1.0/float(self.size)
        else:
            self.temp_delta = 0.1

        # build linear temp list
        if temp_list:
            self.temp_list = temp_list
        else:
            self.load_temp_table()

        if __debug__: logging.info("Temp delta: %f" % self.temp_delta)
        if __debug__: logging.info("Temp list: [%s]" % ",".join(map(str, self.temp_list)))

        settings.RANK = self.rank
        settings.SIZE = self.size
        settings.SIZE_VECTOR = range(self.size)

        self.world_list = dict(zip(range(self.size), range(self.size)))
        if __debug__: logging.info('world [%s]' % ",".join(map(str, self.world_list)))

    def get_temperature(self, rank = None):
        if not rank:
            rank = self.rank

        if len(self.temp_list) == 0 and (self.min_temp - self.max_temp) == 0:
            return self.min_temp
        temp_range = self.temp_list[rank]

        return temp_range

    def load_temp_table(self):
        temp_file = settings.ROOT_DIR + "/%s_%s_temperatures_%d.cfg" % (settings.EXPERIMENT_ROOT_NAME, settings.EXPERIMENT_NAME, self.size)
        self._stepT = 0
        try:
            temp_lines = open(temp_file, 'r').readlines()
            if temp_lines == []:
                raise Exception()

            self._stepT = len(temp_lines)
            self.temp_list = map(float, temp_lines[-1].split(';'))
        except:
            if self.size > 1:
                delta_temp = self.max_temp - self.min_temp
                delta_temp /= float((self.size - 1))
            else:
                delta_temp = 0.1

            self.temp_list = [ self.min_temp + i*delta_temp for i in range(self.size) ]
            self._stepT = 1

            if self.rank == 0:
                temp_cfg = open(temp_file, "w+")
                logging.info("Save init temp table to file %s" % temp_cfg)
                temp_cfg.writelines(";".join(map(str, self.temp_list)))
                temp_cfg.writelines("\n")
                temp_cfg.close()

        logging.info("Read temp_list: %s" % str(self.temp_list))
        logging.info("StepT: %d" % self._stepT)

    def calculate_round_trip(self):
        correct = True
        if self.rank == 0:

            for i in range(self.size):
                if self.round_trip[i] < (self._stepT + 2):
                    logging.info('correct _false %s' % str(self.round_trip))
                    correct = False
                    break

            if correct:
                tmpUpDown = np.zeros(self.size)
                tmpTemp = np.zeros(self.size)

                for i in range(self.size):
                    tmpUpDown[i] = \
                        self.up_histogram[i]/(self.up_histogram[i]+ self.down_histogram[i])

                    tmpTemp[i] = self.temp_list[i]

                tmpTemp = self.generate_new_temp_table(tmpTemp, tmpUpDown)
                temp_file = settings.ROOT_DIR + "/%s_%s_temperatures_%d.cfg" % (settings.EXPERIMENT_ROOT_NAME, settings.EXPERIMENT_NAME, self.size)

                temp_file_f = open(temp_file, 'a')

                l = ";".join(map(str, tmpTemp))
                logging.info("Round trip: %s" % l)
                temp_file_f.writelines(l)
                temp_file_f.writelines("\n")
                temp_file_f.close()

                for i in range(1, self.size):
                    self.send(settings.MC_MODULE.stop_step, rank = i, tag = settings.TAG+5)

            else:
                # need to do more simulation,
                # change of stop_step on MC object,
                # also send this new value to other replics
                new_mc_stop_step = settings.MC_MODULE.stop_step + settings.MC_STEPS
                settings.MC_MODULE.change_stop_step(new_mc_stop_step)
                for i in range(1, self.size):
                    self.send(new_mc_stop_step, rank = i, tag = settings.TAG+5)
        else:
            status = MPI.Status()
            new_stop_mc = self.recv(0, tag = settings.TAG+5, status = status)
            logging.info("%d: Set new stop mc value to: %d" % (settings.RANK,
                            new_stop_mc))
            settings.MC_MODULE.change_stop_step(new_stop_mc)

    def save_histograms(self):
        """
        Save histograms that are used for calculate FOPT
        """
        if self.rank == 0:
            output_switch = "%sSwitch_%d.hst" % (settings.DIRS['calc'], self._stepT)
            output_updown = "%sUpDown_%d.hst" % (settings.DIRS['calc'], self._stepT)
            output_round_trip = "%sRoundTrip_%d.hst" % (settings.DIRS['calc'], self._stepT)

            output_switch = open(output_switch, "w+")
            output_updown_f = open(output_updown, "w+")
            output_round_trip_f = open(output_round_trip, "w+")
            output_files = [output_switch, output_updown_f, output_round_trip_f]

            tmpTrip = 0.0

            for i in range(self.size):
                if i < self.size-1:
                    l = "%d;%f;%f\n" % (i, self.temp_list[i], self.switch_histogram[i]/float(self.switch_count))
                    output_switch.writelines(l)

                l = "%d;%f;%f\n" % (i, self.temp_list[i], self.up_histogram[i]/float(self.up_histogram[i] + self.down_histogram[i]))
                output_updown_f.writelines(l)

                if self.round_trip[i] != 0:
                    tmpTrip = self.round_trip[i] / 2.0
                    l = "%d;%f;%f;%f\n" % (i, self.temp_list[i], tmpTrip, (settings.MC_MODULE.stop_step - settings.MC_MODULE.start_step)/float(tmpTrip))
                    output_round_trip_f.writelines(l)
                else:
                    l = "%d;%f;0.00;0.00\n" % (i, self.temp_list[i])
                    output_round_trip_f.writelines(l)

            for output in output_files:
                output.close()

        ### save histograms of temperature
        output_temp_histogram = "%sTemp_%d_%d.hst" % (settings.DIRS['calc'], self.rank, self._stepT)
        output_temp_histogram_f = open(output_temp_histogram, "a+")
        for i in range(self.size):
            l = "%d;%d;%f;%f\n" % (self.rank, i, self.temp_list[i], self.temperature_histogram[i]/float(self.switch_count))
        output_temp_histogram_f.writelines("\n")
        output_temp_histogram_f.close()


    def send(self, data, rank, tag):
        """
        Send data to the given rank with tag
        """
        #print 'send', data, rank, tag
        return self.comm.send(data, dest=rank, tag=tag)

    def recv(self, rank, tag, status):
        #print 'recv', rank, tag
        data = self.comm.recv(source = rank, tag = tag, status = status)
        return data

    def _new_eta(self, hist, temp):
        eta = np.zeros(self.size)
        for i in range(self.size - 1):
            line_derivative = (hist[i+1] - hist[i])/(float(temp[i+1] - temp[i]))
            eta[i] = np.sqrt( np.fabs(line_derivative)/(temp[i+1] - temp[i]))

        return eta

    def _integrate(self, x, y):
        return sum([ y[i]*(x[i+1] - x[i]) for i in range(len(x) - 1) ])

    def _get_kt(self, temp, eta, k):
        M = float(k) / float(self.size - 1)
        sum = 0.0
        i = 0
        while sum < M - 0.0001:
            sum += eta[i]*(temp[i+1] - temp[i])
            i += 1

        return temp[i] - (sum - M)/eta[i-1]

    def _normalize_y(self, x, y):
        c = 1.0/self._integrate(x, y)
        return c * y

    def _new_temps(self, temp, eta):
        new_temp = [0.0] * self.size
        new_temp[0] = temp[0]
        for k in range(1, self.size - 1):
            new_temp[k] = self._get_kt(temp, eta, k)
        new_temp[self.size - 1] = temp[self.size - 1]
        return new_temp

    def generate_new_temp_table(self, temp, hist):
        eta = self._new_eta(hist, temp)
        eta = self._normalize_y(temp, eta)
        temp = self._new_temps(temp, eta)

        return temp

    def conformation_exchange(self, current_replica):
        """
        Get conformation
        @tuple current_replica: (total_energy, conformations, rank)
        """
        message_tag = 7
        current_replica.append(self.rank)

        if self.rank == 0: ## master node, gather all information about conformations
            time_start = time.time()
            index = 1


            datas = {
                current_replica[0]: current_replica
            }

            while index < self.size:
                status = MPI.Status()
                energy, chain_list, source = self.recv(rank = MPI.ANY_SOURCE, tag = settings.TAG + message_tag, status = status)
                source = status.Get_source()
                print "Get from: ", source
                if source != 0:
                    datas[energy] = (energy, chain_list, source)
                    index += 1

            # check which of conformation is in lowest energy state
            min_energy = min(datas.keys())
            conformation = datas[min_energy]

            print "Conformations on replicas:"
            for e,c in datas.iteritems():
                print "%d: Energy: %f" % (c[2], e)

            print "Found conformation with energy: %f on replica: %d" % (min_energy, conformation[2])

            # send to other
            # tuple: (total_energy, conformations, rank)
            for i in range(1, self.size):
                self.send(conformation, rank = i, tag = settings.TAG + message_tag)

            return conformation[1]

        else:
            self.send(current_replica, rank = 0, tag = settings.TAG + message_tag)

            ##
            status = MPI.Status()
            energy, conformation, source = self.recv(rank = 0, tag = settings.TAG + message_tag, status = status)

            print "%d: Get conformation with energy: %f" % (self.rank, energy)
            return conformation

