"""
Main MC class for Monte Carlo simulation
Author: Jakub Krajniak <jkrajniak@gmail.com>

"""

from lib import Box, Chain, Lattice, Polymer, tools
import kyotocabinet as kb
import logging
import numpy as np
import settings
import sys
import time

class MC(object):
    start_step = None
    stop_step = None
    current_step = None
    box = None

    current_state = None

    stop_simulation = False

    hdfTable = None

    # # STATES
    START = 0
    STOP = 100

    def __init__(self, steps, obj=None, start_step=0):
        self.start_step = start_step
        self.current_step = start_step
        self.stop_step = steps
        if isinstance(obj, Box):
            # # when box is loaded from file
            self.box = obj
            Lattice.box_size = obj.box_size
            settings.BOX_SIZE = obj.box_size
        elif isinstance(obj, Chain):
            self.box = Box(settings.BOX_SIZE, localT=settings.LOCAL_T)
            self.box.add_chain(obj)
        elif isinstance(obj, list):
            self.box = Box(settings.BOX_SIZE, localT=settings.LOCAL_T)
            for c in obj:
                self.box.add_chain(c)

        # # load some specific settings
        Lattice.NEIGHBOURS_POSITION = tools.load_file(settings.ROOT_DIR + "neighbours.dump", {})

        if self.box.localT in settings.EXP_TABLES:
            settings.EXP_TABLE = settings.EXP_TABLES[self.box.localT]
        else:
            settings.EXP_TABLE = tools.load_file(settings.ROOT_DIR + \
                settings.FILES_TEMPLATE['exp'] % self.box.localT, {})
            settings.EXP_TABLES[self.box.localT] = settings.EXP_TABLE

    def prepare_db(self):
        """
        Prepare database to store values
        """
        db_filename = settings.BASE_DIR + "/" + "database.kch"
        self.db = kb.DB()
        if not self.db.open(db_filename, kb.DB.OWRITER | kb.DB.OCREATE):
            print >> sys.stderr, "open error:" + str(self.db.error())
            sys.exit(1)

        print "DB:", db_filename

    def change_stop_step(self, stop_step):
        """
        Set new value of stop step, need for FOPT alghoritm
        """
        self.stop_step = stop_step
        logging.info("New value of stop_step: %d" % stop_step)
        if self.stop_simulation:
            logging.info("Set False to stop simulation")
            self.stop_simulation = False

    def change_current_step(self, current_step):
        """
        Set new value of current step, need for prevent dead lock on MPI
        """
        # if current_step != self.current_step:
        self.current_step = current_step
        # logging.debug("New value of current_step: %d" % self.current_step)

    def set_stop_simulation(self):
        self.stop_simulation = True
        logging.info("Set stop simulation")

    def run(self):
        """
        Run Monte Carlo simulation
        """
        self.start_step = self.start_step + settings.DELTA_STEP

        print "="*10
        print "Local T:", self.box._localT
        print "BOX size:", self.box.box_size
        print "Lattice size:", Lattice.box_size, "/", self.box.lattice.box_size
        print "Number of chains:", len(self.box.chain_list)
        print "Rank:", settings.RANK, "/", settings.SIZE
        print "Logging file:", settings.DIRS['log'] + settings.EXPERIMENT_NAME + ".log"
        print "BASE_DIR:", settings.BASE_DIR
        print "Dirs:"
        for name, dir_name in settings.DIRS.items():
            print "  ", name, ":", dir_name
        print "=" * 10
        print "SETTINGS:"
        print "\tMC steps:", settings.MC_STEPS
        print "\t\tMC start_step:", self.start_step
        print "\t\tMC stop_step:", self.stop_step
        print "\tFREQ_SWAP:", settings.FREQ_SWAP
        print "\tFREQ_COLLECT", settings.FREQ_COLLECT
        print "\tPERFORM_SWAP:", settings.PERFORM_SWAP
        print "\tFREQ_DATA:", settings.FREQ_DATA
        print "\tFREQ_CALC:", settings.FREQ_CALC
        print "\tFREQ_SNAPSHOT:", settings.FREQ_SNAPSHOT
        print "\tFREQ_ATHERMAL:", settings.FREQ_ATHERMAL
        print "\tFREQ_COOLING:", settings.FREQ_COOLING
        print "\tPT:", settings.PT
        print "\tPT_MODULE:", settings.PT_MODULE
        print "\tFOPT:", settings.FOPT
        print "="*10
        print "Chains"
        m_set = set()
        for c in range(len(self.box.chain_list)):
            print "\t%d: %s" % (c, "".join([ m.name for m in self.box.chain_list[c].chain ]))
            for m in self.box.chain_list[c].chain:
                m_set.add((m.name, str(m.interaction_table)))
        print "\tMonomers:"
        for m in m_set:
            print "\t  %s: %s" % (m[0], m[1])
        print "="*10

        # # mc step
        step_time = 0.0

        # self.prepare_hdf_table()
        self.prepare_db()

        settings.FREQ_COLLECT = settings.FREQ_COLLECT + settings.FREQ_ATHERMAL + settings.FREQ_COOLING
        settings.FREQ_COLLECT2 = settings.FREQ_ATHERMAL + settings.FREQ_COOLING  # # second is used in parallel tempering where we want
                                                                                 # # to switch between temperatures after athermal
                                                                                 # # but we want to not collect first part of data

        # ## start simulations
        simulation_time = time.time()

        xrange_chains = [None] * self.box.number_of_chains
        for c in self.box.chain_list:
            xrange_chains[c.idx] = xrange(c._total_length / 8)

        # # MOVES
        ch = Polymer.Chain
        break_move = [ch.PULL_MOVE, ch.PULL_MOVE, ch.B_SNAKE, ch.E_SNAKE]
        # break_move = []
        m = (ch.B_ROTATE, ch.E_ROTATE, ch.S_ROTATE, ch.B_SNAKE, ch.E_SNAKE)
        m_len = len(m)
        b_len = len(break_move)

        # # variables
        acceptance_ratio = dict(zip(ch.TYPE_OF_MOVE, [0] * len(ch.TYPE_OF_MOVE)))
        reject_ratio = acceptance_ratio.copy()
        total_acceptance = 0
        total_reject = 0

        # step file
        step_file = open(settings.DIRS['calc'] + "step", "w+")

        # # MC SWEEP while loop
        while self.current_step < self.stop_step and not self.stop_simulation:
            # print settings.PT_MODULE.rank if settings.PT else 0, 'step', step

            step = self.current_step

            if settings.TIME: start_time = time.time()
            move_name = ""

            settings.CURRENT_STEP = self.current_step

            # # step environment changes like athermal heating and linear cooling
            if (step < settings.FREQ_ATHERMAL or settings.LOCAL_T < 0.0) and not self.box.athermal_state:
                logging.info("Launch athermal mode for %d" % settings.FREQ_ATHERMAL)
                self.box.athermal_state = True
            elif step >= settings.FREQ_ATHERMAL and self.box.athermal_state:
                if self.box.athermal_state:
                    delta_time = time.time() - simulation_time
                    logging.info("Switch off athermal state on %d (time: %f)" % (step, delta_time))
                # end of athermal state, if we use parallel tempering, try to find lowest
                # enery configuration in replicas and run every replica from it

                self.box.athermal_state = False

            move = False
            accept = 0
            reject = 0

            for rr_chain in range(self.box.number_of_chains):

                if self.box.number_of_chains == 1:
                    chain = self.box.chain_list[0]
                else:
                    chain_idx = np.random.randint(0, self.box.number_of_chains)
                    chain = self.box.chain_list[chain_idx]

                move_time = 0.0

#                for rr_chain in xrange_chains[chain.idx]:
                for rr_chain in [1]:
                    if settings.TIME: start_move = time.time()

                    if settings.TIME: m_time = time.time()

                    # random move
                    move_idx = np.random.randint(0, b_len)
                    type_of_move = break_move[move_idx]

                    move = chain.move(type_of_move)

                    if settings.TIME: m_time = time.time() - m_time
                    if settings.TIME: move_time += m_time
                    if settings.TIME:
                        stop_move = time.time()
                        logging.info("Move %s: %f" % (chain.move_name, stop_move - start_move))

                    if move:
                        acceptance_ratio[chain.move_name] += 1
                        total_acceptance += 1
                        accept += 1
                    else:
                        total_reject += 1
                        reject_ratio[chain.move_name] += 1
                        reject += 1

                    if settings.TIME: move_name = chain.move_name
                    if settings.DEBUG: print chain.move_name, move

            if settings.TIME: print step, time.time() - start_time, 'acc', accept, 'rej', reject, 'move_time_avg', move_time

            # # performe some saves
            rare_save = False
            if step <= settings.FREQ_COLLECT: rare_save = (step % settings.RARE_SAVE == 0)

            if settings.TIME: chain_move_time = time.time()

            if (step % 10000) == 0:
                step_file.write("%d\n" % step)
                step_file.flush()

            # save some data
            if (step >= settings.FREQ_COLLECT and (step % settings.FREQ_SNAPSHOT) == 0) or rare_save:
                self.db["jmol_%d" % step ] = self.box.jmol_snapshot_string()


            if (step >= settings.FREQ_COLLECT and (step % settings.FREQ_CALC) == 0):
                for chain in self.box.chain_list:
                    chain.calculate_values()
                    val = [step, self.box._localT] + chain.last_values
                    self.db["values_%d" % step ] = val

            # ## swap
            if settings.PERFORM_SWAP:
                if step >= settings.FREQ_COLLECT2 and settings.PT and ((step + 1) % settings.FREQ_SWAP) == 0:
                    self.box.swap()
                    data_file = settings.DIRS['calc'] + settings.FILES_TEMPLATE['pt_hist'] % (settings.RANK)
                    tools.save_file(data_file, self.box.temperature_history)

            if settings.TIME:
                stop_time = time.time()
                delta_move_time = chain_move_time - start_time
                logging.info("%s;%f;%d" % (move_name, delta_move_time, move))

            # # END of simulation block
            # # increment step counter
            step += 1
            self.current_step += 1

            if step == self.stop_step:  # # last step of simulation, we check feedback
                if settings.PERFORM_SWAP and settings.PT and settings.FOPT:
                    logging.info("Perform Feedback optimization")
                    settings.PT_MODULE.calculate_round_trip()

            ##################################

        simulation_time = time.time() - simulation_time

        # data_file_name = settings.FILES['data'] % (step)
        # tools.save_file(data_file_name, self.box)
        # jmol_file = settings.FILES['model'] % (step)
        # self.box.jmol_snapshot(jmol_file)

        self.db['jmol_%d' % step ] = self.box.jmol_snapshot_string()

        # self.values_tbl.flush()

        # # save remaning values
        # for k,v in values_to_save.iteritems():
        #    f = values_files[k]
        #    f.writelines(v)
        #    f.close()

        # # save histograms
        if settings.PERFORM_SWAP and settings.PT:
            settings.PT_MODULE.save_histograms()

        print "="*5, "Rank", settings.RANK, "of", settings.SIZE, "="*5
        print "Total time:", simulation_time
        print "Rank:", settings.RANK
        print "Steps:", self.stop_step
        print "Chains:", len(self.box.chain_list)
        for chain in self.box.chain_list:
            print "\tChain: %s" % chain
            print "\t  Min Energy:", chain.min_energy
            print "\t  Min Rg:", chain.min_rg
            for i in range(len(Chain.TYPE_OF_MOVE)):
                print "\t  +%s: %d" % (chain.TYPE_OF_MOVE[i], chain.MOVE_STATS[i])
        print "Total accepted: ", total_acceptance
        print "Total reject: ", total_reject
        print "Acceptance move:"
        for k, v in acceptance_ratio.iteritems():
            print " +%s: %s" % (str(k), str(v))
        print "Reject move:"
        for k, v in reject_ratio.iteritems():
            print " +%s: %s" % (str(k), str(v))

        print "Temperature: ", self.box._localT
        print "="*10
