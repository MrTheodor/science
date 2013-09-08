# Jakub Krajniak <jkrajniak@gmail.com>
# Define experiments

from lib import Polymer, tools
from lib.mpi import PT
from lib.Monomer import TypeC, TypeO

from lib.Experiment import Experiment

import logging
import settings


class Experiment4(Experiment):
    """
    DOI 10.1002/pssb.200880252
    """
    box_size = (64, 64, 64)
    name = "exp4"
    length = 64
    epsilon = 1
    chain_seq = None
    global_tag = 4
    freq = { 'pt': 200, 
            'calc': 100, 
            'snapshot': 1000, 
            'data': 1000000, 
            'collect': 0,
            'mc_steps': 1000000
            }
    
    def init(self):
        self.epsilon = 1
        
        self.prepare_chain()

        settings.BASE_DIR = self.get_base_dir()
        settings.EXPERIMENTAL_ROOT_DIR = settings.BASE_DIR
        if settings.CREATE_DIR_STRUCTURE: tools.create_dir(settings.BASE_DIR)

        settings.MC_STEPS = 10000000
        settings.COLLECT = 0

    def prepare_chain(self):
        length = self.length
        co_block = length / 2
        one_block = length / 4
        from lib.Monomer import TypeA, TypeB
        self.chain_seq = Polymer.Chain([ TypeA() if (x % co_block) < one_block else TypeB() for x in range(length)])
        return self.chain_seq

    def athermal_setup(self):
        """
        Athermal experiment 1
        """
        
        self.freq['calc'] = 10000
        self.freq['data'] = 1000000
        self.freq['snapshot'] = 1000

        self.set_collect(**self.freq)
        
        PT.min_temp = settings.MIN_T = -1.0
        PT.max_temp = settings.MAX_T = -1.0
        
        self.pt = PT()
        
        name = "ath"
        self.current_exp_name = name
        self.create_dir_structure("%s_%d_%d_%d" % ((name,) + self.box_size), settings.EXPERIMENT_COPY)

        ## Configure temperature
        settings.LOCAL_T = -1.0
        logging.info("Local T: %f" % settings.LOCAL_T)
        
        self.update_settings_files()

        return (self.current_exp_name, self.chain_seq)

    def const_setup(self):
        """
        Experiment using parallel tempering, with linear temperature function
        """
        name = "const_t_1_14"
        self._setup(name, 1, 14, 2)
        self.set_collect(
            **{'mc_steps': 1000000, 
               'collect': 0, 
               'calc': 1, 
               'data': 1000000, 
               'snapshot': 10000, 
               'pt': False}) 
        
        return (name, self.chain_seq)

    def mpit_setup(self):
        """
        """
        settings.PERFORM_SWAP = True

        name = "mpi_t_1_14"
        self._setup(name, 1, 14, 2, swap = True)
        self.set_collect(**self.freq) 

        return (name, self.chain_seq)

    def mpit2_setup(self):
        settings.COLLECT = 0
        settings.PERFORM_SWAP = True
        name = "mpi2_t_1_14"
        self._setup(name, 1, 14, 20, swap = True)
        self.set_collect(**self.freq)

        return (name, self.chain_seq)
    
class Experiment5(Experiment4):
    """
    """
    box_size = (64, 64, 64)
    name = "exp5"
    length = 64
    global_tag = 5
    freq = {
            'mc_steps': 1000000,
            'pt': 200, 
            'calc': 100, 
            'snapshot': 1000, 
            'data': 1000000, 
            'collect': 0
    }

    def athermal_setup(self):
        
        return super(Experiment5, self).athermal_setup()

    def const_setup(self):
        """
        Experiment using parallel tempering, with linear temperature function
        """
        settings.COLLECT = 0
        settings.PERFORM_SWAP = False
        
        name = "const_t_1_14"
        self._setup(name, 1, 14, 2)
        self.set_collect(**self.freq) 
        
        return (name, self.chain_seq)

    def mpi_setup(self):
        settings.COLLECT = 0
        settings.PERFORM_SWAP = True
        name = "mpi_t_0.1_1.8"
        self._setup(name, 0.1, 1.8, 88, True)
        self.set_collect(**self.freq)
        
        return (name, self.chain_seq)

    def mpi2_setup(self):
        settings.COLLECT = 0
        settings.PERFORM_SWAP = True
        name = "mpi_t_2_4.5"
        self._setup(name, 2.0, 4.5, 99, True)
        self.set_collect(**self.freq)
        
        return (name, self.chain_seq)


class Experiment6(Experiment5):
    box_size = (16, 16, 16)
    name = "exp6"
    length = 64
    global_tag = 7
    freq = {
            'mc_steps': 1000000,
            'pt': 200, 
            'calc': 10, 
            'snapshot': 1000, 
            'data': 1000000, 
            'collect': 0
    }

   

class Experiment7(Experiment5):
    box_size = (32, 32, 32)
    name = "exp7"
    length = 64
    global_tag = 71
    freq = {
            'mc_steps': 1000000,
            'pt': 200, 
            'calc': 1, 
            'snapshot': 1, 
            'data': 1, 
            'collect': 0,
            'athermal': 0
    }


class Experiment8(Experiment5):
    box_size = (128, 128, 128)
    name = "exp8"
    length = 240
    global_tag = 80
    freq = {
            'mc_steps': 5000000,
            'pt': 100, 
            'calc': 100, 
            'snapshot': 1000, 
            'data': 1000000, 
            'collect': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ TypeO() if (x % co_block) < one_block else TypeC() for x in range(length) ])

        return self.chain_seq 
    
    def mpi1_setup(self):
        """
        """
        settings.PERFORM_SWAP = True
        self.temp_list = [1.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0]
        name = "mpi_t_1_100"
        self._setup(name, 1, 100, 2, swap = True)
        self.set_collect(**self.freq) 

        return (name, self.chain_seq)
    
    def mpi4_setup(self):
        """
        """
        self.freq['mc_steps'] = 5000000
        self.freq['collect'] = 0
        self.freq['pt'] = 100
        self.freq['data'] = 1000000
        self.freq['calc'] = 100
        self.freq['snapshot'] = 10000

        settings.PERFORM_SWAP = True

        name = "mpi4_t_1_30"
        self._setup(name, 1, 100, 2, swap = True)
        self.set_collect(**self.freq) 

        return (name, self.chain_seq)


    def mpi2_setup(self):
        settings.COLLECT = 0
        settings.PERFORM_SWAP = True
        name = "mpi2_t_0.5_4"
        self._setup(name, 0.5, 4, 20, swap = True)
        self.set_collect(**self.freq)

        return (name, self.chain_seq)
    
    def mpi3_setup(self):
        """
        """
        settings.PERFORM_SWAP = True

        name = "mpi3_t_14_30"
        self._setup(name, 14, 30, 30, swap = True)
        self.set_collect(**self.freq) 

        return (name, self.chain_seq)

class Experiment9(Experiment7):
    name = "exp9_valid"
    box_size = (32, 32, 32)
    length = 64
    global_tag = 112

    freq = {
        'mc_steps': 2500000,
        'collect': 500000,
        'snapshot': 1000,
        'calc': 100,
        'data': 1000000,
        'pt': 100
    }

    def prepare_chain(self):
        length = self.length
        co_block = 16
        one_block = 8
        from lib.Monomer import TypeA, TypeB

        self.chain_seq = Polymer.Chain([ TypeA() if (x % co_block) < one_block else TypeB() for x in range(length)])
        return self.chain_seq

class Experiment10(Experiment9):
    name = "exp10"
    box_size = (32, 32, 32)
    length = 64
    global_tag = 113

    freq = {
        'mc_steps': 500000,
        'collect': 0,
        'snapshot': 1000,
        'calc': 100,
        'data': 1000000,
        'pt': 100
    }
    def prepare_chain(self):
        length = self.length
        co_block = 16
        one_block = 8
        from lib.MonomerExp9 import TypeA, TypeB

        self.chain_seq = Polymer.Chain([ TypeA() if (x % co_block) < one_block else TypeB() for x in range(length)])
        return self.chain_seq


class HPExperiment(Experiment):
    name = "hp64"
    protein_sequence = "HHHHHHHHHHHHPHPHPPHHPPHHPPHPPHHPPHHPPHPPHHPPHHPPHPHPHHHHHHHHHHHH"
    box_size = (32,32,32)
    length = 64
    global_tag = 15

    freq = {
        'mc_steps': 1000000,
        'collect': 0,
        'snapshot': 1000,
        'calc': 10,
        'data': 1000000,
        'pt': 100
    }

    def prepare_chain(self):
        m = __import__('lib')
        monomer_sequence = []
        for s in self.protein_sequence:
            m_obj = getattr(m, 'Type%s' % s)
            monomer_sequence.append(m_obj())
        self.chain_seq = Polymer.Chain(monomer_sequence, name="hp64")
        print self.name, "Protein sequence", self.protein_sequence
        return self.chain_seq

    def init(self):
        self.epsilon = 1

        self.prepare_chain()

        settings.BASE_DIR = self.get_base_dir()
        settings.EXPERIMENTAL_ROOT_DIR = settings.BASE_DIR
        if settings.CREATE_DIR_STRUCTURE: tools.create_dir(settings.BASE_DIR)

    def athermal_setup(self):
        name = "ath"
        self._setup(name, -1.0, -1.0, 10)
        self.set_collect(**self.freq)

        return (name, self.chain_seq)

    def mpi1_setup(self):
        name = "mpi1"
        self._setup(name, 1, 30, 10, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def mpi2_setup(self):
        name = "mpi1_30_100"
        self._setup(name, 30, 100, 11, swap = True)
        self.freq['mc_steps'] = 5000000
        self.freq['calc'] = 100
        self.freq['data'] = 1000000
        self.freq['snapshot'] = 5000
        self.freq['pt'] = 1000
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)
    
    def mpi3_setup(self):
        name = "mpi1_15_45"
        self._setup(name, 15.0, 45.0, 11, swap = True)
        self.freq['mc_steps'] = 5000000
        self.freq['calc'] = 100
        self.freq['data'] = 1000000
        self.freq['snapshot'] = 5000
        self.freq['pt'] = 1000
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def mpi1130_setup(self):
        name = "mpi1_11_30"
        self._setup(name, 11.0, 30.0, 15, swap = True)
        self.freq['mc_steps'] = 2500000
        self.freq['calc'] = 100
        self.freq['data'] = 1000000
        self.freq['snapshot'] = 1000
        self.freq['pt'] = 500
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)





class ProteinExperiment1(Experiment):
    name = "3MTS"
    protein_sequence = 'GEVEYLCDYKKIREQEYYLVKWRGYPDSESTWEPRQNLKCVRILKQFHKDLERELLRRHHRSKT'
    box_size = (32,32,32)
    length = 64
    global_tag = 10
    chain_seq = None


    def prepare_chain(self):
        import bio
        sequence = bio.seq1_3(self.protein_sequence)
        monomer_sequence = []
        m = __import__('lib')
        for s in sequence:
            m_obj = getattr(m, 'Type%s' % s)
            monomer_sequence.append(m_obj())
        self.chain_seq = Polymer.Chain(monomer_sequence)
        print self.name,"Protein sequence", self.protein_sequence
        return self.chain_seq

    def init(self):
        self.epsilon = 1
        self.prepare_chain()
        settings.BASE_DIR = self.get_base_dir()
        settings.EXPERIMENTAL_ROOT_DIR = settings.BASE_DIR
        if settings.CREATE_DIR_STRUCTURE: tools.create_dir(settings.BASE_DIR)

        settings.MC_STEPS = 10000000
        settings.COLLECT = 0
    
    def const_setup(self):
        settings.COLLECT = settings.MC_STEPS / 2
        settings.PERFORM_SWAP = False
        
        name = "const_t_1_14"
        self._setup(name, 1, 14, 2)
        self.set_collect(settings.MC_STEPS, settings.COLLECT, 100, 1000, 100, False) 
        
        return (name, self.chain_seq)


    def athermal_setup(self):
        """
        Athermal
        """

        settings.MC_STEPS = 1000000
        settings.FREQ_CALC = settings.FREQ_SNAPSHOT = settings.FREQ_DATA = 10000
        name = "ath"
        self._setup(name, -1, -1, 11)
        return (name, self.chain_seq)
