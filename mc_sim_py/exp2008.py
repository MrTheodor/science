# 
# Author: Jakub Krajniak <jkrajniak@gmail.com>
# Define experiments
#

from lib import Polymer, tools
from lib.Experiment import Experiment

from lib.Monomer import TypeA, TypeB


import settings

class Experiment2008(Experiment):
    name = "exp2008"
    box_size = (64, 64, 64)
    length = 64
    global_tag = 113

    freq = {
        'mc_steps': 10000000,
        'collect': 0,
        'snapshot': 100,
        'calc': 10,
        'data': 100,
        'pt': 200,
        'athermal': 1000000,
        'cooling': 0
    }

    def init(self):
        self.epsilon = 1
        
        self.prepare_chain()

        settings.BASE_DIR = self.get_base_dir()
        settings.EXPERIMENTAL_ROOT_DIR = settings.BASE_DIR
        if settings.CREATE_DIR_STRUCTURE: tools.create_dir(settings.BASE_DIR)

    def prepare_chain(self):
        length = self.length
        co_block = length / 4
        one_block = length / 8
        self.chain_seq = Polymer.Chain([ TypeA() if (x % co_block ) < one_block else TypeB() for x in range(length) ])

        return self.chain_seq

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.2, 12, 100, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)
