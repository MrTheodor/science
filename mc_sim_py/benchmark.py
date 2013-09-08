# Jakub Krajniak <jkrajniak@gmail.com>
# Define experiments

from lib import Polymer, tools
from lib.Experiment import Experiment
from lib.MonomerExp8e import TypeC as TypeC05
from lib.MonomerExp8e import TypeO as TypeO05
from lib.Monomer import TypeC, TypeO
import lib.MonomerExp8e as MLib

import settings


class Experiment9(Experiment):
    name = "exp510"
    box_size = (200, 200, 200)
    length = 100
    global_tag = 112

    freq = {
        'mc_steps': 1000,
        'collect': 1,
        'snapshot': 1,
        'calc': 100,
        'data': 100,
        'pt': 1,
        'athermal': 0
    }

    def init(self):
        self.epsilon = 1
        
        self.prepare_chain()

        settings.BASE_DIR = self.get_base_dir()
        settings.EXPERIMENTAL_ROOT_DIR = settings.BASE_DIR
        if settings.CREATE_DIR_STRUCTURE: tools.create_dir(settings.BASE_DIR)

    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 5

        self.chain_seq = Polymer.Chain([ TypeO() for x in range(length)])
        return self.chain_seq

    def c21_setup(self):
        name = "c21"
        #self.temp_list = [ 2.266667, 2.9, 3.2166665, 3.533333, 4.800000, 4.9575, 5.115, 5.2725, 5.43, 5.745, 6.066667, 6.7, 7.333333, 8.600000, 11.13, 12.40 ]
        self._setup(name, 2.27, 20.0, 360, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

class Experiment901(Experiment9):
    name = "exp510_01"
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 5

        self.chain_seq = Polymer.Chain([ TypeO() if (x % co_block) < one_block else TypeC() for x in range(length)])
        return self.chain_seq

