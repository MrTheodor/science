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
    name = "multichain1"
    box_size = (240, 240, 240)
    length = 240
    global_tag = 112

    freq = {
        'mc_steps': 40000000,
        'collect': 15000000,
        'snapshot': 1000,
        'calc': 100,
        'data': 1,
        'pt': 100,
        'swap': 200,
        'athermal': 100000
    }

    def init(self):
        self.epsilon = 1

        self.prepare_chain()

        settings.BASE_DIR = self.get_base_dir()
        settings.EXPERIMENTAL_ROOT_DIR = settings.BASE_DIR
        if settings.CREATE_DIR_STRUCTURE: tools.create_dir(settings.BASE_DIR)

    def const_setup(self):
        """
        Experiment using parallel tempering, with linear temperature function
        """
        settings.COLLECT = 0
        settings.PERFORM_SWAP = False

        name = "const"
        self._setup(name, 1, 14, 2)
        self.set_collect(**self.freq)

        return (name, self.chain_seq)

    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ TypeO() if (x % co_block) < one_block else TypeC() for x in range(length)])
        return self.chain_seq

    def c21_setup(self):
        name = "c21"
        self.temp_list = [2.266667, 2.9, 3.2166665, 3.533333, 4.800000, 4.9575, 5.115, 5.2725, 5.43, 5.745, 6.066667, 6.7, 7.333333, 8.600000, 11.13, 12.40]
        self._setup(name, 2.27, 8.96, 360, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def fopt_setup(self):
        name = "fopt_1"
        self.freq['mc_steps'] = 1000
        self.freq['collect'] = 0
        self.freq['athermal'] = 0
        self.freq['data'] = 5000000
        self.freq['calc'] = 5000000
        self.freq['snapshot'] = 5000000
        self._setup(name, 1.0, 10.0, 320, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True
        settings.FOPT = True
        return (name, self.chain_seq)

class Experiment901M(Experiment9):
    name = "multichain901M"
    length = 64
    box_size = (64, 64, 64)
    freq = {
        'mc_steps': 1000000,
        'collect': 100,
        'snapshot': 1000,
        'calc': 100,
        'data': 1,
        'pt': 200,
        'athermal': 500000
    }

    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = [ Polymer.Chain([ TypeO() if (x % co_block) < one_block else TypeC() for x in range(length)]) for x in range(20) ]
        return self.chain_seq


