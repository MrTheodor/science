# 
# Author: Jakub Krajniak <jkrajniak@gmail.com>
# Define experiments
#

from lib import Polymer, tools
from lib.Experiment import Experiment

import lib.monomers.MonomerExpHP as MonomerExpHP

import lib.monomers.MonomerExpHP36 as MonomerExpHP36

import bio

import settings

# Sequence from Table II
# Phys. Rev. E 84 031934 (2011) DOI: 10.1103/PhysRevE.84.031934
SEQUENCE_DB = {
    "3d1":  "PPHHHHHPPPHHPPPPPHHPPPHPPPPPPHPHPPPHPPHPPHPPPPPHPPPPHHPHHPPHPPHP",
    "3d2":  "PPHPHPPHPPHHHPHHHHPPHHHPPPPHPHPPPHPHPPPHPHPPPPPHPHPPHPHPPPHPPHPP",
    "3d3":  "HPHHPPHHPHPPPPPHHHPHHHHPPHPPHPHHPPPHPHPPHHHPHHPHPPPPPHHHHHHHHPPP",
    "3d4":  "HPPHHPPHPPHPHPPHPPPPHPPPPPPHPHPHHHPPHPHPPPHPHPPHHPPHPPHPPHPHHHPH",
    "3d5":  "HPPPHHPPHPHPPPHPPPHPHHPPPHHPHPHHPHPPHPPPHPPHPHHHPPHPPHPPHHHPHHHH",
    "3d6":  "HPPHHPHHHHPPPPPPHHPPHPPPPHHPPPHPPHPHHPHPPPPHHPPPPHPPPPPHPPPPHPHH",
    "3d7":  "PPPPHPPPHPPPHHHHPHHPPPPPHPPHPHHPHPHPPPPPHPPPPPPPPPPHHHHPPPPHHPPH",
    "3d8":  "PPPHHHPPHPHPPHPPHHPPPHPPHPPHHPHPPPHPPPPPPPHPHHHPHHHHHPPHHPPPHPPH",
    "3d9":  "HPPHPPHHHPPPPHPHPPPHPHHPHHHHHPPPPHPHPHPPPPHPHPPPHHPHPPPPHPPHHPHP",
    "3d10": "PPHPPHPPHHHPPPHPHPPHPPHPPPPPPHPPHHHPPHPPHPPHPHPPPPPPHHHPPPPPHPHP"
}

class ExperimentHP(Experiment):
    name = "hp_3d1"
    box_size = (64, 64, 64)
    length = 64
    global_tag = 113

    freq = {
        'mc_steps': 5000000,
        'collect': 100,
        'snapshot': 100,
        'calc': 10,
        'data': 100,
        'pt': 200,
        'athermal': 10000,
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
        chain_seq = bio.prepare_sequence(SEQUENCE_DB['3d1'], MonomerExpHP)
        self.chain_seq = Polymer.Chain(chain_seq)
        return self.chain_seq

    def c1_setup(self):
        name = "c1"
        self._setup(name, 0.5, 10.0, 100, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c2_setup(self):
        name = "c2"
        self._setup(name, 0.01, 5.0, 200, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False
        
        return (name, self.chain_seq)


    def m1_setup(self):
        name = "m1"
        self._setup(name, 0.5, 5.0, 101, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

class ExperimentHP36(ExperimentHP):
	def prepare_chain(self):
		length = self.length
		chain_seq = bio.prepare_sequence(SEQUENCE_DB['3d1'], MonomerExpHP36)
		self.chain_seq = Polymer.Chain(chain_seq)
		return self.chain_seq 

if __name__ == "__main__":
    print ExperimentHP
