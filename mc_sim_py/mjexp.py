# 
# Author: Jakub Krajniak <jkrajniak@gmail.com>
# Define experiments
#

from lib import Polymer, tools
from lib.Experiment import Experiment

from lib.monomers import MonomersMJ

import bio

import settings

SEQUENCE_DB = {
    'seq1': ("L"*10+"K"*10)*12,
    'seq2': ("L"*10+"K"*10)*12,
    'seq3': ("L"*10+"N"*10)*12, # Lys + Asn
    'seq4': ("G"*10 + "Y"*10) * 12 # Gly (hydrophobic) - Tyr (hydrophylic)
}

class ExperimentMJa(Experiment):
    name = "mj240a"
    box_size = (240, 240, 240)
    length = 240
    global_tag = 113
    seq = 'seq1'

    freq = {
        'mc_steps': 40000000,
        'collect': 10000000,
        'snapshot': 1000,
        'calc': 100,
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
        chain_seq = bio.prepare_sequence(bio.seq1_3(SEQUENCE_DB[self.seq]), MonomersMJ)
        self.chain_seq = Polymer.Chain(chain_seq)
        return self.chain_seq

    def c1_setup(self):
        name = "c1"
        self._setup(name, 0.5, 40.0, 100, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c2_setup(self):
        name = "c2"
        self._setup(name, 20.0, 40.0, 200, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c3_setup(self):
        name = "c3"
        self.temp_list = [50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0]
        self._setup(name, 50.0, 2000.0, 200, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

	return (name, self.chain_seq)



class ExperimentMJb(Experiment):
    name = "mj240b"
    box_size = (240, 240, 240)
    length = 240
    global_tag = 118
    seq = 'seq3'

    freq = {
        'mc_steps': 40000000,
        'collect': 10000000,
        'snapshot': 1000,
        'calc': 100,
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
        chain_seq = bio.prepare_sequence(bio.seq1_3(SEQUENCE_DB[self.seq]), MonomersMJ)
        self.chain_seq = Polymer.Chain(chain_seq)
        return self.chain_seq

    def c1_setup(self):
        name = "c1"
        self._setup(name, 0.5, 40.0, 100, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

class ExperimentMJc(ExperimentMJb):
    name = "mjGlyTyr"
    box_size = (240, 240, 240)
    length = 240
    global_tag = 119
    seq = 'seq4'



class Experiment3MTS(Experiment):
    name = "3mts"
    box_size = (64, 64, 64)
    length = 64
    global_tag = 118
    seq = 'GEVEYLCDYKKIREQEYYLVKWRGYPDSESTWEPRQNLKCVRILKQFHKDLERELLRRHHRSKT'

    freq = {
        'mc_steps': 40000000,
        'collect': 10000000,
        'snapshot': 1000,
        'calc': 100,
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
        chain_seq = bio.prepare_sequence(bio.seq1_3(self.seq), MonomersMJ)
        self.chain_seq = Polymer.Chain(chain_seq)
        return self.chain_seq

    def c1_setup(self):
        name = "c1"
        self._setup(name, 0.5, 40.0, 100, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c2_setup(self):
        name = "c2"
        self._setup(name, 40, 100.0, 100, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

