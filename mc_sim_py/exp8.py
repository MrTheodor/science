# Jakub Krajniak <jkrajniak@gmail.com>
# Define experiments
# Article from 2011

from lib import Polymer, tools
from lib.Monomer import TypeC, TypeO

from lib.MonomerExp8e import TypeC as TypeC05
from lib.MonomerExp8e import TypeO as TypeO05

import lib.MonomerExp8e as MLib

from lib.Experiment import Experiment

import settings


class Experiment8(Experiment):
    box_size = (128, 128, 128)
    name = "exp8"
    length = 240
    global_tag = 80
    freq = {
            'mc_steps': 1000000,
            'pt': 100, 
            'calc': 100, 
            'snapshot': 100, 
            'data': 100, 
            'collect': 0,
            'athermal': 0
    }

    def init(self):
        self.epsilon = 1
        
        self.prepare_chain()

        settings.BASE_DIR = self.get_base_dir()
        settings.EXPERIMENTAL_ROOT_DIR = settings.BASE_DIR
        if settings.CREATE_DIR_STRUCTURE: tools.create_dir(settings.BASE_DIR)

        settings.MC_STEPS = self.freq['mc_steps']
        settings.COLLECT = self.freq['collect']
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ TypeO() if (x % co_block) < one_block else TypeC() for x in range(length) ])

        return self.chain_seq 
    
    def c1_setup(self):
        self.freq['mc_steps'] = 5000000
        #self.temp_list = [0.1, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0]
        name = "c_t_0.1_2.0"
        self._setup(name, 0.1, 2.0, 3, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)
    
    def c2_setup(self):
        name = "c_t_2.0_4.0"
        self._setup(name, 2.0, 4.0, 5, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

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
        self.freq['data'] = 10000
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

class Experiment8b(Experiment8):
    box_size = (240, 240, 240)
    name = "exp8b"
    freq = {
            'mc_steps': 2500000,
            'pt': 100,
            'calc': 10,
            'snapshot': 10,
            'data': 10,
            'collect': 1000,
            'athermal': 100,
            'cooling': 100
    }

    def m1_setup(self):
        self.temp_list = [0.01, 0.1, 1, 10.0]
        name = "m1"
        self._setup(name, 0.01, 2.0, 30, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def m2_setup(self):
        name = "m2"
        self._setup(name, 0.5, 1.0, 150, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)
 
class Experiment8c(Experiment8b):
    box_size = (240, 240, 240)
    name = "exp8c" 
    freq = {
        'mc_steps': 50000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 100,
        'collect': 2000,
        'athermal': 100,
        'cooling': 10
    }
    
    def m4_setup(self):
        settings.PERFORM_SWAP = True
        name = "m4"
        self._setup(name, 0.1, 2.0, 100, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True
        
        return (name, self.chain_seq)

    def m5_setup(self):
        settings.PERFORM_SWAP = True
        name = "m5"
        self._setup(name, 0.5, 2.0, 150, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def m6_setup(self):
        settings.PERFORM_SWAP = True
        name = "m6"
        self._setup(name, 0.70, 0.80, 160, swap = True)
        self.set_collect(**self.freq)

        return (name, self.chain_seq)

    def m8_setup(self):
        settings.PERFORM_SWAP = True
        name = "m8"
        self.temp_list = [0.96, 0.9750, 0.98, 0.99, 0.97]
        self._setup(name, 0.96, 1.0, 160, swap = True)
        self.set_collect(**self.freq)

        return (name, self.chain_seq)

    def m7_setup(self):
        name = "m7"
        #self.temp_list = [3.0, 6.0, 7.0, 10.0]
        self.temp_list = tools.get_list('3.0;3.72716356307;3.8919300086;4.03248996656;4.45029462142;6.32773712706;7.4997732507;8.22677940102;8.76063872225;9.15603842045;9.49544613272;10.0')
        self._setup(name, 3.0, 10.0, 210, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def m9_setup(self):
        name = "m9"
        self._setup(name, 1.8, 3.0, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def m10_setup(self):
        name = "m10"
        self.temp_list = tools.get_list('0.5;0.682184978051;0.822181459851;0.933994546054;1.45231087169;1.84849481007;2.0892410826;2.18968649677;2.37318753997;2.57764955147;2.78018156317;3.0')
        self._setup(name, 0.5, 3.0, 220, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def m13_setup(self):
        name = "m13"
        self.temp_list = tools.get_list('2.0;2.12560282997;2.28684203686;2.54750669268;2.92985687666;3.37317931562;3.63022478269;3.74963892554;3.80798716489;4.16663724091;4.43233034752;4.89163854788;5.23490450547;5.50684521542;5.85620322367;6.0')
        self._setup(name, 2.0, 6.0, 240, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def c16_setup(self):
        name = "c16"
        self.temp_list = tools.get_list('2.0;2.12560282997;2.28684203686;2.54750669268;2.92985687666;3.37317931562;3.63022478269;3.74963892554;3.80798716489;4.16663724091;4.43233034752;4.89163854788;5.23490450547;5.50684521542;5.85620322367;6.0;7.0;8.0;10.0;15.0')
        self._setup(name, 2.0, 6.0, 240, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c17_setup(self):
        name = "c17"
        self._setup(name, 3.74963892554, 3.74963892554, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)


    def m12_setup(self):
        name = "m12"
        self._setup(name, 2.78, 3.1, 230, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def pt1_setup(self):
        name = "pt1"
        self.freq['mc_steps'] = 10000
        self.freq['collect'] = 0
        self.freq['athermal'] = 0
        self.freq['data'] = 5000000
        self.freq['calc'] = 5000000
        self.freq['snapshot'] = 5000000
        self._setup(name, 0.5, 3.0, 300, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True
        settings.FOPT = True

        return (name, self.chain_seq)
 
    def pt2_setup(self):
        name = "pt2"
        self.freq['mc_steps'] = 10000
        self.freq['collect'] = 0
        self.freq['athermal'] = 0
        self.freq['data'] = 5000000
        self.freq['calc'] = 5000000
        self.freq['snapshot'] = 5000000
        self._setup(name, 3.0, 10.0, 320, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True
        settings.FOPT = True

        return (name, self.chain_seq)

    def pt3_setup(self):
        name = "pt3"
        
        self.freq['mc_steps'] = 1000
        self.freq['collect'] = 0
        self.freq['athermal'] = 0
        self.freq['data'] = 5000000
        self.freq['calc'] = 5000000
        self.freq['snapshot'] = 5000000
        self._setup(name, 2.0, 10.0, 330, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True
        settings.FOPT = True

        return (name, self.chain_seq)

        

    def c6_setup(self):
        settings.PERFORM_SWAP = False
        name = "c6"
        self._setup(name, 0.70, 0.80, 160, swap = False)
        settings.PERFORM_SWAP = False
        self.set_collect(**self.freq)

        return (name, self.chain_seq)

    def c7_setup(self):
        settings.PERFORM_SWAP = False
        name = "c7"
        self._setup(name, 0.8, 2.0, 150, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c8_setup(self):
        settings.PERFORM_SWAP = False
        name = "c8"
        self.temp_list = [ 0.6, 1.0, 2.0, 10.0 ]
        self._setup(name, 0.8, 2.0, 190, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c9_setup(self):
        name = "c9"
        self._setup(name, 0.6, 1.1, 200, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c10_setup(self):
        name = "c10"
        self._setup(name, 0.45, 0.55, 210, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c11_setup(self):
        name = "c11"
        self._setup(name, 3.0, 3.0, 210, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c12_setup(self):
        name = "c12"
        self.temp_list = [0.9750, 1.0, 0.850, 0.90 ]
        self._setup(name, 3.0, 3.0, 210, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)


    def c13_setup(self):
        name = "c13"
        self._setup(name, 0.8750, 1.0, 210, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c14_setup(self):
        name = "c14"
        self._setup(name, 3.0, 6.0, 220, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c15_setup(self):
        name = "c15"
        self._setup(name, 0.5, 3.0, 230, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c20_setup(self):
        name = "c20"
        self.temp_list = [ 2.266667, 2.90, 3.533333, 4.16, 4.80, 5.43, 6.066667, 6.7, 7.333333, 7.96 ]
        self._setup(name, 2.27, 7.96, 350, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c21_setup(self):
        name = "c21"
        self.temp_list = [ 2.266667, 2.9, 3.2166665, 3.533333, 4.800000, 4.9575, 5.115, 5.2725, 5.43, 5.745, 6.066667, 6.7, 7.333333, 8.600000, 11.13, 12.40 ]
        self._setup(name, 2.27, 8.96, 360, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)


class Experiment8d(Experiment8c):
    box_size = (240, 240, 240)
    name = "exp8d"
    freq = {
        'mc_steps': 10000000,
        'pt': 200,
        'calc': 100,
        'snapshot': 100,
        'data': 1000,
        'collect': 1000000,
        'athermal': 1000000,
        'cooling': 0
    }
 
class Experiment8e(Experiment8c):
    box_size = (240, 240, 240)
    name = "exp8e"
    length = 120
    freq = {
        'mc_steps': 2000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 100,
        'collect': 5000,
        'athermal': 10000,
        'cooling': 0
    }

class Experiment8g(Experiment8c):
    """
    For epsilon = 0.5
    """

    box_size = (240, 240, 240)
    name = "exp8g0_5"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 100,
        'collect': 20000000,
        'athermal': 5000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ TypeO05() if (x % co_block) < one_block else TypeC05() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def m1_setup(self):
        name = "m1"
        self._setup(name, 1.0, 20.0, 259, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def c2_setup(self):
        name = "c2"
        self._setup(name, 1.0, 3.0, 291, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c3_setup(self):
        name = "c3"
        self._setup(name, 1.0, 7.0, 270, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False
        return (name, self.chain_seq)

    def c4_setup(self):
        name = "c4"
        self._setup(name, 5.0, 8.0, 280, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def m16_setup(self):
        name = "m16"
        self.temp_list = tools.get_list('2.0;2.93869324227;3.07217225203;3.30758457056;3.94316528763;4.29872503277;4.52461320786;5.58411715289;6.17418580166;6.79175962589;7.56188084477;8.53862159296;8.94672367604;9.45417046881;9.85630466154;10.0')
        self._setup(name, 2.0, 10.0, 240, swap = True)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = True

        return (name, self.chain_seq)

    def c20_setup(self):
        name = "c20"
        self.temp_list = [ 2.27, 2.90, 3.53, 4.16, 4.80, 5.43, 6.07, 6.7, 7.33, 7.96 ]
        self._setup(name, 2.27, 7.96, 350, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

    def c21_setup(self):
        name = "c21"
        self.temp_list = [ 2.266667, 3.533333, 4.800000, 6.066667, 7.333333, 8.600000 ]
        self._setup(name, 2.27, 8.96, 360, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)




class Experiment8f(Experiment8g):
    box_size = (240, 240, 240)
    name = "exp8f0_1"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 1000,
        'collect': 10000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ TypeO() if (x % co_block) < one_block else TypeC() for x in range(length) ])

        return self.chain_seq

class Experiment801(Experiment8g):
    box_size = (240, 240, 240)
    name = "exp8g0_1"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 1000,
        'data': 1000,
        'collect': 10000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ TypeO() if (x % co_block) < one_block else TypeC() for x in range(length) ])

        return self.chain_seq
    


class Experiment480a(Experiment8g):
    box_size = (240, 240, 240)
    name = "exp480a_05"
    length = 480
    freq = {
        'mc_steps': 10000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 1000,
        'collect': 1000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ TypeO05() if (x % co_block) < one_block else TypeC05() for x in range(length) ])

        return self.chain_seq


class Experiment802(Experiment8c):
    """
    For epsilon = 0.2
    """

    box_size = (240, 240, 240)
    name = "exp8g0_2"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 1000,
        'data': 100,
        'collect': 10000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO02() if (x % co_block) < one_block else MLib.TypeC02() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

class Experiment803(Experiment8c):
    """
    For epsilon = 0.3
    """

    box_size = (240, 240, 240)
    name = "exp8g0_3"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 1000,
        'data': 100,
        'collect': 10000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO03() if (x % co_block) < one_block else MLib.TypeC03() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

class Experiment804(Experiment8c):
    """
    For epsilon = 0.4
    """

    box_size = (240, 240, 240)
    name = "exp8g0_4"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 1000,
        'data': 100,
        'collect': 20000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO04() if (x % co_block) < one_block else MLib.TypeC04() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
	self.temp_list = [ 2.27, 3.53, 4.80, 6.07, 7.33, 8.60, 9.87, 11.13, 12.40, 13.67, 14.93, 16.20, 17.47, 18.73, 20.0, 21.0]
        self._setup(name, 2.27, 21.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

class Experiment805(Experiment8c):
    """
    For epsilon = 0.5
    """

    box_size = (240, 240, 240)
    name = "exp8g0_5"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 1000,
        'data': 100,
        'collect': 20000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO02() if (x % co_block) < one_block else MLib.TypeC02() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)



class Experiment806(Experiment8c):
    """
    For epsilon = 0.6
    """

    box_size = (240, 240, 240)
    name = "exp8g0_6"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 100,
        'collect': 10000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO06() if (x % co_block) < one_block else MLib.TypeC06() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

class Experiment807(Experiment8c):
    """
    For epsilon = 0.7
    """

    box_size = (240, 240, 240)
    name = "exp8g0_7"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 100,
        'collect': 10000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO07() if (x % co_block) < one_block else MLib.TypeC07() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

class Experiment808(Experiment8c):
    """
    For epsilon = 0.8
    """

    box_size = (240, 240, 240)
    name = "exp8g0_8"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 100,
        'collect': 10000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO08() if (x % co_block) < one_block else MLib.TypeC08() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

class Experiment809(Experiment8c):
    """
    For epsilon = 0.9
    """

    box_size = (240, 240, 240)
    name = "exp8g0_9"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 100,
        'collect': 20000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO09() if (x % co_block) < one_block else MLib.TypeC09() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)


class Experiment810(Experiment8c):
    """
    For epsilon = 1.0
    """

    box_size = (240, 240, 240)
    name = "exp8g0_10"
    length = 240
    freq = {
        'mc_steps': 40000000,
        'pt': 100,
        'calc': 100,
        'snapshot': 100,
        'data': 100,
        'collect': 20000000,
        'athermal': 1000000,
        'cooling': 0
    }
    
    def prepare_chain(self):
        length = self.length
        co_block = 20
        one_block = 10

        self.chain_seq = Polymer.Chain([ MLib.TypeO10() if (x % co_block) < one_block else MLib.TypeC10() for x in range(length) ])

        return self.chain_seq 

    def c1_setup(self):
        name = "c1"
        self._setup(name, 1.0, 20.0, 250, swap = False)
        self.set_collect(**self.freq)
        settings.PERFORM_SWAP = False

        return (name, self.chain_seq)

