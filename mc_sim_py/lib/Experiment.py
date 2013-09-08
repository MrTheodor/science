import os
import logging

import settings

from lib import Polymer, tools
from lib.Box import Box, Lattice
from lib.mpi import PT
from lib.Monomer import TypeA, TypeB, TypeA2, TypeB2
from lib import Monomer


class Experiment(object):
    current_exp_name = None
    box_size = None
    temp_list = None
    
    def __init__(self):
        
        settings.BOX_SIZE = self.box_size

        Lattice.box_size = self.box_size
        self.init()

    def get_base_dir(self, extra = None):
        """
        Build base structure
        """
        
        exp_dir_name = "%s_%d" % (self.name, settings.RANK)
        dir = "%s/%s%s/%s" % (os.getcwdu(), settings.EXPERIMENT_DIR, exp_dir_name, extra if extra else "")
        ## sometimes we need a copy of experiment, additional one will be put in _1, _2 dirs
        return dir

    def create_dir_structure(self, name, index = None):
        """
        @name string: name of dir above base_dir
        @make_copy int: index of dir, if exists create max(index)+1
        Create dir structure
        """
        if index != None:
            base_dir = self.get_base_dir()
            if settings.CREATE_DIR_STRUCTURE:
                while os.path.exists(base_dir + (name + "_%d" % index)):
                    max_index = max([ int(p.split('_')[-1]) for p in os.listdir(base_dir) if name in p ])
                    index = max_index + 1
             
            name_with_index = (name + "_%d") % index
            settings.BASE_DIR = self.get_base_dir(name_with_index)
        else:
            settings.BASE_DIR = self.get_base_dir(name)

        for name,dir in settings.DIRS.iteritems():
            settings.DIRS[name] = settings.BASE_DIR + "/" + settings.DIRS[name]
            if settings.CREATE_DIR_STRUCTURE: tools.create_dir(settings.DIRS[name])

    def update_settings_files(self):
        """
        Update FILES dict in settings
        """
        for name, dir in settings.FILES_TEMPLATE.iteritems():
            if settings.DIRS.has_key(name) and settings.FILES_TEMPLATE.has_key(name):
                settings.FILES[name] = settings.DIRS[name] + settings.FILES_TEMPLATE[name]
    
    def _setup(self, name, min_t, max_t, tag, swap = False):
        """
        Setup experiment
        
        @string name: name of experiment
        @float min_t: min temperature
        @float max_t: max temperature
        @int tag: tag for mpi communication
        @bool swap: perform swap?

        """
        settings.EXPERIMENT_NAME = name
        settings.EXPERIMENT_ROOT_NAME = self.name

        settings.PERFORM_SWAP = swap
        PT.min_temp = settings.MIN_T = min_t
        PT.max_temp = settings.MAX_T = max_t
        if self.temp_list is None:
            self.pt = PT()
        else:
            self.pt = PT(self.temp_list)

        if self.pt.size > 0:
            settings.PT = True

        settings.LOCAL_T = self.pt.get_temperature()
        logging.info("Local T: %f" % settings.LOCAL_T)
        
        settings.TAG = self.global_tag + tag
        settings.DELTA_STEP = 0
        settings.PT_MODULE = self.pt
        self.create_dir_structure("%s_%d_%d_%d" % ((name, ) + self.box_size), settings.EXPERIMENT_COPY)
        
        self.update_settings_files()
        self.current_exp_name = name
        
        temp_file = "%stemperature" % settings.EXPERIMENTAL_ROOT_DIR
        if not os.path.exists(temp_file):
            temp_file_mode = "w+"
        else:
            temp_file_mode = "a+"

        open(temp_file, temp_file_mode).writelines("%s;%f\n" % (name, settings.LOCAL_T))

    def set_collect(self, mc_steps, collect, calc, data = None, snapshot = None, pt = 200, athermal=500000, cooling = 0):
        
        if collect is False: 
            collect = mc_steps*2
        elif collect is None: 
            collect = self.freq['collect']

        if calc is False: 
            calc = mc_steps*2
        elif calc is None:
            calc = self.freq['calc']

        if data is False: 
            data = mc_steps*2
        elif data is None:
            data = self.freq['data']

        if snapshot is False: 
            snapshot = mc_steps*2
        elif snapshot is None:
            snapshot = self.freq['snapshot']

        if pt is False: 
            pt = mc_steps * 2
        elif pt is None:
            pt = self.freq['pt']

        if athermal is False:
            athermal = 0
        elif athermal is None:
            athermal = self.freq['athermal']

        if cooling is False:
            cooling = 0
        elif cooling is None:
            cooling = self.freq['cooling']
        
        settings.MC_STEPS = mc_steps
        settings.FREQ_COLLECT = collect
        settings.FREQ_DATA = data
        settings.FREQ_SNAPSHOT = snapshot
        settings.FREQ_CALC = calc
        settings.FREQ_SWAP = pt
        settings.FREQ_ATHERMAL = athermal
        settings.FREQ_COOLING = cooling
    
