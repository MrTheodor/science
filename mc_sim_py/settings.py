'''
Settings for simulation
@author: Jakub Krajniak <jkrajniak@gmail.com>

It is a global settings file, for each experiment some values can be modified
inside experiment class (especialy FREQ values)
'''

import os
DEBUG = False
PROFILER = False
TIME = False
LOGGING = False

## Simulation settings
MIN_T = None
MAX_T = None
LOCAL_T = None
EPSILON = 1
BOX_SIZE = None

MC_STEPS = 100000

EXPERIMENT_NAME = ""
EXPERIMENT_ROOT_NAME = ""
EXPERIMENT_COPY = 0

## dirs
CREATE_DIR_STRUCTURE = True
EXPERIMENT_DIR = "experiment/"
BASE_DIR = os.getcwdu()
ROOT_DIR = os.getcwdu() + "/" + EXPERIMENT_DIR
EXPERIMENTAL_ROOT_DIR = None

TEMP_TABLE = ""

DIRS = {
     "cache": "cache/",
     "model": "models/",
     "log": "logs/",
     "data": "data/",
     "calc": "calc/"
}


FILES = {
     'log_file': 'logging.log',
     'lattice' : 'fcc.dump',
     'temp_cfg': 'temperatures.cfg'
}

FILES_TEMPLATE = {
    "exp": "exp_%f.dump",
    "data": "box_%d.dump",
    "calc": "calc_%s.csv",
    "model": "jmol_%d.jmol",
    "pt_hist": "pt_temp_%d.data"
}

## cache table for exponential values in given local_t
EXP_TABLE = {}
EXP_TABLES = {}

## settings
MC_MODULE = None
BOX_MODULE = None
LATTICE_MODULE = None
PT_MODULE = None

MPI = True
PT = False
FOPT = False
RANK = 0
TAG = 0
SIZE = 1
SIZE_VECTOR = [0]
PERFORM_SWAP = False
CONFORMATION_EXCHANGE = False

SYSTEM_ID = 0

CURRENT_STEP = 0
DELTA_STEP = 0

## calculation frequency
FREQ_COLLECT = int(MC_STEPS/2.0)
FREQ_CALC = 100                       # calculate values every nth step
FREQ_SNAPSHOT = 1000                  # do xyz snapshot every nth step
FREQ_DATA = 1000                      # dump configuration every nth step
FREQ_SWAP = 100                       # perform swap every nth step
FREQ_GC = 1000                        # perform garbage collection
FREQ_ATHERMAL = 500000                # perform athermal state for nth steps
FREQ_COOLING = 100                    # perform linear cooling for nth steps

FREQ_COLLECT = FREQ_COLLECT + FREQ_ATHERMAL + FREQ_COOLING
FREQ_COLLECT2 = FREQ_ATHERMAL + FREQ_COOLING  # # second is used in parallel tempering where we want
                                                                         # # to switch between temperatures after athermal
                                                                         # # but we want to not collect first part of data

RARE_SAVE = 50000                     # save data in some rare steps

JMOL_SCALE = 1.0 ### scale for jmol   # scale jmol xyz snapshot

K_CONST = 1.3806504e-23
