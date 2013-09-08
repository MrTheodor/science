"""
Main file for starting experiment
Author: Jakub Krajniak <jkrajniak@gmail.com>
"""
#import pyximport; pyximport.install()

import datetime
import logging
import sys
import os
from mpi4py import MPI

import settings
from lib import tools
from lib.Lattice import Lattice
from lib.Box import Box
from main import MC
import hotshot

import signal

## experiment name defince
## package:experiment_class.experiment_setup
## Exp8:Experiment8.c1 -> exp8.py class Experiment8 method c1_setup

#### logging issue
if "." in sys.argv[1]:
    tmp_ex_name = sys.argv[1].split(":")
    ex_module = tmp_ex_name[0].lower()
    ex_name = tmp_ex_name[1].split(".")
else:
    ex_name = sys.argv[1].split(':')
    ex_module = 'experiment'

try:
    job_id = str(os.environ['JOB_ID'])
except KeyError:
    job_id = ''

try: job_name = str(os.environ['JOB_NAME'])
except KeyError: job_name = ''

rank_idx = str(MPI.COMM_WORLD.Get_rank())
rank_idx = "%s_%s_%s_%s_%s" % (sys.argv[1].replace(':','_'), job_name, job_id, rank_idx, datetime.datetime.today().strftime("%Y-%m-%dT%H%M"))

x = logging.getLogger()
x.setLevel(logging.DEBUG)
h = logging.FileHandler(rank_idx + "_log.log")
x.addHandler(h)
###
# pre-settings
## command to run
try:
    cmd = sys.argv[2]
except:
    cmd = ""

valid_cmd = ['clean', 'run', 'start', 'start-box-file', 'continue']
if cmd not in valid_cmd:
    print "Available commands:"
    for v in valid_cmd:
        print "  - %s" % v
    sys.exit(1)

logging.debug('Run with command: %s' % cmd)

if not cmd in ["run", "start-box-file"]:
    settings.CREATE_DIR_STRUCTURE = False
if cmd in ['continue', 'start']:
    settings.EXPERIMENT_COPY = int(sys.argv[3])
    settings.CREATE_DIR_STRUCTURE = True


## run settings for experiment, it will set something in settings.py
m = __import__(ex_module) 

experiment_object = getattr(m, ex_name[0])()

#### PROFILE
PROFILER_MC_FILE = "profile_mc"
PROFILER_SETUP_FILE = "profile_setup"
if settings.PROFILER:
    prof_mc = hotshot.Profile(PROFILER_MC_FILE)
    prof_setup = hotshot.Profile(PROFILER_SETUP_FILE)


try:
    experiment_setup = getattr(experiment_object, ex_name[1] + "_setup")
    if settings.PROFILER:
        (name, chain) = prof_setup.runcall(experiment_setup)
    else:
        (name, chain) = experiment_setup()
    
    settings.TEMP_TABLE = name + "_temperatures.cfg"
except KeyError, e:
    print "Can't run method", ex_name[1]+"_setup"
    print "Available methods:"
    for key,val in experiment_object.__class__.__dict__.iteritems():
        print "   ", key
    print e
    sys.exit(1)

## handle system signals

def handle_12(signum, frame):
    logging.shutdown()
    tools.save_file(settings.ROOT_DIR + "neighbours.dump", Lattice.NEIGHBOURS_POSITION)
    for k,v in settings.EXP_TABLES.iteritems():
        tools.save_file(settings.ROOT_DIR + settings.FILES_TEMPLATE['exp'] % k, settings.EXP_TABLES[k])

signal.signal(12, handle_12)


try:
    mc = None
    if cmd == 'clean':
        for dir in settings.DIRS.values():
            tools.clean_dir(dir)
    elif cmd == 'run':
        
        if settings.PROFILER:
            print "RUN in PROFILER MODE"
            print "MC Steps: 100"
            print "PROFILER_MC_FILE:", PROFILER_MC_FILE
            print "PROFILER_SETUP_FILE:", PROFILER_SETUP_FILE

            mc = MC(steps = 100, obj = chain)
            settings.MC_MODULE = mc
            try:
                temp = float(sys.argv[-1])
                print "Set temp:", temp
                settings.LOCAL_T = temp
                mc.box._localT = temp
            except:
                pass

            prof_mc.runcall(mc.run)
            prof_mc.close()
        else:
            temp = settings.LOCAL_T
            try:
                temp = float(sys.argv[-1])
                print "Set temp:", temp
                settings.LOCAL_T = temp
            except:
                pass

            mc = MC(steps = settings.MC_STEPS, obj = chain)
            if not temp:
                mc.box._localT = temp
            settings.MC_MODULE = mc
            mc.run()
        
    elif cmd == 'start': ## start from some configuration
        print "start start_step temp box_file" 
        start_step = int(sys.argv[4])
        temp = float(sys.argv[5])
        if "jmol" in sys.argv[6]:
            box = Box(settings.BOX_SIZE, localT = temp)
            box.add_chain(chain)
            box.load_jmol_conf(sys.argv[6])
        else:
            box = tools.load_file(sys.argv[6])
        box.localT = temp
        mc = MC(steps = settings.MC_STEPS, obj = box)
        settings.MC_MODULE = mc
        mc.run()
    elif cmd == "start-box-file":
        box_file = sys.argv[3]
        try:
            start_step = int(sys.argv[4])
            print "Start step: ", start_step
        except:
            start_step = 0

        box = tools.load_file(box_file)
        box.localT = settings.LOCAL_T
        mc = MC(steps = settings.MC_STEPS+start_step, obj = box, start_step = start_step)
        settings.MC_MODULE = mc
        if settings.PROFILER:
            import hotshot
            prof = hotshot.Profile(PROFILER_MC_FILE)
            prof.runcall(mc.run)
            prof.close()
        else:
            mc.run()
    elif cmd == 'continue': # continue from some configuration
        #start_step = int(sys.argv[3])
        #temp = float(sys.argv[4])
        max_step = max([ int(o.replace('.dump','').split('_')[-1]) for o in os.listdir(settings.DIRS['data']) ])
        logging.info("Start from step: %d" % max_step)
        box_file_name = settings.DIRS['data'] + settings.FILES_TEMPLATE['data'] % (max_step)
        logging.debug("Box file name: %s" % box_file_name)
        box = tools.load_file(box_file_name)
        if not box:
            raise Exception("Box not found !")
        mc = MC(steps = settings.MC_STEPS, obj = box, start_step = max_step)
        settings.MC_MODULE = mc
        mc.run()


    logging.shutdown()
    for k,v in settings.EXP_TABLES.iteritems():
        tools.save_file(settings.ROOT_DIR + settings.FILES_TEMPLATE['exp'] % k, settings.EXP_TABLES[k])
    #tools.save_file(settings.ROOT_DIR + "neighbours.dump", Lattice.NEIGHBOURS_POSITION)
    #tools.save_file(settings.ROOT_DIR + settings.FILES_TEMPLATE['exp'] % mc.box.localT, settings.EXP_TABLE)
except KeyboardInterrupt: ## stop on Ctrl+C
    logging.shutdown()
    for k,v in settings.EXP_TABLES.iteritems():
        tools.save_file(settings.ROOT_DIR + settings.FILES_TEMPLATE['exp'] % k, settings.EXP_TABLES[k])
    #tools.save_file(settings.ROOT_DIR + "neighbours.dump", Lattice.NEIGHBOURS_POSITION)
    #tools.save_file(settings.ROOT_DIR + settings.FILES_TEMPLATE['exp'] % mc.box.localT, settings.EXP_TABLE)
