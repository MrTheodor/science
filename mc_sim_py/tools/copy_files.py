import os

import argparse as arg
import glob

import shutil


parser = arg.ArgumentParser(description="Process files")
parser.add_argument('--config-file', help="Config file", dest="cfile")
parser.add_argument('-r', help='Root dir where experiment is stored', type=str)
parser.add_argument('-e', help='Experiment name', default='', type=str)
parser.add_argument('-s', help="Setup name, eg. m7", default="", type=list)
parser.add_argument('-o', help="Output dir", default="", type=str)
parser.add_argument('-new', help="Output file with new files", default=None)

result = parser.parse_args()
ROOT_DIR = result.r
EXP_NAME = exp_name = result.e
SETUP_NAME = result.s
OUTPUT_DIR = result.o
try:
    config_file = result.cfile
    execfile(config_file)
except:
    print "Can't read config file"

exp_name = EXP_NAME

########## print config
print "exp_name:", exp_name
print "setup_name:", SETUP_NAME
print "config_file:", config_file
#########################

tmp_s = ROOT_DIR.split('/')
csv_basename = exp_name

exp_dirs = [ "%s/%s" % (ROOT_DIR, d) for d in os.listdir(ROOT_DIR) if os.path.isdir(ROOT_DIR + "/" + d) and exp_name in d ]

print "Found dirs"
print "\t\n".join(exp_dirs)

dirs = []
files = []
for d in exp_dirs:
    for s in SETUP_NAME:
        dd = [ "%s/%s" % (d, o) for o in os.listdir(d) if os.path.isdir(d + "/" + o) and s in o]
        dirs.extend(dd)

for d in dirs:
    for p in glob.glob(d + "/calc/*calc_*.csv"):
        if os.path.exists(p):
            files.append(p)
        else:
            files_list = [ d + "/" + f for f in os.listdir(d) if 'calc' in f and not os.path.isdir(d+"/"+f) ]
            if files_list is not []:
                files.extend(files_list)

print "Create dir", OUTPUT_DIR
try:
    os.mkdir(OUTPUT_DIR)
except:
    print OUTPUT_DIR, "exists"

output_new_file = None
if result.new:
    output_new_file = open(result.new,"w+")

for f in files:
    output_file = "_".join(f.replace(ROOT_DIR + "/", "").split("/"))
    output_path = OUTPUT_DIR + "/" + output_file
    if not os.path.exists(output_path):
        print f, "copy to", output_path
        shutil.copy(f, output_path)
        if output_new_file:
            output_new_file.writelines(output_path + "\n")
    else:
        print f, "skip"
