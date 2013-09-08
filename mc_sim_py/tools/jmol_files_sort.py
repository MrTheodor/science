import os
import sys
import shutil

import argparse as arg

from mpi4py import MPI

parser = arg.ArgumentParser(description="Sort jmol files for calculate structure factor")
parser.add_argument('-r', help="Experiment root dir")
parser.add_argument('-e', help="Experiment name")
parser.add_argument('-d', help='Destination root folder')
parser.add_argument('-c', help='Config file')
parser.add_argument('-eq', help='Equilibrium step', default=0, type=int)
parser.add_argument('-s', help='Autocorrelation step', default=1, type=int)

result = parser.parse_args()

root_dir = result.r
exp_name = result.e
dest_dir = result.d
from_step = result.eq
s = result.s

try:
    config_file = result.c
    execfile(config_file)
except:
    print "Can't read config file", config_file

print "Root dir", root_dir
print "Exp name", exp_name
print "Dest dir", dest_dir
print 'Equilibrium step', from_step
print 'AC step', s
print "Config file", config_file


if not os.path.isdir(dest_dir) and not os.path.exists(dest_dir):
    os.makedirs(dest_dir)

##
## Get list of calc_.csv files
##

exp_dirs = [ root_dir + "/" + d for d in os.listdir(root_dir) if os.path.isdir(root_dir + "/" + d) and exp_name in d ]
dirs = []
files = []
for d in exp_dirs:
    dd = [ d+"/"+o for o in os.listdir(d) if os.path.isdir(d + "/" + o) ]
    dirs.extend(dd)

for d in dirs:
    if os.path.exists(d + "/calc/calc_.csv"):
        files.append((d, d + "/calc/calc_.csv"))
    else:
        files_list = [ (d, d + "/" + f) for f in os.listdir(d) if 'calc' in f and not os.path.isdir(d+"/"+f) ]
        if files_list is not []:
            files.extend(files_list)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print "Total exp dirs", len(files)

if size > 0:
    files_len = len(files)
    chunk_len = files_len / size
    print "Rank", rank, "Size", size, "From", chunk_len*rank, "to", chunk_len*rank + chunk_len
    if rank == size -1:
        files = files[chunk_len*rank:]
    else:
        files = files[chunk_len*rank:chunk_len*rank + chunk_len]

temp_set = set()

for dir, calc_file in files:
    data = open(calc_file).readlines()[1:]
    for line in data:
        sl = line.split(";")
        temp = float(sl[1])
        temp_set.add(temp)
        step = float(sl[0])
        if step < from_step:
            continue
        if step % s != 0:
            continue
        
        tmp_exp_name = dir.split("/")[-1]
        
        jmol_path = "%s/models/jmol_%d.jmol" % (dir, step)
        output_dir = "%s/%f" % (dest_dir, temp)
        output_path = "%s/%s_jmol_%d.xyz" % (output_dir, tmp_exp_name, step)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        try:
            if not os.path.exists(output_path):
                shutil.copy(jmol_path, output_path)
                print "Copy", jmol_path,"to", output_path
            else:
                print "File exists", output_path
        except:
            print "File not found", jmol_path

print "Temp set", temp_set
