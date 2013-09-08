import os
import sys
import numpy as np
from scipy import constants as const

import argparse as arg

import glob

from lib.tools import f_avg, f_avg2, f_cv2, avg_dict, save_dict_to_csv

parser = arg.ArgumentParser(description="Process files")
parser.add_argument('--config-file', help="Config file", dest="cfile")
parser.add_argument('-file', help="Input file", type=str, dest="file")
parser.add_argument('-rdir', help="Result dir", type=str, default=".")
parser.add_argument('-f', help="From temp", default=0.0, type=float)
parser.add_argument('-t', help="To: temp", default=None, type=float)
parser.add_argument('-ac', help="Autocorrelation step", default=1, type=int)
parser.add_argument('-fs', help="From step", default=0, type=int)
parser.add_argument('-N', help="Number of monomers", type=int)
parser.add_argument('--exclude', action='append', help="Exclude temperature", default=[], type=list)
parser.add_argument('-setup', default=[''], type=list)
parser.add_argument('-prefix', default='', type=str)
parser.add_argument('-double', default=False, action='store_true')
parser.add_argument('-e2e', default=False, action='store_true')
parser.add_argument('-c', default=False, action='store_true')

result = parser.parse_args()
parser = result

FROM_TEMP = result.f
TO_TEMP = result.t
EXCLUDED_TEMP = result.exclude
NUMBER_OF_MONOMERS = result.N
AC_STEP = result.ac
FROM_STEP = result.fs
RESULT_DIR = result.rdir
SETUP_NAME = result.setup
PREFIX = result.prefix
FILES = {}

try:
    config_file = result.cfile
    execfile(config_file)
except:
    print "Can't read config file"

AC_STEP = AC_STEP if AC_STEP > 0 else 1


########## print config
print "config_file:", config_file
print "FROM_TEMP:", FROM_TEMP
print "TO_TEMP:", TO_TEMP
print "FROM_STEP:", FROM_STEP
print "EXCLUDED_TEMP:", EXCLUDED_TEMP
print "AUTOCORRELATION:", AC_STEP
print "NUMBER_OF_MONOMERS:", NUMBER_OF_MONOMERS
print "FILES:", FILES
#########################

### process files
data_e = {}
t_data_e = {}
data_rg = {}
t_data_rg = {}
t_data_rg_others = {}
data_e2e = {}
t_data_e2e = {}

dirs = []
files = []


print "Found files"
f = result.file

valid_file = False

for s in SETUP_NAME:
    if s in f:
        valid_file = True
        break

if not valid_file:
    sys.exit(1)

csv_basename = f.replace("/", "_")

if FILES:
    fk = FILES.keys()
    for x in fk:
        if x in f:
            FROM_STEP = FILES[x]['from_step']
            print "FROM_STEP for file %s: %d" % (f, FROM_STEP)

try:
    data = np.loadtxt(f, delimiter=';', skiprows=0)
    temp_list = set()
    for line in data:
        #sl = line.split(';')
        temp = float(line[1])
        temp_list.add(temp)
        step = float(line[0])

        if (temp <= FROM_TEMP) or \
            (TO_TEMP and temp > TO_TEMP) or \
            (EXCLUDED_TEMP and temp in EXCLUDED_TEMP):
            continue

        if not (step % AC_STEP  == 0) or step < FROM_STEP:
            continue

        e = float(line[-1])
        e2e = float(line[-2])**2
        rg = float(line[2])

        try:    
            if step in t_data_e[temp]:
                t_data_e[temp][step].append(e / float(2))
            else:
                t_data_e[temp][step] = [e / float(2)]
        except Exception, exception:
            t_data_e[temp] = {step: [e / float(2)]}
        
        try:
            if step in t_data_rg[temp]:
                t_data_rg[temp][step].append(rg)
            else:
                t_data_rg[temp][step] = [rg]
        except:
            t_data_rg[temp] = {step: [rg]}

        try:
            if step in t_data_e2e[temp]:
                t_data_e2e[temp][step].append(e2e)
            else:
                t_data_e2e[temp][step] = [e2e]
        except:
            t_data_e2e[temp] = {step: [e2e]}

    print "Process:", f, ",".join(map(str, temp_list))
    del data
except Exception, e:
    print "Fail to process:", f, e


error_e = []
error_rg = []


data_e = avg_dict(t_data_e)
data_rg = avg_dict(t_data_rg)
data_e2e = avg_dict(t_data_e2e)

for k,v in data_e.iteritems():
    del v[v.index(max(v))]
    del v[v.index(min(v))]
    data_e[k] = v
    error_e.append((k, np.std(v)/np.math.sqrt(len(v))))

for k,v in data_rg.iteritems():
    del v[v.index(max(v))]
    del v[v.index(min(v))]
    data_rg[k] = v
    error_rg.append((k, np.std(v)/np.math.sqrt(len(v))))

temp = data_e.keys()
temp.sort()

e = []
e2 = []
rg = []
e2e = []

for k,v in data_e.iteritems():
    print "E", k, len(v)
    e.append((k, f_avg(v)))
    e2.append((k, f_avg2(v)))

for k,v in data_rg.iteritems():
    rg.append((k, f_avg(v)))
    e2e.append((k, f_avg(data_e2e[k])))

## sort
e.sort(key = lambda x: x[0])
e2.sort(key = lambda x: x[0])
rg.sort(key = lambda x: x[0])
e2e.sort(key = lambda x: x[0])

error_e.sort(key = lambda x: x[0])
error_rg.sort(key = lambda x: x[0])

error_e = map(lambda x: x[1], error_e)
error_rg = map(lambda x: x[1], error_rg)

#print "\n".join(map(str, error_e))
#print "\n".join(map(str, error_rg))

e = map(lambda x: x[1], e)
e2 = map(lambda x: x[1], e2)
rg = map(lambda x: x[1], rg)
e2e = map(lambda x: x[1], e2e)
cv = f_cv2(data_e, NUMBER_OF_MONOMERS)

#cv = map(lambda x: x, f_cv(e, e2, temp))

try:
    os.mkdir(RESULT_DIR)
except:
    pass


e_all = RESULT_DIR + "/" + PREFIX + "e_all.csv"
rg_all = RESULT_DIR + "/" + PREFIX + "rg_all.csv"
cv_all = RESULT_DIR + "/" + PREFIX + "cv_all.csv"

if not parser.c:
    csv_e = open(e_all, "a+")
    csv_rg = open(rg_all, "a+")
    csv_cv = open(cv_all, "a+")
else:
    csv_e = sys.stdout
    csv_rg = sys.stdout
    csv_cv = sys.stdout
print "E"
csv_e.writelines([ "%f;%f\n" % (temp[x], e[x]) for x in range(len(temp)) ])
if not parser.c: csv_e.close()
print "RG"
csv_rg.writelines([ "%f;%f\n" % (temp[x], rg[x]) for x in range(len(temp)) ])
if not parser.c: csv_rg.close()
print "CV"
csv_cv.writelines([ "%f;%f\n" % (temp[x], cv[x]) for x in range(len(temp)) ])
if not parser.c: csv_cv.close()

print "Save:", e_all, cv_all, rg_all
