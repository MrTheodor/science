import os
import sys
import numpy as np
from scipy import constants as const

import argparse as arg

import glob

from lib.tools import f_avg, f_avg2, f_cv2, avg_dict, save_dict_to_csv

parser = arg.ArgumentParser(description="Process files")
parser.add_argument('--config-file', help="Config file", dest="cfile")
parser.add_argument('-r', help='Root dir where experiment is stored', type=str, default='')
parser.add_argument('-er', help='Root dir where files are stored', type=str, default='')
parser.add_argument('-e', help='Experiment name', default='', type=str)
parser.add_argument('-s', help="Setup name, eg. m7", default="", type=list)
parser.add_argument('-f', help="From temp", default=0.0, type=float)
parser.add_argument('-t', help="To: temp", default=None, type=float)
parser.add_argument('-ac', help="Autocorrelation step", default=1, type=int)
parser.add_argument('-fs', help="From step", default=0, type=int)
parser.add_argument('-N', help="Number of monomers", type=int)
parser.add_argument('--exclude', action='append', help="Exclude temperature", default=[], type=list)
parser.add_argument('--csvplot', help="Create only cvs files for plot", default=False, action='store_true')
parser.add_argument('--csv', help="Save partial csv files for each temperature", default=False, action='store_true')
parser.add_argument('--csv-dir', help="Where to save partial csv files", default="", dest="csv_dir")

result = parser.parse_args()
root_dir = result.r
ex_root_dir = result.er
EXP_NAME = exp_name = result.e
SETUP_NAME = result.s
FROM_TEMP = result.f
TO_TEMP = result.t
EXCLUDED_TEMP = result.exclude
NUMBER_OF_MONOMERS = result.N
CREATE_CSV = result.csvplot
PARTIAL_CSV = result.csv
CSV_DIR = result.csv_dir
AC_STEP = result.ac
FROM_STEP = result.fs
try:
    config_file = result.cfile
    execfile(config_file)
except:
    print "Can't read config file"

if not CREATE_CSV:
    from matplotlib import pyplot as plt

exp_name = EXP_NAME

########## print config
print "root_dir:", root_dir
print "exp_name:", exp_name
print "setup_name:", SETUP_NAME
print "config_file:", config_file
print "FROM_TEMP:", FROM_TEMP
print "TO_TEMP:", TO_TEMP
print "FROM_STEP:", FROM_STEP
print "EXCLUDED_TEMP:", EXCLUDED_TEMP
print "AUTOCORRELATION:", AC_STEP
print "NUMBER_OF_MONOMERS:", NUMBER_OF_MONOMERS
print "CREATE_CSV:", CREATE_CSV
print "PARTIAL_CSV:", PARTIAL_CSV
#########################

### process files
data_e = {}
t_data_e = {}
data_rg = {}
t_data_rg = {}
t_data_rg_others = {}
data_e2e = {}
t_data_e2e = {}

tmp_s = root_dir.split('/')
csv_basename = exp_name
dirs = []
files = []


if not ex_root_dir:
    exp_dirs = [ "%s/%s" % (root_dir, d) for d in os.listdir(root_dir) if os.path.isdir(root_dir + "/" + d) and exp_name in d ]

    print "Found dirs"
    print "\t\n".join(exp_dirs)

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
else:
    for p in glob.glob(ex_root_dir + "/*.csv"):
        if os.path.exists(p):
            files.append(p)
        else:
            files_list = [ ex_root_dir + "/" + f for f in os.listdir(ex_root_dir) if 'calc' in f and not os.path.isdir(d + "/" + f) ]
            if files_list is not []:
                files.extend(files_list)

print "Found files"
print "\n".join(files)

for f in files:
    try:
        data = np.loadtxt(f, delimiter=';', skiprows=0)
#    data = open(f).readlines()[1:]
    ## half of data
    #data = data[len(data)/2:]
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
    except:
        print "Fail to process:", f


error_e = []
error_rg = []


## save pure csv to files temp_exp_name.csv
if PARTIAL_CSV:
    if not os.path.isdir(CSV_DIR):
       os.makedirs(CSV_DIR)

    csv_rg_file = "rg_" + exp_name + ".csv"
    csv_e_file = "e_" + exp_name + ".csv"
    temps = t_data_rg.keys()
    temps.sort()
    for temp in temps:
        rg_file = "%s/%f_%s" % (CSV_DIR, temp, csv_rg_file)
        e_file = "%s/%f_%s" % (CSV_DIR, temp, csv_e_file)
        save_dict_to_csv(t_data_rg[temp], rg_file)
        save_dict_to_csv(t_data_e[temp], e_file)


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

csv_e = open("e_" + csv_basename + ".csv", "w+")
csv_rg = open("rg_" + csv_basename + ".csv", "w+")
csv_cv = open("cv_" + csv_basename + ".csv", "w+")

if CREATE_CSV:
    csv_e.writelines([ "%f;%f\n" % (temp[x], e[x]) for x in range(len(temp)) ])
    csv_e.close()
    csv_rg.writelines([ "%f;%f\n" % (temp[x], rg[x]) for x in range(len(temp)) ])
    csv_rg.close()
    csv_cv.writelines([ "%f;%f\n" % (temp[x], cv[x]) for x in range(len(temp)) ])
    csv_cv.close()


#for idx in range(len(temp)):
#    print temp[idx], e[idx], error_e[idx], rg[idx], error_rg[idx]

#print "\n".join(map(str, temp[5:]))

#temp_ticks = [1.0, 3.0, 4.0, 5.0, 7.0, 15.0, 30.0]
temp_ticks = temp

if not CREATE_CSV:
    fg1 = plt.figure(1, figsize=(8.5, 5))
    fg1.subplots_adjust(hspace=0.5)
    plt.subplot(311)
    plt.xticks(temp_ticks)
    plt.title('E*/n vs T*')
    plt.xlabel('T*')
    plt.ylabel('E*/n')
#plt.xscale('log')
    plt.grid(b=True)
    #plt.errorbar(temp, map(lambda x: float(x)/NUMBER_OF_MONOMERS, e), yerr=error_e, figure=fg1)
#plt.errorbar(temp, e, yerr=error_e, figure=fg1)
    plt.plot(temp, map(lambda x: float(x), e), 'x-', figure=fg1)


#fg2 = plt.figure(2, figsize=(8.5, 5))
    plt.subplot(312)
    plt.xticks(temp_ticks)
    plt.title('Rg^2 vs T*')
    plt.xlabel('T*')
    plt.ylabel('Rg^2')
#plt.xscale('log')
    plt.grid(b=True)
    plt.errorbar(temp, rg, yerr=error_rg, figure=fg1)
#plt.plot(temp, map(lambda x: float(x), rg), 'x-', figure=fg2)

#fg3 = plt.figure(3, figsize=(8.5, 5))
    plt.subplot(313)
    plt.xticks(temp_ticks)
    plt.xlabel('T*')
    plt.title('Cv vs T*')
    plt.ylabel('Cv')
#plt.xscale('log')
    plt.grid(b=True)
#temp = [ t for t in temp if t <= 15.0 ]
    plt.plot(temp, cv[0:len(temp)], 'x-', figure=fg1)

    plt.show()

#plt.show()
#fg1.savefig('e_err_graph.png')
#fg2.savefig('rg_err_graph.png')
#fg3.savefig('cv_err_graph.png')

