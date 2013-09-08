## compute autocorrelation
## Jakub Krajniak <jkrajniak@gmail.com>
##

import sys
import numpy as np
from matplotlib import pyplot as plt

try:
    file_name = sys.argv[1]
    eq_steps = int(sys.argv[2])
    time = int(sys.argv[3])
    column = int(sys.argv[4])
except:
    print "Generate autocorrelation"
    print "file_name, mc_steps, ac_time, column, output file"
    print "\tmc_steps - where is the point of equilibrium"
    print "\tac_time - autocorrelation measure time"
    print "\tcolumn - which column need to be measure"
    output_file = None

try:
    output_file = sys.argv[5]
except:
    output_file = None

data = [ float(x.replace(",",".").split(';')[column]) for x in open(file_name).readlines() if int(x.replace(",",".").split(";")[0]) >= eq_steps ]


l_data = len(data)

var_data = np.var(data)
mean_data = np.mean(data)

ret = []

if output_file:
    out_f = open(output_file, 'w+')
for t in range(time):
    norm = 1.0/(l_data - t)
    
    tmp = norm*sum([ (data[x] - mean_data)*(data[x+t] - mean_data) for x in range(l_data - t)])
    ret.append(tmp)

ret = np.array(ret)
ret = ret / ret[0]
plt.plot(ret)
plt.show()
if output_file:
    out_f.writelines([ "%f\n" % x for x in ret ])
