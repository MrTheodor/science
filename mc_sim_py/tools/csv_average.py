"""
Average every rwo
"""

import sys
import numpy as np

print "file_in col"

file_in = sys.argv[1]
col = int(sys.argv[2])

f_data = [ map(float, x.split(";")) for x in open(file_in, "r+").readlines() ]
avg_data = [ [x[0], np.average(x[col:])] for x in f_data ]

out = open(file_in + "_avg", "w+")

