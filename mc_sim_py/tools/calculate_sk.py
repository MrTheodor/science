# Calculate structure factor
# Based on average over time configurations

import os
import sys

import numpy as np
from matplotlib import pyplot as plt

import itertools

import cPickle

print "xyz_dir, output, box_size, segment"

xyz_dir = sys.argv[1]
output = sys.argv[2]
box_size = float(sys.argv[3])
segment = sys.argv[4]

files = [ "%s/%s" % (xyz_dir, o) for o in os.listdir(xyz_dir) ]

k_min = 2*np.pi/box_size
#k_max = 2*np.pi/np.math.sqrt(2)
k_max = 3.0

## read files
data_files = []
for f in files:
    d = open(f).readlines()[2:]
    data_files.append([ np.array(map(float, x.split(" ")[1:4])) for x in d if x.startswith(segment) ])

##
na = 32.0

k_range = np.arange(k_min, k_max, 0.5)
print "len k", len(k_range)
k_vec = [ np.array(x) for x in itertools.product(k_range, k_range, k_range) ]
print "len k_vec", len(k_vec)

##
Sk = []
for data_idx in range(len(data_files[0:1000])):
    sk = []
    data = data_files[data_idx]
    print "Process",files[data_idx]
    for k in k_vec:
        v1 = 0.0
        v2 = 0.0
        for rm in data:
            v1 += np.cos(k.dot(rm))
            v2 += np.sin(k.dot(rm))
        val = (v1**2) + (v2**2)
        sk.append(val)

    Sk.append(sk)

Sk = np.average(Sk, 0)

cPickle.dump(obj=Sk, file=open('sk_out', 'wb'), protocol=2)
