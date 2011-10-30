##
## Author: Jakub Krajniak
## 
## Compute autocorrelation function
## 
import sys

data = [ int(d) for d in open(sys.argv[1]) ]
data_norm = [ d/float(len(data)) for d in data ]

data_avg = sum(data_norm)/float(len(data_norm))

time = range(16)

chi = []
K = len(data_norm)

for t in time:
    print t
    normalization_factor = 1 / (K - t)
    for s in range(K-t)

