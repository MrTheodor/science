import os
import sys
from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np

rc('text', usetex=True)
rc('font', family='serif')

print "Plot (x,y) some column"
print "file.csv column"

import argparse as arg

parser = arg.ArgumentParser(description="Plot column")
parser.add_argument('file', help='File')
parser.add_argument('-f', help='File', dest='file')
parser.add_argument('-n', help='Column', dest='column', type=int)
parser.add_argument('-multi', help='Multi variable', action='store_true')
parser.add_argument('-d', help='Delimeter', type=str, default=';')
parser.add_argument('-xticks', help='Use first column as xticks', action='store_true')
parser.add_argument('-yticks', help='Use first column as yticks', action='store_true')
parser.add_argument('-title', help='Title', default='')
parser.add_argument('-xlabel', default='')
parser.add_argument('-ylabel', default='')
parser.add_argument('-legend', default='')
parser.add_argument('-s', default='-x', type=str)
parser.add_argument('-yscale', default=1.0, type=float)
parser.add_argument('-png', default='', type=str)

parser = parser.parse_args()

file = parser.file
column = parser.column
multi = parser.multi
delimiter = parser.d

if "," in file:
    files = file.split(",")
else:
    files = [file]

data = []

for f in files:
    data.append(np.loadtxt(f, delimiter=delimiter))

#lines.sort(0)

def process_lines(lines):
    print "Multi"
    out = {}
    for k,l in lines:
        if k in out:
            out[k].append(l)
        else:
            out[k] = [l]


    x_val = []
    y_val = []
    
    x_val = out.keys()
    x_val.sort()

    for k in x_val:
        y_val.append(np.average(out[k]))

    return (x_val, y_val)

    #max_idx = y_val.index(max(y_val))
    #del y_val[max_idx]
    #del x_val[max_idx]
    #max_idx = y_val.index(max(y_val))
    #del y_val[max_idx]
    #del x_val[max_idx]

plot_list = []

for d in data:
    lines = d[:,[0, column]]
    x_val, y_val = process_lines(lines)
    
    if parser.yscale != 1.0:
        y_val = map(lambda x: x / parser.yscale, y_val)

    if parser.xticks:
        plt.xticks(x_val, ["%.2f" % x for x in x_val ])
    if parser.yticks:
        plt.yticks(y_val, ["%.2f" % y for y in y_val ])
    
    plot_list.append(plt.plot(x_val, y_val, parser.s))


plt.title(parser.title)
plt.xlabel(parser.xlabel)
plt.ylabel(parser.ylabel)

if parser.legend and ";" in parser.legend:
    l = parser.legend.split(';')
    plt.legend(plot_list, l)

if parser.png:
    plt.savefig(parser.png)
else:
    plt.show()

#open(file+"_out", "w+").writelines(map(lambda x: "%s\n" % str(x), lines))
