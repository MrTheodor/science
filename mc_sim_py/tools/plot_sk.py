import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
import os

import argparse as arg

rc('text', usetex=True)
rc('font', family='serif')

parser = arg.ArgumentParser(description="Plot SK graph")
parser.add_argument('-f', default='')
parser.add_argument('-legend', default='')
parser.add_argument('-ymove', default=10.0, type=float)
parser.add_argument('-d', default='.')
parser.add_argument('-title', default='')
parser.add_argument('-on_line', default=False, action='store_true')
parser.add_argument('-legend_prefix', default='')
parser.add_argument('-png', default='', type=str)
parser.add_argument('-xrange', default='', type=str)
parser.add_argument('-arrow', default='')

parser = parser.parse_args()

if parser.d != ".":
    files = os.listdir(parser.d)
else:
    files = [parser.f]
    parser.d = ""

files.sort()

data_list = []

for f in files:
    print f
    dir = parser.d + "/" if parser.d else ""
    if os.path.isfile(dir + f):
        data_list.append(np.loadtxt(dir + f))

x = data_list[0][:,0]
min_x = min(x)
max_x = max(x)
if parser.xrange:
    max_x = float(parser.xrange)
    plt.xlim(min(x), max_x)

plt.xlim( (min_x, max_x) )

legend = parser.legend.split(";") if ';' in parser.legend else []

#ll = [ "%s%s" % (parser.legend_prefix, x) for x in legend ]

#legend = ll

plot_list = []
scale_idx = 0

if legend:
    data_list = data_list[0:len(data_list)]

for d in data_list:
    y = d[:,1] * np.math.pow(10, scale_idx)
    plot_list.append(plt.semilogy(x, y))
    scale_idx += 1
    print scale_idx

plt.title(parser.title)

plt.xlabel('k')
plt.ylabel('S(k)')

idx = 0
if legend:
    if parser.on_line:
        for l in legend:
            plt.annotate(l, xy= ((max_x - min_x)/2, np.math.pow(11, idx) + 0.5))
            idx += 1
    else:
        plt.legend(plot_list[0:len(legend)], legend, 'upper right')

#if parser.arrow:
#    arrow_position = [ tuple([ float(z) for z in x.split(",")]) for x in parser.arrow.split(";") ]
#    for pos in arrow_position:



if parser.png:
    plt.savefig(parser.png)
else:
    plt.show()
