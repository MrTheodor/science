import numpy as np
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from scipy.interpolate import interp1d

import os

from matplotlib import rc
from matplotlib import pyplot as plt

import argparse as arg


rc('text', usetex=True)
rc('font', family='serif')

parser = arg.ArgumentParser(description="Plot SK graph")
parser.add_argument('-f', default='')
parser.add_argument('-legend', default='')
parser.add_argument('-ymove', default=10.0, type=float)
parser.add_argument('-title', default='')
parser.add_argument('-on_line', default=False, action='store_true')
parser.add_argument('-legend_prefix', default='')
parser.add_argument('-png', default='', type=str)

parser = parser.parse_args()

if "," in parser.f:
    files = parser.f.split(",")

files.sort()

data_list = []
for f in files:
    print f
    if os.path.isfile(f):
        data_list.append(np.loadtxt(f, delimiter=";", skiprows=1))

x = data_list[0][:,0]

legend = parser.legend.split(";") if ";" in parser.legend else []
titles = parser.title.split(";")

plot_list = []
plt_idx = 1
xi = np.linspace(min(x), max(x), 100)

for d in data_list:
    plt.subplot(len(data_list)/2, 2, plt_idx)
    plt.xticks(x, [ "%.2f" % k for k in x ])
    plt.xlim( min(x), max(x) )
    plt.ylim(0.0, 1.0)
    for z in range(1,d.shape[1]):
        #yi = InterpolatedUnivariateSpline(x, d[:,z])(xi)
        y = d[:,z]
        f = interp1d(x, y, kind='cubic')
        plt.plot(xi, f(xi))
    plt_idx += 1

plt.show()
