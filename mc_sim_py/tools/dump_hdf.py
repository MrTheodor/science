import sys
import tables
import numpy as np

import argparse as arg

parser = arg.ArgumentParser(description='dump')
parser.add_argument('-f', help='File')
parser.add_argument('-jmol', action='store_true', default=False)
parser.add_argument('-node', default='')
parser.add_argument('-csv', action='store_true', default=False)
parser.add_argument('-o')

parser = parser.parse_args()


h5file = tables.openFile(parser.f, 'r')
dataset = h5file.getNode(parser.node)

if parser.csv:
    output_file = open(parser.o, 'w+')
    output_file.writelines("%s\n" % ";".join(dataset.description._v_names))
    data = dataset.read()
    output_file.writelines([ "%s\n" % ";".join(map(str, x )) for x in data ])
    output_file.close()
    print "Wrote %d lines" % len(data)
    
