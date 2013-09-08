import sys
from lib import tools
import cPickle, gzip
import os

box_dir = sys.argv[1]
out_file = sys.argv[2]

out_file = open(out_file, 'w+')

for box_file in os.listdir(box_dir):
    box = tools.load_file(box_dir + "/" +box_file)
    
    for c in box.chain_list:
        out_file.writelines(["%s\n" % ";".join(map(str, c.calculate_rdf()))])

out_file.close()

