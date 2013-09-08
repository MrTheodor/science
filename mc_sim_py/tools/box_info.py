from lib.tools import load_file
from lib.Box import Box

import sys

print "Basic information in Box file"
print "<box_file_name>"


f_box = sys.argv[1]
box = load_file(f_box)

print "Box size", box.box_size
print "Number of chains", len(box.chain_list)
print "Local T", box.localT
