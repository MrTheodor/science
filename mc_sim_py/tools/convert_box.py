from lib import Box, Chain, tools
from lib import tools

import os
import sys

try:
    box_file = sys.argv[2]
    xyz_file = sys.argv[3]
    
    cmd = sys.argv[1]
except:
    print "Convert box to xyz and vs."
    print "box2xyz|xyz2box box_file xyz_file"
    print "For xyz2box also: (box_size) and <localT>"


if cmd is "xyz2box":
    box_size = eval(sys.argv[4])
    localT = float(sys.argv[5])
    xyz = open(os.path.realpath(xyz_file)).readlines()[2:]


elif cmd is "box2xyz":
    b = tools.load_file(box_file)
    b.jmol_snapshot(xyz_file)
