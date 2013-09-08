import sys
import itertools

print "jmol_file box_size solvent_symbol"

input_jmol = sys.argv[1]
box_size = int(sys.argv[2])
symbol = sys.argv[3]

r = range(0, box_size)
coordinate = set([ k for k in itertools.product(r,r,r) if sum(k) % 2 == 0 ])
jmol_coordinate = set(map(lambda x: tuple(map(float, x.split(" ")[1:])), open(input_jmol, "r").readlines()[2:]))

diff = coordinate.difference(jmol_coordinate)

output_jmol = open(input_jmol, 'a+')

for cr in diff:
    line = "%s %f %f %f\n" % (symbol, cr[0], cr[1], cr[2])
    output_jmol.write(line)

output_jmol.close(0)
