import sys
import numpy as np

input_file = open(sys.argv[1], "r")

in_lines = input_file.readlines()
last_symbol = ''

return_dict = {}

table = []
tmp = []
for l in in_lines:
    try:
        value = float(l)
        tmp.append(value)
    except:
        key = l.replace("\n", "")
        if key != "":
           if key not in return_dict:
               return_dict[key] = []
               if tmp != []:
                   table.append(tmp[:])
                   tmp = [key]
               else:
                   tmp = [key]
            
           last_symbol = key

table.append(tmp[:])


keys = []
for l in table:
    keys.append(l[0])

keys.insert(0, '')
table.insert(0, keys)

tmp_table = []
for x in range(1, len(table)):
    tmp_table.append(table[x][1:])
tab = np.array(tmp_table)

ta = []

for row in range(len(tab)):
    tt = []
    tt.extend(list(tab[:,row][0:row]))
    tt.extend(list(tab[row][row:]))
    ta.append(tt)

np.savetxt('output_mj', ta, delimiter=';', fmt="%.4f")

for l in ta:
    row = zip(keys[1:], l)
    print ",".join(map(lambda x: "'%s': %.4f" % x, row))
