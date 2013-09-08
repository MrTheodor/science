#!/usr/bin/python
import argparse as arg
import re
import kyotocabinet as kb
import sys

parser = arg.ArgumentParser(description="Dump KCH")
parser.add_argument('f', help='File')
parser.add_argument('-f', help='File')
parser.add_argument('-delimiter', help='Delimiter', default=';')
parser.add_argument('-o', help='Output file')
parser.add_argument('-csv', help='Output as csv', action='store_true')
parser.add_argument('-jmol', help='Output as jmol', action='store_true')
parser.add_argument('-d',  help='Output dir for files', default='.')
parser.add_argument('-k', help='Key prefix', default='')
parser.add_argument('-regex', help='Regexp prefix', default='')
parser.add_argument('-get', help='Get value', default='')
parser.add_argument('-c', help='Output on console', action="store_true")
parser.add_argument('-repair', help='Repair DB', action='store_true')

parser = parser.parse_args()

delimiter = parser.delimiter

db = kb.DB()
print "Open database %d" % db.open(parser.f, db.OREADER)

if parser.k != '':
    values = db.match_prefix(parser.k)
    print "match values: %d" % len(values)

elif parser.get != '':
    if parser.jmol:
        values = [parser.get]
    else:
        values = [db[parser.get]]
elif parser.regex != '':
    values = db.match_regex(parser.regex)

if parser.csv:
    output_file = sys.stdout if parser.c else open(parser.o, "w+") 
    idx = 0

    for v in values:
        line = "%s\n" % delimiter.join(map(str, eval(db[v])))
        output_file.writelines(line)
        idx += 1
    output_file.close()

    print "Stored %d values" % idx

elif parser.jmol:
    idx = 0

    for v in values:

        if parser.o:
            output_file = open(parser.o, "w+")
        elif parser.d:
            output_file = open(parser.d + "/" + v + ".jmol", "w+")

        t = eval(db[v])
        output_file.write(str(len(t)) + "\n\n");

        for z in t:
            line = "%s\n" % " ".join(map(str, z))
            output_file.writelines(line)

    	output_file.close()
        print "wrote", output_file

elif parser.repair:
	l = db.match_prefix('')
	print "Values: %d" % len(l)
	db.close()
    	output_file.close()
else:
    values.sort(key = lambda x: re.match(".*(\d+).*",x).groups()[0])
    print values[-1]

