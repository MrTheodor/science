import re
import sys
import argparse as arg
from lib.Lattice import Lattice

import numpy as np

from __init__ import DummyMonomer


parser = arg.ArgumentParser(description="Calculate clusters")
parser.add_argument('-f', help="JMOL file", dest='file')
parser.add_argument('-tar', help='TAR file', dest='tar', action='store_true')
parser.add_argument('-scale', help='Rescale coordinates', default = -1, type=float)

parser.add_argument('-fs', help='For TAR/KC, from step', default=0, type=int)
parser.add_argument('-es', help='For TAR/KC, to step', default=0, type=int)
parser.add_argument('-fss', help='For TAR, step increment', default=100, type=int)
parser.add_argument('-output_file', type=str)

parser.add_argument('-box', help='Box size', required=True)
parser.add_argument('-jmol', help='Generate JMOL file', action='store_true')
parser.add_argument('-nr', help='Return number of clusters', action='store_true')
parser.add_argument('-step', help='Return step;number of clusters', action='store_true')
parser.add_argument('-re' )
parser.add_argument('-interface', help='Calculate interface', action='store_true')

parser.add_argument('-kch', help='Calculate from KyotoCabinet database', default="")
parser.add_argument('-d', help='JMOL dir', default='')

parser.add_argument('-jmol_list', help='JMOL list file')

parser.add_argument('-avg', action='store_true')

parser = parser.parse_args()

file = parser.file
scale = parser.scale if parser.scale > -1 else 1
box_size = map(int, parser.box.split(","))

Lattice.box_size = box_size

avg_dict = {}
total_count = 0

def process(file, step=None, output_file=None):

    if isinstance(file, str):
        f = open(file)
    else:
        f = file
    try:
        line = f.readlines()
        nr = int(line[0].replace("\n", ""))
    except:
        print line
        sys.exit(1)
    
    data = np.loadtxt(file, dtype=str, skiprows=1)
    
    #data = f.readlines()[1:]

    point_monomer = {}
    monomer_list = []
    position_set = set()

    if parser.scale == -1 and "." in data[0]:
        print "Atom position are float, please set -scale parameter to rescale them!"
        sys.exit(1)

    for line in data:
        #tmp = line.replace("\n","").split(" ")
        type = line[0]
        ll = map(float, line[1:])
        
        x = int(ll[0] * scale)
        y = int(ll[1] * scale)
        z = int(ll[2] * scale)

        position = (x, y, z)

        m = DummyMonomer(position, type, Lattice.get_neighbours(tuple(position)))

        point_monomer[position] = m
        monomer_list.append(m)
        position_set.add(position)

    return point_monomer, monomer_list, position_set

def create_valid_lists(list):
    point_monomer = {}
    monomer_list = []
    position_set = set()

    for line in list:
        type = line[0]
        ll = np.array(line[1:])*scale

        position = tuple(ll)

        m = DummyMonomer(position, type, Lattice.get_neighbours(tuple(position)))

        point_monomer[position] = m
        monomer_list.append(m)
        position_set.add(position)

    return point_monomer, monomer_list, position_set

def process_data(point_monomer, monomer_list, position_set, output_file=None, step=None):

    surface_bulk = {}
    clusters_list = {}

    interface_monomers = set()
    interface_monomers_count = 0

## build of neighbour list
    for m in monomer_list:
        m.nb_list = m.neighbours_list.intersection(position_set)
        m.nb_types = [ point_monomer[x].type for x in m.nb_list ]
        m.nb_count = len(m.nb_list)
        clusters_list[m.type] = []
        
        ## count how many neighbours is in given type
        if m.nb_types.count(m.type) != len(m.nb_types):
            m.on_interface = True
            if m not in interface_monomers:
                interface_monomers_count += 1
            interface_monomers.add(m)

## assign monomers to cluster
    for m in monomer_list:
        c_list = clusters_list[m.type]
        found_cluster = False
        for cluster in c_list:
            if cluster.intersection(m.nb_list):
                cluster.add(m.position)
                found_cluster = True

        if found_cluster == False:
            c_list.append(set([m.position]))

#for c in clusters_list.values():
#    print len(c)
#    print sum([ len(x) for x in c ])

## join clusters together

    taken = None

    def dfs(node, index):
        taken[index] = True
        ret = node
        for i, item in enumerate(l):
            if not taken[i] and not ret.isdisjoint(item):
                ret.update(dfs(item, i))

        return ret

    def merge_all(k):
        ret = []
        for i, node in enumerate(k):
            if not taken[i]:
                ret.append(list(dfs(node,i )))

        return ret
              

    for t in clusters_list:
        l = clusters_list[t]
        taken = [False]*len(l)
        clusters_list[t] = merge_all(l)

    if parser.jmol:
        for t, c in clusters_list.iteritems():
            for x in c:
                print len(x)
                print 
                for s in x:
                    print t,s[0]/2.0,s[1]/2.0,s[2]/2.0

    elif parser.nr:
        cluster_count = 0
        for t, c in clusters_list.iteritems():
            cluster_count += len(c)

        if output_file:
            output_file.write("%d\n" % cluster_count)
        print cluster_count

    elif parser.step:
        cluster_count = 0
        for t, c in clusters_list.iteritems():
            cluster_count += len(c)

        if parser.avg:
            if cluster_count in avg_dict:
                avg_dict[cluster_count] += 1
            else:
                avg_dict[cluster_count] = 0

            total_count += 1
        
        if step:
            if output_file:
                output_file.write("%d;%d;%d\n" % (step, cluster_count, interface_monomers_count))
            else:
                print "%d;%d;%d" % (step, cluster_count, interface_monomers_count)
        else:
            if output_file:
                output_file.write("%s;%d\n" % (re.match(r'.*jmol_(\d+)\.jmol', file).groups()[0], cluster_count))
            else:
                print "%s;%d" % (re.match(r'.*jmol_(\d+)\.jmol', file).groups()[0], cluster_count)


    elif parser.interface:
        inner_contacts = 0
        outer_contacts = 0
        cluster_count = 0
        for t,c in clusters_list.iteritems():
            cluster_count += len(c)

        print "%d;%d" % (cluster_count, interface_monomers_count)

if parser.tar:
    pattern = re.compile(r".*_(\d+).*$")
    import tarfile
    tarFile = tarfile.open(file)
    members = tarFile.getmembers()
    output_file = open(parser.output_file, "w+")
    if parser.fs > 0:
        step = parser.fs
        mm = 0
        for m in members:
           mt = int(re.match(pattern, m.path).groups()[0])
           if mt > step:
              f = tarFile.extractfile(m)
              point_monomer, monomer_list, position_set = process(f)
              process_data(point_monomer, monomer_list, position_set, step=mt, output_file = output_file)
           mm += 1
elif parser.file:
    point_monomer, monomer_list, position_set = process(parser.file)
    process_data(point_monomer, monomer_list, position_set)

elif parser.kch != "":
    import kyotocabinet as kb
    db = kb.DB()
    if not db.open(parser.kch, db.OREADER|db.ONOLOCK):
        print "ERROR"
        sys.exit(1)
    else:
        print "Open", parser.kch
    
    output_file = open(parser.output_file, "w+")

    jmols = [ (x, "jmol_%d" % x) for x in range(parser.fs, parser.es, parser.fss) ]
    for step, jmol in jmols:
            data = eval(db[jmol])
            if data:
                point_monomer, monomer_list, position_set = create_valid_lists(data)
                process_data(point_monomer, monomer_list, position_set, output_file = output_file, step = step)
                output_file.flush()
                #print "Process file", jmol

    output_file.close()
elif parser.d != "":
    import os
    files = os.listdir(parser.d)

    file_list = []
    idx = 0

    output_file = open(parser.output_file, "w+")
    jmols = [ (x, "%s/jmol_%d.jmol" % (parser.d, x)) for x in range(parser.fs, parser.es, parser.fss) ]
    idx = len(jmols)
    print "%d files" % idx
    for step, jmol in jmols:
        try:
            point_monomer, monomer_list, position_set = process(jmol)
            process_data(point_monomer, monomer_list, position_set, output_file = output_file, step = step)
            output_file.flush()
        except:
            continue
elif parser.jmol_list != "" and parser.kch != "":
    jmols = [ x.replace("\n", "") for x in open(parser.jmol_list, "r").readlines() ]
    jmols = [ (int(x.split("_")[1]), x) for x in jmols ]

    for step, jmol in jmols:
        try:
            data = eval(db[jmol])
            if data:
                point_monomer, monomer_list, position_set = create_valid_lists(data)
                process_data(point_monomer, monomer_list, position_set, output_file = output_file, step = step)
                output_file.flush()
        except:
            pass

else:
    point_monomer, monomer_list, position_set = create_valid_lists(file)
    process_data(point_monomer, monomer_list, position_set)
