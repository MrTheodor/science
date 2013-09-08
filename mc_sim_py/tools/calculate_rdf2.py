import sys
sys.path += ['/home/teodor/Documents/UAM/Studia magisterskie/Praca magisterska/Workspace/intraglobular/src']


from lib import tools
import cPickle, gzip
import os


from mpi4py import MPI

print "Calculate for set of experiments"
print "root_dir exp_name, setup_name, mc_steps, debug"

exp_name = sys.argv[2]
root_dir = sys.argv[1]
setup_name= sys.argv[3]
mc_steps = int(sys.argv[4])
debug = eval(sys.argv[5]) if len(sys.argv) == 6 else False

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print "Rank:", rank, "size", size

dirs = [ o for o in os.listdir(root_dir) if exp_name in o ]
d = len(dirs)/size

mpi_dir = dirs[d*rank:d*rank+d]

for dir in mpi_dir:
    if debug: print "read dir", dir
    box_dir = root_dir + "/" + dir + "/" + setup_name + "/data"
    out_file = open("%s_%s.rdf" % (dir, setup_name), "w+")
    file_list = [ o for o in os.listdir(box_dir) if int(o.replace('.dump','').split('_')[1]) >= mc_steps ]
    for box_file in file_list:
        if debug: print "  read file",box_file
        box = tools.load_file(box_dir + "/" +box_file)
    
        for c in box.chain_list:
            out_file.writelines(["%s\n" % ";".join(map(str, c.calculate_rdf2()))])

    out_file.close()

