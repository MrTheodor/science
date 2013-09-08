import os
import sys
from mpi4py import MPI

root_dir = sys.argv[1]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dirs = [ os.path.realpath(root_dir + "/" + o) for o in os.listdir(root_dir) ]

dirs_len = len(dirs)
chunk_len = dirs_len / size

print "="*10
print "Total dirs", dirs_len
print "Chunk len", chunk_len
print "Rank", rank
print "Size", size
print "="*10

if rank == size - 1:
    dirs = dirs[chunk_len*rank:]
else:
    dirs =  dirs[chunk_len*rank: chunk_len*rank + chunk_len]

for d in dirs:
    os.chdir(d)
    try:
        print "Run calc_rdf in", d
        os.system("./calc_rdf")
    except:
        print "Problem in", d
