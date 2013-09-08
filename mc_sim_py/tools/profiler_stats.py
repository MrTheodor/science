import hotshot.stats
import sys

file_name = sys.argv[1]

stats = hotshot.stats.load(file_name)
stats.dump_stats(file_name+"_stats")
