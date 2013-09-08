import sys
import kyotocabinet as kb

db = kb.DB()

db.open(sys.argv[1], db.OREADER)

m = db.match_prefix(sys.argv[2])
m.sort(key=lambda x: x.split('_')[1])
print m[-1]
