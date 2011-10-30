import sys
import itertools

def bcc(n):
    product = itertools.product(range(n), range(n), range(n))
    return [ c for c in product if (c[0] + c[1] + c[2]) % 2 != 0 ]

def fcc(n):
    product = itertools.product(range(n), range(n), range(n))
    
    return [ c for c in product if (c[0] + c[1] + c[2]) % 2 == 0 ]

def sc(n):
    product = itertools.product(range(n), range(n), range(n))

    return list(product)
    


type = sys.argv[1]
n = int(sys.argv[2])
file = open(sys.argv[3], "w+")
scale = float(sys.argv[4]) if len(sys.argv) > 4 else 1

select = {"sc": sc, "bcc": bcc, "fcc": fcc }

data = select[type](n)

file.writelines("%d\n\n" % len(data))
file.writelines([ "C %s\n" % (" ".join(map(str, map(lambda x: x*scale, d)))) for d in data ])
file.close()

