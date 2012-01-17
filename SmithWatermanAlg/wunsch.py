# Implementation of Needleman/Wunsch
# 
# Jakub Krajniak <jkrajniak at gmail.com>
# 
# 

import sys
import os
import numpy as np

## AA functional groups
SEQ1 = "acdefghiklmnopqrstuvwy-x" # all AA 
SEQ1_HP = "acfgilmpuvw" # hydrop
SEQ1_P = "dehknqrsty" # polar
SEQ1_ACID = "cdet" # acid AA
SEQ1_BASIC = "hlr" # basic AA
SEQ1_AROMATIC = "fhwy" # aromatic AA
SEQ1_ALIPHATIC = "ilv" # aliphatic AA

F_GROUPS = [SEQ1_HP, SEQ1_P, SEQ1_ACID, SEQ1_BASIC, SEQ1_AROMATIC, SEQ1_ALIPHATIC]

BLAST_DIR="blast_matrices"

try:
    seq1 = sys.argv[1].lower()
    seq2 = sys.argv[2].lower()
except:
    print "Sequence alignment using Needleman/Wunsch technique"
    print "Usage:"
    print "\tseq1 seq2 matrix"
    print ""
    
    sys.exit(1)


m = len(seq1)
n = len(seq2)

## const
MATCH = 2
GROUP_MATCH = 1 ## if residue is in the same functional group
GAP = -2 ## penalty
MISMATCH = -1
##

H = np.zeros((m+1, n+1))
trace = np.zeros((m+1, n+1))
bifurcation = np.zeros((m+1, n+1))

def match(s1, s2):
    if s1 == s2:
        return MATCH
    else:
        #for group in F_GROUPS:
        #    if s1 in group and s2 in group:
        #        return GROUP_MATCH
        return MISMATCH

def print_matrix(row, col, matrix):
    print "      "," ".join(map(lambda x: str(x).center(6), row))
    for lidx in range(len(matrix)):
        print col[lidx-1]," ".join(map(lambda x: str(x).center(6), matrix[lidx]))

max_point = (0,0)
max_score = 0

for i in range(1, m+1):
    for j in range(1, n+1):
        s1 = H[i-1][j-1] + match(seq1[i-1], seq2[j-1])
        s2 = H[i][j-1] + GAP ## insert 
        s3 = H[i-1][j] + GAP ## delete
        v = [s1, s2, s3]
        max_v = max(v)
        H[i][j] = max_v
        trace[i][j] = v.index(max_v)
        #bifurcation[i][j] = set([ v.index(max_v), v[1:].index(max_v), v[2:].index(max_v)  ])
        bifurcation[i][j] = max([ v.count(max_v) for x in v ])
        if H[i][j] >= max_score:
            max_point = (i, j)
            max_score = H[i][j]

print "="*10
print "H matrix"
H = H.transpose()
print_matrix(seq1, seq2,H)
print "-"*10

print "Trace matrix"
print "0 - diagonal, 1.0 - top, 2.0 - left"
trace = trace.transpose()
print_matrix(seq1, seq2, trace)
print "-"*10

print "Bifurcation"
print_matrix(seq1, seq2, bifurcation.transpose())

## trace
i = m
j = n
H = H.transpose()
alignmentA = ''
alignmentB = ''

scoring = 0
path = []

trace = trace.transpose()

while (i > 0 and j > 0):
    t = trace[i][j]
    s1 = seq1[i-1]
    s2 = seq2[j-1]
    path.append((i,j))
    if t == 0.0:
        alignmentA = s1 + alignmentA
        alignmentB = s2 + alignmentB
        i -= 1
        j -= 1
    elif t == 2.0: ## delete in seq2
        alignmentA = s1 + alignmentA
        alignmentB = '-' + alignmentB
        i -= 1
    else: ## inser in seq2
        alignmentA = '-' + alignmentA
        alignmentB = s2 + alignmentB
        j -= 1
    scoring += H[i][j]

while i > 0:
    alignmentA = s1 + alignmentA
    alignmentB = '-' +  alignmentB
    i -= 1

while j > 0:
    alignmentA = '-' + alignmentA
    alignmentB = seq2[j] + alignmentB
    j -= 1

alignmentA = list(alignmentA)
alignmentB = list(alignmentB)

#alignmentA.reverse()
#alignmentB.reverse()

sco = []

for idx in range(len(alignmentA)):
    s1 = alignmentA[idx]
    s2 = alignmentB[idx]
    for f in F_GROUPS:
        if s1 in f and s2 in f:
            sco.append('+')
            break
        else:
            sco.append('-')
            break

print "A:", "".join(alignmentA)
print "B:", "".join(alignmentB)
print 'M:', "".join(sco)
print "Scoring:", scoring
print "Path:", path
