# Implementation of Smith Waterman alghoritm for sequence local alignment
# 
# Jakub Krajniak <jkrajniak at gmail.com>
# 
# Based on http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#cite_note-Altshul1986-1
#
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
    matrix_group = sys.argv[3]
except:
    print "Sequence alignment using Smith-Waterman algorithm"
    print "Usage:"
    print "\tseq1 seq2 matrix"
    print "\tmatrix: simple|BLAST"
    print ""
    print "BLAST:"
    print "|".join(os.listdir(BLAST_DIR))

    sys.exit(1)


m = len(seq1)
n = len(seq2)

## const
MATCH = 2
GROUP_MATCH = MATCH*2 ## if residue is in the same functional group
MISMATCH = -4 ## penalty
##

def match(s1, s2):
    if s1 == s2:
        return MATCH
    else:
        for group in F_GROUPS:
            if s1 in group and s2 in group:
                return GROUP_MATCH
        return MISMATCH

H = np.zeros((m+1, n+1))
trace = np.zeros((m+1, n+1))

max_score = 0
max_point = (0,0) ## starting point of backtracking

for i in range(1, m+1):
    for j in range(1, m+1):
        s1 = H[i-1][j] + MISMATCH # up
        s2 = H[i][j-1] + MISMATCH # down
        s3 = H[i-j][j-1] + match(seq1[i-1], seq2[j-1]) # diag
        v = [0, s1, s2, s3]
        H[i][j] = max(v)
        trace[i][j] = v.index(H[i][j])
        if H[i][j] >= max_score:
            max_point = (i, j)
            max_score = H[i][j]

alseq1, alseq2 = [],[]

## Backtracking, follow path stored in :trace
i,j = max_point
while trace[i][j] != 0: ## indicates end of path
    point = trace[i][j]
    if point == 3: ## diagonal
        alseq1.append(seq1[i-1])
        alseq2.append(seq2[j-1])
        i, j = i-1, j-1
    elif point == 2: ## down
        alseq1.append('-')
        alseq2.append(seq2[j-1])
        j -= 1
    elif point == 1: ## up
        alseq1.append(seq1[i-1])
        alseq2.append('-')
        i -= 1

alseq1.reverse()
alseq2.reverse()

print "seq1:", seq1
print "seq2:", seq2
print "="*10
print "seq1:", "".join(alseq1)
print "seq2:", "".join(alseq2)

