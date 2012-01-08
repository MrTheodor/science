import sys
import os

SEQ1 = "acdefghiklmnopqrstuvwy-x"
SEQ1_HP = "acfgilmpuvw"
SEQ1_P = "dehknqrsty"
SEQ1_ACID = "cdet"
SEQ1_BASIC = "hlr"
SEQ1_AROMATIC = "fhwy"
SEQ1_ALIPHATIC = "ilv"

F_GROUPS = [SEQ1_HP, SEQ1_P, SEQ1_ACID, SEQ1_BASIC, SEQ1_AROMATIC, SEQ1_ALIPHATIC]

seq1 = sys.argv[1].lower()
seq2 = sys.argv[2].lower()

m = len(seq1)
n = len(seq2)

## const
MATCH = 2
GROUP_MATCH = MATCH*2 ## if residue is in the same functional group
MISMATCH = -1
GAP = 0
##

