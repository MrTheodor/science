SEQ1 = "acdefghiklmnopqrstuvwy-x"
SEQ1_HP = "acfgilmpuvw"
SEQ1_P = "dehknqrsty"
SEQ1_ACID = "cdet"
SEQ1_BASIC = "hlr"
SEQ1_AROMATIC = "fhwy"
SEQ1_ALIPHATIC = "ilv"

SEQ3 = ['Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile', 'Lys',
        'Leu', 'Met', 'Asn', 'Pyl', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr',
        'Sec', 'Val', 'Trp', 'Tyr', '---', 'UNK']

seq1to3 = dict(zip(SEQ1, SEQ3))
seq3to1 = dict(zip(SEQ3, SEQ1))

def chunk(l, n):
    """
    Return chunks of string
    """
    return (l[i:i+n] for i in xrange(0, len(l), n)) 

def seq1_3(seq1):
    seq1 = seq1.lower()
    return [ seq1to3[x] if seq1to3.has_key(x) else seq1to3['x'] for x in seq1 ]

def seq3_1(seq3):
    if not isinstance(seq3, list):
        if len(seq3) % 3: raise Exception('It is not valid 3-letter sequence')
        seq3 = [ x.title() for x in chunk(seq3, 3) ]

    return [ seq3to1[x] if seq3to1.has_key(x) else seq3to1['UNK'] for x in seq3 ]


def prepare_sequence(sequence_code, library):
    """
    Return valid Monomer object sequence, base on one letere sequence code
    
    @string sequence_code: one letter sequence
    @Monomer library: module with Monomers class that are used for build sequence

    return @list
    """
    monomer_sequence = []
    for s in sequence_code:
        m_obj = getattr(library, 'Type%s' % s)
        monomer_sequence.append(m_obj())

    return monomer_sequence

class PDBFile(object):
    """
    Object represents PDB file
    """
    pass
