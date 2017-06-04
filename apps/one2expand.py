'''
    ============================================================================
                GreenBox - Script Collection for Metagenomics
    ============================================================================
    Convert 1 letter coded DNA sequence into a set of sequences permitting only
    standard letters 'A', 'C', 'G', 'T' in FASTA format. "name" plus counter 
    will be used as an accession number.
    
    ============================================================================
    Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
    ============================================================================

'''

import sys
import copy

d = {'B': ['C', 'G', 'T'], \
    'D': ['A', 'G', 'T'], \
    'H': ['A', 'C', 'T'], \
    'I': ['A','C','G','T'], \
    'K': ['G', 'T'], \
    'M': ['A', 'C'], \
    'N': ['A', 'C', 'G', 'T'], \
    'R': ['A', 'G'], \
    'S': ['C', 'G'], \
    'V': ['A', 'C', 'G'], \
    'W': ['A', 'T'], \
    'Y': ['C', 'T']}

def is_converted(seq):
    for key in d.keys():
        if seq.find(key) > -1:
            return False
    return True

def next_idx(seq):
    key_set = set(d.keys())
    for i in range(len(seq)):
        if seq[i] in key_set:
            return i
    return -1

def expand(seq):
    l = [seq]
    while is_converted(l[0]) is False:
        idx = next_idx(l[0])
        print idx
        l_aux = []
        for s in l:
            for letter in d[l[0][idx]]:
                s_new = copy.copy(s[:idx] + letter + s[idx+1:])
                l_aux.append(s_new)
        l = l_aux
    print l[:10]
    print(len(l))
    print(len(set(l)))
    return l

def writeFASTA(seqs, headers, filename):
    with open(filename, 'w') as f:
        for seq, header in zip(seqs, headers):
            f.write('>' + header + '\n' + seq + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage python one2four <name> <sequence> <outputfile>")
    else:
        seqs = expand(sys.argv[2])
        writeFASTA(seqs, [sys.argv[1] + str(i) for i in range(1, len(seqs)+1)], sys.argv[3])
