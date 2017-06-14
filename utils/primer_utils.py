'''
    ============================================================================
                GreenBox - Script Collection for Metagenomics
    ============================================================================

    Set of helper functions to predict chemical primer pair properties.

    ============================================================================
    Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
    ============================================================================
'''

import math
import logging


'''
    B = C or G or T
    D = A or G or T
    H = A or C or T
    K = G or T
    M = A or C
    N = A or C or G or T
    R = A or G
    S = C or G
    V = A or C or G
    W = A or T
    Y = C or T
'''
one_letter_encode = {1: 'A', 2: 'C', 4: 'G', 8: 'T', \
                    3: 'M', 5: 'R', 9: 'W', 6: 'S', 10: 'Y', 12: 'K', \
                    7: 'V', 11: 'H', 13: 'D', 14: 'B', 15: 'N'}

# Wallace rule to compute primer melting temperature
def primer_melt_wallace(primer):
    return 2*(sum([1 for nucl in primer if nucl in ['A', 'T']])) + 4*(sum([1 for nucl in primer if nucl in ['C', 'G']]))

# salt-adjusted method to compute primer melting temperature
# input primer:string sequence, Na:float molar Natrium ion concentration
def primer_melt_salt(primer, Na):
    primer = primer.lower()
    return 100.5 + 41.0*(sum([1 for nucl in primer if nucl in ['C', 'G']]))/len(primer) - 820.0/len(primer) + 16.6*math.log10(Na)

# variation measure for a block of aligned sequences in terms of number of columns
# that have more than one nt. A low score indicate higher conservation along all
# sequences.
# input sequences:[[]], left:int, right:int start and exclusive end position
def variation_score(aligned_sequences, left, right):
    num_cols = right - left
    transposed = [[aseq.seq[i] for aseq in aligned_sequences] for i in range(left, right)]
    #print transposed
    score = sum([len(set(col)) for col in transposed]) - num_cols
    #print score
    return score

# check if all target sequences satisfy melting temperature range and do not differ too much
def filter_melt(aligned_sequences, left, right, cfg):
    melt = [primer_melt_wallace(aseq.seq[left:right]) for aseq in aligned_sequences]
    max_melt, min_melt = max(melt), min(melt)
    if min_melt >= cfg.var['min_melt_temp'] and max_melt <= cfg.var['max_melt_temp'] and max_melt - min_melt <= cfg.var['max_melt_diff']:
        return True, min_melt, max_melt
    return False, min_melt, max_melt

# GC content in the primer should be between 40-60%
def filter_GC_content(aligned_sequences, pos, length):
    GC_min, GC_max = int(.4*length), int(.6*length)
    logging.debug('GC_min, max = [{}, {}]'.format(GC_min, GC_max))
    GC_cnts = [len([1 for nt in aseq.seq[pos:pos + length] if nt.upper() in ['C', 'G']]) for aseq in aligned_sequences]
    GC_cnts.sort()
    logging.debug('GC_ctns = ' + str(GC_cnts))
    if GC_cnts[0] < GC_min or GC_cnts[-1] > GC_max:
        return False, GC_cnts[0], GC_cnts[-1]
    return True, GC_cnts[0], GC_cnts[-1]

# one-letter encoding for set of aligned sequences, no gaps
def compress_helper(aligned_sequences, pos, length, bin_codes):
    seq_x = ''
    for i in range(length):
        code = reduce(lambda x, y: x|y, [bin_codes[aseq.seq[i]] for aseq in aligned_sequences])
        seq_x += one_letter_encode[code]
    return seq_x

def compress(aligned_sequences, pos, length):
    return compress_helper(aligned_sequences, pos, length, {'A': 1, 'C': 2, 'G': 4, 'T': 8})

def complement_compress(aligned_sequences, pos, length):
    return compress_helper(aligned_sequences, pos, length, {'A': 8, 'C': 4, 'G': 2, 'T': 1})[::-1]

if __name__ == '__main__':
    sequences = ['aaaaaa', 'aaaaaa', 'aaaaaa']
    variation_score(sequences, 2, 4)
