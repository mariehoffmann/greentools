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


one_letter_decode = "\
    B = C or G or T\n\
    D = A or G or T\n\
    H = A or C or T\n\
    K = G or T\n\
    M = A or C\n\
    N = A or C or G or T\n\
    R = A or G\n\
    S = C or G\n\
    V = A or C or G\n\
    W = A or T\n\
    Y = C or T"


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
# having different nucleotides. A low score indicates high conservation.
# input sequences:[[]], pos:int, offset:int
def ambiguous_positions(aligned_sequences, pos, offset):
    transposed = [[aseq.seq[i] for aseq in aligned_sequences] for i in range(pos, pos+offset)]
    #print transposed
    score = sum([1 if len(set(col)) > 1 else 0 for col in transposed])
    #print score
    return score

# check if all target sequences satisfy melting temperature range and do not differ too much
def filter_melt(aligned_sequences, pos, offset, cfg):
    melt = [primer_melt_wallace(aseq.seq[pos:pos+offset]) for aseq in aligned_sequences]
    max_melt, min_melt = max(melt), min(melt)
    if min_melt >= cfg.var['melt_temp'][0] and max_melt <= cfg.var['melt_temp'][1] and max_melt - min_melt <= cfg.var['max_melt_diff']:
        return True, min_melt, max_melt
    return False, min_melt, max_melt

# GC content in the primer should be between 40-60%
def filter_GC_content(aligned_sequences, pos, offset, cfg):
    GC_min, GC_max = int(cfg.var['gc_content'][0]*offset), int(cfg.var['gc_content'][1]*offset)
    logging.debug('GC_min, max = [{}, {}]'.format(GC_min, GC_max))
    GC_cnts = [len([1 for nt in aseq.seq[pos:pos+offset] if nt.upper() in ['C', 'G']]) for aseq in aligned_sequences]
    GC_cnts.sort()
    logging.debug('GC_ctns = ' + str(GC_cnts))
    if GC_cnts[0] < GC_min or GC_cnts[-1] > GC_max:
        return False, GC_cnts[0], GC_cnts[-1]
    return True, GC_cnts[0], GC_cnts[-1]

# check for GC at 3' end, DNA sense/'+': 5' to 3', antisense/'-': 3' to 5', should be <= 3 in last 5 bps
def filter_GC_clamp(sequence, sense='+'):
    if sense == '+':
        gc = len([1 for nt in sequence[-5:] if nt in ['C', 'G']])
    else:
        gc = len([1 for nt in sequence[:5] if nt in ['C', 'G']])
    if gc > 3:
        return False, gc
    return True, gc

# check for 2ndary structure hairpin, may only be present at 3' end with a delta(G) = -2 kcal/mol, or internally with a delta(G) of -3 kcal/mol
# TODO: upper limit for loop length is disrespected currently
def filter_hairpin(seq, cfg):
    n, min_loop_len = len(seq), int(cfg.var['hairpin_loop_len'][0])
    palindrome_len_rng = range(3, len(seq)/2 - min_loop_len + 1)
    seq_ci = complement(seq)[::-1] # inverted, complemented sequence
    for i in range(len(seq) - 2*palindrome_len_rng[0] - min_loop_len):
        for m in palindrome_len_rng:
            for i_inv in range(n - 2*m - min_loop_len):
                if seq[i:i+m] == seq_ci[i_inv:i_inv+m]:
                    #print seq[i:i+m], ' == ', seq_ci[i_inv:i_inv+m]
                    return True
    return False

# translate into complementary string without reversing
def complement(dna_sequence):
    m = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([m[nt] for nt in dna_sequence])


# one-letter encoding for set of aligned sequences, no gaps
def compress_helper(aligned_sequences, pos, length, bin_codes):
    seq_x = ''
    for i in range(pos, pos+length):
        code = reduce(lambda x, y: x|y, [bin_codes[aseq.seq[i]] for aseq in aligned_sequences])
        seq_x += one_letter_encode[code]
    return seq_x

def compress(aligned_sequences, pos, length):
    return compress_helper(aligned_sequences, pos, length, {'A': 1, 'C': 2, 'G': 4, 'T': 8})

def complement_compress(aligned_sequences, pos, length):
    return compress_helper(aligned_sequences, pos, length, {'A': 8, 'C': 4, 'G': 2, 'T': 1})[::-1]



if __name__ == '__main__':
    class Config(object):
        def __init__(self):
            self.var = {'hairpin_loop_range': [3.0, 8.0]}

    cfg = Config()
    print(filter_hairpin('AAACCTTTGGAA', cfg))
    print(filter_hairpin('AAACCTTTGGTA', cfg))
