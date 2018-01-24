'''
    ============================================================================
        GreenBox - Script Collection for Metagenomics
    ============================================================================
        Task: barcode search (brcd) with low complex primers

        Compute good primer pairs for a set of aligned sequences such that species resolution
        is optimal. In other words, search a region that serves as a marker and for which
        sound primers can be found. For setting path to aligned sequence file, primer pair
        melting temperatures, etc. use brcd.cfg.
        To produce aligned sequences you can use PAGAN:
            wasabiapp.org/software/pagan/phylogenetic_multiple_alignment

        Call:   python brcd_search.py <path_to/brcd_search.cfg>

        Output: csv file with forward (fw) and reverse (rv) primer start positions and their
                sequences and a score based on how well the design constraints match. Format:
                [fw_pos, fw_len, fw_primer_seq, fw_pos, fw_len, fw_primer_seq, score]

        Algorithm:
                1. Identify regions that allow a separation corresponding to
                given labels, i.e. markers.
                2. Search for more highly conserved blocks surrounding the markers.

    ============================================================================
        Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
    ============================================================================
'''

import sys
import os
import logging

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'utils'))
import config_reader as config
import primer_utils

theta = 0.18 # procentual difference between 2 items in same cluster

# storage object for aligned sequences
class AlignedSequence(object):
    def __init__(self, _seq, _label, _acc, _header):
        self.seq = _seq
        self.label = _label
        self.acc = _acc     # unique
        self.header = _header
        self.seq_print = self.seq if len(_seq) <= 70 else self.seq[:50] + ' ... ' + self.seq[-20:]

    def write(self):
        print self.acc + ', ' +  self.label + ', ' + self.seq_print

    def write_site(self, site_obj):
        print 'accession\tlabel\tsubsequence\n' + self.acc + '\t' +  self.label + '\t' + self.seq[site.pos:site.pos+site.len]



# store potential primer site by start index, length, 'conservation' score, melting temperature range
class Site(object):
    def __init__(self, _pos, _len, _score, _melt_rng):
        self.pos, self.len, self.score, self.melt_range = _pos, _len, _score, _melt_rng
    def write(self):
        print 'start_pos\tlength\tscore\tmelt_range\n' + '\t'.join([self.pos, self.len, self.score, str(self.melt_range)])

# test for gap freedom in transcript block
def is_gapfree(transcripts, i, j):
    if cfg.var['gap_symbol'] in set([symbol for transcript in transcripts for symbol in transcript.seq[i:j]]):
        return False
    return True

# read aligned sequences and return dict {accession: AlignedSequence}
def aligned_sequence_reader(fasta_aligned):
    print cfg.__dict__
    d = {}  # {header: seq_obj}
    with open(fasta_aligned, 'r') as f:
        acc = None
        label = None
        seq = ''
        header = ''
        for line in f.readlines():
            line = line.strip()
            if line.startswith('>'):
                if acc is not None: # write previously parsed entry
                    d[acc] = AlignedSequence(seq, label, acc, header)
                # reset header, sequence, label, accession
                header = line[1:]
                seq = ''
                acc = line[1:line.find('_')]
                for l in cfg.var['labels']:
                    if line.find(l) > -1:
                        label = l
                        print "label = ", label
                        break
            else:
                seq += line.strip()
        # write last parsed entry
        if seq != '':
            d[acc] = AlignedSequence(seq, label, acc, header)

        keys = sorted(d.keys())
        print 'num keys = ' + str(len(keys))
        print '0..199'
        print '\n'.join([key + ': ' + d[key].seq[:200] for key in keys])
        print '200..399'
        print '\n'.join([key + ': ' + d[key].seq[200:400] for key in keys])
        print '400..599'
        print '\n'.join([key + ': ' + d[key].seq[400:600] for key in keys])
        print '600..799'
        print '\n'.join([key + ': ' + d[key].seq[600:800] for key in keys])
    print d.keys()
    return d.values()

# compute best blocks of fulfilling primer constraints (see brcd.cfg)
def find_primer_sites(aligned_sequences):
    n = len(aligned_sequences[0].seq) # total sequence length
    pr_min, pr_max = cfg.var['min_primer_len'], cfg.var['max_primer_len']
    print type(pr_min)
    print 'primer_len_interval = [' + str(pr_min) + ', ' + str(pr_max) + ']'
    sites = []
    i = 0
    while i < n - pr_min:
        print 'i = ', i, ', min block is gap-free: ', is_gapfree(aligned_sequences, i, i+pr_min)
        # block has to be gap-free
        # TODO: skip more than one nt
        while i < n - pr_min and is_gapfree(aligned_sequences, i, i+pr_min) == False:
            gaps_right_most = [aseq.seq[i:i+pr_min].rfind('-') for aseq in aligned_sequences]
            print 'i = ', i, ', min block is not gap-free: ', is_gapfree(aligned_sequences, i, i+pr_min), ', skip by ',max(gaps_right_most) + 1
            i += max(gaps_right_most) + 1
        if i >= n - pr_min:
            break
        for j in range(pr_min, pr_max):
            if i + j >= n:
                break
            # test for gaps beyond min primer length
            if j >= pr_min:
                if is_gapfree(aligned_sequences, i+pr_min, i+pr_max) == False:
                    i += pr_min-1
                    break
            print 'potential site at [', i, ', ', i+j, ']'
            score = primer_utils.variation_score(aligned_sequences, i, j)

            check_melt, min_melt, max_melt = primer_utils.filter_melt(aligned_sequences, i, i+j, cfg)
            print 'score = ', score, ', check_melt = ', check_melt, ', min_melt = ', min_melt, ', max_melt = ', max_melt
            if check_melt == True:
                sites.append(Site(i-1, j, score, (min_melt, max_melt)))
        i += 1

    # filter top k
    '''
    min_k = sorted([site[1] for score in sites], reverse=True)[top_k-1]
    print "min k: ", min_k
    sites = [Site(site[0], site[1]) for site in sites if site[1] >= min_k]
    for idx, site in enumerate(sites):
        print "\nblock #" + str(idx) + ':'
        for val in d.values():
            print val.seq[site.pos:site.pos+primer_len+1]
    #sys.exit(0)
    '''

    for site in sites:
        site.write()
        for seq_obj in aligned_sequences:
            seq_obj.write_site(site)
    return sites

# 2. compute transcript clusters for all primer combinations
# distance as number of mismatches
def distance(s1, s2):
    return sum([1 if a != b else 0 for (a,b) in zip(s1, s2)])

# distance between 2 clusters is smallest distance between 2 individuals
def single_linkage(c1, c2, left, right):
    seq_set1 = set([t_obj.seq[left:right] for t_obj in c1])
    seq_set2 = set([t_obj.seq[left:right] for t_obj in c2])
    return min([distance(s1, s2) for s1 in seq_set1 for s2 in seq_set2])

# seq_obj_list: [AlignedSequence], posx: transcript start/end position
def agglomerative_clustering(seq_obj_list, left, right):
    threshold = int(theta * len(seq_obj_list[0].seq[left:right]))
    clusters = [[seq_obj] for seq_obj in seq_obj_list]
    min_dist_pair = [(single_linkage(c1, c2, left, right), (i,j+i+1)) for i, c1 in enumerate(clusters[:-1]) for j, c2 in enumerate(clusters[i+1:])]
    min_dist_pair.sort(key=lambda item: item[0])
    while min_dist_pair[0][0] <= threshold and len(clusters) > 1:
        i, j = min_dist_pair[0][1][0], min_dist_pair[0][1][1]
        # merge cluster i and j
        cluster_merged = clusters[i] + clusters[j]
        del clusters[j]
        del clusters[i]
        clusters.append(cluster_merged)
        if len(clusters) < 2:
            return clusters
        # compute new smallest distance
        min_dist_pair = [(single_linkage(c1, c2, left, right), (i,j+i+1)) for i, c1 in enumerate(clusters[:-1]) for j, c2 in enumerate(clusters[i+1:])]
        min_dist_pair.sort(key=lambda item: item[0])
    return clusters

# classify transcripts for all site combinations and select if discrimination is in same manner like the target one
def combine_sites(sites):
    cluster_select = []

    for i, fw in enumerate(sites[:-1]):
        for rv in sites[i+1:]:
            if rv[0] - 20 - fw[0] in range(cfg.min_amplicon_len, cfg.max_amplicon_len + 1) and is_gapfree(d.values(), fw[0]+primer_len, rv[0]) is True:
                left, right = fw[0] + 20, rv[0]
                #transcripts = [Transcript(d[key]val[0][pos_start:pos_end], val[1], val[2]) for key in d.keys()]
                clusters = agglomerative_clustering(d.values(), left, right)
                label_identities = [0 if len(set([transcript.label for transcript in cluster])) == 1 else 1 for cluster in clusters]
                if sum(label_identities) <= cfg.max_num_cluster_mixed and len(clusters) in range(cfg.min_num_cluster, cfg.max_num_cluster+1):  #  0.1*len(label_identities)
                    print 'discr comb found!'
                    cluster_select.append((fw, rv))
                    for cluster in clusters:
                        seqs = set()
                        print 'number of classes: ', len(cluster)
                        for transcript in cluster:
                            transcript.write_block(left, right)
                        print '-------------------------------'
    print 'discriminative clusters could be produced with ', cluster_select
    return cluster_select

def write_output(clusters):
    pass

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Wrong number of arguments: python brcd.py <path_to/brcd.cfg>"
        sys.exit(0)
    logging.basicConfig(filename='brcd.log', level=logging.DEBUG)
    global cfg
    cfg = config.Config(sys.argv[1])
    print cfg.var['gap_symbol']

    aligned_sequences = aligned_sequence_reader(cfg.var['fasta_aligned'])
    # 1st filter: select all primer sites fulfilling melting temperature range
    sites = find_primer_sites(aligned_sequences)
    clusters = combine_sites(sites)
