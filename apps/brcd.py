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

        Call:   python brcd.py <path_to/brcd.cfg>

        Output: csv file with forward (fw) and reverse (rv) primer start positions and their
                sequences and a score based on how well the design constraints match. Format:
                [fw_pos, fw_len, fw_primer_seq, fw_pos, fw_len, fw_primer_seq, score]

        Algorithm:
                1. Search of conservative blocks, i.e. low number of different basepairs per
                column.
                2. For each block combination with sound transcript lengths an agglomerative
                clustering is performed.
                3. Reject or accept solution: since the labels are known, their correspondance
                to the resulting clusters is checked. In a best case, only identically labeled
                individuals end up in the same cluster.

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

cfg = None
log_file = 'brcd.log'

GAP_SYMBOL = '-'

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

    def write_site(self, site, flag_w_header=True):
        if flag_w_header == True:
            print 'accession\tlabel\t\tsubsequence'
        print self.acc + '\t' +  self.label + '\t' + self.seq[site.pos:site.pos+site.len]

    def write_transcript(self, site1, site2, flag_w_header=True):
        if flag_w_header == True:
            print 'accession\tlabel\t\ttranscript_len\tsubsequence'
        left, right = site1.pos+site1.len, site2.pos
        seq_display = self.seq[left:right]
        if len(seq_display) < 80:
            print '{}\t{}\t{}\t\t[{}:{}]{}'.format(self.acc, self.label, right-left, left, right, self.seq[left:right])
        else:
            offset = 25
            print '{}\t{}\t{}\t\t[{}:{}...{}:{}] {} ... {}'.format(self.acc, self.label, right-left, left, left+offset, right-offset, right, self.seq[left:left+offset], self.seq[right-offset:right])

# store potential primer site by start index, length, 'conservation' score, melting temperature range
class Site(object):
    def __init__(self, _pos, _len, _score, _melt_rng):
        self.pos, self.len, self.score, self.melt_range = _pos, _len, _score, _melt_rng
    def write(self):
        print 'start_pos\tlength\tscore\tmelt_range\n' + '{}\t{}\t{}\t{}'.format(self.pos, self.len, self.score, self.melt_range)

# test for gap freedom in transcript block
def is_gapfree(transcripts, i, j):
    #print set([symbol for transcript in transcripts for symbol in transcript.seq[i:j]])
    if GAP_SYMBOL in set([symbol for transcript in transcripts for symbol in transcript.seq[i:j]]):
        #print 'is gap-free: False'
        return False
    #print 'is gap-free: True'
    return True

# read aligned sequences and return dict {accession: AlignedSequence}
def aligned_sequence_reader(fasta_aligned):
    print 'cfg in as_reader: ', cfg
    print 'cfg.labels = ', cfg.var['labels']
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
    pr_min, pr_max = int(cfg.var['min_primer_len']), int(cfg.var['max_primer_len'])
    print 'primer_len_interval = [' + str(pr_min) + ', ' + str(pr_max) + ']'
    sites = []
    i = 0
    while i < n - pr_min:
        print 'i = ', i, ', min block is gap-free: ', is_gapfree(aligned_sequences, i, i+pr_min)
        # block has to be gap-free
        # TODO: skip more than one nt
        while i < n - pr_min and is_gapfree(aligned_sequences, i, i+pr_min) == False:
            gaps_right_most = [aseq.seq[i:i+pr_min].rfind('-') for aseq in aligned_sequences]
            #print 'i = ', i, ', min block is not gap-free: ', is_gapfree(aligned_sequences, i, i+pr_min), ', skip by ',max(gaps_right_most) + 1
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
            logging.debug('potential site at [{}, {}]'.format(i, i+j))
            score = primer_utils.variation_score(aligned_sequences, i, i+j)

            check_melt, min_melt, max_melt = primer_utils.filter_melt(aligned_sequences, i, i+j, cfg)
            print 'score = ', score, ', check_melt = ', check_melt, ', min_melt = ', min_melt, ', max_melt = ', max_melt
            check_GC, GC_min, GC_max = primer_utils.filter_GC_content(aligned_sequences, i-1, j)
            logging.debug('GC content check: {}, range [{}:{}]'.format(check_GC, float(GC_min)/j, float(GC_max)/j))
            if check_melt == True: # and check_GC == True:
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
        aligned_sequences[0].write_site(site)
        for seq_obj in aligned_sequences[1:]:
            seq_obj.write_site(site, False)
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
    # similarity in terms of differing nt columns
    threshold = int(cfg.var['theta'] * (right-left))
    logging.debug('threshold for clustering = ' + str(threshold) + ', theta = {}, lenseq = {}'.format(cfg.var['theta'], right-left))
    clusters = [[seq_obj] for seq_obj in seq_obj_list]
    min_dist_pair = [(single_linkage(c1, c2, left, right), (i,j+i+1)) for i, c1 in enumerate(clusters[:-1]) for j, c2 in enumerate(clusters[i+1:])]
    min_dist_pair.sort(key=lambda item: item[0])
    #print 'min_dist_pair = ', min_dist_pair[0]
    logging.debug('threshold = {}'.format(threshold))
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

    #sys.exit(0)
    return clusters

# classify transcripts for all site combinations and select if discrimination is in same manner like the target one
def combine_sites(aligned_sequences, sites):
    cluster_select = []
    k = 0
    output = '#fw_pos\tfw_len\tfw_seq\tfw_melt\trv_pos\trv_len\trv_seq\trv_melt\ttranscript_len\n'
    sites.sort(key=lambda site: (site.pos, site.len))
    for i, fw in enumerate(sites[:-1]):
        i2 = i+1
        while i2 < len(sites):
            rv = sites[i2]
            left, right = fw.pos + fw.len, rv.pos
            if right - left in range(int(cfg.var['min_PCR_prod_len']), int(cfg.var['max_PCR_prod_len'])+1) and is_gapfree(aligned_sequences, left, right) is True:
                logging.debug('combine (left,right) = ({},{}), diff = {}'.format(left, right, right-left))
                k += 1
                #transcripts = [Transcript(d[key]val[0][pos_start:pos_end], val[1], val[2]) for key in d.keys()]
                #print "current block from " + str(left) + " to " + str(right)
                #for transcript in d.values():
                #    print transcript.seq[left:right]
                clusters = agglomerative_clustering(aligned_sequences, left, right)
                #print '\n##################################'
                label_identities = [0 if len(set([transcript.label for transcript in cluster])) == 1 else 1 for cluster in clusters]
                logging.debug('label identity in clusters: ' + str(label_identities))
                #print "sum(label_identities) < 2: ", sum(label_identities) < 2
                #print "len(clusters) in range(NUM_CLUSTER_MIN, NUM_CLUSTER_MAX+1) is True", len(clusters) in range(NUM_CLUSTER_MIN, NUM_CLUSTER_MAX+1)
                logging.debug('len(clusters) = {}'.format(len(clusters)))
                logging.debug("sum(label_identities) <= cfg.var['max_num_clusters_nunique'] is {}".format(sum(label_identities) <= cfg.var['max_num_clusters_nunique']))
                logging.debug("len(clusters) in range(int(cfg.var['min_num_cluster']), int(cfg.var['max_num_cluster'])+1) == True is {}".format(len(clusters) in range(int(cfg.var['min_num_cluster']), int(cfg.var['max_num_cluster'])+1)))
                logging.debug("cluster_size range is [{}:{}]".format(int(cfg.var['min_num_cluster']), int(cfg.var['max_num_cluster']+1)))
                if sum(label_identities) <= cfg.var['max_num_clusters_nunique'] and len(clusters) in range(int(cfg.var['min_num_cluster']), int(cfg.var['max_num_cluster'])+1):  #  0.1*len(label_identities)
                    print 'discr comb found!'
                    cluster_select.append((fw, rv))
                    print '####### NEW CLASSIFICATION ########'
                    print 'number of classes: ', len(clusters)

                    for cluster in clusters:
                        cluster[0].write_transcript(fw, rv)
                        for transcript in cluster[1:]:
                            transcript.write_transcript(fw, rv, False)
                        print '-------------------------------'
                        check_GC, min_GC, max_GC = primer_utils.filter_GC_content(aligned_sequences, fw.pos, fw.len)
                        logging.debug("check_GC, min_GC, max_GC = {}, [{}:{}]".format(check_GC, float(min_GC)/fw.len, float(max_GC)/fw.len))
                    #    output = '#fw_pos\tfw_len\tfw_seq\tfw_melt\trv_pos\trv_len\trv_seq\trv_melt\ttranscript_len\n'
                    # debug: should be gap free here
                    fw_seq = primer_utils.compress(aligned_sequences, fw.pos, fw.len)
                    rv_seq = primer_utils.complement_compress(aligned_sequences, rv.pos, rv.len)

                    output += '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(\
                    fw.pos, fw.len, fw_seq, primer_utils.primer_melt_wallace(fw.seq), \
                    rv.pos, rv.len, rv_seq, primer_utils.primer_melt_wallace(rv.seq), \
                    right-left)
                    #sys.exit(0)
            # skip 2nd index if next site is only extension of previous one and did not produce correct
            # TODO: check this logic for correctness
            if len(cluster_select) == 0 or cluster_select[-1][1].pos != rv.pos:
                while i2 < len(sites)-1 and sites[i2].pos == sites[i2+1].pos:
                    i2 += 1
                #i2 -= 1
                #if k >= 500:
                #    sys.exit(0)
            i2 += 1
    return cluster_select

# return all sites with scores smaller or equal to k-th score in sorted list
def filter_topk(sites, k):
    topk_score = sorted(list(set([site.score for site in sites])))[k]
    print 'topk_score = ', topk_score

    return [site for site in sites if site.score <= topk_score]

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Wrong number of arguments: python brcd.py <path_to/brcd.cfg>"
        sys.exit(0)
    os.remove(log_file)
    logging.basicConfig(filename=log_file, file_mode='w', level=logging.DEBUG)
    #logger = logging.getLogger()
    global cfg
    cfg = config.Config(sys.argv[1])
    print 'cfg = ', cfg
    aligned_sequences = aligned_sequence_reader(cfg.var['fasta_aligned'])
    # 1st filter: select all primer sites fulfilling melting temperature range
    sites = find_primer_sites(aligned_sequences)
    logging.info('Number of single primer sites found: ' + str(len(sites)))
    print 'Number of single primer sites found: {}'.format(len(sites))
    k = 1
    if len(sites) >= 1024:
        sites = filter_topk(sites, k)
    logging.debug('Best score: ' + str(sites[0].score))
    logging.info('Number of sites filtered for {} best scores: {}'.format(k, len(sites)))
    #sys.exit(0)
    clusters = combine_sites(aligned_sequences, sites)
