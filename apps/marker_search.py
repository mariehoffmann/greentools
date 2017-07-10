'''
    ============================================================================
        GreenBox - Script Collection for Metagenomics
    ============================================================================
        Task: marker search with low complex primers

        Compute good primer pairs for a set of aligned sequences such that species resolution
        is optimal. In other words, search a region that serves as a marker and for which
        sound primers can be found. For setting path to aligned sequence file, primer pair
        melting temperatures, etc. use marker_search.cfg.
        To produce aligned sequences you can use PAGAN:
            wasabiapp.org/software/pagan/phylogenetic_multiple_alignment

        Call:   python marker_search.py <path_to/marker_search.cfg>

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
log = None
log_file = 'marker_search.log'

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
    def __init__(self, _pos, _len, _score, _melt_rng=None, _gc_content=None):
        self.pos, self.len, self.score = _pos, _len, _score
        self.melt_range = _melt_rng
        self.gc_content = _gc_content

    def write(self):
        print 'start_pos\tlength\tscore\tmelt_range\n' + '{}\t{}\t{}\t{}'.format(self.pos, self.len, self.score, self.melt_range)

class Log(object):
    def __init__(self):
        self.cnt = []  # count for potential sites [(count, filter_name)]
    # count:int, filter:str
    def add(self, count, filter_name):
        self.cnt.append((count, filter_name))

# test for gap freedom in transcript block in
def is_gapfree(aligned_sequences, i, offset):
    if GAP_SYMBOL in set([symbol for aseq in aligned_sequences for symbol in aseq.seq[i:i+offset]]):
        return False
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
        keys_srt = sorted(d.keys())
        block_len = 200  # num characters per row for terminal output
        for block_id in range(len(d[keys_srt[0]].seq)/block_len + 1):
            print '{}..{}'.format(block_id*block_len, (block_id+1)*block_len-1)
            print '\n'.join([key + ': ' + d[key].seq[block_id*block_len:(block_id+1)*block_len] for key in keys_srt])
            print ' '*(len(key)+2) + primer_utils.column_matches(d.values(), block_id*block_len, block_len, GAP_SYMBOL)
    return d.values()

# apply primer site filter constraints and store in site obj
# i: start position, j: offset
def apply_site_filter(aligned_sequences, pos, offset, site):
    if is_gapfree(aligned_sequences, pos, offset) == False:
        logging.DEBUG('Site at {}:{} is not gap free'. format(pos, pos+offset))
        print 'site not gap free'
        return False
    #print [aseq.seq[pos:pos+offset] for aseq in aligned_sequences]
    # check melting temperature in range
    check_melt, min_melt, max_melt = primer_utils.filter_melt(aligned_sequences, pos, offset, cfg)
    # check GC content in range
    check_GC, GC_min, GC_max = primer_utils.filter_GC_content(aligned_sequences, pos, offset, cfg)
    logging.debug('GC content check: {}, range [{}:{}]'.format(check_GC, float(GC_min)/offset, float(GC_max)/offset))
    site.melt_range = (min_melt, max_melt)
    site.gc_content = (float(GC_min*100/offset)/100, float(GC_max*100/offset)/100)
    # check for GC clamps in all sequences at target region
    check_clamp, GC_cnt = reduce(lambda c1, c2: c1 and c2, [primer_utils.filter_GC_clamp(aseq.seq[pos:pos+offset],'+') for aseq in aligned_sequences])
    # check for hairpins
    check_hairpin = reduce(lambda c1, c2: c1 and c2, [primer_utils.filter_hairpin(aseq.seq[pos:pos+offset], cfg) for aseq in aligned_sequences])
    # check for not more max_num_ambig_pos ambiguous columns
    check_ambig_pos = primer_utils.ambiguous_positions(aligned_sequences, pos, offset) <= cfg.var['max_num_ambig_pos']
    if pos in [213, 216, 222, 223, 321, 322]:
        logging.debug('apply_site_filter, pos = {}, check_GC: {}, check_melt: {} ([{}:{}]), check_clamp: {}, check_hairpin: {}, check_ambig_pos: {}'.\
        format(pos, check_GC, check_melt, site.melt_range[0], site.melt_range[1], check_clamp, check_hairpin, check_ambig_pos))
    return reduce(lambda c1, c2: c1 and c2, [check_GC, check_melt, check_clamp, check_ambig_pos]), site  # check_melt, check_clamp, check_hairpin

def apply_site_combined_filter(aligned_sequences, fw, rv):
    ts_pos = fw.pos + fw.len    # transcript start position
    ts_offset = rv.pos - ts_pos    # transcript length
    check_PCR_prod_len = ts_offset in range(int(cfg.var['PCR_prod_len'][0]), int(cfg.var['PCR_prod_len'][1])+1)
    check_PCR_prod_gapfree = is_gapfree(aligned_sequences, ts_pos, ts_offset)
    check_melt_diff = max(fw.melt_range[1], rv.melt_range[1])-min(fw.melt_range[0], rv.melt_range[0]) <= cfg.var['max_melt_diff']
    logging.debug('apply combined filter, fw.pos = {}, rv.pos = {}. check_PCR_prod_len: {}, check_PCR_prod_gapfree: {}, check_melt_diff: {}'.\
    format(fw.pos, rv.pos, check_PCR_prod_len, check_PCR_prod_gapfree, check_melt_diff))
    return reduce(lambda c1, c2: c1 and c2, [check_PCR_prod_len, check_PCR_prod_gapfree, check_melt_diff])

# compute best blocks of fulfilling primer constraints (see marker_search.cfg)
def find_primer_sites(aligned_sequences):
    n = len(aligned_sequences[0].seq) # total sequence length
    pr_min, pr_max = int(cfg.var['opt_primer_len'][0]), int(cfg.var['opt_primer_len'][1])
    sites = []
    i = 0
    while i < n - pr_min:
        # block has to be gap-free
        while i < n - pr_min and is_gapfree(aligned_sequences, i, pr_min) == False:
            #print '{}:{} is not gap free: {}'.format(i, i+pr_min, is_gapfree(aligned_sequences, i, pr_min))
            gaps_right_most = [aseq.seq[i:i+pr_min].rfind('-') for aseq in aligned_sequences]
            #print gaps_right_most
            #print 'i = ', i, ', min block is not gap-free: ', is_gapfree(aligned_sequences, i, i+pr_min), ', skip by ',max(gaps_right_most) + 1
            skip = pr_min+1 if max(gaps_right_most) == -1 else max(gaps_right_most)+1
            i += skip
        if i >= n - pr_min:
            break
        for j in range(pr_min, pr_max):
            if i + j >= n:
                break
            # test for gaps beyond min primer length
            if j >= pr_min:
                if is_gapfree(aligned_sequences, i+pr_min, pr_max-pr_min) == False:
                    i += pr_min-1
                    break
            logging.debug('potential site at [{}, {}]'.format(i, i+j))
            # TODO: score will later include other stats like distance of ambigous positions towards end
            score = primer_utils.ambiguous_positions(aligned_sequences, i, j)
            site = Site(i, j, score, None)
            check, site = apply_site_filter(aligned_sequences, i, j, site)

            if check == True:
                sites.append(site)
        i += 1
        #print 'i = ', i

    log.add(len(sites), 'apply_site_filter')
    # filter for top k scoring conserved sites
    if len(sites) == 0:
        print 'no sites found!'
        sys.exit(0)
    sites.sort(key=lambda s: s.score)
    print 'best and worst scores: ', sites[0].score, sites[-1].score
    len_sites_before = len(sites)
    best_scores = sorted(list(set([site.score for site in sites])))

    if len(best_scores) > cfg.var['top_k_conserved']:
        kth_best_score = best_scores[int(cfg.var['top_k_conserved'])-1]

        idx = len(sites)
        for idx in range(len(sites)):
            if sites[idx].score > kth_best_score:
                break
        del sites[idx:]
    log.add(len(sites), 'top_k_conserved')
    print  'filter for top k conserved regions reduced number of sites from ', len_sites_before, ' to ', len(sites) #logging.info('')
    '''
    for site in sites:
        site.write()
        aligned_sequences[0].write_site(site)
        for seq_obj in aligned_sequences[1:]:
            seq_obj.write_site(site, False)
    '''
    return sites

# 2. compute transcript clusters for all primer combinations
# distance as number of mismatches
def distance(s1, s2):
    return sum([1 if a != b else 0 for (a,b) in zip(s1, s2)])

# distance between 2 clusters is smallest distance between 2 individuals
def single_linkage(c1, c2, pos, offset):
    seq_set1 = set([t_obj.seq[pos:pos+offset] for t_obj in c1])
    seq_set2 = set([t_obj.seq[pos:pos+offset] for t_obj in c2])
    return min([distance(s1, s2) for s1 in seq_set1 for s2 in seq_set2])

# seq_obj_list: [AlignedSequence], pos: transcript start position, length: transcript length
def agglomerative_clustering(seq_obj_list, pos, offset):
    # similarity in terms of differing nt columns
    threshold = int(cfg.var['theta'] * offset)
    logging.debug('threshold for clustering = ' + str(threshold) + ', theta = {}, lenseq = {}'.format(cfg.var['theta'], offset))
    clusters = [[seq_obj] for seq_obj in seq_obj_list]
    min_dist_pair = [(single_linkage(c1, c2, pos, offset), (i,j+i+1)) for i, c1 in enumerate(clusters[:-1]) for j, c2 in enumerate(clusters[i+1:])]
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
        min_dist_pair = [(single_linkage(c1, c2, pos, offset), (i,j+i+1)) for i, c1 in enumerate(clusters[:-1]) for j, c2 in enumerate(clusters[i+1:])]
        min_dist_pair.sort(key=lambda item: item[0])

    #sys.exit(0)
    return clusters

# classify transcripts for all site combinations and select if discrimination is in same manner like the target one
def combine_sites(aligned_sequences, sites):
    cluster_select = []
    k = 0
    output = primer_utils.one_letter_decode + '\n#fw_pos\tfw_len\tfw_seq\t\t\tfw_melt\
    \tGC_content\trv_pos\trv_len\trv_seq\t\t\trv_melt\tGC_content\ttranscript_len\
    hairpin\n'
    sites.sort(key=lambda site: (site.pos, site.len))
    print 'site start positions = ', [site.pos for site in sites]
    for i, fw in enumerate(sites[:-1]):
        i2 = i+1
        while i2 < len(sites):
            rv = sites[i2]
            ts_pos = fw.pos + fw.len    # transcript start position
            ts_offset = rv.pos - ts_pos    # transcript length
            #print 'try (left,right) = ({},{}), diff = {}'.format(ts_pos, ts_pos+ts_offset, ts_offset)
            check = apply_site_combined_filter(aligned_sequences, fw, rv)
            if check == True:
                #print 'product in range'
                logging.debug('combine (left,right) = ({},{}), diff = {}'.format(ts_pos, ts_pos+ts_offset, ts_offset))
                k += 1
                #transcripts = [Transcript(d[key]val[0][pos_start:pos_end], val[1], val[2]) for key in d.keys()]
                #print "current block from " + str(left) + " to " + str(right)
                #for transcript in d.values():
                #    print transcript.seq[left:right]
                clusters = agglomerative_clustering(aligned_sequences, ts_pos, ts_offset)
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

                    #    output = '#fw_pos\tfw_len\tfw_seq\tfw_melt\trv_pos\trv_len\trv_seq\trv_melt\ttranscript_len\n'
                    # debug: should be gap free here
                    for s in aligned_sequences:
                        print s.seq[fw.pos:fw.pos+fw.len]
                    fw_seq = primer_utils.compress(aligned_sequences, fw.pos, fw.len)
                    rv_seq = primer_utils.complement_compress(aligned_sequences, rv.pos, rv.len)
                    has_hairpin = "Yes" if reduce(lambda c1, c2: c1 and c2, [primer_utils.filter_hairpin(aseq.seq[fw.pos+fw.len:rv.pos], cfg) for aseq in aligned_sequences]) == True else "No"
                    output += '{}\t{}\t{}\t[{},{}]\t\t[{},{}]\t{}\t{}\t{}\t[{},{}]\t[{},{}]\t{}\t{}\n'.format(\
                    fw.pos, fw.len, fw_seq, fw.melt_range[0], fw.melt_range[1], fw.gc_content[0], fw.gc_content[1], \
                    rv.pos, rv.len, rv_seq, rv.melt_range[0], rv.melt_range[1], rv.gc_content[0], rv.gc_content[1], \
                    ts_offset, has_hairpin)
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
    print output
    return cluster_select

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Wrong number of arguments: python marker_search.py <path_to/marker_search.cfg>"
        sys.exit(0)
    os.remove(log_file)
    logging.basicConfig(filename=log_file, file_mode='w', level=logging.DEBUG)
    global cfg
    cfg = config.Config(sys.argv[1])
    print 'cfg = ', cfg
    aligned_sequences = aligned_sequence_reader(cfg.var['fasta_aligned'])
    # 1st filter: select all primer sites fulfilling melting temperature range
    global log
    log = Log()
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
