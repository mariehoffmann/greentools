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
import primer_utils as prutils
import config

labels = set(['Kellicottia', 'Keratella', 'Synchaeta', 'Trichocerca'])
primer_len = 20
top_k = 5
transcript_range = [20, 500]
theta = 0.18 # procentual difference between 2 items in same cluster
#fasta_aligned = 'Rotifer_all_NCBI_Pagan_realignment.fa'
fasta_aligned = 'Rotifer_few_NCBI_Pagan_realignment.fa'
NUM_CLUSTER_MIN = 3
NUM_CLUSTER_MAX = 5
MAX_NUM_CLUSTER_DIVERS = 0
GAP_SYMBOL = '-'
max_melt_diff=15
max_melt_temp=77.0
min_melt_temp=55.0

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

    def write_block(self, i, j):
        print self.acc + ', ' +  self.label + ', ' + self.seq[i:j]

# store potential primer site by index and 'conservation' score
class Site(object):
    self __init__(self, _pos, _score):
        self.pos, self.score = _pos, _score

# test for gap freedom in transcript block
def is_gapfree(transcripts, i, j):
    if GAP_SYMBOL in set([symbol for transcript in transcripts for symbol in transcript.seq[i:j]]) is True:
        return False
    return True

# read aligned sequences and return dict {accession: AlignedSequence}
def aligned_sequence_reader(fasta_aligned):
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
                for l in labels:
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
def get_primer_sites(aligned_sequences):
    seq_len = len(d.values()[0].seq)
    sites = []
    i = 0
    while i < seq_len-primer_len:
        # block has to be gap-free
        while i <  seq_len-primer_len and is_gapfree(d.values(), i, i+primer_len) is False:
            i += 1
        if i >= seq_len - primer_len:
            break
        for align_seq in d.keys():
            if align_seq[i:i+primer_len+1].find('-') > -1:
                i += align_seq[i:i+primer_len+1].find('-')
        score = 0
        for j in range(primer_len):
            # abort if one sequence has gap
            for val in d.values():
                #print val.seq
                seq = val.seq
                if seq[i] == '-':
                    score = -1
            if score < 0:
                break;
            score += 6 - len(set([val.seq[i+j] for val in d.values()]))
        i += 1
        if score < 0:
            continue
        sites.append((i-1, score))

    # filter top k
    min_k = sorted([site[1] for score in sites], reverse=True)[top_k-1]
    print "min k: ", min_k
    sites = [Site(site[0], site[1]) for site in sites if site[1] >= min_k]
    for idx, site in enumerate(sites):
        print "\nblock #" + str(idx) + ':'
        for val in d.values():
            print val.seq[site.pos:site.pos+primer_len+1]
    #sys.exit(0)
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
    #for seqobj in seq_obj_list:
    #    print seqobj.seq[left:right]
    threshold = int(theta * len(seq_obj_list[0].seq[left:right]))
    clusters = [[seq_obj] for seq_obj in seq_obj_list]
    min_dist_pair = [(single_linkage(c1, c2, left, right), (i,j+i+1)) for i, c1 in enumerate(clusters[:-1]) for j, c2 in enumerate(clusters[i+1:])]
    min_dist_pair.sort(key=lambda item: item[0])
    #print 'min_dist_pair = ', min_dist_pair[0]
    while min_dist_pair[0][0] <= threshold and len(clusters) > 1:
        i, j = min_dist_pair[0][1][0], min_dist_pair[0][1][1]
        # merge cluster i and j
        #print "clusters before merge: ", len(clusters), ', ', clusters
        #print "i + j = ", clusters[i] + clusters[j]
        #clusters = clusters[:i] + [clusters[i] + clusters[j]] + clusters[i+1:j] + clusters[j+1:]
        #print 'i = ', i, ', j = ', j
        cluster_merged = clusters[i] + clusters[j]
        del clusters[j]
        del clusters[i]
        clusters.append(cluster_merged)
        #print "clusters after merge: ", len(clusters), ', ', clusters
        if len(clusters) < 2:
            return clusters
        # compute new smallest distance
        min_dist_pair = [(single_linkage(c1, c2, left, right), (i,j+i+1)) for i, c1 in enumerate(clusters[:-1]) for j, c2 in enumerate(clusters[i+1:])]
        min_dist_pair.sort(key=lambda item: item[0])
        #print "new min dist pair = ", min_dist_pair, ' with i,j = ', min_dist_pair[0][1][0],  min_dist_pair[0][1][1]
    return clusters

# classify transcripts for all site combinations and select if discrimination is in same manner like the target one
def combine_sites(sites):
    cluster_select = []

    for i, fw in enumerate(sites[:-1]):
        for rv in sites[i+1:]:
            #print 'fw = ', fw
            #print 'rv = ', rv
            if rv[0] - 20 - fw[0] in range(transcript_range[0], transcript_range[1]+1) and is_gapfree(d.values(), fw[0]+primer_len, rv[0]) is True:
                left, right = fw[0] + 20, rv[0]
                #print "pos_start, end = ", pos_start, pos_end
                #transcripts = [Transcript(d[key]val[0][pos_start:pos_end], val[1], val[2]) for key in d.keys()]
                #print "current block from " + str(left) + " to " + str(right)
                #for transcript in d.values():
                #    print transcript.seq[left:right]
                clusters = agglomerative_clustering(d.values(), left, right)
                #print '\n##################################'
                label_identities = [0 if len(set([transcript.label for transcript in cluster])) == 1 else 1 for cluster in clusters]
                #print 'label identity in clusters: ', label_identities
                #print "sum(label_identities) < 2: ", sum(label_identities) < 2
                #print "len(clusters) in range(NUM_CLUSTER_MIN, NUM_CLUSTER_MAX+1) is True", len(clusters) in range(NUM_CLUSTER_MIN, NUM_CLUSTER_MAX+1)
                #print 'len(clusters) = ', len(clusters)
                if sum(label_identities) <= MAX_NUM_CLUSTER_DIVERS and len(clusters) in range(NUM_CLUSTER_MIN, NUM_CLUSTER_MAX+1):  #  0.1*len(label_identities)
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
    cfg = config.Config(argv[1])
    aligned_sequences = aligned_sequence_reader(cfg)
    sites = get_primer_sites(aligned_sequences)
    clusters = combine_sites(sites)
