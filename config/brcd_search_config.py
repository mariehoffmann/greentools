## configuration file for primer design w.r.t. label distinction, i.e.
## an unsupervised clustering should result in groups containing the same label
# aligned sequence file in fasta format, nucleotides in upper case
fasta_aligned = '<path_to_aligned_sequence_file.fa>'
# gap symbol
gap_symbol = -
# similarity overlap
similarity_overlap = 10
# tail length
tail_len = 8
# maximal primer length
max_primer_len = 25
# minimal primer length
min_primer_len = 18
# maximal melting temperature
max_melt_temp = 77.0
# minimal melting temperature
min_melt_temp = 55.0
# critical melting temperature
critical_melt_temp = 58.0
# minimal melting temperature with ambiguity positions
min_melt_temp_with_ambig_pos = 58.0
# minimal number of 3'-end matches
min_num_3end_matches = 2
# maximum number of ambiguity positions
max_num_ambig_pos = 4
# minimum length with three or more ambiguity positions
min_len_mult_ambig_pos = 25
# minimum length with two ambiguity positions
min_len_2_ambig_pos = 21
# critical ambiguity position distance from 3'-end
crit_ambig_dist_from_3 = 5
# maximum diversity per nucleotide position
max_divers_per_pos = 4
# maximum number of ambiguity positions with maximum diversity
max_num_ambig_pos_max_divers = 1
# minimum number of nucleotides per nucleotide position
min_num_nucl_per_pos = 2
# primer concentration in nanomolar
primer_concent_nanomolar = 250
# salt concentration in molar
salt_concent_molar = 0.05
# maximum melting temperature difference
max_melt_diff = 15
# minimum PCR product length
min_amplicon_len = 20
# maximum PCR product length
max_amplicon_len = 500
# optimal primer length dispensation with no ambiguity positions
# conserved region length
conserved_region_len = 2
# conservation window length
conservation_window_len = 80
# minimum percentage of conserved nucleotides within conservation window
min_conserved_in_window = 0.9
# introns in sequences
# At least one primer must be based on three sequences
#### Settings for classification (agglomerative clustering) ####
# labels are uniquely matching keywords or regular expressions in fasta header lines
labels = ['<label1>', '<label2>', ...]
# similarity threshold, i.e. relative number of varying columns
theta= 0.18
# Minimal number of clusters resulting from transcript based classification
min_num_cluster = 3
# Maximal number of clusters resulting from transcript based classification
max_num_cluster = 5
# Maximal number of final clusters with mixed items w.r.t. their labels
max_num_clusters_mixed = 2
