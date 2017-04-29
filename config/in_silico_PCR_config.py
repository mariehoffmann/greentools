## Constraints for in-silico PCR
COVERAGE_RATE=float(0.8)
MISMATCH=int(3)
# minimal distance between forward and reverse ending
DISTANCE_MIN=int(30)
# maximal distance
DISTANCE_MAX=int(700)
## IO
# File (fasta) containing forward primers with unique IDs
PRIMER_FW=<path_to_forward_primer(s).fsa>
# File (fasta) containing reverse primers with unique IDs
PRIMER_RV=<path_to_reverse_primer(s).fsa>
# Target file with row format "<species_name>;<other>;[<acc_nums>]"
FILE_ACC=['<path_to_accession_file1.csv>', '<path_to_accession_file2.csv>', ...]
# read only row in range [FILE_TAXIDS_ROW1, FILE_TAXIDS_ROW2], maxint indicates all rows
FILE_TAXIDS_ROW1=int(0)
FILE_TAXIDS_ROW2=int(sys.maxint)
# Output files are stored under PATH_OUTPUT/results.csv and PATH_OUTPUT/no_hits.csv
PATH_OUTPUT=<path_to_output_folder>
# if result file exists append or delete
APPEND=bool(False)
PATH_TMP=<path_to_tmp_folder>
# database name given when executing makeblastdb
DB_NAME=['<database_alias1>', '<database_alias2>', ...]
