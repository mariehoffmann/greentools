'''
    ============================================================================
                GreenBox - Script Collection for Metagenomics
    ============================================================================
    in-silico_PCR simulates the success of a PCR given sets of forward and reverse
    primer sequences by testing for matches against a blast database built on the
    target species genomes.
    Blast queries will have the form
        blastn -db non_human -query $Primer{Fw|Rev} -seqidlist $AccNums -task blastn -outfmt 6
    and all pairwise results for forward and reverse primer sequences are tested for
    suitability. Criteria for a successful in-silico PCR, e.g. less than 3 mismatches,
    matched pair distance in range [150, 400], etc., are set in the configuration
    file ../config/in_silico_PCR_config.py.
    Target sequences for each species have to be provided as a list of accession numbers.
    Use the script name2acc.py to receive all accession numbers as a csv file for
    a species given its name and the fasta file from which the blast db was built.

    This script requires:
        1. BLAST: local install of blast command line tools and execution of makeblastdb
        with option -parse_seqids in order to receive sequence ids. Example: use
        the extracted fasta file nt available ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/
        as a source to create a database named 'non_human':
            makeblastdb -in nt -out non_human -parse_seqids -dbtype nucl
        2. Linux/Unix tools: echo, sed, wc
    ============================================================================
    Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
    ============================================================================

'''

import os
import sys
import subprocess
import re
import logging
from datetime import datetime

LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}

# outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = 0,1,2,3,4,5,6,7,8,9,10,11
eps = 1e-6

class Primer(object):
    def __init__(self, primer_id_seq):
        self.id = primer_id_seq[0]
        self.seq = primer_id_seq[1]

    def __str__(self):
        return "({:s}, {:s})".format(self.id, self.seq)

class Config(object):
    def __init__(self, config_file):
        self.settings = {'COVERAGE_RATE': None, 'MISMATCH': 0, 'DISTANCE_MIN': None,
        'DISTANCE_MAX': None, 'PATH_PRIMERS': None, 'PRIMER_FW': None, 'PRIMER_RV': None,
        'FILE_ACCS': None, 'PATH_OUTPUT': None, 'DB_NAME': None, 'DIR_TMP': None,
        'FROM_LINE': 0, 'TO_LINE': None, 'APPEND': False}
        eval_rx = re.compile("(int|float|bool)\(\d*\.?(\d+|True|False)\)")
        with open(config_file, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                line = line.rstrip().split('=')
                if line[0] in self.settings.keys():
                    if eval_rx.match(line[1]) is not None:
                        self.settings.update({line[0]: eval(line[1])})
                    else:
                        self.settings.update({line[0]: line[1]})
                else:
                    logging.warning("Unknown config var: " + line[0])

        self.primer_fw_list, self.primer_rv_list = None, None
        self.taxid2accs = {}  # {(name, taxid): [acc_nums]}
        self.result_file = os.path.join(self.settings["PATH_OUTPUT"], "results.csv")
        self.no_hits_file = os.path.join(self.settings["PATH_OUTPUT"], "no_hits.csv")
        self.output_details = [qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore]
        # delete result files from previous runs
        print self.settings['APPEND']
        if self.settings['APPEND'] is False:
            logging.info("delete old files")
            if os.path.isfile(self.result_file):
                os.system("rm " + self.result_file)
            if os.path.isfile(self.no_hits_file):
                os.system("rm " + self.no_hits_file)
            os.system('echo "#species;taxid;qseqid_fw;qseqid_rv;sseqid;pident_fw;pident_rv;evalue_fw;evalue_rv;sstart_fw;send_fw;sstart_rv;send_rv;dist" > {:s}'.format(self.result_file))
        # create tmp directory if not existing
        if not os.path.isdir(self.settings["DIR_TMP"]):
            os.makedirs(self.settings["DIR_TMP"])

# convert non-string attributes from blastn result string
# [qseqid:str sseqid:str pident:float length:int mismatch:int gapopen:int qstart:int qend:int sstart:int send:int evalue:float bitscore:float]
def blast_fmt6(query):
    return [query[qseqid], query[sseqid], float(query[pident]), int(query[length]),
    int(query[mismatch]), int(query[gapopen]), int(query[qstart]), int(query[qend]),
    int(query[sstart]), int(query[send]), float(query[evalue]), float(query[bitscore])]

## connect two query results from blastn with output format 6 as a list

class Match(object):
    def __init__(self, query_fw, query_rv):
        self.fw, self.rv = blast_fmt6(query_fw), blast_fmt6(query_rv)

    def __str__(self):
        return "forward match: " + '\t'.join([str(attr) for attr in self.fw]) + \
        "\nreverse match: " + '\t'.join([str(attr) for attr in self.rv])


taxids_dict = {}    # {taxid: species_name}, for later sorting
query_results = {}  # {taxid: QueryResult}

# primers  [[fws], [rvs]]
def query_all(config):
    wc_res = subprocess.check_output(["wc", "-l", config.settings['FILE_ACCS']]).strip().split(' ')[0]
    for row in range(config.settings['FROM_LINE'], min(config.settings['TO_LINE'], int(wc_res))):
        logging.info("STATUS: row " + str(row))
        success, name_taxid_accs = query(config, row)
        if success is False:
            with open(config.no_hits_file, 'a') as fhdlr:
                fhdlr.write(';'.join([name_taxid_accs[0], name_taxid_accs[1]]) + '\n')

# BLAST query for all forward (reverse, resp.) primers restricted to given accession numbers
def query(config, row):
    # get accession numbers from taxid file and write them row-wise in tmp file
    name_taxid_accs = subprocess.check_output(["sed", "-n", "{:d}p".format(row+1), config.settings['FILE_ACCS']]).rstrip().split(';')

    logging.debug(name_taxid_accs[:50])
    acc_file = os.path.join(config.settings['DIR_TMP'], name_taxid_accs[0])
    accs = '\n'.join(name_taxid_accs[2:min(10, len(name_taxid_accs))])
    logging.debug("write tmp acc file: " + acc_file)
    logging.debug("accs = " + accs)
    with open(acc_file, 'w') as f:
        f.write(accs)

    #blastn -db non_human -query $Primer -seqidlist $Acanthoceras -task blastn -outfmt "6
    blast_fw = subprocess.check_output(["blastn", "-db", config.settings["DB_NAME"], "-query",
    config.settings["PRIMER_FW"], "-seqidlist", acc_file, "-task", "blastn", "-outfmt", "6"])
    blast_fw = [line.split('\t') for line in blast_fw.strip().split('\n')]

    logging.debug("blast fw: " + str(blast_fw))

    blast_rv = subprocess.check_output(["blastn", "-db", config.settings["DB_NAME"], "-query",
    config.settings["PRIMER_RV"], "-seqidlist", acc_file, "-task", "blastn", "-outfmt", "6"])
    blast_rv = [line.split('\t') for line in blast_rv.strip().split('\n')]

    logging.debug("blast rv: " + str(blast_rv))

    # no success if one of the queries returned empty string
    if len(blast_fw[0]) == 1 or len(blast_rv[0]) == 1:
        return False, name_taxid_accs
    # check all pairwise forward-reverse matches if they meet constraints:
    # identical query sequence, mismatches not exceeding, distance in range, etc.
    matches = []
    for fw in blast_fw:
        for rv in blast_rv:
            if fw[sseqid] == rv[sseqid] \
            and int(fw[mismatch]) <= config.settings["MISMATCH"] \
            and int(rv[mismatch]) <= config.settings["MISMATCH"]:
                if int(fw[sstart]) < int(rv[sstart]) and int(rv[sstart]) - int(fw[send]) in range(config.settings["DISTANCE_MIN"], config.settings["DISTANCE_MAX"]):
                    matches.append(Match(fw, rv))
                    logging.debug("append match")
                elif int(rv[sstart]) < int(fw[sstart]) and int(fw[sstart]) - int(rv[send]) in range(config.settings["DISTANCE_MIN"], config.settings["DISTANCE_MAX"]):
                    logging.debug("Unhandled case: rv[sstart] < fw[sstart] for query " + fw[qseqid])
    # delete tmp accession number file
    # todo
    logging.debug("Delete tmp accession file: " + acc_file)
    os.remove(acc_file)

    logging.debug("num matches = " + str(len(matches)))
    if len(matches) == 0:
        return False, name_taxid_accs
    # select best matches w.r.t. evalue of fw and rv sequences
    matches.sort(key=lambda m: m.fw[evalue] + m.rv[evalue])
    for match in matches:
        logging.debug(str(match))
        logging.debug(match.fw[evalue] + match.rv[evalue])
    min_evalue = matches[0].fw[evalue] + matches[0].rv[evalue]
    matches = filter(lambda m:  m.fw[evalue] + m.rv[evalue] <= min_evalue + eps, matches)

    # write intermediate results to csv file
    # species;taxid;qseqid_fw;qseqid_rv;sseqid;pident_fw;pident_rv;evalue_fw;evalue_rv;sstart_fw;sstart_rv;dist
    logging.debug(config.result_file)
    with open(config.result_file, 'a') as fhdlr:
        logging.info("Write match results for " + name_taxid_accs[0])
        for match in matches:
            fhdlr.write(';'.join([name_taxid_accs[0], name_taxid_accs[1], match.fw[qseqid],
            match.rv[qseqid], match.fw[sseqid], str(match.fw[pident]), str(match.rv[pident]),
            str(match.fw[evalue]), str(match.rv[evalue]), str(match.fw[sstart]), str(match.fw[send]),
            str(match.rv[sstart]), str(match.rv[send]), str(match.rv[sstart] - match.fw[send])]) + '\n')
        return True, None
    return False, None

if __name__ == "__main__":
    if len(sys.argv) <  2:
        print("Usage: python in-silico_PCR.py <path_to_config.py> [<logging_level>]")
        sys.exit(-1)
    level_name = 'info'
    if len(sys.argv) == 3:
        level_name = sys.argv[1]
    level = LEVELS.get(level_name, logging.NOTSET)
    logging.basicConfig(level=level)

    config = Config(sys.argv[1])
    query_all(config)
