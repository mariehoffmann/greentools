'''
    ============================================================================
                GreenBox - Script Collection for Metagenomics
    ============================================================================
    Retrieve accession numbers by name from a fasta file with header line format
            >accession_number description
    Names are given (lexicographically sorted) line by line in a csv file
    (delimiter ';'). First column is the name to be searched
    for, additional columns, if present, will be copied for the output.
    Format: "name; col2; col3; ...; acc1; acc2; ..." if accession number could
    be retrieved, else names are printed std out.
    For every name pattern found in a title line, the accession numbers are
    collected. This script is optimized under the assumption that the name
    pattern file is much smaller than the fasta source file. Additional speedup
    is yielded by extracting first only the header lines from the fasta file and
    using them as input source.
    This script requires: grep

    ============================================================================
    Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
    ============================================================================

'''

import sys
import re
import math
import os
import subprocess
import time

# output csv species, taxid, acc1, acc2, ...
def parse_fasta(fasta_file, name_file, output_file, rng=(0,None)):
    # dictionary of names to be searched and a regex
    name_dict = {}  # {name: rest_line}
    with open(name_file, 'r') as nf:
        for line in nf.readlines():
            line = line.strip().split(';')
            name_dict[line[0]] = ';'.join(line[1:])

    names = sorted(name_dict.keys())
    if rng[1] is None:
        rng = (rng[0], len(names))
    acc_dict = {key: [] for key in name_dict.keys()}
    output_file_no_acc = output_file[:-4] + "_no_acc.csv"
    #subprocess.call(["rm", output_file])
    for i, name in enumerate(names):
        if i >= rng[1]:
            break
        if i < rng[0]:
            continue
        print("parsing fasta file for " + name + " ...")
        # grep return codes: 0 success, 1 no match, 2 error
        try:
            #cmd = 'grep -i "(?=%s)(?=[rRNA|16S|18S])" %s' % (name, fasta_file) #cmd = 'grep -i "%s" %s' %(name, fasta_file)
            cmd = "awk '/%s/ && /[18S|16S|rRNA]/' %s" % (name, fasta_file)
            
            print(cmd)
            lines_name = subprocess.check_output(cmd, shell=True)
            print("grep result: " + str(lines_name))
            acc_nums = []
            for line in lines_name.strip().split('\n'):
                acc_nums.append(line[1:].split(' ')[0])
            acc_dict.update({name: acc_nums})
            with open(output_file, 'a') as of:
                of.write(';'.join([name, name_dict[name]] + acc_dict[name]) + "\n")
        except subprocess.CalledProcessError as e:
            print("No match for " + name)
            if e.returncode > 1:
                raise
            with open(output_file_no_acc, 'a') as f:
                f.write(name + '\n')

if __name__ == "__main__":
    if len(sys.argv) not in range(4, 7):
        print "Usage: python name2acc.py <src_fasta|headers_only> <name_file> <output_file.csv> [<from_line> [<to_line>]]"
        sys.exit(-1)
    start_t = time.time()
    if len(sys.argv) == 4:
        parse_fasta(fasta_file=sys.argv[1], name_file=sys.argv[2], output_file=sys.argv[3])
    elif len(sys.argv) == 5:
        parse_fasta(fasta_file=sys.argv[1], name_file=sys.argv[2], output_file=sys.argv[3], rng=[int(sys.argv[4]), None])
    else:
        parse_fasta(fasta_file=sys.argv[1], name_file=sys.argv[2], output_file=sys.argv[3], rng=[int(sys.argv[4]), int(sys.argv[5])])
    print("exec time %s seconds" % (time.time() - start_t))
