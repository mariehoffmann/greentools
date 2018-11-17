'''
    Convert a fastq file into a fasta file by removing the quality score line. 
    The output filename will be the same.
    
    Author: Marie Hoffmann marie.hoffmann[at]fu-berlin.de
'''

import os
import sys
import re
import datetime
import subprocess

def run(fastq_file):
    output_file = re.sub(r"(fq|fastq)$", "fa", fastq_file)
    with open(fastq_file, 'r') as f, open(output_file, 'w') as fw:
        for i, line in enumerate(f.readlines()):
            if ((i+1) % 4) == 0 or ((i+2) % 4 == 0):
                continue
            if line.startswith('@'):
                line = '>' + line[1:]
            fw.write(line)

    print("output written to {}".format(output_file))

if __name__ == "__main__":
  if len(sys.argv) == 2:
    run(sys.argv[1])
  else:
    print("Usage: python fq2fa.py <fastq_file>")
