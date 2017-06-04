'''
    Author: Marie Hoffmann ozymandiaz147@googlemail.com
    Script for UNIX/LINUX to perform PEAR analysis for 18s datasets.
    Call $python PEAR_script.py <indir> <outdir>
    Give absolute paths of input and output directories. A logfile will store the PEAR command line output.
'''

import os
import sys
import re
import datetime

# regex to catch forward files, backward read files are named same way ending with "2"
rx = re.compile("demultipl__barcode(.+)?_1.fastq")

def callPEAR(indir, outdir):
  logfname = os.path.join(outdir, "PEAR.log")
  print("Logfile written to ", logfname)
  with open(logfname, 'w') as logfile:
    logfile.write("Analysis time {}\n\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
  
  for root, dirs, filenames in os.walk(indir):
    for f in filenames:
      mo = rx.match(f)
      if mo is not None:
        fname_fw = mo.group(0)
        fname_bw = os.path.join(indir, fname_fw.replace("1.fastq", "2.fastq"))
        fname_out = os.path.join(outdir, fname_fw.replace("_1.fastq", ""))
        fname_fw = os.path.join(indir, fname_fw)
        command = "PEAR -f {} -r {} -o {} >> {}".format(fname_fw, fname_bw, fname_out, logfname)
        os.system(command)

if __name__ == "__main__":
  if len(sys.argv) == 3:
    callPEAR(sys.argv[1], sys.argv[2])
  else:
    print("Usage: python PEAR_script.py <indir> <outdir>")

