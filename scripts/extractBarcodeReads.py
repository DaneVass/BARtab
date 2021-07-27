""" 
Dane Vassiliadis. 
Peter MacCallum Cancer Centre.

extractBarcodeReads.py

From a fastq or BAM input, extract reads containing specified barcode pattern.
Return a fasta/fastq file containing the extracted barcode sequences.

From a fastq input, filter out reads that meet the following criteria:
- Exact match to the specified 5' and/or 3' constant regions if given
- [OPTIONAL] A correct sample index
- Average phred quality >= 30 across barcode
- No ambiguous N residues
- Barcode length at least 30 [to do - fix logic for minumum barcode length required]

"""

import re
import sys
import time
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from datetime import datetime
import gzip
import shutil
from optparse import OptionParser

startTime = datetime.now()

# paths for testing
#path = "/Users/vassiliadisdane/Dropbox/PostDoc_PeterMac/Projects/Dawson_Lab/Katie/barcoding_experiment/scripts/test.fastq"

usage = "USAGE: extractBarcodeReads.py --input [input file] --output [output file] --pattern [barcode pattern] [OPTIONS]"

# Define command line options
parser = OptionParser()
parser.add_option('-i', '--input', type = 'string', dest = 'input', help='Path to input fastq file')
parser.add_option('-o', '--output', type = 'string', dest = 'outfile', help='Path to output fastq file')
parser.add_option('-f', '--format', type = 'string', dest = 'format', help='Input file format. Currently supports "fastq" or "bam"')
parser.add_option('-u', '--upconstant', type = 'string', dest = 'upconstant', help="5' constant region")
parser.add_option('-d', '--downconstant', type = 'string', dest = 'downconstant', help="3' constant region")
parser.add_option('-x', '--index', type = 'string', dest = 'index_file', help='Path to file containing sequencing index information for input fastq file')
parser.add_option('-p', '--pattern', type = 'string', dest = 'pattern', help='Regex pattern to search for in reads')
parser.add_option('-p', '--pattern', type = 'string', dest = 'pattern', help='Regex pattern to search for in reads')


(options,args) = parser.parse_args()

# parse input
if options.input is not None:
    try:
        os.path.isfile(options.input)
        path = options.input
        print("Input file set to ", os.path.basename(path))
    except:
        print("No input file defined. Exiting")
        print(usage)
        sys.exit(1)

# parse format
if options.format is not None:
    try:
        os.path.isfile(options.input)
        path = options.input
        print("Input file set to ", os.path.basename(path))
    except:
        print("No input file defined. Exiting")
        print(usage)
        sys.exit(1)
else:
    
    print("no file format given. ")

# parse output
if options.outfile is not None:
    print("Writing output to ", options.outfile)
    outfile_name = options.outfile
else:
    print("No output file defined. Exiting")
    print(usage)
    sys.exit(1)

# parse constant regions
if options.upconstant is not None:
    upstream_constant = options.upconstant
    print("5' constant region: ", upstream_constant.upper())
else:
    upstream_constant = "CGGATCctgaccatgtacgattgacta"
    print("no 5' constant region given, setting constant region to: ", upstream_constant.upper())

# parse index file
if options.index_file is not None:
    try:
        os.path.isfile(options.index_file)
    except:
        print("No valid index file specified. Exiting")
        print(usage)
        sys.exit(1)
    print("Index file: ", options.index_file)
    indexes = options.index_file
    
else:
    print("No index file given. Continuing without")

# parse pattern
if options.pattern is not None:
    print("Regex will search for follwoing pattern: " + pattern)
    regex = pattern
    barcode_re = re.compile(regex)
else:
    regex = "([ATCG][ATCG][GC][AT][GC][ATCG][ATCG][AT][GC][AT]){2,6}"
    print("No regex pattern given. Reverting to SPLINTR pattern" + regex)
    barcode_re = re.compile(regex)
    
# Parsing function


