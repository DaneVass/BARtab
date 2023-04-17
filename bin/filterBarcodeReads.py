#!/usr/bin/env python3

""" 
Dane Vassiliadis. 
Peter MacCallum Cancer Centre.

filterBarcodeReads.py

From a FASTQ file input, filter reads meeting desired quality thresholds

From a fastq input, filter out reads that meet the following criteria:
- [OPTIONAL] A correct sample index
- Average phred quality >= 30 across barcode
- No ambiguous N residues

"""

import re
import sys
import time
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
#from Bio.Alphabet import IUPAC
from datetime import datetime
import gzip
import shutil
from optparse import OptionParser

startTime = datetime.now()

usage = "USAGE: filterBarcodeReads.py --input [input fastq] --output [output fastq] [options] "

# Define command line options
parser = OptionParser()
parser.add_option('-i', '--input', type = 'string', dest = 'input', help='Path to input fastq file')
parser.add_option('-o', '--output', type = 'string', dest = 'outfile', help='Path to output fastq file')
parser.add_option('-x', '--index', type = 'string', dest = 'index_file', help='[OPTIONAL] Path to file containing sequencing index information for input fastq file')
parser.add_option('-p', '--minqual', type = 'string', dest = 'pattern', help='Minimum average phred quality across read [30]')

(options,args) = parser.parse_args()

# make sure input string is a file
if options.input is not None:
    try:
        os.path.isfile(options.input)
        path = options.input
        print("Input file set to ", os.path.basename(path))
    except:
        print("No input file defined. Exiting")
        print(usage)
        sys.exit(1)

# make sure we are working with fastq files
try:
    os.path.isfile(options.input)
    path = options.input
    print("Input file set to ", os.path.basename(path))
except:
    print("Input file does not seem to be a fastq file. Exiting")
    print(usage)
    sys.exit(1)

# check an output location is given and its location exists
if options.outfile is not None:
    print("Writing output to ", options.outfile)
    outfile_name = options.outfile
else:
    print("No output file defined. Exiting")
    print(usage)
    sys.exit(1)

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

# parse quality parameters
if options.minqual is not None:
    print("Minimum average phred quality score required to pass: " + options.minqual)
    qualityCutoff = options.minqual
else:
    print("Minimum average phred quality score required to pass: " + options.minqual)
    qualityCutoff = 30


#------------------
# Filtering function
#------------------

def filter_reads(fastq_generator, qualityCutoff):
    # collect counts of total and filtered reads
    total_count = 0
    filtered_count = 0
    filtered_N = 0
    filtered_qual = 0
        
    # main parsing block
    for record in fastq_generator:
        total_count += 1

        # check if read contains any N residues
        read = str(record.seq)
            
        # barcodes cannot contain N residues
        if "N" not in str(read.seq):
            # average quality of read must be above quality threshold

            avg_qual = sum(read.letter_annotations["phred_quality"]) / len(read)
            if avg_qual >= qualityCutoff:
                yield read
            else:
                filtered_qual += 1
        else:
            filtered_N += 1

    print("filterBarcodeReads.py")
    print("")
    print("Total reads parsed: ", total_count)
    print("Number of reads passing filters: ", filtered_count)
    print("Number of reads removed - N residues: ", filtered_N)
    print("Number of reads removed - Quality: ", filtered_qual)
    print("Percentage of reads kept: ", round((filtered_count / total_count) * 100,2))
    
        

#---------------------
# Run Filtering script
#---------------------
def parse_fastq(file, outfile_name, index, minqual):

    # get sample name from filename
    filename = os.path.basename(file)
    filename = filename.split("_")[0]
    filename = filename.split(".")[0]
    print("")
    print("Parsing", filename)
    
    # setup generator to parse fastq
    with gzip.open(file, 'rt') as fq_handle:
        generator = SeqIO.parse(fq_handle, "fastq")
        
        # parse and write filtered and trimmed reads to file
        with open(outfile_name, "w") as out_handle:
            SeqIO.write(filter_reads(generator, indexes, minqual), out_handle, "fastq")


    # gzip the outfile, can probably do this in one step?
    with open(outfile_name, 'rb') as f_in:
        with gzip.open(outfile_name + '.gz' , 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

parse_fastq(path, outfile_name, indexes, qualityCutoff)

# finish and print runtime
print("Script runtime: ", datetime.now() - startTime)
