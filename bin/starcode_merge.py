#!/usr/bin/env python3

"""
Henrietta Holze
Peter MacCallum Cancer Centre

starcode_merge.py

Merges multiple starcode output files into one.

Input: file paths of starcode output, output file path

Run e.g. "starcode_merge.py test1_starcode.txt test2_starcode.txt all_counts_combined.txt"
"""

import sys

files = sys.argv[1:len(sys.argv)-1]
outfile = sys.argv[-1]

# read every file with barcode and readcount
# increase readcount of barcode or add new barcode
barcode_counts = {}
for file in files:
    with open(file) as f:
        rows = [line.strip().split('\t') for line in f]
        for row in rows:
            if row[0] not in barcode_counts:
                barcode_counts[row[0]] = 0
            barcode_counts[row[0]] += int(row[1])

# write dict to file
with open(outfile, 'w') as f:
    for key, value in barcode_counts.items():
        f.write(f"{key}\t{value}\n")
