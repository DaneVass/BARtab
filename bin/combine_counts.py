#!/usr/bin/env python3

import pandas as pd
import os
import sys

files = sys.argv[1:len(sys.argv)-1]
outfile = sys.argv[-1]

print("reading count tables")

samples = []
for count_file in files:
    # sample name is everything before last _
    sample_name = "_".join(os.path.basename(count_file).split("_")[:-1])
    # should be faster if merging on index
    samples.append(pd.read_csv(count_file, sep="\t", names=[sample_name], index_col=[0]))

print("merging count tables")

if len(samples) > 1:
    # in case of multiple samples, do an iterative outer join
    combined_counts = samples[0]
    for i in range(1, len(samples)):
        combined_counts = combined_counts.merge(samples[i], how="outer", left_index=True, right_index=True)
    combined_counts = combined_counts.fillna(0).astype(int).reset_index(names="Barcode")
else:
    combined_counts = samples[0].reset_index(names="Barcode")

print("writing counts to file")
combined_counts.to_csv(outfile, sep="\t", index=False)

print("combine counts complete!")
