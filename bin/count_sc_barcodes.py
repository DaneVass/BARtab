#!/usr/bin/env python3
import sys
import pandas as pd

in_file = sys.argv[1]
out_file = sys.argv[2]

print("\n\n############ Count UMI ############\n")

# read output of umi-tools count_tab
umitools_demultiplexed = pd.read_csv(in_file, sep="\t")

# fix weird string encoding of cell barcodes
umitools_demultiplexed["cell"] = umitools_demultiplexed["cell"].str.replace("b", "").str.replace("\'", "")

# count # UMIs supporting barcode in cell
# umitools count_tab does not have --per_contig option like umitools count. Therefore need to sum. 
umitools_demultiplexed_counts = umitools_demultiplexed.groupby(["cell", "gene"]).sum().reset_index()
print(f"Counted {umitools_demultiplexed_counts['gene'].unique().shape[0]} barcodes in {umitools_demultiplexed_counts['cell'].unique().shape[0]} cells")

# write to file
umitools_demultiplexed_counts[["gene", "cell", "count"]].to_csv(out_file, index=False, sep="\t")
