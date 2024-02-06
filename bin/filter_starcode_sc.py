#!/usr/bin/env python3

import pandas as pd
# import os
import sys

infile = sys.argv[1]
# infile = "/dawson_genomics/Projects/pfizer_lung_cancer/SR03_04/scripts/work/fb/73fe746ce7d69ae46b7239a3cc0aef/SR03-04-Pool-1_S1_unmapped_starcode.tsv"
outfile = sys.argv[2]
# length of cell barcode and UMI together
cb_umi_length = int(sys.argv[3])
# cb_umi_length = 28
# length of only the cell barcode
cb_length = int(sys.argv[4])
# cb_length = 16

print("############ Clean PCR chimerism ############\n")

# read starcode output into dataframe
starcode_output = pd.read_csv(infile, sep="\t", header=None, names=["sequence", "count"])

# split sequences into CB, UMI and lineage barcode
starcode_output["barcode"] = starcode_output["sequence"].str[cb_umi_length:]
starcode_output['CB_UMI'] = starcode_output['sequence'].str[0:cb_umi_length]
starcode_output["UMI"] = starcode_output["CB_UMI"].str[cb_length:cb_umi_length]
starcode_output["CB"] = starcode_output["CB_UMI"].str[0:cb_length]

print(f"Parsed {starcode_output.shape[0]} cell ID-UMI-barcode combinations ({starcode_output['CB'].unique().shape[0]} cells, {starcode_output['barcode'].unique().shape[0]} barcodes)")

# only keep CB_UMI barcode combination with highest UMI count
idx = starcode_output.groupby(["CB", "UMI"])['count'].transform(max) == starcode_output['count']
starcode_output_max = starcode_output[idx]
print(f"Removed {starcode_output.shape[0] - starcode_output_max.shape[0]} cell ID-UMI-barcode combinations with maximum read count")
print(f"Kept {starcode_output_max.shape[0]} ({starcode_output_max['CB'].unique().shape[0]} cells, {starcode_output_max['barcode'].unique().shape[0]} barcodes)")

# remove ties by removing all duplicated CB_UMI values
starcode_output_max_no_ties = starcode_output_max.drop_duplicates(subset=["CB", "UMI"], keep=False)
print(f"Removed {starcode_output_max.shape[0] - starcode_output_max_no_ties.shape[0]} cell ID-UMI-barcode combinations with read count ties")
print(f"Kept {starcode_output_max_no_ties.shape[0]} ({starcode_output_max_no_ties['CB'].unique().shape[0]} cells, {starcode_output_max_no_ties['barcode'].unique().shape[0]} barcodes)")

# write file to pass to umi-tools count_tab
# create read name that contains a unique index, UMI and cell barcode
starcode_output_max_no_ties = starcode_output_max_no_ties.reset_index()
# it is essential to order by barcode before umi-tools count_tab
starcode_output_max_no_ties = starcode_output_max_no_ties.sort_values("barcode")
starcode_output_max_no_ties["read"] = starcode_output_max_no_ties["index"].astype(str) + "_" + starcode_output_max_no_ties["UMI"] + "_" + starcode_output_max_no_ties["CB"]

starcode_output_max_no_ties[["read", "barcode"]].to_csv(outfile, index=False, sep="\t", header=False)
