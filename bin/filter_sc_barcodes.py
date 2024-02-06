#!/usr/bin/env python3
import sys
import pandas as pd

in_file = sys.argv[1]
out_file = sys.argv[2]

# in_file = "SR03-04-Pool-2_S2.mapped.sam"
# out_file = "SR03-04-Pool-2_S2_barcodes.tsv"

print("############ Clean PCR chimerism ############\n")

# only use read header and barcode column
mapping_output = pd.read_csv(in_file, sep="\t", header=None, usecols=[0,2], comment='@', names=["readname", "barcode"])
print(f"Parsed {mapping_output.shape[0]} reads")

# get cell barcode and UMI from header
mapping_output[["read", "CB", "UMI"]] = mapping_output["readname"].str.split("_", expand=True)

# only keep CB_UMI barcode combination with highest UMI count
mapping_output = mapping_output.value_counts(["barcode", "CB", "UMI"]).reset_index()
print(f"Counted {mapping_output.shape[0]} cell ID-UMI-barcode combinations ({mapping_output['CB'].unique().shape[0]} cells, {mapping_output['barcode'].unique().shape[0]} barcodes)")

idx = mapping_output.groupby(["CB", "UMI"])['count'].transform(max) == mapping_output['count']
mapping_output_max = mapping_output[idx]
print(f"Removed {mapping_output.shape[0] - mapping_output_max.shape[0]} cell ID-UMI-barcode combinations with maximum read count")
print(f"Kept {mapping_output_max.shape[0]} ({mapping_output_max['CB'].unique().shape[0]} cells, {mapping_output_max['barcode'].unique().shape[0]} barcodes)")

# remove ties by removing all duplicated CB_UMI values
mapping_output_max_no_ties = mapping_output_max.drop_duplicates(subset=["CB", "UMI"], keep=False)
print(f"Removed {mapping_output_max.shape[0] - mapping_output_max_no_ties.shape[0]} cell ID-UMI-barcode combinations with read count ties")
print(f"Kept {mapping_output_max_no_ties.shape[0]} ({mapping_output_max_no_ties['CB'].unique().shape[0]} cells, {mapping_output_max_no_ties['barcode'].unique().shape[0]} barcodes)")

# write file to pass to umi-tools count_tab
# create read name that contains a unique index, UMI and cell barcode
mapping_output_max_no_ties = mapping_output_max_no_ties.reset_index()
# it is essential to order by barcode before umi-tools count_tab
mapping_output_max_no_ties = mapping_output_max_no_ties.sort_values("barcode")
mapping_output_max_no_ties["read"] = mapping_output_max_no_ties["index"].astype(str) + "_" + mapping_output_max_no_ties["UMI"] + "_" + mapping_output_max_no_ties["CB"]

mapping_output_max_no_ties[["read", "barcode"]].to_csv(out_file, index=False, sep="\t", header=False)
