#!/usr/bin/env python3
import sys
import pandas as pd

in_file = sys.argv[1]
out_file = sys.argv[2]
# input is either sam file with alignments or output of starcode_umi
input_type = sys.argv[3]
if input_type == "sam":
    # cell barcode UMI pattern
    delimiter = sys.argv[4]
if input_type == "starcode_umi":
    # cell barcode UMI pattern
    cb_umi_pattern = sys.argv[4]
    delimiter = "_"

print("############ Clean PCR chimerism ############\n")

######### parse input files

if input_type == "starcode_umi":
    # infile = "/dawson_genomics/Projects/pfizer_lung_cancer/SR03_04/scripts/work/fb/73fe746ce7d69ae46b7239a3cc0aef/SR03-04-Pool-1_S1_unmapped_starcode.tsv"
    # cb_umi_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNNNN"

    # length of cell barcode and UMI together
    cb_umi_length = len(cb_umi_pattern)
    # length of only the cell barcode
    cb_length = len(cb_umi_pattern.strip("N"))

    print("Parsing starcode_umi output\n")
    # read starcode output into dataframe
    starcode_output = pd.read_csv(in_file, sep="\t", header=None, names=["sequence", "count"])

    # split sequences into CB, UMI and lineage barcode
    starcode_output["barcode"] = starcode_output["sequence"].str[cb_umi_length:]
    starcode_output['CB_UMI'] = starcode_output['sequence'].str[0:cb_umi_length]
    starcode_output["UMI"] = starcode_output["CB_UMI"].str[cb_length:cb_umi_length]
    starcode_output["CB"] = starcode_output["CB_UMI"].str[0:cb_length]

    output = starcode_output[["CB", "UMI", "barcode", "count"]]

if input_type == "sam":
    # in_file = "SR03-04-Pool-2_S2.mapped.sam"
    # out_file = "SR03-04-Pool-2_S2_barcodes.tsv"
    print("Parsing sam file\n")

    # only use read header and barcode column
    mapping_output = pd.read_csv(in_file, sep="\t", header=None, usecols=[0,2], comment='@', names=["readname", "barcode"])
    print(f"Parsed {mapping_output.shape[0]} reads")

    # get cell barcode and UMI from header
    mapping_output[["read", "CB", "UMI"]] = mapping_output["readname"].str.split(delimiter, expand=True)
    mapping_output = mapping_output.value_counts(["barcode", "CB", "UMI"]).reset_index()

    output = mapping_output[["CB", "UMI", "barcode", "count"]]


# only keep CB_UMI barcode combination with highest UMI count
print(f"Counted {output.shape[0]} cell ID-UMI-barcode combinations ({output['CB'].unique().shape[0]} cells, {output['barcode'].unique().shape[0]} barcodes)")

idx = output.groupby(["CB", "UMI"])['count'].transform(max) == output['count']
output_max = output[idx]
print(f"Removed {output.shape[0] - output_max.shape[0]} cell ID-UMI-barcode combinations with maximum read count")
print(f"Kept {output_max.shape[0]} ({output_max['CB'].unique().shape[0]} cells, {output_max['barcode'].unique().shape[0]} barcodes)")

# remove ties by removing all duplicated CB_UMI values
output_max_no_ties = output_max.drop_duplicates(subset=["CB", "UMI"], keep=False)
print(f"Removed {output_max.shape[0] - output_max_no_ties.shape[0]} cell ID-UMI-barcode combinations with read count ties")
print(f"Kept {output_max_no_ties.shape[0]} ({output_max_no_ties['CB'].unique().shape[0]} cells, {output_max_no_ties['barcode'].unique().shape[0]} barcodes)")

# write file to pass to umi-tools count_tab
# create read name that contains a unique index, UMI and cell barcode
output_max_no_ties = output_max_no_ties.reset_index()
# it is essential to order by barcode before umi-tools count_tab
output_max_no_ties = output_max_no_ties.sort_values("barcode")
output_max_no_ties["read"] = output_max_no_ties["index"].astype(str) + delimiter + output_max_no_ties["UMI"] + delimiter + output_max_no_ties["CB"]

output_max_no_ties[["read", "barcode"]].to_csv(out_file, index=False, sep="\t", header=False)
