#!/usr/bin/env python3
import sys
import pandas as pd

in_file = sys.argv[1]
out_file = sys.argv[2]

sam = pd.read_csv(in_file, sep="\t", header=None, usecols=[0,2], comment='@', names=["readname", "barcode"])
sam[["read", "CID", "MID"]] = sam["readname"].str.split("\|", expand=True)

sam = sam[["CID", "MID", "barcode"]].drop_duplicates().reset_index(drop=True)

counts = sam.groupby(by=["CID", "barcode"]).count().sort_values("MID", ascending=False)

# create same format as umitools count
counts = counts.reset_index()[["barcode", "CID", "MID"]]
counts = counts.rename(columns={'barcode': 'gene', "CID": "cell", "MID": "count"})
                   
counts.to_csv(out_file, sep="\t", index=False)
