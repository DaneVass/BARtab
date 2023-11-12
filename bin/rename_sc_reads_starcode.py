#!/usr/bin/env python3
import sys
import gzip

in_file = sys.argv[1]
out_file = sys.argv[2]

# run with e.g. rename_sc_reads_starcode /researchers/henrietta.holze/splintr_tools/BARtab/test/test_out/sc/trimmed_reads/test_sc.trimmed.fastq.gz /researchers/henrietta.holze/splintr_tools/BARtab/test/test_out/sc/trimmed_reads/test_sc.trimmed.fasta

i = 0
read_chunk = []

new_fasta = []

with gzip.open(in_file, 'r') as f:
	for line in f:
		i += 1
		read_chunk.append(line.strip().decode())
		if i % 4 == 0:
			read_chunk[1] = "".join(["".join(read_chunk[0].split()[0].split("_")[-2:]), read_chunk[1]])
			fasta_header = ">" + read_chunk[0][1:]
			# write read to fasta file
			new_fasta.extend([fasta_header, read_chunk[1]])
			read_chunk = []

with open(out_file, 'w') as f:
	f.write("\n".join(new_fasta))    
