# BARtab 1.4.0

- Remove PCR chimerism in single-cell data: for each cell barcode - UMI combination, only keep lineage barcode with the most reads; remove all ambiguous ties
- Replace fastx_toolkit with fastp for filtering reads based on quality
- Add option for filtering reads: fastp low complexity threshold
- Filter R2 in single-cell data with fastp
- Update conda environment and Docker image with fastp dependency
- Run starcode_umi, not starcode, when clustering unmapped reads from single-cell data
- Update and improve documentation

# BARtab 1.3.1

- Clustering unmapped reads
    - always save unmapped reads to fastq file with bowtie
    - lower resource usage for starcode on unmapped reads
    - trim reads before running starcode on unmapped reads (if only upstream or only downstream adapter trimmed)

- add NEWS.md to track changes
- fix resource labels in slurm and lsf config
