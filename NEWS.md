# BARtab 1.4.1
- filter clustered barcodes from single-cell data: for each cell barcode - UMI combination, only keep lineage barcode with the most reads; remove all ambiguous ties.

# BARtab 1.3.1

- Clustering unmapped reads
    - always save unmapped reads to file
    - lower resource usage for starcode on unmapped reads
    - trim reads before running starcode on unmapped reads (if only upstream or only downstream adapter trimmed)

- add NEWS.md to track changes
- fix resource labels in slurm and lsf config
