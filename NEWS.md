# BARtab 1.4.0

- Replace fastx_toolkit with fastp for filtering reads based on quality
- Update conda environment and Docker image with fastp dependency
- run starcode_umi, not starcode, when clustering unmapped reads from single-cell data. 
- filter clustered barcodes from single-cell data: for each cell barcode - UMI combination, only keep lineage barcode with the most reads; remove all ambiguous ties.

# BARtab 1.3.1

- Clustering unmapped reads
    - always save unmapped reads to file
    - lower resource usage for starcode on unmapped reads
    - trim reads before running starcode on unmapped reads (if only upstream or only downstream adapter trimmed)

- add NEWS.md to track changes
- fix resource labels in slurm and lsf config
