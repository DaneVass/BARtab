#!/bin/bash

BARtab="/home/dvassiliadis/my_scripts/nextflow_pipelines/BARtab/BARtab.nf" # FULL path to BARtab.nf
indir="" # FULL path to folder containing input fastq files
outdir="" # FULL path to desired output directory
ref="" # FULL path to reference genome fasta file
upconstant="CGATTGACTA" #SPLINTR 1st gen upstream constant region # upstream constant 
downconstant="TGCTAATGCG" #SPLINTR 1st gen downstream constant region # downstream constant
alnmismatches=1 # allowed mismatches during bowtie alignment
merge=false # merge paired end reads before extracting barcodes?
minqual=20 # minimum quality required across barcode
pctqual=80 # percent of barcode bases that must reach quality threshold
constants="up" # constant regions to use for barcode extraction
constantmismatches=0.1 # proportion of mismatches allowed in constant regions
#email="" # user email address (pipeline updates will be sent here)

nextflow run $BARtab --indir ${indir} \
    --outdir ${outdir} \
    --ref ${ref} \
    --upconstant ${upconstant} \
    --downconstant ${downconstant} \
    --alnmismatches ${alnmismatches} \
    --minqual ${minqual} \
    --pctqual ${pctqual} \
    --constants ${constants} \
    --constantmismatches ${constantmismatches} \
    -profile slurm # slurm based execution

