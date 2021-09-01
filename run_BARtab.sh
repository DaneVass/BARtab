#!/bin/bash

BARtab="/home/dvassiliadis/my_scripts/nextflow_pipelines/BARtab/BARtab.nf" # path to BARtab.nf

indir = "test/dat" # path to folder containing input fastq files
outdir = "test/test_out" # path to desired output directory
index = "$PWD/test/ref/mCHERRY_barcode" # path to bowtie index for reference genome (include index prefix)
ref = "test/ref/mCHERRY_barcode_reference_library.fasta" # path to reference genome fasta file
upconstant = "CGATTGACTA" // SPLINTR 1st gen upstream constant region # upstream constant 
downconstant = "TGCTAATGCG" // SPLINTR 1st gen downstream constant region # downstream constant
alnmismatches = 1 # allowed mismatches during bowtie alignment
merge = false # merge paired end reads before extracting barcodes?
minqual = 20 # minimum quality required across barcode
pctqual = 80 # percent of barcode bases that must reach quality threshold
constants = "up" # constant regions to use for barcode extraction
constantmismatches = 0.1 # proportion of mismatches allowed in constant regions
#email = "" # user email address (pipeline updates will be sent here)

nextflow run $BARtab --indir ${indir} \
    --outdir ${outdir} \
    --index ${index} \
    --ref ${ref} \
    --upconstant ${upconstant} \
    --downconstant ${downconstant} \
    --alnmismatches ${alnmismatches} \
    --minqual ${minqual} \
    --pctqual ${pctqual} \
    --constants = ${constants} \
    --constantmismatches ${constantmismatches}

