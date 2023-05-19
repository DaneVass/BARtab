#!/usr/bin/env bash

out=$1 
echo "Nextflow version:" >> $out
nextflow -v >> $out

echo "------------------------------------" >> $out
echo "fastQC version:" >> $out
fastqc --version >> $out

echo "------------------------------------" >> $out
echo "FLASh version:" >> $out
flash -v >> $out

echo "------------------------------------" >> $out
echo "fastx-toolkit - fastq_quality_filter version:" >> $out
fastq_quality_filter -h >> $out

echo "------------------------------------" >> $out
echo "bowtie -version" >> $out
bowtie --version >> $out

echo "------------------------------------" >> $out
echo "cutadapt version:" >> $out 
cutadapt --version >> $out

echo "------------------------------------" >> $out
echo "starcode version:" >> $out
starcode --version &>> $out

echo "------------------------------------" >> $out
echo "R version:" >> $out
R --version >> $out

echo "------------------------------------" >> $out
echo "Python version:" >> $out
python --version >> $out

echo "------------------------------------" >> $out
echo "samtools version:" >> $out
samtools --version >> $out

echo "------------------------------------" >> $out
echo "UMI-tools version:" >> $out
umi_tools --version >> $out

echo "------------------------------------" >> $out
echo "MultiQC version:" >> $out
multiqc --version >> $out
