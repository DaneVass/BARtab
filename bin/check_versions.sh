#!/usr/bin/env bash

out=$1 
echo "Nextflow version:" >> $out
nextflow -v >> $out
  
echo "------------------" >> $out
echo "fastQC version:" >> $out
fastqc --version >> $out
  
echo "------------------" >> $out
echo "FLASh version:" >> $out
flash -v >> $out
  
echo "------------------" >> $out
echo "fastx-toolkit - fastq_quality_filter version:" >> $out
fastq_quality_filter -h >> $out
  
echo "------------------" >> $out
echo "bowtie -version" >> $out
bowtie --version >> $out
  
echo "------------------" >> $out
echo "cutadapt version:" >> $out 
cutadapt --version >> $out
  
echo "------------------" >> $out
echo "starcode version:" >> $out
sc=$(starcode --version)
echo $sc >> $out
  
echo "------------------" >> $out
echo "R version:" >> $out
R --version >> $out
