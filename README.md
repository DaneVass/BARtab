# BARtab
A Nextflow pipeline to tabulate synthetic barcode counts from NGS data

  Usage: nextflow run BARtab.nf --input <input dir> 
                                --output <output dir> 
                                --index <path>/<prefix> 
                                --contrasts CONTRASTS.csv 
                                -profile local
                                --help

    Required arguments:
      --indir                        Directory containing raw *.fastq.gz files
      --index                        Path to the bowtie2 index for the sgRNA library. Include prefix.

    Filtering arguments:
      --minqual                      Minimum PHRED quality across read.
    
    Trimming arguments:
      --error                        Proportion of mismatches allowed in constant regions.

    Mapping arguments:
      --mismatches                   Number of allowed mismatches during reference mapping.

    Optional arguments:
      --contrasts                    CSV file detailing the comparisons to test [contrasts.csv]
      -profile                       Configuration profile to use. Can use multiple (comma separated)
                                            Available: local, singularity, slurm
      --outdir                       Output directory to place output [./]
      --threads                      Number of CPUs to use [4]
      --help                         Print this help statement.


    Profiles:
      local                          local execution
      singularity                    local execution with singularity
      slurm                          SLURM execution 

    Author:
      Dane Vassiliadis (dane.vassiliadis@petermac.org)
    
## Stages:
- Check raw data quality using `fastqc`
- [OPTIONAL] merge paired end reads using `FLASh`
- Quality filter reads using `fastx-toolkit`
- filter barcode reads and trim 5' and/or 3' constant regions using `cutadapt`
- align to reference barcode library using `bowtie`
- [OPTIONAL] if no reference library, derive consensus barcode repertoire using `starcode`
- count number of reads aligning per barcode using `samtools`
- merge counts files for multiple samples
- report metrics for individual samples

## Dependiencies
* [Nextflow](https://bitbucket.org/snakemake/snakemake)
* [R > 3.3.0](https://www.r-project.org/)
    * The [tidyverse package](https://www.tidyverse.org/)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)
* [Starcode](https://github.com/gui11aume/starcode)
* [FLASh](http://ccb.jhu.edu/software/FLASH/)
* [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [bowtie1](http://bowtie-bio.sourceforge.net/index.shtml)
* [samtools](http://www.htslib.org/)

## Installing the pipeline
1. `conda create -f BARtab_environment.yaml`
2. `conda activate BARtab`
3. `nextflow run BARtab.nf --help`

## Running the pipeline

1. Run the test datasets - `nextflow run BARtab.nf --indir test`
2. Symbolic link fastq files into required directory. Or specify input and output directories.
3. Run the pipeline - `nextflow run BARtab.nf -i <input_dir> -o <output_dir> <other args>`


