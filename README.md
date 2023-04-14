# BARtab
A Nextflow pipeline to tabulate synthetic barcode counts from NGS data

```
    Usage: nextflow run BARtab.nf --input <input dir> 
                                  --output <output dir> 
                                  --index <path>/<prefix> 
                                  --contrasts CONTRASTS.csv 
                                  -profile local
                                  --help

    Required arguments:
      --indir                     Directory containing raw *.fastq.gz files
      --index                     Path to the bowtie2 index for the sgRNA library. Include prefix.

    Filtering arguments:
      --minqual                   Minimum PHRED quality across read.
    
    Trimming arguments:
      --error                     Proportion of mismatches allowed in constant regions.

    Mapping arguments:
      --mismatches                Number of allowed mismatches during reference mapping.

    Optional arguments:
      --contrasts                 CSV file detailing the comparisons to test [contrasts.csv]
      -profile                    Configuration profile to use. Can use multiple (comma separated)
                                         Available: local, singularity, slurm
      --outdir                    Output directory to place output [./]
      --threads                   Number of CPUs to use [4]
      --help                      Print this help statement.


    Profiles:
      local                       local execution
      slurm                       SLURM execution 

    Author:
      Dane Vassiliadis (dane.vassiliadis@petermac.org)
```
## Stages:
- Check raw data quality using `fastqc`
- [OPTIONAL] Merge paired end reads using `FLASh`
- Quality filter reads using `fastx-toolkit`
- Filter barcode reads and trim 5' and/or 3' constant regions using `cutadapt`
- Align to reference barcode library using `bowtie`
- [OPTIONAL] If no reference library, derive consensus barcode repertoire using `starcode`
- Count number of reads aligning per barcode using `samtools`
- Merge counts files for multiple samples
- Report metrics for individual samples

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

## Installing the pipeline (Conda method)
1. Install Nextflow using the instructions found here: https://www.nextflow.io/docs/latest/getstarted.html
2. Install miniconda using the instructions found here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html 
3. Install BARtab dependencies by running `conda env create -f environment.yaml`
4. Activate BARtab environment by running `conda activate BARtab_env`
5. Print BARtab help message `nextflow run BARtab.nf --help`

## Running the pipeline
1. Run the test datasets using `nextflow run BARtab.nf --indir test/dat/test_SE`
2. Setup a directory containing input fastq files.
3. Run the pipeline - `nextflow run BARtab.nf -i <input_dir> -o <output_dir> <other args>`
4. Alternatively, configure run parameters in `run_BARtab.sh` and run the pipeline using `bash run_BARtab.sh`
5. If running the pipeline on a remote HPC it is recommended to run the pipeline within a tmux session as follows
   1. Initialise a new tmux session called "bartab" by running `tmux new -s bartab`
   2. Re activate the BARtab environment within the tmux session by running `conda activate BARtab_env`
   3. Run the pipeline via steps 3. or 4. above.
   4. exit the tmux session using ^B then D (Mac) or Ctrl+B then D Windows/UNIX
   5. Continue working on something else
   6. re-enter the tmux session by running `tmux a -t bartab`


