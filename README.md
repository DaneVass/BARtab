# BARtab
A Nextflow pipeline to tabulate synthetic barcode counts from NGS data

```
  Usage: nextflow run BARtab.nf --indir <input dir> 
                                --outdir <output dir> 
                                --ref <path/to/reference/fasta> 
                                --mode <single-bulk | paired-bulk | single-cell> 
                                -profile local
                                --help

    Input arguments:
      --input                    Directory containing input *.fastq.gz files. Must contain R1 and R2 if running in mode paired-bulk or single-cell.
                                        For single-cell mode, a BAM file can be provided instead (see --bam)
      --ref                      Path to a reference fasta file for the barcode / sgRNA library.
                                        If null, reference-free workflow will be used for single-bulk and paired-bulk modes.
      --mode                     Workflow to run. <single-bulk, paired-bulk, single-cell> [default = 'single-bulk']

    Read merging arguments:
      --mergeoverlap             Length of overlap required to merge paired-end reads [default = 10]

    Filtering arguments:
      --minqual                  Minimum PHRED quality per base [default = 20]
      --pctqual                  Percentage of bases within a read that must meet --minqual [default = 80]

    Trimming arguments:
      --constants                Which constant regions flanking barcode to search for in reads <up, down, both> [default = 'up']
      --upconstant               Sequence of upstream constant region [default = 'CGATTGACTA'] // SPLINTR 1st gen upstream constant region
      --downconstant             Sequence of downstream constant region [default = 'TGCTAATGCG'] // SPLINTR 1st gen downstream constant region
      --constantmismatches       Proportion of mismatched bases allowed in constant regions [default = 0.1]
      --min_readlength           Minimum read length [default = 15]

    Mapping arguments:
      --alnmismatches            Number of allowed mismatches during reference mapping [default = 1]

    Sincle-cell arguments:
      --bam                      Path to BAM file output of Cell Ranger, containing reads that do not map to the reference genome. Only permitted in single-cell mode
      --cellnumber               Number of cells expected in sample, only when no BAM provided [default = 5000]
      --umi_dist                 Hamming distance between UMIs to be collapsed during counting [default = 1]

    Optional arguments:
      -profile                   Configuration profile to use. Can use multiple (comma separated) [default = 'local']
                                        Available: local, singularity, slurm
      --outdir                   Output directory to place output [default = './']
      --threads                  Number of CPUs to use [default = 4]
      --email                    Direct output messages to this address [default = '']
      --help                     Print this help statement.

    Profiles:
      local                      local execution
      singularity                use singularity container
      slurm                      SLURM execution 

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
1. Run the test datasets using `nextflow run BARtab.nf --indir test/dat`
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


