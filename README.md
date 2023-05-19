# BARtab
A Nextflow pipeline to tabulate synthetic barcode counts from NGS data

```
  Usage: nextflow run danevas/bartab --indir <input dir> 
                                --outdir <output dir> 
                                --ref <path/to/reference/fasta> 
                                --mode <single-bulk | paired-bulk | single-cell>

    Input arguments:
      --input                    Directory containing input *.fastq.gz files. Must contain R1 and R2 if running in mode paired-bulk or single-cell.
                                        For single-cell mode, a BAM file can be provided instead (see --bam)
      --ref                      Path to a reference fasta file for the barcode / sgRNA library.
                                        If null, reference-free workflow will be used for single-bulk and paired-bulk modes.
      --mode                     Workflow to run. <single-bulk,paired-bulk,single-cell>

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
      -profile                   Configuration profile to use. Can use multiple (comma separated)
                                        Available: conda, singularity, docker, slurm
      --outdir                   Output directory to place output [default = './']
      --threads                  Number of CPUs to use [default = 4]
      --email                    Direct output messages to this address [default = '']
      --help                     Print this help statement.

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
See [citations](../CITATIONS.md)
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
1. Install Nextflow using the instructions found [here](https://www.nextflow.io/docs/latest/getstarted.html) (and [here](https://www.nextflow.io/blog/2021/nextflow-developer-environment.html))
    ```
    # download the executable
    curl get.nextflow.io | bash
    # move the nextflow file to a directory accessible by your $PATH variable
    sudo mv nextflow /usr/local/bin
    ```

2. Try out the pipeline   
    (this will automatically [pull](https://www.nextflow.io/docs/latest/sharing.html#pulling-or-updating-a-project) the pipeline, usually into `~/.nextflow/assets/`)
    ```
    nextflow run danevas/bartab --help
    ```
    Alternatively, [clone](https://www.nextflow.io/docs/latest/sharing.html#cloning-a-project-into-a-folder) the pipeline into a directory of your choice first with
    ```
    nextflow clone danevas/bartab target_dir/
    ```

2. Install dependencies
    ### Conda
    1. Install miniconda using the instructions found here: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html 
    3. It is recommended to use mamba to create the conda environment `conda install -c conda-forge mamba`
    4. Install BARtab dependencies by running `mamba env create -f environment.yaml` or `conda env create -f environment.yaml`
    5. Run the pipeline with `nextflow run danevas/bartab -profile conda [options]` 
    
    The location of the conda environment is specified in `conf/conda.config`.

    ### Docker
    
    ### Singularity

## Running the pipeline
Print the help message with `nextflow run danevas/bartab --help`

Run any of the test datasets using `nextflow run danevas/bartab -profile <test_SE,test_PE, test_SE_ref_free,test_sc,test_sc_bam>,<conda,docker,singularity>`

Run the pipeline with your own data
- Single-end bulk workflow: `nextflow run danevas/bartab --mode single-bulk --indir <input_dir> --outdir <output_dir> --ref <reference> [options]`
- Paired-end bulk workflow: `nextflow run danevas/bartab --mode paired-bulk --indir <input_dir> --outdir <output_dir> --ref <reference> [options]`
- Reference-free single-end bulk workflow: `nextflow run danevas/bartab --mode single-bulk --indir <input_dir> --outdir <output_dir> [options]`
- Single-cell workflow with fastq input: `nextflow run danevas/bartab --mode single-cell --indir <input_dir> --outdir <output_dir> --ref <reference> [options]`
- Single-cell workflow with BAM file input: `nextflow run danevas/bartab --mode single-cell --bam <input_dir> --outdir <output_dir> --ref <reference> [options]`

Use `-w` to specify the location of the work directory and `-resume` when only parts of the input have changed or only a subset of process has to be re-run. 