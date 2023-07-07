# BARtab
A Nextflow pipeline to tabulate synthetic barcode counts from NGS data

```
  Usage: nextflow run danevas/bartab --indir <input dir>
                                     --outdir <output dir>
                                     --ref <path/to/reference/fasta>
                                     --mode <single-bulk | paired-bulk | single-cell>

    Input arguments:
      --indir                    Directory containing input *.fastq.gz files. Must contain R1 and R2 if running in mode paired-bulk or single-cell.
                                        For single-cell mode, directory can contain BAM files.
      --input_type               Input file type, either fastq or bam, only relevant for single-cell mode [default = fastq]
      --ref                      Path to a reference fasta file for the barcode / sgRNA library.
                                        If null, reference-free workflow will be used for single-bulk and paired-bulk modes.
      --mode                     Workflow to run. <single-bulk, paired-bulk, single-cell>

    Read merging arguments:
      --mergeoverlap             Length of overlap required to merge paired-end reads [default = 10]

    Filtering arguments:
      --minqual                  Minimum PHRED quality per base [default = 20]
      --pctqual                  Percentage of bases within a read that must meet --minqual [default = 80]

    Trimming arguments:
      --constants                Which constant regions flanking barcode to search for in reads: up, down or both. "all" runs all 3 modes and combines the results. 
                                 Single-cell mode always runs with "all". <up, down, both, all> [default = 'up']
      --upconstant               Sequence of upstream constant region [default = 'CGATTGACTA'] // SPLINTR 1st gen upstream constant region
      --downconstant             Sequence of downstream constant region [default = 'TGCTAATGCG'] // SPLINTR 1st gen downstream constant region
      --constantmismatches       Proportion of mismatched bases allowed in constant regions [default = 0.1]
      --min_readlength           Minimum read length [default = 20]

    Mapping arguments:
      --alnmismatches            Number of allowed mismatches during reference mapping [default = 2]

    Sincle-cell arguments:
      --cb_umi_pattern           Cell barcode and UMI pattern on read 1, required for fastq input. N = UMI position, C = cell barcode position [defauls = CCCCCCCCCCCCCCCCNNNNNNNNNNNN]
      --cellnumber               Number of cells expected in sample, only when no BAM provided [default = 5000]
      --umi_dist                 Hamming distance between UMIs to be collapsed during counting [default = 1]
      --umi_count_filter         Minimum number of UMIs per barcode per cell [default = 1]
      --umi_fraction_filter      Minimum fraction of UMIs per barcode per cell compared to dominant barcode in cell (barcode supported by most UMIs) [default = 0.3]

    Resources:
      --max_cpus                 Maximum number of CPUs [default = 6]
      --max_memory               Maximum memory [default = "14.GB"]
      --max_time                 Maximum time [default = "40.h"]

    Optional arguments:
      -profile                   Configuration profile to use. Can use multiple (comma separated)
                                        Available: conda, singularity, docker, slurm
      --outdir                   Output directory to place output [default = './']
      --email                    Direct output messages to this address [default = '']
      --help                     Print this help statement.

    Author:
      Dane Vassiliadis (dane.vassiliadis@petermac.org)
      Henrietta Holze (henrietta.holze@petermac.org)
```

## Pipeline summary 

The pipeline can extract barcode counts from bulk or single-cell RNA-seq data. 
For bulk RNA-seq data, paired-end or single-end fastq files can be provided. BARtab can perform reference-free barcode extraction or perform alignment to a reference. 
Single-cell data can be provided as either BAM files containing reads that do not map to the reference or fastq files.

### Bulk workflow

The bulk workflow is executed with mode `single-bulk` and `paired-bulk` for single-end or paired-end reads, respectively. 

- Check raw data quality using `fastqc` [FASTQC](#fastqc)
- [Paired-end] Merge paired end reads using `FLASh` [MERGE_READS](#merge_reads)
- Quality filter reads using `fastx-toolkit` [FILTER_READS](#filter_reads)
- Filter barcode reads and trim 5' and/or 3' constant regions using `cutadapt` [CUTADAPT_READS](#cutadapt_reads)
- [Reference-based] Align to reference barcode library using `bowtie` [BUILD_BOWTIE_INDEX](#build_bowtie_index), [BOWTIE_ALIGN](#bowtie_align)
- [Reference-based optional] Filter alignments for sequences mapping to either end of a barcode [FILTER_ALIGNMENTS](#filter_alignments)
- [Reference-based] Count number of reads aligning per barcode using `samtools` [SAMTOOLS](#samtools), [GET_BARCODE_COUNTS](#get_barcode_counts)
- [Reference-free] If no reference library, derive consensus barcode repertoire using `starcode` [STARCODE](#starcode)
- Merge counts files for multiple samples [COMBINE_BARCODE_COUNTS](#combine_barcode_counts)
- Report metrics for individual samples [MULTIQC](#multiqc)

### Single-cell workflow
The single-cell workflow either expects fastq files or a BAM files as input. 

Fastq files must match the regex `*_R{1,2}*.{fastq,fq}.gz`.

Alternatively, if raw data was already processed with Cell Ranger, BAM files can be used as input. 
This way, cell calling and UMI extraction can be skipped. 
Unmapped reads can be extracted from the BAM file produced by cellranger with  
`samtools view -b -f 4 <sample_id>/outs/possorted_genome_bam.bam > <sample_id>_unmapped_reads.bam`.  
All BAM files can then be symlinked to an input directory and the parameter `input_type` set to `bam`.

- [fastq] Check raw data quality using `fastqc` [FASTQC](#fastqc)
- [fastq] Extraction of cell barcodes and UMIs using `umi-tools` [UMITOOLS_WHITELIST](#umitools_whitelist), [UMITOOLS_EXTRACT](#umitools_extract)
- [BAM] Filter reads containing cell barcode and UMI and convert to fastq using `samtools` [PROCESS_BAM](#process_bam)
- Filter barcode reads and trim 5' and/or 3' constant regions using `cutadapt` [CUTADAPT_READS](#cutadapt_reads)
- Align to reference barcode library using `bowtie` [BUILD_BOWTIE_INDEX](#build_bowtie_index), [BOWTIE_ALIGN](#bowtie_align)
- [Optional] Filter alignments for sequences mapping to either end of a barcode [FILTER_ALIGNMENTS](#filter_alignments)
- Extract barcode counts using `umi-tools` [SAMTOOLS](#samtools), [UMITOOLS_COUNT](#umitools_count)
- Filter and tabulate barcodes per cell and produce QC plots [PARSE_BARCODES_SC](#parse_barcodes_sc)


## Dependiencies
See [citations](../CITATIONS.md)
* [Nextflow](https://www.nextflow.io/)
* [R](https://www.r-project.org/)
    * The [tidyverse package](https://www.tidyverse.org/)
* [Python](https://www.python.org/)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FLASh](http://ccb.jhu.edu/software/FLASH/)
* [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)
* [bowtie1](http://bowtie-bio.sourceforge.net/index.shtml)
* [samtools](http://www.htslib.org/)
* [Starcode](https://github.com/gui11aume/starcode)
* [MultiQC](https://multiqc.info/)
* [umi-tools](https://github.com/CGATOxford/UMI-tools)
* [parallel](https://www.gnu.org/software/parallel/)

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
    4. Install BARtab dependencies by running `mamba env create -f environment.yaml` (or `conda env create -f environment.yaml`)
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
- Single-cell workflow with BAM file input: `nextflow run danevas/bartab --mode single-cell --indir <input_dir> --input_type bam --outdir <output_dir> --ref <reference> [options]`

Use `-w` to specify the location of the work directory and `-resume` when only parts of the input have changed or only a subset of process has to be re-run. 

It is recommended to have a look at the log files!

## Module Descriptions

### SOFTWARE_CHECK

Software check is always performed as first module.  
Output files:
- `reports/software_check.txt`: Report of all software versions.

### FASTQC
QC of fastq files is performed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

Output files:
- `qc/<sample_id>/<sample_id>.html`: html report

### MERGE_READS
If running in mode paired-bulk, forward and reverse reads are merged using [FLASh](http://ccb.jhu.edu/software/FLASH/).  
The minimum overlap of reads can be specified with the parameter `mergeoverlap` (default 10 bases).

Output files:
- `merged_reads/<sample_id>/<sample_id>.extendedFrags.fastq.gz`: merged reads  
- `merged_reads/<sample_id>/<sample_id>.notCombined_<1,2>.fastq.gz`: reads that could not be merged 
- `merged_reads/<sample_id>/<sample_id>.flash.log`: log
- `merged_reads/<sample_id>/<sample_id>.hist`: Numeric histogram of merged read lengths.
- `merged_reads/<sample_id>/<sample_id>.histogram`: Visual histogram of merged read lengths.

### FILTER_READS

Reads are quality filtered using [fastx-toolkit](http://hannonlab.cshl.edu/fastx_toolkit/) `fastq_quality_filter` command.  

The minimum quality score to keep can be specified with the parameter `minqual`. 
The minimum percent of bases that must have `minqual` quality can be specified with the parameter `pctqual`.

Output files:
- `filtered_reads/<sample_id>.filtered.fastq.gz`: filtered reads
- `filtered_reads/<sample_id>.filter.log`: log

### CUTADAPT_READS

Adapter sequences are trimmed and reads are filtered for length and N bases using [cutadapt](https://cutadapt.readthedocs.io/en/stable/).

Constants can be specified with the parameters `upconstant` and `downconstant`.  
In bulk mode, reads can be filtered for containing either upconstant `up`, downconstant `down` or both `both` with the parameter `constants`.  
In single-cell mode or when `contstants` is set to `all`, reads are filtered in all three ways. Fastq and log files are merged.

The minimum read length can be specified with `min_readlength` (default 20).  
If a constant barcode length is set with `barcode_length`, this is set as maximum sequence length. 
For `both`, only sequences matching exactly `barcode_length` will be retained.  
The fraction of mismatches in the constant region can be specified with `constantmismatches` (default 0.1).

Output files: 
- `trimmed_reads/<sample_id>.trimmed.fastq`: filtered and trimmed reads
- `trimmed_reads/<sample_id>.cutadapt.log`: log

### BUILD_BOWTIE_INDEX

If a reference is provided, it is indexed using [bowtie1](http://bowtie-bio.sourceforge.net/index.shtml).

### BOWTIE_ALIGN
If a reference is provided, trimmed and filtered reads are aligned to the indexed reference using [bowtie1](http://bowtie-bio.sourceforge.net/index.shtml).

`--norc` is specified, bowtie will not attempt to align against the reverse-complement reference strand. 
Only non-ambiguous alignments are reported with the flags `-a --best --strata -m1`. 
Sequences that map with the same number of mismatches to multiple barcodes will be discarded. 
The number of allowed mismatches can be specified with the parameter `alnmismatches` (default 1).

Output files:
- `mapped_reads/<sample_id>.mapped.sam`: Aligned reads
- `mapped_reads/<sample_id>.bowtie.log`: log

### FILTER_ALIGNMENTS

If the barcodes have a consistent length specified with `barcode_length`, alignments to the middle of a barcode sequence are filtered out.
Alignments that start at the first position or end at the last are retained.  
This ensures confidence in barcodes detected with short mapping sequences (`min_readlength`).

Output files:
- `mapped_reads/<sample_id>.mapped_filtered.sam`: Aligned and filtered reads

### SAMTOOLS

The SAM file of aligned barcode reads is sorted, indexed and compressed using [samtools](http://www.htslib.org/).

Output files:
- `mapped_reads/<sample_id>.mapped.bam`
- `mapped_reads/<sample_id>.mapped.bam.bai`

### GET_BARCODE_COUNTS

Barcode counts for each sample are extracted with [samtools](http://www.htslib.org/) `indexstats`. 

Output files:
- `counts/<sample_id>_rawcounts.txt`: tsv containing barcode and count

### STARCODE

If no reference is provided, the consensus barcode repertoire is derived using [starcode](https://github.com/gui11aume/starcode).  
Starcode clusters the filtered and trimmed barcode sequences based on their Levenshtein distance. The maximum distance by default is `min(8, 2 + [median seq length]/30)`. 

Output files:
- `starcode/<sample_id>_starcode.tsv`: barcode counts with sequence of centroid of each barcode cluster and read count

### COMBINE_BARCODE_COUNTS

Barcode counts of all samples are combined into one table with an outer join.

Output files:
- `counts/all_counts_combined.tsv`: table of barcodes and counts for each sample

### MULTIQC

MultiQC aggregates 

### UMITOOLS_WHITELIST

Cell barcodes are identified in R1 using [umi-tools whitelist](https://umi-tools.readthedocs.io/en/latest/reference/whitelist.html).

The expected number of cells should be specified with the parameter `cellnumber`.  
This should be approximately the number of cells loaded. The command is only utilized to extract cell barcodes and not to perform cell calling. 
Cell calling should be done with tools such as [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). 
Barcodes identified in droplets that do not contain cells or doublets will be removed when merging the barcode counts table with e.g. QC'd Seurat object.  
If the number of cells loaded differs a lot between samples, they must be processed separately with adjusted `cellnumber` values. 

Output files:
- `extract/<sample_id>_whitelist.tsv`: whitelisted cell barcodes and counts
- `extract/<sample_id>_whitelist.log`: log

### UMITOOLS_EXTRACT

Reads that contain cell barcode and UMI are extracted using [umi-tools extract](https://umi-tools.readthedocs.io/en/latest/reference/extract.html).

Output files:
- `extract/<sample_id>_R2_extracted.fastq`: reads that contain cell barcode and UMI, both added to the read name
- `extract/<sample_id>_exctract.log`: log

### PROCESS_BAM

Reads are filtered for flags CB and UB to obtain reads that contain a cell barcode and UMI.
At a later step (for efficiency), cell ID and UMI are added to the read headers with the module RENAME_READS.

Output files:
- `process_bam/<sample_id>_R2.fastq.gz`: reads containing cell barcode and UMI
- `process_bam/<sample_id>.filtered.bam`: reads containing cell barcode and UMI

### UMITOOLS_COUNT

Trimmed, filetered and aligned barcodes are counted using [umi-tools count](https://umi-tools.readthedocs.io/en/latest/reference/count.html#).

The Hamming distance between UMIs to be collapsed within cells during counting can be specified with parameter `umi_dist` (default 1). Collapsing barcodes can lower the number of UMIs supporting each barcode. 

Output files:
- `counts/<sample_id>.counts.tsv`: barcode counts with columns barcode, cell barcode and deduplicated UMI count
- `counts/<sample_id>_counts.log`: log

### PARSE_BARCODES_SC

Since multiple barcodes can be detected in a cell, the counts table generated by `umi-tools count` needs to be aggregated. 
This allows the results to be merged into the metadata of a single-cell object such as a Seurat or AnnData object. 

Barcodes can be filtered based on the number of supporting UMIs (`--umi_count_filter`) and by the number of UMIs in comparison to the dominant barcode per cell (`--umi_fraction_filter`). 
E.g. if barcode a has 5 supporting UMIs in a cell and a second barcode with 2 supporting UMIs and `umi_fraction_filter` set to 0.3, `5 / 2 = 0.25 < 0.3`, so the second barcode will be discarded.  

Barcodes and UMIs are semicolon-separated if multiple barcodes were detected per cell.

Output files:
- `counts/<sample_id>_cell-barcode-anno.tsv`: aggregated barcode counts per cell with cell barcode as row index and barcode and UMI count as columns
- `counts/<sample_id>_barcodes_per_cell.pdf`: QC plot, number of detected barcode per cell
- `counts/<sample_id>_UMIs_per_bc.pdf`: QC plot, UMIs supporting the most frequent barcode per cell
- `counts/<sample_id>_avg_sequence_length.pdf`: QC plot, average mapped sequence length per barcode
- `counts/<sample_id>_barcodes_per_cell_filtered.pdf`: QC plot, number of detected barcode per cell
- `counts/<sample_id>_UMIs_per_bc_filtered.pdf`: QC plot, UMIs supporting the most frequent barcode per cell
- `counts/<sample_id>_avg_sequence_length_filtered.pdf`: QC plot, average mapped sequence length per barcode
- `counts/<sample_id>_avg_sequence_length.tsv`: average mapped sequence length per barcode as table
