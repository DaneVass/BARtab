/**
BARtab: A Nextflow pipeline to tabulate synthetic barcode counts from NGS data.

Author: Dane Vassiliadis
Affiliation: PeterMac Cancer Epigenetics Laboratory

Processes:
- fastqc
- [OPTIONAL] merge paired end reads using FLASh
- Quality filter reads using Fastx-toolkit
- filter barcode reads and trim 5' and/or 3' constant regions using Cutadapt
- align to reference barcode library using Bowtie
- [OPTIONAL] if no reference library, derive consensus barcode repertoire using Starcode
- count number of reads aligning per barcode using samtools
- merge counts files for multiple samples
- report metrics for individual samples

**/

def helpMessage() {
  log.info ""
  log.info """
*************************************************
  -- BARtab - Tabulate Barcode Counts in NGS --
*************************************************

  Description: A Nextflow pipeline to tabulate synthetic barcode counts from NGS data.

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
    
    """
    log.info ""
 }

//--------------------------------------------------------------------------------------
// Pipeline Config
//--------------------------------------------------------------------------------------

// config file
configfile: "envs/BARtab_config.yaml"

/*
// Set defaults
//process 1 file unless --n used at run time, e.g. -n 32
/params.indir = "./test_SE"
params.outdir = "$PWD"
params.index = "/Users/vassiliadisdane/Dropbox/PostDoc_PeterMac/Projects/Barcoding_project/barcode_reference_libs/mCHERRY_reference_library/bowtie/mCHERRY_barcode"
params.ref = "/Users/vassiliadisdane/Dropbox/PostDoc_PeterMac/Projects/Barcoding_project/barcode_reference_libs/mCHERRY_reference_library/mCHERRY_barcode_reference_library.fasta"
params.threads = 2
params.help = false
params.upconstant = "CGATTGACTA" // SPLINTR 1st gen upstream constant region 
params.downconstant = "TGCTAATGCG" // SPLINTR 1st gen downstream constant region 
params.alnmismatches = 1
params.merge = false
params.minqual = 20
params.pctqual = 100
params.constants = "up"
params.constantmismatches = 0.1
params.email = ""

*/

// Show help message
 if (params.help) {
     helpMessage()
     exit 0
 }

// setup run info for logging
log.info ""
log.info """
▀█████████▄     ▄████████    ▄████████     ███        ▄████████ ▀█████████▄  
  ███    ███   ███    ███   ███    ███ ▀█████████▄   ███    ███   ███    ███ 
  ███    ███   ███    ███   ███    ███    ▀███▀▀██   ███    ███   ███    ███ 
 ▄███▄▄▄██▀    ███    ███  ▄███▄▄▄▄██▀     ███   ▀   ███    ███  ▄███▄▄▄██▀  
▀▀███▀▀▀██▄  ▀███████████ ▀▀███▀▀▀▀▀       ███     ▀███████████ ▀▀███▀▀▀██▄  
  ███    ██▄   ███    ███ ▀███████████     ███       ███    ███   ███    ██▄ 
  ███    ███   ███    ███   ███    ███     ███       ███    ███   ███    ███ 
▄█████████▀    ███    █▀    ███    ███    ▄████▀     ███    █▀  ▄█████████▀  
                            ███    ███                                       
"""
// https://www.coolgenerator.com/ascii-text-generator Delta Corps Priest 1
log.info ""
log.info " ---------------------- Tabulate Barcode Counts in NGS ----------------------"
log.info ""

log.info "    Run parameters: "
log.info " ======================"
log.info " input directory          : ${params.indir}"
log.info " output directory         : ${params.outdir}"
log.info " reference index          : ${params.index}"
log.info " reference fasta          : ${params.ref}"
log.info " upstream constant        : ${params.upconstant}"
log.info " downstream constant      : ${params.downconstant}"
log.info " constants to use         : ${params.constants}"
log.info " constant mismatches      : ${params.constantmismatches}"
log.info " alignment mismatches     : ${params.alnmismatches}"
log.info " CPU threads              : ${params.threads}"
log.info " Minimum PHRED quality    : ${params.minqual}"
log.info " Quality percentage       : ${params.pctqual}"
log.info " ======================"
log.info ""

// check and report on software versions used in the pipeline run
/* process software_check {
  label 'software_check'

  publishDir params.outdir

  output:
    "${params.outdir}/software_check.txt"

  script:
  """
  
  echo "Nextflow version:" 
  nextflow -v 

  echo "fastQC version:" 
  fastqc --version 

  echo "FLASh version:" 
  flash -v 

  echo "fastx-toolkit - fastq_quality_filter version:" 
  fastq_quality_filter -h
  
  echo "bowtie -version" 
  bowtie --version 
  
  echo "cutadapt version:" 
  cutadapt --version 

  echo "starcode version:"
  starcode -v

  echo "R version:" 
  R --version 
   
  """
} */

//--------------------------------------------------------------------------------------
// Main Pipeline 
//--------------------------------------------------------------------------------------

//process 1 file unless --n used at run time, e.g. -n 32
params.n = 1 

// 00 Process inputs

// assume single end reads by default
// if --merge == true assume paired end reads and merge prior to filtering

if(params.merge) {
    println "Assuming paired-end reads"
    println ""
    Channel
      .fromFilePairs( "${params.indir}/*_R{1,2}.{fastq,fq}.gz" )
      .ifEmpty { error "Cannot find any *_R{1,2}.{fastq,fq}.gz files in: ${params.reads}" }
      .into { readsChannel; qcChannel }
}
else {
    println "Assuming single-end reads"
    println ""
    Channel
      .fromPath( "${params.indir}/*.{fastq,fq}.gz" )
      .map { file -> tuple( file.baseName.replaceAll(/\.fastq|\.fq/, ''), file ) }
      .ifEmpty { error "Cannot find any *.{fastq,fq}.gz files in: ${params.indir}" }
      .into { readsChannel; qcChannel }
}

//readsChannel.view { "file: $it" }
//qcChannel.view { "file: $it" }

// 01_FastQC
process fastqc{
  tag "FastQC on $sample_id"
  publishDir "${params.outdir}/qc/$sample_id", mode: 'copy'
  
  input:
    tuple val(sample_id), path(reads) from qcChannel

  output:
    tuple file('*.html'), file('*.zip') into ch_out_fastqc
  
  script:

  """
  fastqc --threads ${params.threads} ${reads}
  """
}

 
//02_Read_merging
if(params.merge) {
    process merge_reads{
      tag "FLASh on $sample_id"
      publishDir "${params.outdir}/merged_reads/$sample_id", mode: 'symlink'

      input:
        tuple val(sample_id), path(reads) from readsChannel

      output:
        set val(sample_id), file("${sample_id}.extendedFrags.fastq.gz") into mergedReadsChannel
        file "${sample_id}.notCombined_1.fastq.gz"
        file "${sample_id}.notCombined_2.fastq.gz"
        file "${sample_id}.hist"
        file "${sample_id}.histogram"

      script: 
      """
      flash -z --min-overlap=10 -t ${params.threads} --output-prefix=${sample_id} ${reads}
      """
    }
//mergedReadsChannel.view { "file: $it" }
}

// 03_gunzip_files
if(params.merge) {
  process gunzip_reads_PE{
    tag "gunzip -f on $sample_id"
  
    input:
      tuple val(sample_id), path(reads) from mergedReadsChannel

    output:
      set val(sample_id), file("${sample_id}.extendedFrags.fastq") into prefilterReadsChannel
  
    script:
    """
    gunzip -f ${reads}
    """
    }
  }
else {
  process gunzip_reads_SE{
    tag "gunzip -f on $sample_id"

    input:
      tuple val(sample_id), path(reads) from readsChannel

    output:
      set val(sample_id), file("*.fastq") into prefilterReadsChannel
  
    script:
    """
    gunzip -f ${reads}
    """
  }
}

// 04_Filter_reads_on_quality
if(params.merge) {
  process filter_reads_PE{
    tag "fastx-toolkit fastq_filter_quality on $sample_id"
    publishDir "${params.outdir}/filtered_reads/", mode: 'symlink'

    input:
      tuple val(sample_id), path(reads) from prefilterReadsChannel

    output:
      set val(sample_id), file("${sample_id}.filtered.fastq.gz") into filteredReadsChannel
  
    script:
    """
    fastq_quality_filter -z -v -p ${params.pctqual} -q ${params.minqual} -i ${reads} -o "${sample_id}.filtered.fastq.gz"
    """
  }
  //filteredReadsChannel.view { "file: $it" }
}
else {
  process filter_reads_SE{
    tag "fastx-toolkit fastq_filter_quality on $sample_id"
    publishDir "${params.outdir}/filtered_reads/", mode: 'symlink'

    input:
      tuple val(sample_id), path(reads) from prefilterReadsChannel

    output:
      set val(sample_id), file("${sample_id}.filtered.fastq.gz") into filteredReadsChannel
  
    script:
    """
    fastq_quality_filter -z -v -p ${params.pctqual} -q ${params.minqual} -i ${reads} -o "${sample_id}.filtered.fastq.gz"
    """ 
  }
  //filteredReadsChannel.view { "file: $it" }
}

// 05_cutadapt // use cutadapt to filter for length
process cutadapt_reads{
  tag "cutadapt on $sample_id"
  publishDir "${params.outdir}/trimmed_reads/", mode: 'symlink'

  input:
    tuple val(sample_id), path(reads) from filteredReadsChannel

  output:
    set val(sample_id), file("${sample_id}.trimmed.fastq") into trimmedReadsChannel
  
  script:

  if(params.merge)
    """
      cutadapt -g "${params.upconstant}...${params.downconstant}" ${reads} > ${sample_id}.trimmed.fastq
    """
  
  if(params.merge && params.constants == "up")
    """
      cutadapt -g "${params.upconstant}" ${reads} > ${sample_id}.trimmed.fastq
    """
  
  if(params.constants == "up")
    """
      cutadapt -g "${params.upconstant}" ${reads} > ${sample_id}.trimmed.fastq
    """
  
}

// 06_align_barcodes
process align_barcodes{
  tag "bowtie on $sample_id"
  publishDir "${params.outdir}/mapped_reads/", mode: 'symlink'

  input:
    tuple val(sample_id), file(reads) from trimmedReadsChannel

  output:
    set val(sample_id), "${sample_id}.mapped.bam" into mappedReadsChannel
    file "${sample_id}.unmapped.bam"
    file "${sample_id}.mapped.bam.bai"

  script:
  """
  bowtie -v ${params.alnmismatches} --norc -t -p ${params.threads} --sam --un ${sample_id}.unmapped.bam ${params.index} ${reads} | samtools view -Sb - | samtools sort - > ${sample_id}.mapped.bam
  samtools index ${sample_id}.mapped.bam ${sample_id}.mapped.bam.bai
  
  """
}

// 07_get_barcode_counts
process get_barcode_counts{
  tag "samtools idxstats on $sample_id"
  publishDir "${params.outdir}/counts/", mode: 'symlink'

  input:
    tuple val(sample_id), file(reads) from mappedReadsChannel

  output:
    set val(sample_id), file("${sample_id}_rawcounts.txt") into rawCountsChannel
  
  shell:
  """
  samtools idxstats ${reads} | cut -f1,3 > ${sample_id}_rawcounts.txt
  """
}

// 08_combine_barcode_counts
process combine_barcode_counts{
  publishDir "${params.outdir}/counts/", mode: 'copy'

  input:
  file(counts) from rawCountsChannel.collect()

  output:
    file "all_counts_combined.txt"
  
  script: 
  """
  Rscript $PWD/scripts/combine_counts.R $PWD all_counts_combined.txt
  """
}

/*
// 09_multiqc
process multiqc {

  tag { 'all' }

  publishDir "${params.outdir}/counts/", mode: 'copy', overwrite: 'true'

  input:
  file (fastqc: 'fastqc/*') from ch_out_fastqc.collect()
  //file (align: 'align/*') from alignResults.collect()
  //file (featurecounts: 'featureCounts/*') from featureCountsResults.collect()

  output:
  file "*multiqc_report.html" into multiqc_report

  script:
  """
  export LC_ALL=C.UTF-8
  export LANG=C.UTF-8
  multiqc -f -x *.run .
  """
}

/*
// Mail notification

if (params.email == "yourmail@yourdomain" || params.email == "") { 
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: params.email, subject: "BARtab execution", body: msg,  attach: "${outputMultiQC}/multiqc_report.html")
    }
}

workflow.onComplete {
    println "Pipeline BIOCORE@CRG vectorQC completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
*/

// Print completion messages
workflow.onComplete {
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    NC='\033[0m'
    
    
    log.info ""
    log.info " ---------------------- BARtab Pipeline has finished ----------------------"
    log.info ""
    log.info "Status:   " + (workflow.success ? "${GREEN}SUCCESS${NC}" : "${RED}ERROR${NC}")
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Pipeline runtime: ${workflow.duration}\n"
}