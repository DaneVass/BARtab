/**
BARtab: A Nextflow pipeline to tabulate synthetic barcode counts from NGS data.

Author: Dane Vassiliadis, Henrietta Holze
Affiliation: Peter MacCallum Cancer Centre, Melbourne, Australia

**/

nextflow.enable.dsl = 2


//--------------------------------------------------------------------------------------
// Help message
//--------------------------------------------------------------------------------------

// https://www.coolgenerator.com/ascii-text-generator Delta Corps Priest 1

logo = """
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

def helpMessage() {
  log.info logo + """

---------------------- Tabulate Barcode Counts in NGS data ----------------------
                               Version = 1.5.0

  Usage: nextflow run danevass/bartab --indir <input dir>
                                      --outdir <output dir>
                                      --ref <path/to/reference/fasta>
                                      --mode <single-bulk | paired-bulk | single-cell>

    Input/output arguments:
      --indir                    Directory containing input *.fastq.gz files. 
                                    Must contain R1 and R2 if running in mode paired-bulk or single-cell.
                                    For single-cell mode, directory can contain BAM files.
      --input_type               Input file type, either fastq or bam, only relevant for single-cell mode [default = fastq]
      --ref                      Path to a reference fasta file for the barcode / sgRNA library.
                                    If null, reference-free workflow will be used for single-bulk and paired-bulk modes.
      --mode                     Workflow to run. <single-bulk, paired-bulk, single-cell>
      --outdir                   Output directory to place output [default = './']

    Read merging arguments:
      --mergeoverlap             Length of overlap required to merge paired-end reads [default = 10]

    Filtering arguments:
      --minqual                  Minimum PHRED quality per base [default = 20]
      --pctqual                  Percentage of bases within a read that must meet --minqual [default = 80]
      --complexity_threshold     Complexity filter [default = 0]
                                    Minimum percentage of bases that are different from their next base (base[i] != base[i+1])

    Trimming arguments:
      --constants                Which constant regions flanking barcode to search for in reads: up, down or both. 
                                 "all" runs all 3 modes and combines the results. <up, down, both, all> [default = 'up']
      --upconstant               Sequence of upstream constant region [default = 'CGATTGACTA'] // SPLINTR 1st gen upstream constant region
      --downconstant             Sequence of downstream constant region [default = 'TGCTAATGCG'] // SPLINTR 1st gen downstream constant region
      --up_coverage              Number of bases of the upstream constant that must be covered by the sequence [default = 3]
      --down_coverage            Number of bases of the downstream constant that must be covered by the sequence [default = 3]
      --constantmismatches       Proportion of mismatched bases allowed in constant regions [default = 0.1]
      --min_readlength           Minimum length of barcode sequence [default = 20]
      --barcode_length           Optional. Length of barcode if it is the same for all barcodes. 
                                    If constant regions are trimmed on both ends, reads are filtered for this length. 
                                    If either constant region is trimmed, this is the maximum sequence length. 
                                    If barcode_length is set, alignments to the middle of a barcode sequence are filtered out.

    Mapping arguments:
      --alnmismatches            Number of allowed mismatches during reference mapping [default = 2]
      --barcode_length           (see trimming arguments)
      --cluster_unmapped         Cluster unmapped reads with starcode [default = false]

    Reference-free arguments:
      --cluster_distance         Defines the Levenshtein distance for clustering lineage barcodes [default = 3].
      --cluster_ratio            Cluster ratio for message passing clustering. 
                                    A cluster of barcode sequences can absorb a smaller one only if it is at least x times bigger [default = 3].

    Sincle-cell arguments:
      --cb_umi_pattern           Cell barcode and UMI pattern on read 1, required for fastq input. 
                                    N = UMI position, C = cell barcode position [default = CCCCCCCCCCCCCCCCNNNNNNNNNNNN]
      --cellnumber               Number of cells expected in sample, only required when fastq provided. whitelist_indir and cellnumber are mutually exclusive
      --whitelist_indir          Directory that contains a cell ID whitelist for each sample <sample_id>_whitelist.tsv
      --umi_dist                 Hamming distance between UMIs to be collapsed during counting [default = 1]
      --umi_count_filter         Minimum number of UMIs per barcode per cell [default = 1]
      --umi_fraction_filter      Minimum fraction of UMIs per barcode per cell compared to dominant barcode in cell 
                                    (barcode supported by most UMIs) [default = 0.3]
      --pipeline                 To specify if input fastq files were created by SAW pipeline

    Resources:
      --max_cpus                 Maximum number of CPUs [default = 6]
      --max_memory               Maximum memory [default = "14.GB"]
      --max_time                 Maximum time [default = "40.h"]

    Optional arguments:
      -profile                   Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, singularity, docker, slurm, lsf
      --email                    Direct output messages to this address [default = '']
      --help                     Print this help statement.

    Author:
      Dane Vassiliadis (dane.vassiliadis@petermac.org)
      Henrietta Holze (henrietta.holze@petermac.org)
    """
}

//--------------------------------------------------------------------------------------
// Preflight checks
//--------------------------------------------------------------------------------------

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

if (!params.mode) {
  error "Error: please set parameter --mode <single-bulk,paired-bulk,single-cell>."
}
if (!["single-bulk", "paired-bulk", "single-cell"].contains(params.mode)) {
  error "Error: please set parameter --mode <single-bulk,paired-bulk,single-cell>."
}
if (params.input_type != "fastq" && params.input_type != "bam") {
  error "Error: please choose a valid value for --input_type <fastq,bam>."
}
if (params.mode != "single-cell" && params.input_type == "bam") {
  error "Error: bulk workflows do not accept BAM file input."
}
if (!params.indir) {
  error "Error: please provide the location of input files via the parameter indir."
}
if (!params.outdir) {
  error "Error: please specify location of output directory via parameter outdir."
}
if (!["up", "down", "both", "all"].contains(params.constants)) {
  error "Error: unsupported value for parameter constants. Choose either up, down, both or all (default up)."
}
if (params.constants == "both" && params.barcode_length && params.min_readlength) {
  println "Warning: min_readlength=${params.min_readlength} will be ignored because barcode_length=${params.barcode_length} and constants=${params.constants}. Reads will be filtered to match the exact barcode length."
}
if (params.mode == "single-cell" && params.input_type == "fastq" && params.pipeline != "saw" && !params.whitelist_indir && !params.cellnumber) {
  error "Error: Please provide either a whitelist or the expected number of cells for cell ID and UMI extraction."
}

//--------------------------------------------------------------------------------------
// Pipeline Config
//--------------------------------------------------------------------------------------

// setup run info for logging
log.info ""
log.info logo
// https://www.coolgenerator.com/ascii-text-generator Delta Corps Priest 1
log.info ""
log.info " ---------------------- Tabulate Barcode Counts in NGS ----------------------"
log.info "                               Version = 1.5.0 "
log.info ""
log.info "      Run parameters: "
log.info " ========================"
  log.info " Mode                     : ${params.mode}"
  log.info " Input directory          : ${params.indir}"
  log.info " Input type               : ${params.input_type}"
if (params.whitelist_indir) {
  log.info " Whitelist directory      : ${params.whitelist_indir}"
}
  log.info " Output directory         : ${params.outdir}"
if (params.ref) {
  log.info " Reference fasta          : ${params.ref}"
  log.info " Cluster unmapped         : ${params.cluster_unmapped}"
}
if (params.mode == "paired-bulk") {
  log.info " Merge overlap            : ${params.mergeoverlap}"
}
if (params.mode != "single-cell" || (params.input_type == "fastq" && params.pipeline != "saw")) {
  log.info " Minimum PHRED quality    : ${params.minqual}"
  log.info " Quality percentage       : ${params.pctqual}%"
  log.info " Complexity threshold     : ${params.complexity_threshold}%"
}
  log.info " Upstream constant        : ${params.upconstant}"
  log.info " Downstream constant      : ${params.downconstant}"
  log.info " Upstream coverage        : ${params.up_coverage}"
  log.info " Downstream coverage      : ${params.down_coverage}"
  log.info " Constants to use         : ${params.constants}"
  log.info " Constant mismatches      : ${params.constantmismatches}"
  log.info " Min. barcode read length : ${params.min_readlength}"
if (params.barcode_length) {
  log.info " Barcode length           : ${params.barcode_length}"
}
if (params.ref || params.cluster_unmapped) {
  log.info " Alignment mismatches     : ${params.alnmismatches}"
} else {
  log.info " Cluster distance         : ${params.cluster_distance}"
  log.info " Cluster ratio            : ${params.cluster_ratio}"
}
if (params.mode == "single-cell" && params.pipeline != "saw") {
  log.info " UMI distance             : ${params.umi_dist}"
  log.info " UMI count filter         : ${params.umi_count_filter}"
  log.info " UMI fraction filter      : ${params.umi_fraction_filter}"
}
if (params.mode == "single-cell" && params.input_type == "fastq" && params.pipeline != "saw") {
  log.info " Cell barcode UMI pattern : ${params.cb_umi_pattern}"
}
if (params.mode == "single-cell" && params.input_type == "fastq" && !params.whitelist_indir && params.pipeline != "saw") {
  log.info " Cell number              : ${params.cellnumber}"
}


  log.info " Email                    : ${params.email}"
  log.info " ========================"
  log.info ""



//--------------------------------------------------------------------------------------
// Named workflow for pipeline
//--------------------------------------------------------------------------------------

include { SINGLE_CELL } from './workflows/single_cell'
include { BULK } from './workflows/bulk'

workflow {
  if (params.mode == "single-bulk") {

    println "Running single-end bulk workflow"
    println ""

    BULK ()
  }
  else if (params.mode == "paired-bulk") {
    
    println "Running paired-end bulk workflow"
    println ""

    BULK ()
  }
  else if (params.mode == "single-cell") {
    println "Running single-cell workflow"
    println ""

    SINGLE_CELL ()
  }
}

//--------------------------------------------------------------------------------------
// Post processing
//--------------------------------------------------------------------------------------

// Mail notification

if (!params.email) { 
    log.info '\n'
}
else {
    log.info "\n"
    log.info "Sending runtime report to ${params.email}\n"

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
    
    sendMail(to: params.email, subject: "BARtab execution report", body: msg,  attach: "${params.outdir}/multiqc_report.html")
    }
}

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
