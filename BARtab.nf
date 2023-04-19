/**
BARtab: A Nextflow pipeline to tabulate synthetic barcode counts from NGS data.

Author: Dane Vassiliadis
Affiliation: PeterMac Cancer Epigenetics Laboratory

Processes:
- fastqc on raw reads
- [OPTIONAL] merge paired end reads using FLASh
- Quality filter reads using Fastx-toolkit
- filter barcode reads and trim 5' and/or 3' constant regions using Cutadapt
- generate bowtie index using reference fasta
- align to reference using Bowtie
- count number of reads aligning per barcode using samtools
- merge counts files for multiple samples

- [OPTIONAL] if no reference library given, derive consensus barcode repertoire using Starcode

- report metrics for individual samples

**/

nextflow.enable.dsl = 2


//--------------------------------------------------------------------------------------
// Help message
//--------------------------------------------------------------------------------------

// https://www.coolgenerator.com/ascii-text-generator Delta Corps Priest 1

def helpMessage() {
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

---------------------- Tabulate Barcode Counts in NGS ----------------------

  Usage: nextflow run BARtab.nf --indir <input dir> 
                                --outdir <output dir> 
                                --ref <path/to/reference/fasta> 
                                --mode <single-bulk | paired-bulk | single-cell> 
                                -profile local
                                --help

    Required arguments:
      --input                    Directory containing input *.fastq.gz files
      --ref                      Path to a reference fasta file for the barcode / sgRNA library.
                                        If null, reference-free workflow will be used for single-bulk and paired-bulk modes
      --mode                     Workflow to run. <single-bulk, paired-bulk, single-cell> [default = 'single-bulk']

    Read merging arguments:
      --merge                    Boolean. Merge overlapping reads? [default = FALSE]
      --mergeoverlap             Length of overlap required to merge paired-end reads [default = 10]

    Filtering arguments:
      --minqual                  Minimum PHRED quality per base [default = 20]
      --pctqual                  Percentage of bases within a read that must meet --minqual [default = 80]
      --upconstant               Sequence of upstream constant region [default = 'CGATTGACTA'] // SPLINTR 1st gen upstream constant region
      --downconstant             Sequence of downstream constant region [default = 'TGCTAATGCG'] // SPLINTR 1st gen downstream constant region
      --constants                Which constant regions flanking barcode to search for in reads <up, down, both> [default = 'up']

    Trimming arguments:
      --constantmismatches       Proportion of mismatched bases allowed in constant regions [default = 0.1]

    Mapping arguments:
      --alnmismatches            Number of allowed mismatches during reference mapping [default = 1]

    Optional arguments:
      -profile                   Configuration profile to use. Can use multiple (comma separated) [default = 'local']
                                        Available: local, singularity, slurm
      --outdir                   Output directory to place output [default = './']
      --threads                  Number of CPUs to use [default = 4]
      --email                    Direct output messages to this address [default = '']
      --help                     Print this help statement.

    Modes:
      single-bulk                single-end bulk workflow
      paired-bulk                paired-end bulk workflow
      single-cell                paired-end single-cell annotation workflow

    Profiles:
      local                      local execution
      singularity                use singularity container
      slurm                      SLURM execution 

    Author:
      Dane Vassiliadis (dane.vassiliadis@petermac.org)
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

// if --merge == true throw error because single end reads should not be merged
if (params.merge && params.mode == "single-bulk") {
  log.info("mode has been set to 'single-bulk' but params.merge is TRUE. Exiting.")
  exit 0
}

// if --merge == false throw error because paired end reads should be merged
if (!params.merge && params.mode == "paired-bulk") {    
  log.info("mode has been set to "paired-bulk" but params.merge is FALSE."+e,e)
  log.info("BARtab does not currently support paired-end bulk mode without merging. Exiting"+e,e)
  exit 0
}

//--------------------------------------------------------------------------------------
// Pipeline Config
//--------------------------------------------------------------------------------------

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
log.info "                               Version = 1.1.0 "
log.info ""


log.info "      Run parameters: "
log.info " ========================"
log.info " input directory          : ${params.indir}"
log.info " output directory         : ${params.outdir}"
log.info " reference fasta          : ${params.ref}"
log.info " upstream constant        : ${params.upconstant}"
log.info " downstream constant      : ${params.downconstant}"
log.info " constants to use         : ${params.constants}"
log.info " constant mismatches      : ${params.constantmismatches}"
log.info " alignment mismatches     : ${params.alnmismatches}"
log.info " CPU threads              : ${params.threads}"
log.info " Minimum PHRED quality    : ${params.minqual}"
log.info " Quality percentage       : ${params.pctqual}"
log.info " Email                    : ${params.email}"
log.info " Merge paired-end reads   : ${params.merge}"
log.info " Workflow                 : ${params.mode}"
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

if (params.email == "yourmail@yourdomain" || params.email == "") { 
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
