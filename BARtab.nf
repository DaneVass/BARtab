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
// Pipeline Config
//--------------------------------------------------------------------------------------

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

if (params.mode == "single-bulk") {
  include { single_bulk } from './workflows/single-bulk'
  println "Running single-end bulk workflow"
  println ""
    
  workflow {
    single_bulk ()
  }
}

if (params.mode == "paired-bulk") {
  //include { paired-bulk } from './workflows/paired-bulk'
  println "Running paired-end bulk workflow"
  println ""

  //workflow paired-bulk {
  //  paired-bulk ()
  //}
}

if (params.mode == "single-cell") {
  //include { single-cell } from './workflows/single-cell'
  println "Running single-cell workflow"
  println ""

  //workflow single-cell {
  //  single-cell ()
  //}
}

