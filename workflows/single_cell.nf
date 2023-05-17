include { SOFTWARE_CHECK } from '../modules/local/software_check'
include { FASTQC } from '../modules/local/fastqc'
// include { GUNZIP_READS } from '../modules/local/gunzip_reads'
// include { FILTER_READS } from '../modules/local/filter_reads'
include { UMITOOLS_WHITELIST } from '../modules/local/umitools_whitelist'
include { UMITOOLS_EXTRACT } from '../modules/local/umitools_extract'
include { PROCESS_CR } from '../modules/local/process_cr'
include { CUTADAPT_READS } from '../modules/local/cutadapt_reads'
include { BUILD_BOWTIE_INDEX } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN } from '../modules/local/bowtie_align'
include { RENAME_READS } from '../modules/local/rename_reads'
include { SAMTOOLS } from '../modules/local/samtools'
include { UMITOOLS_COUNT } from '../modules/local/umitools_count'
include { PARSE_BARCODES_SC } from '../modules/local/parse_barcodes_sc'
include { MULTIQC } from '../modules/local/multiqc'

workflow SINGLE_CELL {

    main:

        if (params.bam) {
            readsChannel = Channel.fromPath( "${params.bam}" )
                // creates the sample name
                .map { file -> tuple( file.baseName.replaceAll(/\.bam/, ''), file ) }
                .ifEmpty { error "Cannot find file ${params.bam}" }
        } else {
            readsChannel = Channel.fromFilePairs( "${params.indir}/*_R{1,2}*.{fastq,fq}.gz" )
                .ifEmpty { error "Cannot find any *_R{1,2}.{fastq,fq}.gz files in: ${params.indir}" }
        }
        readsChannel.view { "file: $it" }

        reference = file(params.ref)

        params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"

        multiqcConfig = Channel.fromPath(params.multiqc_config, checkIfExists: true)

        output = Channel.fromPath( params.outdir, type: 'dir', relative: true)

        SOFTWARE_CHECK()

        if (!params.bam) {
            FASTQC(readsChannel)

            // filtering reads for quality

            // extract reads with cell barcode from fastq input
            UMITOOLS_WHITELIST(readsChannel)
            r2_fastq = UMITOOLS_EXTRACT(readsChannel, UMITOOLS_WHITELIST.out)
        }
        else {
            // extract reads with cell barcode and UMI and convert to fastq
            r2_fastq = PROCESS_CR(readsChannel)
        }

        CUTADAPT_READS(r2_fastq)

        bowtie_index = BUILD_BOWTIE_INDEX(reference)
        BOWTIE_ALIGN(bowtie_index, CUTADAPT_READS.out.reads)

        if (params.bam) {
            // add CB and UMI info in header
            mapped_reads = RENAME_READS(BOWTIE_ALIGN.out.mapped_reads, readsChannel)
        } else {
            mapped_reads = BOWTIE_ALIGN.out.mapped_reads
        }
        SAMTOOLS(mapped_reads)

        UMITOOLS_COUNT(SAMTOOLS.out.bam, SAMTOOLS.out.bai)

        PARSE_BARCODES_SC(UMITOOLS_COUNT.out)

        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC(multiqcConfig, output, PARSE_BARCODES_SC.out) 
}
