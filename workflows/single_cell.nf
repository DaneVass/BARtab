// include { SOFTWARE_CHECK } from '../modules/local/software_check'
// include { FASTQC } from '../modules/local/fastqc'
// include { GUNZIP_READS } from '../modules/local/gunzip_reads'
// include { FILTER_READS } from '../modules/local/filter_reads'
include { UMITOOLS_WHITELIST } from '../modules/local/umitools_whitelist'
include { UMITOOLS_EXTRACT } from '../modules/local/umitools_extract'
include { CUTADAPT_READS } from '../modules/local/cutadapt_reads'
include { BUILD_BOWTIE_INDEX } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN } from '../modules/local/bowtie_align'
include { SAMTOOLS } from '../modules/local/samtools'
include { UMITOOLS_COUNT } from '../modules/local/umitools_count'
include { PARSE_BARCODES_SC } from '../modules/local/parse_barcodes_sc'

// include { GET_BARCODE_COUNTS } from '../modules/local/get_barcode_counts'
// include { COMBINE_BARCODE_COUNTS } from '../modules/local/combine_barcode_counts'
// include { MULTIQC } from '../modules/local/multiqc'

workflow SINGLE_CELL {

    main:

        scChannel = Channel
        // TODO I think this should be R1 and R2
        // .fromPath( "${params.indir}/*_R1.{fastq,fq}.gz" )
        // .map { file -> tuple( file.baseName.replaceAll(/\.fastq|\.fq/, ''), file ) }
        // .ifEmpty { error "Cannot find any *.{fastq,fq}.gz files in: ${params.indir}" }
        .fromFilePairs( "${params.indir}/*_R{1,2}.{fastq,fq}.gz" )
        .ifEmpty { error "Cannot find any *_R{1,2}.{fastq,fq}.gz files in: ${params.reads}" }

        scChannel.view { "file: $it" }

        reference = path(params.ref)

        UMITOOLS_WHITELIST(scChannel)
        UMITOOLS_EXTRACT(scChannel, UMITOOLS_WHITELIST.out)

        // TODO cutadapt module needs to be adapted
        CUTADAPT_READS(UMITOOLS_EXTRACT.out)

        bowtie_index = BUILD_BOWTIE_INDEX(reference)
        BOWTIE_ALIGN(bowtie_index, CUTADAPT_READS.out.reads)

        SAMTOOLS(BOWTIE_ALIGN.out.mapped_reads)

        UMITOOLS_COUNT(SAMTOOLS.out.bam, SAMTOOLS.out.bai)

        PARSE_BARCODES_SC(UMITOOLS_COUNT.out)
}
