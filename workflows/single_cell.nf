include { SOFTWARE_CHECK } from '../modules/local/software_check'
include { FASTQC } from '../modules/local/fastqc'
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
include { MULTIQC } from '../modules/local/multiqc'

workflow SINGLE_CELL {

    main:

        readsChannel = Channel.fromFilePairs( "${params.indir}/*_R{1,2}.{fastq,fq}.gz" )
            .ifEmpty { error "Cannot find any *_R{1,2}.{fastq,fq}.gz files in: ${params.indir}" }

        readsChannel.view { "file: $it" }

        reference = file(params.ref)

        params.multiqc_config = "$baseDir/config/multiqc_config.yaml"

        multiqcConfig = Channel.fromPath(params.multiqc_config, checkIfExists: true)

        output = Channel.fromPath( params.outdir, type: 'dir', relative: true)

        // TODO add single-cell tools and starcode?
        SOFTWARE_CHECK()

        FASTQC(readsChannel)

        // filtering

        UMITOOLS_WHITELIST(readsChannel)
        UMITOOLS_EXTRACT(readsChannel, UMITOOLS_WHITELIST.out)

        // TODO cutadapt module needs to be adapted, merging 
        CUTADAPT_READS(UMITOOLS_EXTRACT.out)

        bowtie_index = BUILD_BOWTIE_INDEX(reference)
        BOWTIE_ALIGN(bowtie_index, CUTADAPT_READS.out.reads)

        SAMTOOLS(BOWTIE_ALIGN.out.mapped_reads)

        // cellranger input instead

        UMITOOLS_COUNT(SAMTOOLS.out.bam, SAMTOOLS.out.bai)

        // TODO does not work properly potentially 
        PARSE_BARCODES_SC(UMITOOLS_COUNT.out)

        // pass counts to multiqc so it waits to run until all samples are processed
        // MULTIQC(multiqcConfig, output, PARSE_BARCODES_SC.out) 
}
