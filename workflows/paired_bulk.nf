include { SOFTWARE_CHECK } from '../modules/local/software_check'
include { FASTQC } from '../modules/local/fastqc'
include { MERGE_READS } from '../modules/local/merge_reads'
include { GUNZIP_READS_PE } from '../modules/local/gunzip_reads_pe'
include { FILTER_READS } from '../modules/local/filter_reads'
include { CUTADAPT_READS } from '../modules/local/cutadapt_reads'
include { BUILD_BOWTIE_INDEX } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN } from '../modules/local/bowtie_align'
include { SAMTOOLS } from '../modules/local/samtools'
include { GET_BARCODE_COUNTS } from '../modules/local/get_barcode_counts'
include { COMBINE_BARCODE_COUNTS } from '../modules/local/combine_barcode_counts'
include { MULTIQC } from '../modules/local/multiqc'


workflow PAIRED_BULK {
    
    // take:
    //     readsChannel
    //     reference
    //     multiqcConfig
    //     output

    main:

        readsChannel = Channel
            .fromFilePairs( "${params.indir}/*_R{1,2}.{fastq,fq}.gz" )
            .ifEmpty { error "Cannot find any *_R{1,2}.{fastq,fq}.gz files in: ${params.reads}" }

        readsChannel.view { "file: $it" }

        reference = path(params.ref)

        // 10_multiqc_report
        params.multiqc_config = "$baseDir/config/multiqc_config.yaml"

        Channel
            .fromPath(params.multiqc_config, checkIfExists: true)
            .set { multiqcConfig }

        Channel
            .fromPath( params.outdir, type: 'dir', relative: true)
            .set { output }

        SOFTWARE_CHECK()

        FASTQC(readsChannel)

        MERGE_READS(readsChannel)
        GUNZIP_READS_PE(MERGE_READS.out.merged_reads)
        FILTER_READS(GUNZIP_READS_PE.out)
        CUTADAPT_READS(FILTER_READS.out.reads)

        bowtie_index = BUILD_BOWTIE_INDEX(reference)
        BOWTIE_ALIGN(bowtie_index, CUTADAPT_READS.out.reads)

        SAMTOOLS(BOWTIE_ALIGN.out.mapped_reads)
        GET_BARCODE_COUNTS(SAMTOOLS.out.bam)

        COMBINE_BARCODE_COUNTS(GET_BARCODE_COUNTS.out.collect())

        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC(multiqcConfig, output, COMBINE_BARCODE_COUNTS.out) 
}
