include { SOFTWARE_CHECK } from '../modules/local/software_check'
include { FASTQC } from '../modules/local/fastqc'
include { GUNZIP_READS } from '../modules/local/gunzip_reads'
include { FILTER_READS } from '../modules/local/filter_reads'
include { CUTADAPT_READS } from '../modules/local/cutadapt_reads'
include { STARCODE } from '../modules/local/starcode'
include { BUILD_BOWTIE_INDEX } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN } from '../modules/local/bowtie_align'
include { SAMTOOLS } from '../modules/local/samtools'
include { GET_BARCODE_COUNTS } from '../modules/local/get_barcode_counts'
include { COMBINE_BARCODE_COUNTS } from '../modules/local/combine_barcode_counts'
include { MULTIQC } from '../modules/local/multiqc'


workflow SINGLE_BULK {
    
    // take:
    //     readsChannel
    //     reference
    //     multiqcConfig
    //     output

    main:

        readsChannel = Channel
        .fromPath( ["${params.indir}/*.fq.gz", "${params.indir}/*.fastq.gz"] )
        .map { file -> tuple( file.baseName.replaceAll(/\.fastq|\.fq/, ''), file ) }
        .ifEmpty { error "Cannot find any *.{fastq,fq}.gz files in: ${params.indir}" }

        readsChannel.view { "file: $it" }

        reference = file(params.ref)

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

        GUNZIP_READS(readsChannel)
        FILTER_READS(GUNZIP_READS.out)
        CUTADAPT_READS(FILTER_READS.out.reads)

        // if reference-free, use starcode to cluster barcodes
        // STARCODE(CUTADAPT_READS.out.reads)

        bowtie_index = BUILD_BOWTIE_INDEX(reference)
        BOWTIE_ALIGN(bowtie_index, CUTADAPT_READS.out.reads)

        SAMTOOLS(BOWTIE_ALIGN.out.mapped_reads)
        GET_BARCODE_COUNTS(SAMTOOLS.out.bam)

        COMBINE_BARCODE_COUNTS(GET_BARCODE_COUNTS.out.collect())

        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC(multiqcConfig, output, COMBINE_BARCODE_COUNTS.out) 
}
