include { SOFTWARE_CHECK } from '../modules/local/software_check'
include { FASTQC } from '../modules/local/fastqc'
include { MERGE_READS } from '../modules/local/merge_reads'
include { FILTER_READS } from '../modules/local/filter_reads'
include { CUTADAPT_READS } from '../modules/local/cutadapt_reads'
include { STARCODE } from '../modules/local/starcode'
include { BUILD_BOWTIE_INDEX } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN } from '../modules/local/bowtie_align'
include { SAMTOOLS } from '../modules/local/samtools'
include { GET_BARCODE_COUNTS } from '../modules/local/get_barcode_counts'
include { COMBINE_BARCODE_COUNTS } from '../modules/local/combine_barcode_counts'
include { MULTIQC } from '../modules/local/multiqc'


workflow BULK {
    
    main:

        ///////////////////
        // reading files //
        ///////////////////
        
        if (params.mode == "single-bulk") {
            readsChannel = Channel.fromPath( ["${params.indir}/*.fq.gz", "${params.indir}/*.fastq.gz"] )
                // creates the sample name
                .map { file -> tuple( file.baseName.replaceAll(/\.fastq|\.fq/, ''), file ) }
                .ifEmpty { error "Cannot find any *.{fastq,fq}.gz files in: ${params.indir}" }
        } else if (params.mode == "paired-bulk") {
            readsChannel = Channel.fromFilePairs( "${params.indir}/*_R{1,2}.{fastq,fq}.gz" )
                .ifEmpty { error "Cannot find any *_R{1,2}.{fastq,fq}.gz files in: ${params.indir}" }
        }

        if (params.ref) {
            reference = file(params.ref)
        }

        params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"

        multiqcConfig = Channel.fromPath(params.multiqc_config, checkIfExists: true)

        output = Channel.fromPath( params.outdir, type: 'dir', relative: true)

        ///////////////
        // workflows //
        ///////////////

        SOFTWARE_CHECK()

        FASTQC(readsChannel)

        if (params.mode == "single-bulk") {
            reads = readsChannel
        } else if (params.mode == "paired-bulk") {
            MERGE_READS(readsChannel)
            reads = MERGE_READS.out.merged_reads
        }
        
        FILTER_READS(reads)

        CUTADAPT_READS(FILTER_READS.out.reads)

        if (params.ref) {

            bowtie_index = BUILD_BOWTIE_INDEX(reference)
            BOWTIE_ALIGN(bowtie_index, CUTADAPT_READS.out.reads)

            SAMTOOLS(BOWTIE_ALIGN.out.mapped_reads)
            GET_BARCODE_COUNTS(SAMTOOLS.out)

            combined_reads = COMBINE_BARCODE_COUNTS(GET_BARCODE_COUNTS.out.collect())
        } 
        else {
            // if reference-free, use starcode to cluster barcodes
            STARCODE(CUTADAPT_READS.out.reads)
            combined_reads = COMBINE_BARCODE_COUNTS(STARCODE.out.counts.collect())
        }

        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC(multiqcConfig, output, combined_reads) 
}
