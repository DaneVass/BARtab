include { SOFTWARE_CHECK } from '../modules/local/software_check'
include { FASTQC } from '../modules/local/fastqc'
// include { FILTER_READS } from '../modules/local/filter_reads'
include { UMITOOLS_WHITELIST } from '../modules/local/umitools_whitelist'
include { UMITOOLS_EXTRACT } from '../modules/local/umitools_extract'
include { PROCESS_BAM } from '../modules/local/process_bam'
include { CUTADAPT_READS } from '../modules/local/cutadapt_reads'
include { BUILD_BOWTIE_INDEX } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN } from '../modules/local/bowtie_align'
include { FILTER_ALIGNMENTS } from '../modules/local/filter_alignments'
include { RENAME_READS } from '../modules/local/rename_reads'
include { RENAME_READS_SAW } from '../modules/local/rename_reads_saw'
include { SAMTOOLS } from '../modules/local/samtools'
include { UMITOOLS_COUNT } from '../modules/local/umitools_count'
include { COUNT_BARCODES_SAM } from '../modules/local/count_barcodes_sam'
include { PARSE_BARCODES_SC } from '../modules/local/parse_barcodes_sc'
include { MULTIQC } from '../modules/local/multiqc'

include { INPUT_CHECK } from './input_check'

workflow SINGLE_CELL {

    main:

        // if (params.input_type == "bam") {
        //     readsChannel = Channel.fromPath( "${params.indir}/*.bam" )
        //         // creates the sample name
        //         // baseName already removes file type extension
        //         .map { file -> tuple( file.baseName.replaceAll(/\.bam/, ''), file ) }
        //         .ifEmpty { error "Cannot find any *.bam files in: ${params.indir}" }
        // } else if (params.input_type == "fastq" & params.pipeline == "saw") {
        //     readsChannel = Channel.fromPath( "${params.indir}/*.fq.gz" )
        //         // baseName removes .gz
        //         .map { file -> tuple( file.baseName.replaceAll(/\.fq/, ''), file ) }
        //         .ifEmpty { error "Cannot find any *.fq.gz files in: ${params.indir}" }
        // } else {
        //     readsChannel = Channel.fromFilePairs( "${params.indir}/*_R{1,2}*.{fastq,fq}.gz" )
        //         .ifEmpty { error "Cannot find any *_R{1,2}.{fastq,fq}.gz files in: ${params.indir}" }
        // }

        // if (params.whitelist_indir) {
        //     // read in whitelists if provided
        //     whitelistChannel = Channel.fromPath( "${params.whitelist_indir}/*_whitelist.tsv" )
        //         // creates the sample name, baseName removes .tsv
        //         .map { file -> tuple( file.baseName.replaceAll(/_whitelist$/, ''), file ) }
        //         .ifEmpty { error "Cannot find any *_whitelist.tsv files in: ${params.whitelist_indir}" }
        // }

        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //
        ch_input = file(params.input)
        INPUT_CHECK (
            ch_input
        )

        // INPUT_CHECK.out.view()

        // fastqc only expects read name and samples
        // this is only relevant
        reads_only = INPUT_CHECK.out
            .map{it ->
            [ it[0], it[1] ]}

        reference_only = INPUT_CHECK.out
            .map{it ->
            [ it[2] ]}

        params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"

        multiqcConfig = Channel.fromPath(params.multiqc_config, checkIfExists: true)

        output = Channel.fromPath( params.outdir, type: 'dir', relative: true)

        SOFTWARE_CHECK()

        if (params.input_type == "fastq" & params.pipeline == "saw") {
            r2_fastq = reads_only

        } else 
        if (params.input_type == "fastq") {

            // only pass sample id and read pair to fastqc
            FASTQC(reads_only)

            // Make a tuple with sample, reads and cellnumber
            reads_cellnumber = INPUT_CHECK.out
                .map{it ->
                [ it[0], it[1], it[4] ]}

            reads_whitelist = INPUT_CHECK.out
                .map{it ->
                [ it[0], it[3] ]}

            // will only run if cellnumber is provided
            UMITOOLS_WHITELIST(reads_cellnumber)
            // UMITOOLS_WHITELIST.out.whitelist.view()

            // combine output of whitelist and input by sample id
            // to either take provided whitelist or generated whitelist
            reads_whitelist_joined = reads_whitelist
                .concat(UMITOOLS_WHITELIST.out.whitelist)
                .groupTuple(by: 0)
                // .view()
                .map{it ->
                [ it[0], it[1][0] ? it[1][0] : it[1][1] ]}
                .combine(reads_only, by: 0)
                // .view()

            // extract reads with whitelisted cell barcode from fastq input
            r2_fastq = UMITOOLS_EXTRACT(reads_whitelist_joined).reads

        } else if (params.input_type == "bam") {

            // extract reads with cell barcode and UMI and convert to fastq
            r2_fastq = PROCESS_BAM(reads_only).reads
        }

        // r2_fastq.view()
        trimmed_reads = CUTADAPT_READS(r2_fastq).reads

        trimmed_reads.view()

        if (params.input_type == "fastq" & params.pipeline == "saw") {
            trimmed_reads = RENAME_READS_SAW(CUTADAPT_READS.out.reads)
        }

        // Build bowtie index for all references provided. Reference MUST have UNIQUE FILE NAME. Path is not considered when merging back with samples. 
        // TODO check if index is provided, uniq references without index, build index, merge with samples based on reference
        // This has some loop holes, like if different indexes have been build for the same reference or if an index is specified only for some of the samples of that reference
        BUILD_BOWTIE_INDEX(reference_only.unique())

        // // BUILD_BOWTIE_INDEX.out.view()
        // BUILD_BOWTIE_INDEX.out
        //     .map(it -> it[0].getClass())
        //     .view()

        // combine trimmed reads with which ref they need to be aligned to
        trimmed_reads_with_indexed_ref = trimmed_reads
            .combine(INPUT_CHECK.out, by: 0)
            .view()
            // get reference file name string into position 0 to merge on
            // reference file name, sample name, trimmed reads
            .map{it -> [ it[3].fileName.asType(String), it[0], it[1]]}
            // .view()
            .combine(BUILD_BOWTIE_INDEX.out, by: 0)
            .view()
        
        BOWTIE_ALIGN(trimmed_reads_with_indexed_ref)

        BOWTIE_ALIGN.out.mapped_reads.view()

        // filter alignments if barcode has fixed length
        mapped_reads = params.barcode_length ? FILTER_ALIGNMENTS(BOWTIE_ALIGN.out.mapped_reads) : BOWTIE_ALIGN.out.mapped_reads

        if (params.input_type == "bam") {
            // add CB and UMI info in header
            mapped_reads = RENAME_READS(mapped_reads.combine(INPUT_CHECK.out, by: 0))
        }

        if (params.pipeline == "saw") {
            // count barcodes from sam file
            counts = COUNT_BARCODES_SAM(mapped_reads)
        } else {
            SAMTOOLS(mapped_reads)

            counts = UMITOOLS_COUNT(SAMTOOLS.out).counts
        }
        
        PARSE_BARCODES_SC(counts.combine(mapped_reads, by: 0))

        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC(multiqcConfig, output, PARSE_BARCODES_SC.out.counts)
}
