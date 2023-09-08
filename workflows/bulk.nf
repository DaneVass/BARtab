include { SOFTWARE_CHECK } from '../modules/local/software_check'
include { FASTQC } from '../modules/local/fastqc'
include { MERGE_READS } from '../modules/local/merge_reads'
include { FILTER_READS } from '../modules/local/filter_reads'
include { CUTADAPT_READS } from '../modules/local/cutadapt_reads'
include { STARCODE } from '../modules/local/starcode'
include { BUILD_BOWTIE_INDEX } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN } from '../modules/local/bowtie_align'
include { FILTER_ALIGNMENTS } from '../modules/local/filter_alignments'
include { SAMTOOLS } from '../modules/local/samtools'
include { GET_BARCODE_COUNTS } from '../modules/local/get_barcode_counts'
include { COMBINE_BARCODE_COUNTS } from '../modules/local/combine_barcode_counts'
include { MULTIQC } from '../modules/local/multiqc'

include { INPUT_CHECK } from './input_check'

workflow BULK {
    
    main:

        ///////////////////
        // reading files //
        ///////////////////

        params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"

        multiqcConfig = Channel.fromPath(params.multiqc_config, checkIfExists: true)

        output = Channel.fromPath( params.outdir, type: 'dir', relative: true)

        ch_input = file(params.input)

        ///////////////
        // workflows //
        ///////////////

        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        INPUT_CHECK (
            ch_input
        )

        // split intput into only sample ID and reads
        // creates [ sample, [ read1, read2 ] ]
        reads_only = INPUT_CHECK.out
            .map{it ->
            [ it[0], it[1] ]}

        // create channel with only reference for indexing
        reference_only = INPUT_CHECK.out
            .map{it ->
            [ it[2] ]}

        // Print out all software versions
        SOFTWARE_CHECK()

        // Generate QC report for read quality of each sample
        FASTQC(reads_only)

        // Merge reads if paired-end
        if (params.mode == "single-bulk") {
            reads = reads_only
        } else if (params.mode == "paired-bulk") {
            MERGE_READS(reads_only)
            reads = MERGE_READS.out.merged_reads
        }
        
        // Filter out poor quality reads
        FILTER_READS(reads)

        // Identify reads with constant region and trim it
        trimmed_reads = CUTADAPT_READS(FILTER_READS.out.reads).reads

        // If reference is provided, build reference, align to reference and reads per barcode
        if (params.reference) {
            // see sc workflow for more comments
            // Build bowtie index for all references provided. Reference MUST have UNIQUE FILE NAME. Path is not considered when merging back with samples. 
            BUILD_BOWTIE_INDEX(reference_only.unique())

            // combine trimmed reads with which ref they need to be aligned to
            trimmed_reads_with_indexed_ref = trimmed_reads
                .combine(INPUT_CHECK.out, by: 0)
                // get reference file name string into position 0 to merge on
                // reference file name, sample name, trimmed reads
                .map{it -> [ it[3].fileName.asType(String), it[0], it[1]]}
                .combine(BUILD_BOWTIE_INDEX.out, by: 0)
            
            // Align trimmed reads to reference
            BOWTIE_ALIGN(trimmed_reads_with_indexed_ref)

            // filter alignments if barcode has fixed length
            mapped_reads = params.barcode_length ? FILTER_ALIGNMENTS(BOWTIE_ALIGN.out.mapped_reads) : BOWTIE_ALIGN.out.mapped_reads

            SAMTOOLS(mapped_reads)
            GET_BARCODE_COUNTS(SAMTOOLS.out)

            // Create barcode table with all samples and all detected barcodes
            combined_reads = COMBINE_BARCODE_COUNTS(GET_BARCODE_COUNTS.out.collect())
        } else {
            // If reference-free, use starcode to cluster barcodes
            STARCODE(trimmed_reads)
            combined_reads = COMBINE_BARCODE_COUNTS(STARCODE.out.counts.collect())
        }

        // Aggregate log files with MultiQC across all samples
        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC(multiqcConfig, output, combined_reads) 
}
