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

        if (params.input_type == "fastq" & params.pipeline == "saw") {
            r2_fastq = reads_only

        } else 
        if (params.input_type == "fastq") {

            // Generate QC report for read quality of each sample
            FASTQC(reads_only)

            // Make a tuple with [ sample, [ read1, read2 ], cellnumber ] to pass to whitelist
            reads_cellnumber = INPUT_CHECK.out
                .map{it ->
                [ it[0], it[1], it[4] ]}

            // Make another tuple with [ sample, whitelist ]
            reads_whitelist = INPUT_CHECK.out
                .map{it ->
                [ it[0], it[3] ]}

            // Extract valid cells from sample i.e. create whitelist of cell barcodes
            // will only run if cellnumber is provided
            UMITOOLS_WHITELIST(reads_cellnumber)

            // For each sample, either take provided whitelist or generated whitelist
            // Combine output of whitelist and input by sample id
            // creates [ sample, final_whitelist, [ read1, read2 ] ]
            reads_whitelist_joined = reads_whitelist
                .concat(UMITOOLS_WHITELIST.out.whitelist)
                .groupTuple(by: 0)
                .map{it ->
                [ it[0], it[1][0] ? it[1][0] : it[1][1] ]}
                .combine(reads_only, by: 0)

            // Extract reads with whitelisted cell barcode from fastq input
            // Created [ sample, extracted_read2 ]
            r2_fastq = UMITOOLS_EXTRACT(reads_whitelist_joined).reads

        } else if (params.input_type == "bam") {

            // extract reads with cell barcode and UMI and convert to fastq
            r2_fastq = PROCESS_BAM(reads_only).reads
        }

        // Identify reads with constant region and trim it
        trimmed_reads = CUTADAPT_READS(r2_fastq).reads

        if (params.input_type == "fastq" & params.pipeline == "saw") {
            trimmed_reads = RENAME_READS_SAW(CUTADAPT_READS.out.reads)
        }

        // Build bowtie index for all references provided. Reference MUST have UNIQUE FILE NAME. Path is not considered when merging back with samples. 
        // TODO check if index is provided, uniq references without index, build index, merge with samples based on reference
        // TODO index does not seem to be cached when running with bulk workflow
        // This has some loop holes, like if different indexes have been build for the same reference or if an index is specified only for some of the samples of that reference
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
        
        // Filter barcodes and aggregate barcodes per cell
        PARSE_BARCODES_SC(counts.combine(mapped_reads, by: 0))

        // Aggregate log files with MultiQC across all samples
        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC(multiqcConfig, output, PARSE_BARCODES_SC.out.counts)
}
