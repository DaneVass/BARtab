include { SOFTWARE_CHECK } from '../modules/local/software_check'
include { FASTQC } from '../modules/local/fastqc'
// include { FILTER_READS } from '../modules/local/filter_reads'
include { UMITOOLS_WHITELIST } from '../modules/local/umitools_whitelist'
include { UMITOOLS_EXTRACT } from '../modules/local/umitools_extract'
include { BAM_TO_FASTQ } from '../modules/local/bam_to_fastq'
include { CUTADAPT_READS } from '../modules/local/cutadapt_reads'
include { STARCODE } from '../modules/local/starcode'
include { STARCODE_SC } from '../modules/local/starcode_sc'
include { TRIM_BARCODE_LENGTH } from '../modules/local/trim_barcode_length'
include { BUILD_BOWTIE_INDEX } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN } from '../modules/local/bowtie_align'
include { FILTER_ALIGNMENTS } from '../modules/local/filter_alignments'
include { RENAME_READS_BAM } from '../modules/local/rename_reads_bam'
include { RENAME_READS_SAW } from '../modules/local/rename_reads_saw'
include { SAMTOOLS } from '../modules/local/samtools'
include { UMITOOLS_COUNT } from '../modules/local/umitools_count'
include { COUNT_BARCODES_SAM } from '../modules/local/count_barcodes_sam'
include { PARSE_BARCODES_SC } from '../modules/local/parse_barcodes_sc'
include { MULTIQC } from '../modules/local/multiqc'

workflow SINGLE_CELL {

    main:

        if (params.input_type == "bam") {
            readsChannel = Channel.fromPath( "${params.indir}/*.bam" )
                // creates the sample name
                // baseName already removes file type extension
                .map { file -> tuple( file.baseName.replaceAll(/\.bam/, ''), file ) }
                .ifEmpty { error "Cannot find any *.bam files in: ${params.indir}" }
        } else if (params.input_type == "fastq" & params.pipeline == "saw") {
            readsChannel = Channel.fromPath( "${params.indir}/*.fq.gz" )
                // baseName removes .gz
                .map { file -> tuple( file.baseName.replaceAll(/\.fq/, ''), file ) }
                .ifEmpty { error "Cannot find any *.fq.gz files in: ${params.indir}" }
        } else {
            readsChannel = Channel.fromFilePairs( "${params.indir}/*_R{1,2}*.{fastq,fq}.gz" )
                .ifEmpty { error "Cannot find any *_R{1,2}*.{fastq,fq}.gz files in: ${params.indir}" }
        }

        if (params.whitelist_indir) {
            // read in whitelists if provided
            whitelistChannel = Channel.fromPath( "${params.whitelist_indir}/*_whitelist.tsv" )
                // creates the sample name, baseName removes .tsv
                .map { file -> tuple( file.baseName.replaceAll(/_whitelist$/, ''), file ) }
                .ifEmpty { error "Cannot find any *_whitelist.tsv files in: ${params.whitelist_indir}" }
        }

        if (params.ref) {
            reference = file(params.ref)
        }

        params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"

        multiqcConfig = Channel.fromPath(params.multiqc_config, checkIfExists: true)

        output = Channel.fromPath( params.outdir, type: 'dir', relative: true)

        SOFTWARE_CHECK()

        if (params.input_type == "fastq" & params.pipeline == "saw") {
            r2_fastq = readsChannel

        } else if (params.input_type == "fastq") {
            FASTQC(readsChannel)

            // filtering reads for quality
            // TODO

            // use provided whitelist of cell barcodes (e.g. cellranger) or generate a whitelist with provided number of cells
            whitelist = (params.whitelist_indir ? whitelistChannel : UMITOOLS_WHITELIST(readsChannel).whitelist)

            // extract reads with whitelisted cell barcode from fastq input
            r2_fastq = UMITOOLS_EXTRACT(readsChannel.combine(whitelist, by: 0)).reads

        } else if (params.input_type == "bam") {

            // extract unmapped reads with cell barcode and UMI and convert to fastq
            r2_fastq = BAM_TO_FASTQ(readsChannel).reads
        }

        trimmed_reads = CUTADAPT_READS(r2_fastq).reads

        if (params.input_type == "fastq" & params.pipeline == "saw") {
            trimmed_reads = RENAME_READS_SAW(CUTADAPT_READS.out.reads)
        }

        if (params.ref) {

            bowtie_index = BUILD_BOWTIE_INDEX(reference)
            BOWTIE_ALIGN(bowtie_index, trimmed_reads)

            // cluster unmapped reads
            // "true" indicates that starcode is running on unmapped reads, will indicate this in output file name
            if (params.cluster_unmapped) {
                unmapped_reads = BOWTIE_ALIGN.out.unmapped_reads
                // trim barcodes to same length if only one adapter and stagger (see same in ref-free workflow)
                if (params.constants == "up" | params.constants == "down") {
                    unmapped_reads = TRIM_BARCODE_LENGTH(unmapped_reads).reads
                } else if (params.constants == "all") {
                    // not implemented
                    // i.e. clustering unmapped barcodes driectly from scRNA-seq data (bam files) is not possible. 
                    error "Error: this function has not been implemented. Please contact henrietta.holze[at]petermac.org"
                }
                STARCODE(unmapped_reads, true)
            }

            // filter alignments if barcode has fixed length
            mapped_reads = params.barcode_length ? FILTER_ALIGNMENTS(BOWTIE_ALIGN.out.mapped_reads) : BOWTIE_ALIGN.out.mapped_reads

            if (params.input_type == "bam") {
                // add CB and UMI info in header
                mapped_reads = RENAME_READS_BAM(mapped_reads.combine(readsChannel, by: 0))
            }

            if (params.pipeline == "saw") {
                // count barcodes from sam file
                counts = COUNT_BARCODES_SAM(mapped_reads)
            } else {
                SAMTOOLS(mapped_reads)

                counts = UMITOOLS_COUNT(SAMTOOLS.out).counts
            }

            // pass SAM file from mapping for mapped read length QC figures
            PARSE_BARCODES_SC(counts.combine(mapped_reads, by: 0))
        } else {
            if (params.constants == "both") {
                trimmed_reads = trimmed_reads
            } else if (params.constants == "up" | params.constants == "down") {
                // trim reads to same length (min_readlength) befor running starcode
                // this is only necessary if only one adapter was trimmed and the difference in barcode length is due to a stagger
                // and not sequencing errors (indels)
                TRIM_BARCODE_LENGTH(trimmed_reads)

                trimmed_reads = TRIM_BARCODE_LENGTH.out.reads
            } else {
                // not implemented
                error "Error: this function has not been implemented. Please contact henrietta.holze[at]petermac.org"
            }

            counts = STARCODE_SC(trimmed_reads).counts

            // place holder empty file instead of SAM file from bowtie mapping
            PARSE_BARCODES_SC(counts.combine(Channel.of("$projectDir/assets/NO_FILE")))
        }

        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC(multiqcConfig, output, PARSE_BARCODES_SC.out.counts)
}
