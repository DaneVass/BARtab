include { SOFTWARE_CHECK                                        } from '../modules/local/software_check'
include { FASTQC                                                } from '../modules/local/fastqc'
include { UMITOOLS_WHITELIST                                    } from '../modules/local/umitools_whitelist'
include { UMITOOLS_EXTRACT                                      } from '../modules/local/umitools_extract'
include { FILTER_READS                                          } from '../modules/local/filter_reads'
include { BAM_TO_FASTQ                                          } from '../modules/local/bam_to_fastq'
include { CUTADAPT_READS                                        } from '../modules/local/cutadapt_reads'
include { STARCODE_SC                                           } from '../modules/local/starcode_sc'
include { STARCODE_SC as STARCODE_SC_UNMAPPED                   } from '../modules/local/starcode_sc'
include { TRIM_BARCODE_LENGTH                                   } from '../modules/local/trim_barcode_length'
include { BUILD_BOWTIE_INDEX                                    } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN                                          } from '../modules/local/bowtie_align'
include { FILTER_ALIGNMENTS                                     } from '../modules/local/filter_alignments'
include { RENAME_READS_BAM                                      } from '../modules/local/rename_reads_bam'
include { RENAME_READS_SAW                                      } from '../modules/local/rename_reads_saw'
include { REMOVE_PCR_CHIMERISM                                  } from '../modules/local/remove_pcr_chimerism'
include { REMOVE_PCR_CHIMERISM as REMOVE_PCR_CHIMERISM_UNMAPPED } from '../modules/local/remove_pcr_chimerism'
include { UMITOOLS_COUNT                                        } from '../modules/local/umitools_count'
include { UMITOOLS_COUNT as UMITOOLS_COUNT_UNMAPPED             } from '../modules/local/umitools_count'
include { COUNT_BARCODES_SAM                                    } from '../modules/local/count_barcodes_sam'
include { PARSE_BARCODES_SC                                     } from '../modules/local/parse_barcodes_sc'
include { MULTIQC                                               } from '../modules/local/multiqc'

workflow SINGLE_CELL {

    main:

        ///////////////////
        // reading files //
        ///////////////////

        if ( params.input_type == "bam" ) {
            readsChannel = Channel.fromPath ( "${params.indir}/*.bam" )
                // creates the sample name
                // baseName already removes file type extension
                .map { file -> tuple ( file.baseName.replaceAll( /\.bam/, '' ), file ) }
                .ifEmpty { error "Cannot find any *.bam files in: ${params.indir}" }
        } else if ( params.input_type == "fastq" & params.pipeline == "saw" ) {
            readsChannel = Channel.fromPath ( "${params.indir}/*.fq.gz" )
                // baseName removes .gz
                .map { file -> tuple ( file.baseName.replaceAll( /\.fq/, '' ), file ) }
                .ifEmpty { error "Cannot find any *.fq.gz files in: ${params.indir}" }
        } else {
            readsChannel = Channel.fromFilePairs ( "${params.indir}/*_R{1,2}*.{fastq,fq}.gz" )
                .ifEmpty { error "Cannot find any *_R{1,2}*.{fastq,fq}.gz files in: ${params.indir}" }
        }

        if ( params.whitelist_indir ) {
            // read in whitelists if provided
            whitelistChannel = Channel.fromPath ( "${params.whitelist_indir}/*_whitelist.tsv" )
                // creates the sample name, baseName removes .tsv
                .map { file -> tuple( file.baseName.replaceAll( /_whitelist$/, '' ), file ) }
                .ifEmpty { error "Cannot find any *_whitelist.tsv files in: ${params.whitelist_indir}" }
        }

        if ( params.ref ) {
            reference = file( params.ref )
        }

        params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"
        multiqcConfig         = Channel.fromPath ( params.multiqc_config, checkIfExists: true )
        output                = Channel.fromPath ( params.outdir, type: 'dir', relative: true )

        ///////////////
        // workflows //
        ///////////////

        SOFTWARE_CHECK()

        if ( params.input_type == "fastq" & params.pipeline == "saw" ) {
            r2_fastq = readsChannel

        ///// amplicon data /////
        // extract reads with cell barcodes, filter R2

        } else if ( params.input_type == "fastq" ) {
            FASTQC( readsChannel )

            // use provided whitelist of cell barcodes (e.g. cellranger) or generate a whitelist with provided number of cells
            whitelist = ( params.whitelist_indir ? whitelistChannel : UMITOOLS_WHITELIST ( readsChannel ).whitelist )

            // extract reads with whitelisted cell barcode from fastq input
            r2_fastq = UMITOOLS_EXTRACT ( readsChannel.combine( whitelist, by: 0 ) ).reads

            // filter read2 for quality and complexity
            // filtering before umi-tools whitelist and extract would probably be slightly faster
            // skip quality filtering if minimum phred score is set to 0 and also no complexity threshold
            // fastp by default will remove reads with >5 N bases
            // Cutadapt will filter reads with any N bases in any case --max-n=0
            if ( ( params.minqual == 0 | params.pctqual == 0 ) & params.complexity_threshold == 0 ) {
                filtered_reads = r2_fastq
            } else {
                r2_fastq = FILTER_READS ( r2_fastq ).reads
            }

        } else if (params.input_type == "bam") {

            // extract unmapped reads with cell barcode and UMI and convert to fastq
            r2_fastq = BAM_TO_FASTQ ( readsChannel ).reads
        }

        trimmed_reads = CUTADAPT_READS ( r2_fastq ).reads

        if ( params.input_type == "fastq" & params.pipeline == "saw" ) {
            trimmed_reads = RENAME_READS_SAW ( CUTADAPT_READS.out.reads )
        }

        ////////// reference-based workflow //////////

        if ( params.ref ) {

            BUILD_BOWTIE_INDEX ( reference                             )
            BOWTIE_ALIGN       ( BUILD_BOWTIE_INDEX.out, trimmed_reads )

            // filter alignments if barcode has fixed length
            // barcodes must align to start or end position of the reference, not the middle
            // necessary when extracting very short barcode reads from scRNA-seq data
            mapped_reads = params.barcode_length ? FILTER_ALIGNMENTS ( BOWTIE_ALIGN.out.mapped_reads ) : BOWTIE_ALIGN.out.mapped_reads

            if ( params.input_type == "bam" ) {
                // add CB and UMI info in header
                mapped_reads = RENAME_READS_BAM ( mapped_reads.combine( readsChannel, by: 0 ) )
            }

            if ( params.pipeline == "saw" ) {
                // count barcodes from sam file
                counts = COUNT_BARCODES_SAM ( mapped_reads )
            } else {
                REMOVE_PCR_CHIMERISM    ( mapped_reads, "sam", false                      )
                counts = UMITOOLS_COUNT ( REMOVE_PCR_CHIMERISM.out.barcodes, false, false ).counts
            }

            // pass SAM file from mapping for mapped read length QC figures
            PARSE_BARCODES_SC ( counts.combine ( mapped_reads, by: 0 ) )

            ///// cluster unmapped reads /////

            // cluster unmapped reads
            // "true" indicates that starcode is running on unmapped reads, will indicate this in output file name
            if ( params.cluster_unmapped ) {
                unmapped_reads = BOWTIE_ALIGN.out.unmapped_reads
                // trim barcodes to same length if only one adapter and stagger (see same in ref-free workflow)
                if ( params.constants == "up" | params.constants == "down" ) {
                    unmapped_reads = TRIM_BARCODE_LENGTH ( unmapped_reads ).reads
                } else if ( params.constants == "all" ) {
                    // not implemented
                    // i.e. clustering unmapped barcodes driectly from scRNA-seq data (bam files) is not possible. 
                    error "Error: this function has not been implemented. Please contact henrietta.holze[at]petermac.org"
                }
                STARCODE_SC_UNMAPPED          ( unmapped_reads, true                                    )
                REMOVE_PCR_CHIMERISM_UNMAPPED ( STARCODE_SC_UNMAPPED.out.barcodes, "starcode_umi", true )
                UMITOOLS_COUNT_UNMAPPED       ( REMOVE_PCR_CHIMERISM_UNMAPPED.out.barcodes, true, true  )
            }

        ////////// reference-free workflow //////////

        } else {
            // no trimming required if adapters are trimmed from both ends
            if ( params.constants == "up" | params.constants == "down" ) {
                // trim reads to same length (min_readlength) befor running starcode
                // this is only necessary if only one adapter was trimmed and the difference in barcode length is due to a stagger
                // and not sequencing errors (indels)
                trimmed_reads = TRIM_BARCODE_LENGTH ( trimmed_reads ).reads
            } else if ( params.constants == "all" ) {
                // not implemented
                error "Error: this function has not been implemented. Please contact henrietta.holze[at]petermac.org"
            }

            STARCODE_SC             ( trimmed_reads, false                            )
            REMOVE_PCR_CHIMERISM    ( STARCODE_SC.out.barcodes, "starcode_umi", false )
            counts = UMITOOLS_COUNT ( REMOVE_PCR_CHIMERISM.out.barcodes, false, true  ).counts

            // place holder empty file instead of SAM file from bowtie mapping
            PARSE_BARCODES_SC ( counts.combine ( Channel.of( "$projectDir/assets/NO_FILE" ) ) )
        }

        //// report ////

        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC ( multiqcConfig, output, PARSE_BARCODES_SC.out.counts )
}
