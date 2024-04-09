include { SOFTWARE_CHECK                } from '../modules/local/software_check'
include { FASTQC                        } from '../modules/local/fastqc'
include { MERGE_READS                   } from '../modules/local/merge_reads'
include { FILTER_READS                  } from '../modules/local/filter_reads'
include { CUTADAPT_READS                } from '../modules/local/cutadapt_reads'
include { STARCODE                      } from '../modules/local/starcode'
include { STARCODE as STARCODE_UNMAPPED } from '../modules/local/starcode'
include { TRIM_BARCODE_LENGTH           } from '../modules/local/trim_barcode_length'
include { BUILD_BOWTIE_INDEX            } from '../modules/local/build_bowtie_index'
include { BOWTIE_ALIGN                  } from '../modules/local/bowtie_align'
include { FILTER_ALIGNMENTS             } from '../modules/local/filter_alignments'
include { GET_BARCODE_COUNTS            } from '../modules/local/get_barcode_counts'
include { COMBINE_BARCODE_COUNTS        } from '../modules/local/combine_barcode_counts'
include { MULTIQC                       } from '../modules/local/multiqc'


workflow BULK {
    
    main:

        ///////////////////
        // reading files //
        ///////////////////
        
        if ( params.mode == "single-bulk" ) {
            readsChannel = Channel.fromPath ( ["${params.indir}/*.fq.gz", "${params.indir}/*.fastq.gz"] )
                // creates the sample name
                .map { file -> tuple ( file.baseName.replaceAll( /\.fastq|\.fq/, '' ), file ) }
                .ifEmpty { error "Cannot find any *.{fastq,fq}.gz files in: ${params.indir}" }
        } else if ( params.mode == "paired-bulk" ) {
            readsChannel = Channel.fromFilePairs ( "${params.indir}/*_R{1,2}*.{fastq,fq}.gz" )
                .ifEmpty { error "Cannot find any *_R{1,2}*.{fastq,fq}.gz files in: ${params.indir}" }
        }

        if ( params.ref ) {
            reference = file ( params.ref )
        }

        params.multiqc_config = "$baseDir/assets/multiqc_config.yaml"
        multiqcConfig         = Channel.fromPath ( params.multiqc_config, checkIfExists: true )
        output                = Channel.fromPath ( params.outdir, type: 'dir', relative: true )

        ///////////////
        // workflows //
        ///////////////

        SOFTWARE_CHECK()

        FASTQC( readsChannel )

        if ( params.mode == "single-bulk" ) {
            reads = readsChannel
        } else if ( params.mode == "paired-bulk" ) {
            reads = MERGE_READS ( readsChannel ).merged_reads
        }
        
        // skip quality filtering if minimum phred score is set to 0 and also no complexity threshold
        // fastp by default will remove reads with >5 N bases
        // Cutadapt will filter reads with any N bases in any case --max-n=0
        if ( ( params.minqual == 0 | params.pctqual == 0 ) & params.complexity_threshold == 0 ) {
            filtered_reads = reads
        } else {
            filtered_reads = FILTER_READS ( reads ).reads
        }

        trimmed_reads = CUTADAPT_READS ( filtered_reads ).reads

        ////////// reference-based workflow //////////

        if ( params.ref ) {

            // trim reads to same length (min_readlength)
            if ( params.trim_length == True ) {
                if ( params.constants == "up" | params.constants == "down" ) {
                    trimmed_reads = TRIM_BARCODE_LENGTH ( trimmed_reads ).reads
                } else if ( params.constants == "all" ) {
                    // not implemented
                    error "Error: it is currently not possible to trim barcodes to the same length with constants=all."
                }
            }

            BUILD_BOWTIE_INDEX ( reference                             )
            BOWTIE_ALIGN       ( BUILD_BOWTIE_INDEX.out, trimmed_reads )

            // filter alignments if barcode has fixed length
            // this checks if the barcode aligns to the either 3' or 5' end of the reference and not in the middle (which is not possible if an adapter has been trimmed)
            mapped_reads = params.barcode_length ? FILTER_ALIGNMENTS ( BOWTIE_ALIGN.out.mapped_reads ) : BOWTIE_ALIGN.out.mapped_reads

            GET_BARCODE_COUNTS     ( mapped_reads                     )
            COMBINE_BARCODE_COUNTS ( GET_BARCODE_COUNTS.out.collect() )

            ///// cluster unmapped reads /////

            // "true" indicates that starcode is running on unmapped reads, will indicate this in output file name
            if ( params.cluster_unmapped ) {
                unmapped_reads = BOWTIE_ALIGN.out.unmapped_reads
                // trim barcodes to same length if only one adapter and stagger (see same in ref-free workflow)
                if ( params.constants == "up" | params.constants == "down" ) {
                    unmapped_reads = TRIM_BARCODE_LENGTH( unmapped_reads ).reads
                } else if ( params.constants == "all" ) {
                    // not implemented
                    error "Error: this function has not been implemented. Please contact henrietta.holze[at]petermac.org"
                }
                STARCODE_UNMAPPED ( unmapped_reads, true )
            }
        } 
        
        ////////// reference-free workflow //////////

        else {
            // if reference-free, use starcode to cluster barcodes
            if ( params.constants == "up" | params.constants == "down" ) {
                // trim reads to same length (min_readlength) befor running starcode
                // this is only necessary if only one adapter was trimmed and the difference in barcode length is due to a stagger
                // and not sequencing errors (indels)
                trimmed_reads = TRIM_BARCODE_LENGTH ( trimmed_reads ).reads
            } else if ( params.constants == "all" ) {
                // not implemented
                error "Error: this function has not been implemented. Please contact henrietta.holze[at]petermac.org"
            }

            STARCODE               ( trimmed_reads, false          )
            COMBINE_BARCODE_COUNTS ( STARCODE.out.counts.collect() )
        }

        //// report ////

        // pass counts to multiqc so it waits to run until all samples are processed
        MULTIQC ( multiqcConfig, output, COMBINE_BARCODE_COUNTS.out ) 
}
