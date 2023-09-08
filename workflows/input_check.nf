//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_input_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_input_channel(LinkedHashMap row) {
    
    if (params.reference != false) {
        if (!file(row.reference).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Reference file does not exist!\n${row.reference}"
        }
    }

    if (params.mode == "single-bulk") {
        if (!file(row.fastq_1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        if (params.reference != false) {
            reads_tuple = [ row.sample, file(row.fastq_1), file(row.reference) ]
        } else {
            reads_tuple = [ row.sample, file(row.fastq_1) ]
        }
    } else if (params.mode == "paired-bulk") {
        if (!file(row.fastq_1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        if (params.reference != false) {
            reads_tuple = [ row.sample, [ file(row.fastq_1), file(row.fastq_2) ], file(row.reference) ]
        } else {
            reads_tuple = [ row.sample, file(row.fastq_1), file(row.fastq_2) ]
        }
    } else if (params.mode == "single-cell") {
        if (params.input_type == "fastq") {
            if (params.pipeline == "saw") {
                if (!file(row.fastq_1).exists()) {
                    exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
                }
                reads_tuple = [ row.sample, file(row.fastq_1), file(row.reference) ]
            } else {
                if (!file(row.fastq_1).exists()) {
                    exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
                }
                if (!file(row.fastq_2).exists()) {
                    exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
                }
                if (row.whitelist) {
                    if (!file(row.whitelist).exists()) {
                        exit 1, "ERROR: Please check input samplesheet -> Whitelist file does not exist!\n${row.whitelist}"
                    }
                    whitelist = file(row.whitelist)
                } else {
                    whitelist = false
                }
                reads_tuple = [ row.sample, [ file(row.fastq_1), file(row.fastq_2) ], file(row.reference), whitelist, row.cellnumber]
            }
        } else if (params.input_type == "bam") {
            if (!file(row.bam).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.bam}"
            }
            reads_tuple = [ row.sample, file(row.bam), file(row.reference) ]
        }
    }

    return reads_tuple
}
