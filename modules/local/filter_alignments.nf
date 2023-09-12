// Filters for alignments in the beginning or end of a barcode
// only works for barcodes with fixed length
process FILTER_ALIGNMENTS {
    tag "$sample_id"
    label "process_low"
    publishDir "${params.outdir}/mapped_reads/", mode: 'symlink'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.mapped_filtered.sam"), emit: bam
    
    script:
        """
        cat <(grep '^@' ${reads}) \
            <(grep -v '^@' ${reads} |\
            awk -v bc_len="${params.barcode_length}" -F'\t' '\$4 == 1 || length(\$10) + \$4 == bc_len + 1') \
        > ${sample_id}.mapped_filtered.sam
        """
}