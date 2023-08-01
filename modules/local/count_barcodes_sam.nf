process COUNT_BARCODES_SAM {
    publishDir "${params.outdir}/counts", mode: 'copy'
    label "process_low"
    input:
        tuple val(sample_id), path(sam)

    output:
        tuple val(sample_id), path("${sample_id}.counts.tsv"), emit: counts
        
    script:
    """
    count_barcodes_sam.py ${sam} ${sample_id}.counts.tsv
    """
}