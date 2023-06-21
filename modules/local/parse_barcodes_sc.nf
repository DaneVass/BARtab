process PARSE_BARCODES_SC {
    publishDir "${params.outdir}/counts", mode: 'copy'
    label "process_low"
    input:
        tuple val(sample_id), path(counts)

    output:
        path "${sample_id}_cell-barcode-anno.tsv", emit: counts
        path "${sample_id}_barcodes_per_cell.pdf"
        path "${sample_id}_UMIs_per_bc.pdf"
        
    script:
    """
    parse_barcodes.R ${counts} ${sample_id}
    """
}