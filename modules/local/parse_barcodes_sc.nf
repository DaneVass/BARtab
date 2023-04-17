process PARSE_BARCODES_SC {
    publishDir "${params.outdir}", mode: 'symlink'
    input:
        tuple val(sample_id), path(counts)

    output:
        path "${sample_id}_cell-barcode-anno.tsv"
        path "barcodes_per_cell.pdf"
        path "UMIs_per_bc.pdf"
        
    script:
    """
    parse_barcodes.R $counts $sample_id
    """
}