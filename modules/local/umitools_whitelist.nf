process UMITOOLS_WHITELIST {
    tag "$sample_id"
    label "process_medium"
    publishDir "${params.outdir}/extract", mode: 'symlink'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_whitelist.tsv"), emit: whitelist
        path("${sample_id}_whitelist.log"), emit: log

    script:
        """
        umi_tools whitelist --stdin ${reads[0]} \\
        --bc-pattern=${params.cb_umi_pattern} \\
        --set-cell-number ${params.cellnumber} \\
        -L ${sample_id}_whitelist.log \\
        -S ${sample_id}_whitelist.tsv
        """
}