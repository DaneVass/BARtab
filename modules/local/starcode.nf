process STARCODE {
    tag { "starcode on ${sample_id}" }
    label "process_low"
    publishDir "${params.outdir}/starcode/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_starcode.txt"
    
    script:
    """
    starcode -t 8 $reads > ${sample_id}_starcode.txt
    """
}
