process STARCODE {
    tag { "starcode on ${sample_id}" }
    label "process_medium"
    publishDir "${params.outdir}/starcode/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_starcode.txt"
    
    script:
    """
    starcode -t ${task.cpus} $reads > ${sample_id}_starcode.txt
    """
}
