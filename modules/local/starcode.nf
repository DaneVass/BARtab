process STARCODE {
    tag { "starcode on ${sample_id}" }
    label "process_medium"
    publishDir "${params.outdir}/starcode/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}_starcode.tsv", emit: counts
    path "${sample_id}_starcode.log", emit: log
    
    script:
    """
    gunzip -c $reads > reads.fastq
    starcode -t ${task.cpus} reads.fastq -o ${sample_id}_starcode.tsv &> ${sample_id}_starcode.log
    rm reads.fastq
    """
}
