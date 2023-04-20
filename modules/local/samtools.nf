process SAMTOOLS {
    tag { "samtools on ${sample_id}" }
    label "process_low"
    publishDir "${params.outdir}/mapped_reads/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.mapped.bam"), emit: bam
    path "${sample_id}.mapped.bam.bai", emit: bai
    
    script:
    """
    samtools view -@ ${params.threads} -Sb ${reads} | samtools sort -@ ${params.threads} - > ${sample_id}.mapped.bam
    samtools index -@ ${params.threads} ${sample_id}.mapped.bam ${sample_id}.mapped.bam.bai
    """
}