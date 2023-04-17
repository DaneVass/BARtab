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
    samtools view -Sb ${reads} | samtools sort - > ${sample_id}.mapped.bam
    samtools index ${sample_id}.mapped.bam ${sample_id}.mapped.bam.bai
    """
}