process SAMTOOLS {
    tag "$sample_id"
    label "process_medium"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.mapped.bam"), path("${sample_id}.mapped.bam.bai")
    
    script:
        """
        samtools sort -@ ${task.cpus} ${reads} -o ${sample_id}.mapped.bam
        samtools index -@ ${task.cpus} ${sample_id}.mapped.bam ${sample_id}.mapped.bam.bai
        """
}