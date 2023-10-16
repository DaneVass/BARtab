process STARCODE {
    tag "$sample_id"
    label "process_medium"

    input:
        tuple val(sample_id), path(reads)

    output:
        path "${sample_id}_starcode.tsv", emit: counts
        path "${sample_id}_starcode.log", emit: log
    
    script:
        def cluster_distance = params.cluster_distance ? "-d ${params.cluster_distance}" : ""
        """
        gunzip -c $reads > reads.fastq
        starcode -t ${task.cpus} ${cluster_distance} reads.fastq -o ${sample_id}_starcode.tsv &> ${sample_id}_starcode.log
        rm reads.fastq
        """
}
