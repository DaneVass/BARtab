process STARCODE {
    tag "$sample_id"
    label "process_medium"

    input:
        tuple val(sample_id), path(reads)
        val(cluster_unmapped)

    output:
        path "${sample_id}*_starcode.tsv", emit: counts
        path "${sample_id}*_starcode.log", emit: log
    
    script:
        // if starcode is run on unmapped reads, that should be visible in output file name
        def unmapped = cluster_unmapped ? "_unmapped" : ""
        """
        gunzip -c $reads > reads.fastq
        starcode -t ${task.cpus} -d ${params.cluster_distance} -r ${params.cluster_ratio} reads.fastq -o ${sample_id}${unmapped}_starcode.tsv &> ${sample_id}${unmapped}_starcode.log
        rm reads.fastq
        """
}
