process MERGE_READS {
    tag "$sample_id"
    label "process_medium_bulk"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.extendedFrags.fastq.gz"), emit: merged_reads
        path "${sample_id}.notCombined_1.fastq.gz"
        path "${sample_id}.notCombined_2.fastq.gz"
        path "${sample_id}.hist", emit: hist
        path "${sample_id}.histogram"
        path "${sample_id}.flash.log", emit: log

    script:
        """
        flash -z --min-overlap=${params.mergeoverlap} -t ${task.cpus} --output-prefix=${sample_id} ${reads} 2>&1 | tee "${sample_id}.flash.log"
        """
}
