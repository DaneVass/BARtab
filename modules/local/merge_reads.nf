process MERGE_READS {
    tag "$sample_id"
    label "process_low_bulk"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.merged.fastq.gz"), emit: merged_reads
        path "${sample_id}.notCombined_1.fastq.gz"
        path "${sample_id}.notCombined_2.fastq.gz"
        path "${sample_id}.merge.log", emit: log
        path "${sample_id}.merge.fastp.html"
        path "${sample_id}.merge.fastp.json"

    script:
        """
        fastp --in1 ${reads[0]} --in2 ${reads[1]} \\
            --merge --correction \\
            --disable_trim_poly_g --disable_adapter_trimming --disable_quality_filtering --disable_length_filtering --dont_eval_duplication \\
            --out1 ${sample_id}.notCombined_1.fastq.gz --out2 ${sample_id}.notCombined_2.fastq.gz \\
            --merged_out ${sample_id}.merged.fastq.gz \\
            --overlap_len_require ${params.mergeoverlap} \\
            --thread ${task.cpus} \\
            2> ${sample_id}.merge.log

            mv fastp.json ${sample_id}.merge.fastp.json
            mv fastp.html ${sample_id}.merge.fastp.html
        """
}
