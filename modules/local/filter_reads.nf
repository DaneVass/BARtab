process FILTER_READS{
    tag "$sample_id"
    // only implemented for bulk
    label "process_medium_bulk"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.filtered.fastq.gz"), emit: reads
        path "${sample_id}.filter.log", emit: log
        path "${sample_id}.filter.fastp.html"
        path "${sample_id}.filter.fastp.json"
    
    script:
        // if parameter is set, turn on low complexity filtering and set parameter
        def complexity_threshold = (params.complexity_threshold == 0) ? "" : "--low_complexity_filter --complexity_threshold ${params.complexity_threshold}"
        """
        fastp --in1 ${reads} --out1 ${sample_id}.filtered.fastq.gz \\
            --thread ${task.cpus} \\
            --qualified_quality_phred ${params.minqual} \\
            --unqualified_percent_limit \$(( 100 - ${params.pctqual} )) ${complexity_threshold} \\
            --dont_eval_duplication --disable_adapter_trimming --disable_trim_poly_g --disable_length_filtering \\
            --json ${sample_id}.filter.fastp.json --html ${sample_id}.filter.fastp.html \\
            2> ${sample_id}.filter.log

        # check if any reads passed filtering, else throw error to exclude sample
        if LC_ALL=C gzip -l ${sample_id}.filtered.fastq.gz | awk 'NR==2 {exit(\$2!=0)}'; then
            printf '%s\n' "No reads passed quality filtering. Consider excluding sample ${sample_id}." >&2  # write error message to stderr
            exit 1
        fi
        """
}