process FILTER_READS{
    tag "$sample_id"
    // unfortunately no parallelization implemented
    label "$params.mode" == "single-cell" ? "process_low_sc" : "process_low_bulk"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.filtered.fastq.gz"), emit: reads
        path "${sample_id}.filter.log", emit: log
    
    script:
        // Must write decompressed file, otherwise EOFError if piping gunzip output directly into fastq_quality_filter
        // fastx-toolkit is slow and should be replaced with fastp in the future. 
        """
        gunzip -c ${reads} > ${sample_id}.tmp.fastq
        fastq_quality_filter -i ${sample_id}.tmp.fastq \\
            -z -v -p ${params.pctqual} \\
            -q ${params.minqual} \\
            > ${sample_id}.filtered.fastq.gz \\
            2> ${sample_id}.filter.log

        rm ${sample_id}.tmp.fastq

        # check if any reads passed filtering, else throw error to exclude sample
        if LC_ALL=C gzip -l ${sample_id}.filtered.fastq.gz | awk 'NR==2 {exit(\$2!=0)}'; then
            printf '%s\n' "No reads passed quality filtering. Consider excluding sample ${sample_id}." >&2  # write error message to stderr
            exit 1
        fi
        """
}