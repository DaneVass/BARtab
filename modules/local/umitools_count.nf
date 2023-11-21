process UMITOOLS_COUNT {
    tag "$sample_id"
    label "process_low_sc"

    input:
        // index file needs to be linked to work directory
        tuple val(sample_id), path(bam), path(bai)

    output:
        tuple val(sample_id), path("${sample_id}.counts.tsv"), emit: counts
        path("${sample_id}_count.log"), emit: log
        
    // if pipeline is saw, use different delimiter
    script:
        def delim = (params.pipeline == "saw") ? "|||" : "_"
        """
        umi_tools count \\
        --per-contig --per-cell \\
        --umi-separator='${delim}' \\
        --edit-distance-threshold=${params.umi_dist} \\
        --random-seed=10101 \\
        -I ${bam} \\
        -S ${sample_id}.counts.tsv \\
        -L ${sample_id}_count.log
        """  
}