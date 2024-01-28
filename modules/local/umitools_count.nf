process UMITOOLS_COUNT {
    tag "$sample_id"
    label "process_low_sc"

    input:
        // index file needs to be linked to work directory
        tuple val(sample_id), path(barcodes)
        val(cluster_unmapped)
        val(starcode)

    output:
        tuple val(sample_id), path("${sample_id}*.counts.tsv"), emit: counts
        path("${sample_id}*_count.log"), emit: log

    script:
        def unmapped = cluster_unmapped ? "_unmapped" : ""
        def starcode = starcode ? "_starcode" : ""
        """
        umi_tools count_tab \\
        --per-cell \\
        --edit-distance-threshold=${params.umi_dist} \\
        -I ${barcodes} \\
        -S ${sample_id}${starcode}${unmapped}_umi_demultiplexed.tsv \\
        > ${sample_id}_count.log

        count_sc_barcodes.py ${sample_id}${starcode}${unmapped}_umi_demultiplexed.tsv ${sample_id}${starcode}${unmapped}.counts.tsv >> ${sample_id}${starcode}${unmapped}_count.log
        """
}
