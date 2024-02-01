process STARCODE_SC {
    tag "$sample_id"
    label "process_high_sc"

    input:
        tuple val(sample_id), path(reads)
        val(cluster_unmapped)

    output:
        tuple val(sample_id), path("${sample_id}*_starcode_counts.tsv"), emit: counts
        path "${sample_id}*_starcode.tsv"     // the output of starcode is only interesting when looking at how many reads support a CB-UMI-barcode combination. 
        path "${sample_id}*_starcode.log", emit: log
    
    script:
        // if starcode is run on unmapped reads, that should be visible in output file name
        def unmapped = cluster_unmapped ? "_unmapped" : ""
        """
        # get length of cell barcode
        cb_length=\$(echo ${params.cb_umi_pattern} | tr -cd 'C' | wc -c)
        # get length of cell barcode and UMI
        cb_umi_length=\$(expr length ${params.cb_umi_pattern})

        # prepend cell barcode and UMI to sequence for starcode-umi
        rename_sc_reads_starcode.py ${reads} ${sample_id}_trimmed.fasta

        # cluster CB-UMI and lineage barcodes
        starcode-umi --starcode-path starcode \\
            --umi-len \$cb_umi_length \\
            --umi-d ${params.umi_dist} \\
            --umi-cluster-ratio 1 \\
            --seq-d ${params.cluster_distance} \\
            --seq-cluster-ratio ${params.cluster_ratio} \\
            --seq-trim 0 \\
            --seq-threads ${task.cpus} \\
            ${sample_id}_trimmed.fasta \\
            > ${sample_id}${unmapped}_starcode.tsv \\
            2> ${sample_id}${unmapped}_starcode.log

        # split sequence into cell barcode, UMI and lineage barcode
        # take lineage barcode with most reads for all cell barcode - UMI combinations, remove all ties
        # discard UMI since it is no longer needed# discard UMI since it is no longer needed
        # count cell barcode lineage barcode combinations to get UMI count per barcode per cell
        # reorder columns and insert header

        filter_starcode_sc.py ${sample_id}${unmapped}_starcode.tsv ${sample_id}${unmapped}_starcode_counts.tsv \$cb_umi_length \$cb_length
        """
}
