process STARCODE_SC {
    tag "$sample_id"
    label "process_medium"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_starcode_counts.tsv"), emit: counts
        path "${sample_id}_starcode.log", emit: log
    
    script:
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
            --umi-cluster "s" \\
            --seq-d ${params.cluster_distance} \\
            --seq-cluster-ratio ${params.cluster_ratio} \\
            --seq-trim 0 \\
            ${sample_id}_trimmed.fasta \\
            > ${sample_id}_starcode.tsv \\
            2> ${sample_id}_starcode.log

        # split starcode results into cell barcode and lineage barcode
        # discard UMI since it is no longer needed
        # count cell barcode lineage barcode combinations to get UMI count per barcode per cell
        # reorder columns and insert header

        cut -f1 ${sample_id}_starcode.tsv |\\
            cut -c 1-\$cb_length,\$((cb_umi_length+1))-1000 --output-delimiter \$'\t' |\\
            sort |\\
            uniq -c |\\
            sed 's/^ *//g;s/ /\\t/g' |\\
            awk '{OFS="\\t"; print \$3,\$2,\$1}' |\\
            sed '1 i\\gene\\tcell\\tcount' > ${sample_id}_starcode_counts.tsv
        """
}
