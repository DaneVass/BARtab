process REMOVE_PCR_CHIMERISM {
    tag "$sample_id"
    label "process_low_sc"

    input:
        tuple val(sample_id), path(data)
        val(format)
        val(cluster_unmapped)

    output:
        tuple val(sample_id), path("${sample_id}*_barcodes.tsv"), emit: barcodes
        path("${sample_id}*_pcr_chimerism.log"), emit: log

    script:
        def unmapped = cluster_unmapped ? "_unmapped" : ""
        if(format == "sam") {
            """
            samtools view ${data} > ${sample_id}.sam
            filter_sc_barcodes.py ${sample_id}.sam ${sample_id}_mapped_barcodes.tsv sam > ${sample_id}_pcr_chimerism.log
            """
        }
        else if( format == "starcode_umi" ) {
            """
            filter_sc_barcodes.py ${data} ${sample_id}${unmapped}_starcode_barcodes.tsv starcode_umi ${params.cb_umi_pattern} > ${sample_id}${unmapped}_starcode_pcr_chimerism.log
            """
        }
}