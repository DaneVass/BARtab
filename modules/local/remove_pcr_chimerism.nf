process REMOVE_PCR_CHIMERISM {
    tag "$sample_id"
    label "process_low_sc"

    input:
        // index file needs to be linked to work directory
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}_mapped_barcodes.tsv"), emit: barcodes
        path("${sample_id}_pcr_chimerism.log"), emit: log

    script:
        def delim = (params.pipeline == "saw") ? "|||" : "_"
        """
        samtools view ${bam} > ${sample_id}.sam

        filter_sc_barcodes.py ${sample_id}.sam ${sample_id}_mapped_barcodes.tsv > ${sample_id}_pcr_chimerism.log
        """
}