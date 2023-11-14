process BOWTIE_ALIGN {
    tag "$sample_id"
    label "$params.mode" == "single-cell" ? "process_medium_sc" : "process_medium_bulk"

    input:
        tuple val(refname), path(ref_files)
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.mapped.sam"), emit: mapped_reads
        tuple val(sample_id), path("${sample_id}.unmapped.fastq.gz"), emit: unmapped_reads, optional: true
        path "${sample_id}.bowtie.log", emit: log

    script:
        def unmapped = params.cluster_unmapped ? "--un ${sample_id}.unmapped.fastq" : "--no-unal"
        """
        bowtie \\
        -x ${refname} \\
        -q ${reads} \\
        -p ${task.cpus} \\
        -v ${params.alnmismatches} \\
        --norc \\
        -t ${unmapped} \\
        -a --best --strata -m1 \\
        -S ${sample_id}.mapped.sam \\
        2> ${sample_id}.bowtie.log

        if [ ${params.cluster_unmapped} ]; then
          gzip ${sample_id}.unmapped.fastq
        fi
        """
}