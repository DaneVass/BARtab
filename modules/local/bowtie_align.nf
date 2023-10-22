process BOWTIE_ALIGN {
    tag "$sample_id"
    label "process_medium"

    input:
        tuple val(refname), path(ref_files)
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.mapped.sam"), emit: mapped_reads
        path "${sample_id}.bowtie.log", emit: log

    script:
        """
        bowtie \\
        -x ${refname} \\
        -q ${reads} \\
        -p ${task.cpus} \\
        -v ${params.alnmismatches} \\
        --norc \\
        -t \\
        --no-unal \\
        -a --best --strata -m1 \\
        -S ${sample_id}.mapped.sam \\
        2> ${sample_id}.bowtie.log
        """
}