process BOWTIE_ALIGN {
    tag { "bowtie on ${sample_id}" }
    label "process_medium"
    publishDir "${params.outdir}/mapped_reads/", mode: 'symlink'

    input:
    tuple val(refname), val(sample_id), path(reads), path(ref_files)

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